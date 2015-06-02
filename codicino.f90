      program codicino

      implicit none

!definition of the problem
      integer                 :: nx,ny,system,nsites,ndim,plane
      real*8                  :: t_h,U_h,mu
      integer, dimension(:,:),   allocatable  :: flag_t !nn
      integer, dimension(:),     allocatable  :: site,spin
      real*8, dimension(:,:),     allocatable :: RS

!****************************************************************************
!Hilbert space (nptot particles)
      integer                                :: nptot
      integer                                :: ord !basis number of elements
      integer, dimension(:,:),   allocatable :: phi !basis set
      integer, dimension(:),     allocatable :: occupation
      integer, dimension(:,:),   allocatable :: occ
      integer, dimension(:),     allocatable :: bra,ket
      integer, dimension(:),     allocatable :: flag_u !interaction
!hamiltonian & co (nptot particles)
      real*8, dimension(:,:),     allocatable :: H,U
      real*8, dimension(:),       allocatable :: epsilon

! *****************************************************************************
!Hilbert space (nptot-1 particles)
      integer                                :: nptot2
      integer                                :: ord2 !basis number of elements
      integer, dimension(:,:),   allocatable :: phi2 !basis set
      integer, dimension(:),     allocatable :: occupation2
      integer, dimension(:,:),   allocatable :: occ2
      integer, dimension(:),     allocatable :: bra2,ket2
      integer, dimension(:),     allocatable :: flag_u2 !interaction
!hamiltonian & co (nptot-1 particles)
      real*8, dimension(:,:),     allocatable :: H2,U2
      real*8, dimension(:),       allocatable :: epsilon2

! ****************************************************************************
!Hilbert space (nptot+1 particles)
      integer                                :: nptot3
      integer                                :: ord3 !basis number of elements
      integer, dimension(:,:),   allocatable :: phi3 !basis set
      integer, dimension(:),     allocatable :: occupation3
      integer, dimension(:,:),   allocatable :: occ3
      integer, dimension(:),     allocatable :: bra3,ket3
      integer, dimension(:),     allocatable :: flag_u3 !interaction
!hamiltonian & co (nptot+1 particles)
      real*8, dimension(:,:),     allocatable :: H3,U3
      real*8, dimension(:),       allocatable :: epsilon3


!observables
      integer                 :: nq,ntimes
      real*8                  :: dt
      real*8, dimension(:,:),     allocatable :: Q
      real*8, dimension(:,:),     allocatable :: G,Gp,Rho
      real*8, dimension(:,:,:),   allocatable :: Gd,Gpd,Rhod
      

!flags and labels
      integer  :: i,j,k,l,m,r,s,t,cont,i1,i2,i3,i4,i5,inn,it,idim
      integer  :: alpha,beta,gamma,eta,n,Sz
      real*8   :: time

! working stuff
      real *8   :: eps,pi,delta,fattore_bra,fattore_ket,overlap,ovrlp
      integer, dimension(:),     allocatable :: bra3a,ket3e
      parameter    (eps=1E-10)
      parameter    (pi=dacos(-1.d0))
      complex*16                              :: eiqr,temp
      real*8                                  :: dot

! lapack stuff

      real*8, dimension(:), allocatable       :: work
      integer lwork,info

!     lettura da file dei dati

      write(*,*)
      write(*,*)
      open(14,file='input.dat',status='old')
      read(14,*)ndim
      read(14,*)nx,ny
      read(14,*)nptot
      read(14,*)t_h
      read(14,*)U_h
      read(14,*)nq
      read(14,*)ntimes
      read(14,*)dt
      close(14)

      write(*,*)
      write(*,*)'******************************************************'
      write(*,*)'*    HUBBARD HAMILTONIAN EXACT DIAGONALIZATION       *'
      write(*,*)'*                     VERSION 1.b                    *'
      write(*,*)'******************************************************'

      nsites=nx*ny
      if(ndim.ne.2)then
       write(*,*)'Sorry, the code works only for d = 2'
       stop
      endif
      if(nsites.ne.4.and.nsites.ne.2)then
       write(*,*)'Sorry, the code works only for 2 and 4 sites '
       stop
      endif
      if(nsites.eq.4)then
       if(nx.eq.2.and.ny.eq.2)then
        if(nptot.ne.4)then
         write(*,*)'Sorry, the code works for 4 particles on 4 sites '
        endif
        system=1
        write(*,*)'Square 2 x 2 , periodic boundary conditions '
        t_h=2.d0*t_h !periodic boundary conditions in special 2 x 2
!sites
        allocate(RS(ndim,nptot))
        RS(1,1)=-0.5
        RS(2,1)=0.5
        RS(1,2)=0.5
        RS(2,2)=0.5
        RS(1,3)=-0.5
        RS(2,3)=-0.5
        RS(1,4)=0.5
        RS(2,4)=-0.5
       elseif(nx.eq.1.or.ny.eq.1)then
        if(nptot.ne.4)then
         write(*,*)'Sorry, the code works for 4 particles on 4 sites '
        endif
        system=2
        write(*,*)'Square 4 x 1 , periodic boundary conditions '
        allocate(RS(ndim,nptot))
        RS(1,1)=-1.5
        RS(2,1)=0.0
        RS(1,2)=-0.5
        RS(2,2)=0.0
        RS(1,3)=0.5
        RS(2,3)=0.0
        RS(1,4)=1.5
        RS(2,4)=0.0
       endif
      elseif(nsites.eq.2)then
       if(nptot.ne.2)then
        write(*,*)'Sorry, the code works for 2 particles on 2 sites '
       endif
       system=3
       allocate(RS(ndim,nptot))
       RS(1,1)=-0.5
       RS(2,1)=0.0
       RS(1,2)=0.5
       RS(2,2)=0.0
       t_h=2.d0*t_h !periodic boundary conditions in special 2 x 1
      endif


      write(*,*)'number of sites ',nsites
      write(*,*)'number of particles ',nptot
      plane=2*nx*ny
      write(*,*)'single particle basis, number of elements = ',plane

      allocate(site(plane))
      allocate(spin(plane))
! mapping between labeling of states and site/spin
      j=0
      do i=1,plane,2
       j=j+1
       site(i)=j
       site(i+1)=j
       spin(i)=1
       spin(i+1)=-1
      enddo
      write(*,*)'check '
      write(*,*)site
      write(*,*)spin

!hopping map
      allocate(flag_t(plane,plane))
      flag_t(:,:)=0.d0
      do i=1,plane
       do j=i+1,plane
        call nn(system,site(i),site(j),spin(i),spin(j),inn)
        if(inn.eq.1)then
         write(*,*)'sites ',site(i),site(j)
         write(*,*)'nn ',i,j
         write(*,*)
         flag_t(i,j)=1
         flag_t(j,i)=1
        endif
       enddo
      enddo

! Diagonalization of the nptot-problem
      write(*,*)'Diagonalization of the nptot-problem****************'
      write(*,*)


!dimension of the Hilbert space
      if(system.eq.3)then
       cont=0
       do i1=1,plane
        do i2=i1+1,plane
         Sz=spin(i1)+spin(i2)
         if(Sz.eq.0)cont=cont+1
        enddo
       enddo
      else
       cont=0
       do i1=1,plane
        do i2=i1+1,plane
         do i3=i2+1,plane
          do i4=i3+1,plane
           Sz=spin(i1)+spin(i2)+spin(i3)+spin(i4)
           if(Sz.eq.0)cont=cont+1
          enddo
         enddo
        enddo
       enddo
      endif
      ord=cont !binomial(plane,nptot)
      write(*,*)'slater determinant basis, number of elements = ',ord

      allocate(occupation(nsites))
      allocate(occ(ord,plane))
      allocate(flag_u(ord))
      allocate(phi(ord,nptot))
      allocate(H(ord,ord))
      allocate(bra(nptot))
      allocate(ket(nptot))
      occ(:,:)=0
      flag_u(:)=0
      H(:,:)=0


      open(11,file='nbody_basis.out')
! build up states
      write(11,*)'hilbert space basis: '
      if(system.eq.3)then
       cont=1
       do i1=1,plane
        do i2=i1+1,plane
           occupation(:)=0

           Sz=spin(i1)+spin(i2)
           if(Sz.eq.0)then
            occupation(:)=0

            phi(cont,1)=i1
            phi(cont,2)=i2

            occ(cont,i1)=occ(cont,i1)+1
            occ(cont,i2)=occ(cont,i2)+1
            occupation(site(i1))=occupation(site(i1))+1
            occupation(site(i2))=occupation(site(i2))+1

!          write(*,*)'phi(',cont,') -> ',i1,i2,i3,i4
            do i=1,nsites
             if(occupation(i).gt.1)then
              flag_u(cont)=flag_u(cont)+1
             endif
            enddo

            write(11,*)'*********************************************'
            write(11,*)'determinant ',cont
            write(11,*)'check flag_u '
            write(11,*)'occupation ',occupation
            write(11,*)'orbitals ',i1,i2
            write(11,*)'sites ',site(i1),site(i2)
            write(11,*)'spins ',spin(i1),spin(i2)
            write(11,*)'flag_u ',flag_u(cont)
            write(11,*)'total spin ',Sz
            write(11,*)
            write(11,*)
            write(11,*)'*********************************************'
            write(11,*)
            write(11,*)

            cont=cont+1
           endif
        enddo
       enddo
      else
       cont=1
       do i1=1,plane
        do i2=i1+1,plane
         do i3=i2+1,plane
          do i4=i3+1,plane
           occupation(:)=0

           Sz=spin(i1)+spin(i2)+spin(i3)+spin(i4)
           if(Sz.eq.0)then
            occupation(:)=0
        
            phi(cont,1)=i1
            phi(cont,2)=i2
            phi(cont,3)=i3
            phi(cont,4)=i4
 
            occ(cont,i1)=occ(cont,i1)+1
            occ(cont,i2)=occ(cont,i2)+1
            occ(cont,i3)=occ(cont,i3)+1
            occ(cont,i4)=occ(cont,i4)+1       
            occupation(site(i1))=occupation(site(i1))+1
            occupation(site(i2))=occupation(site(i2))+1
            occupation(site(i3))=occupation(site(i3))+1
            occupation(site(i4))=occupation(site(i4))+1

!          write(*,*)'phi(',cont,') -> ',i1,i2,i3,i4
            do i=1,nsites
             if(occupation(i).gt.1)then
              flag_u(cont)=flag_u(cont)+1
             endif
            enddo

            write(11,*)'*********************************************'
            write(11,*)'determinant ',cont
            write(11,*)'check flag_u '
            write(11,*)'occupation ',occupation
            write(11,*)'orbitals ',i1,i2,i3,i4
            write(11,*)'sites ',site(i1),site(i2),site(i3),site(i4)
            write(11,*)'spins ',spin(i1),spin(i2),spin(i3),spin(i4)
            write(11,*)'flag_u ',flag_u(cont)
            write(11,*)'total spin ',Sz
            write(11,*)
            write(11,*)
            write(11,*)'*********************************************'
            write(11,*)
            write(11,*)

            cont=cont+1
           endif
          enddo
         enddo
        enddo
       enddo
      endif
      close(11)

      write(*,*)'******************************************************'
      write(*,*)'Building hamiltonian matrix for the toy problem ...'
      write(*,*)'Please wait :-) '
      write(*,*)

! Kinetic energy
      do i=1,ord
       do j=1,ord
        do m=1,plane
         do t=1,plane
          if(flag_t(m,t).ne.0)then
           do k=1,nptot
            bra(k)=phi(i,k)
            ket(k)=phi(j,k)
           enddo
           fattore_bra=1.d0
           fattore_ket=1.d0
           call destroy(ket,nptot,t,plane,fattore_ket)
           if(fattore_ket.eq.0)go to 99
           call destroy(bra,nptot,m,plane,fattore_bra)
           if(fattore_bra.eq.0)go to 99
           overlap=fattore_bra*fattore_ket
           do k=1,nptot-1
            overlap=overlap*delta(bra(k),ket(k))
           enddo
           H(i,j)=H(i,j)-t_h*overlap
          endif
 99       continue
         enddo
        enddo
       enddo
      enddo

! Potential energy (diagonal in coordinates representation)
      do i=1,ord
       H(i,i)=H(i,i)+U_h*flag_u(i)
      enddo
         
      write(*,*)'matrix H ready, see h_ij.out for details'
      write(*,*)
      call print_matrix('h_ij.out',H,ord)

      write(*,*)'******************************************************'

      allocate(U(ord,ord))
      allocate(epsilon(ord))
      lwork=3*ord
      allocate(work(lwork))
      do i=1,ord
       do j=1,ord
        U(i,j)=H(i,j)
       enddo
      enddo

      call dsyev('V','U',ord,U,ord,epsilon,work,lwork,info)
      call print_array('epsilon_i.out',epsilon,ord)
      deallocate(work)
      write(*,*)'see epsilon_i.out for eigenvalues of H'
      write(*,*)'see psi_ij.out for eigenvectors of H'
      write(*,*)'******************************************************'

      open(14,file='psi_ij.out')
      write(14,*)'eigenvectors of H'
      write(14,*)'*****************************************************'
      do i=1,ord
       write(14,*)'psi(',i,') = '
       do j=1,ord
        if(abs(U(j,i)).le.eps)then
         write(14,*)'component ',j,0.d0
        else
         write(14,*)'component ',j,U(j,i)
        endif
       enddo
      enddo
      write(14,*)'*****************************************************'
      close(14)

      deallocate(occupation)
      deallocate(flag_u)
      deallocate(H)



! Diagonalization of the nptot+1-problem
      write(*,*)'Diagonalization of the nptot+1-problem***************'
      write(*,*)
      nptot3=nptot+1

!dimension of the Hilbert space
      if(system.eq.3)then
       cont=0
       do i1=1,plane
        do i2=i1+1,plane
         do i3=i2+1,plane
            Sz=spin(i1)+spin(i2)+spin(i3)
            if(Sz.eq.1)cont=cont+1
         enddo
        enddo
       enddo
      else
       cont=0
       do i1=1,plane
        do i2=i1+1,plane
         do i3=i2+1,plane
          do i4=i3+1,plane
           do i5=i4+1,plane
            Sz=spin(i1)+spin(i2)+spin(i3)+spin(i4)+spin(i5)
            if(Sz.eq.1)cont=cont+1
           enddo
          enddo
         enddo
        enddo
       enddo
      endif
      ord3=cont !binomial(plane,nptot)
      write(*,*)'slater determinant basis, number of elements = ',ord3

      allocate(occupation3(nsites))
      allocate(occ3(ord3,plane))
      allocate(flag_u3(ord3))
      allocate(phi3(ord3,nptot3))
      allocate(H3(ord3,ord3))
      allocate(bra3(nptot3))
      allocate(ket3(nptot3))
      occ3(:,:)=0
      flag_u3(:)=0
      H3(:,:)=0


! build up states
      write(*,*)'hilbert space basis: '
      if(system.eq.3)then
       cont=1
       do i1=1,plane
        do i2=i1+1,plane
         do i3=i2+1,plane
            occupation3(:)=0

            Sz=spin(i1)+spin(i2)+spin(i3)
            if(Sz.eq.1)then
             occupation3(:)=0

             phi3(cont,1)=i1
             phi3(cont,2)=i2
             phi3(cont,3)=i3

             occ3(cont,i1)=occ3(cont,i1)+1
             occ3(cont,i2)=occ3(cont,i2)+1
             occ3(cont,i3)=occ3(cont,i3)+1
             occupation3(site(i1))=occupation3(site(i1))+1
             occupation3(site(i2))=occupation3(site(i2))+1
             occupation3(site(i3))=occupation3(site(i3))+1

!          write(*,*)'phi(',cont,') -> ',i1,i2,i3,i4
             do i=1,nsites
              if(occupation3(i).gt.1)then
               flag_u3(cont)=flag_u3(cont)+1
              endif
             enddo

             write(*,*)'*********************************************'
             write(*,*)'determinant ',cont
             write(*,*)'check flag_u '
             write(*,*)'occupation ',occupation3
             write(*,*)'orbitals ',i1,i2,i3,i4,i5
             write(*,*)'sit',site(i1),site(i2),site(i3)
             write(*,*)'sp',spin(i1),spin(i2),spin(i3)
             write(*,*)'flag_u ',flag_u3(cont)
             write(*,*)'total spin ',Sz
             write(*,*)
             write(*,*)
             write(*,*)'*********************************************'
             write(*,*)
             write(*,*)

             cont=cont+1
            endif
         enddo
        enddo
       enddo
      else
       cont=1
       do i1=1,plane
        do i2=i1+1,plane
         do i3=i2+1,plane
          do i4=i3+1,plane
           do i5=i4+1,plane
            occupation3(:)=0

            Sz=spin(i1)+spin(i2)+spin(i3)+spin(i4)+spin(i5)
            if(Sz.eq.1)then
             occupation3(:)=0

             phi3(cont,1)=i1
             phi3(cont,2)=i2
             phi3(cont,3)=i3
             phi3(cont,4)=i4
             phi3(cont,5)=i5

             occ3(cont,i1)=occ3(cont,i1)+1
             occ3(cont,i2)=occ3(cont,i2)+1
             occ3(cont,i3)=occ3(cont,i3)+1
             occ3(cont,i4)=occ3(cont,i4)+1
             occ3(cont,i5)=occ3(cont,i5)+1
             occupation3(site(i1))=occupation3(site(i1))+1
             occupation3(site(i2))=occupation3(site(i2))+1
             occupation3(site(i3))=occupation3(site(i3))+1
             occupation3(site(i4))=occupation3(site(i4))+1
             occupation3(site(i5))=occupation3(site(i5))+1

!          write(*,*)'phi(',cont,') -> ',i1,i2,i3,i4
             do i=1,nsites
              if(occupation3(i).gt.1)then
               flag_u3(cont)=flag_u3(cont)+1
              endif
             enddo

             write(*,*)'*********************************************'
             write(*,*)'determinant ',cont
             write(*,*)'check flag_u '
             write(*,*)'occupation ',occupation3
             write(*,*)'orbitals ',i1,i2,i3,i4,i5
             write(*,*)'sit',site(i1),site(i2),site(i3),site(i4),site(i5)
             write(*,*)'sp ',spin(i1),spin(i2),spin(i3),spin(i4),spin(i5)
             write(*,*)'flag_u ',flag_u3(cont)
             write(*,*)'total spin ',Sz
             write(*,*)
             write(*,*)
             write(*,*)'*********************************************'
             write(*,*)
             write(*,*)

             cont=cont+1
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
      endif

      write(*,*)'******************************************************'
      write(*,*)'Building hamiltonian matrix for the toy problem ...'
      write(*,*)'Please wait :-) '
      write(*,*)

! Kinetic energy
      do i=1,ord3
       do j=1,ord3
        do m=1,plane
         do t=1,plane
          if(flag_t(m,t).ne.0)then
           do k=1,nptot3
            bra3(k)=phi3(i,k)
            ket3(k)=phi3(j,k)
           enddo
           fattore_bra=1.d0
           fattore_ket=1.d0
           call destroy(ket3,nptot3,t,plane,fattore_ket)
           if(fattore_ket.eq.0)go to 97
           call destroy(bra3,nptot3,m,plane,fattore_bra)
           if(fattore_bra.eq.0)go to 97
           overlap=fattore_bra*fattore_ket
           do k=1,nptot3-1
            overlap=overlap*delta(bra3(k),ket3(k))
           enddo
           H3(i,j)=H3(i,j)-t_h*overlap
          endif
 97       continue
         enddo
        enddo
       enddo
      enddo

      

! Potential energy (diagonal in coordinates representation)
      do i=1,ord3
       H3(i,i)=H3(i,i)+U_h*flag_u3(i)
      enddo

      write(*,*)'matrix H ready, see h_ij_n+1.out for details'
      write(*,*)
      call print_matrix('h_ij_n+1.out',H3,ord3)

      write(*,*)'******************************************************'

      allocate(U3(ord3,ord3))
      allocate(epsilon3(ord3))
      lwork=3*ord
      allocate(work(lwork))
      do i=1,ord3
       do j=1,ord3
        U3(i,j)=H3(i,j)
       enddo
      enddo

      call dsyev('V','U',ord3,U3,ord3,epsilon3,work,lwork,info)
      call print_array('epsilon_i_n+1.out',epsilon3,ord3)
      deallocate(work)
      write(*,*)'see epsilon_i_n+1.out for eigenvalues of H'
      write(*,*)'see psi_ij.out for eigenvectors of H'
      write(*,*)'******************************************************'

      open(14,file='psi_ij_n+1.out')
      write(14,*)'eigenvectors of H'
      write(14,*)'*****************************************************'
      do i=1,ord3
       write(14,*)'psi(',i,') = '
       do j=1,ord3
        if(abs(U3(j,i)).le.eps)then
         write(14,*)'component ',j,0.d0
        else
         write(14,*)'component ',j,U3(j,i)
        endif
       enddo
      enddo
      write(14,*)'*****************************************************'
      close(14)

      deallocate(occupation3)
      deallocate(occ3)
      deallocate(flag_u3)
      deallocate(H3)



! Diagonalization of the nptot-1-problem
      write(*,*)'Diagonalization of the nptot-1-problem***************'
      write(*,*)
      nptot2=nptot-1

!dimension of the Hilbert space
      if(system.eq.3)then
       cont=0
       do i1=1,plane
           Sz=spin(i1)
           if(Sz.eq.-1)cont=cont+1         !I destroy one up particle
       enddo
      else
       cont=0
       do i1=1,plane
        do i2=i1+1,plane
         do i3=i2+1,plane
           Sz=spin(i1)+spin(i2)+spin(i3)
           if(Sz.eq.-1)cont=cont+1         !I destroy one up particle
         enddo
        enddo
       enddo
      endif
      ord2=cont !binomial(plane,nptot)
      write(*,*)'slater determinant basis, number of elements = ',ord2


      allocate(occupation2(nsites))
      allocate(occ2(ord2,plane))
      allocate(flag_u2(ord2))
      allocate(phi2(ord2,nptot2))
      allocate(H2(ord2,ord2))
      allocate(bra2(nptot2))
      allocate(ket2(nptot2))
      occ2(:,:)=0
      flag_u2(:)=0
      H2(:,:)=0


! build up states
      write(*,*)'hilbert space basis: '
      if(system.eq.3)then
       cont=1
       do i1=1,plane
            occupation2(:)=0

            Sz=spin(i1)
            if(Sz.eq.-1)then
             occupation2(:)=0

             phi2(cont,1)=i1

             occ2(cont,i1)=occ2(cont,i1)+1
             occupation2(site(i1))=occupation2(site(i1))+1

!          write(*,*)'phi(',cont,') -> ',i1,i2,i3,i4
             do i=1,nsites
              if(occupation2(i).gt.1)then
               flag_u2(cont)=flag_u2(cont)+1
              endif
             enddo

             write(*,*)'*********************************************'
             write(*,*)'determinant ',cont
             write(*,*)'check flag_u '
             write(*,*)'occupation ',occupation3
             write(*,*)'orbitals ',i1
             write(*,*)'sit',site(i1)
             write(*,*)'sp ',spin(i1)
             write(*,*)'flag_u ',flag_u2(cont)
             write(*,*)'total spin ',Sz
             write(*,*)
             write(*,*)
             write(*,*)'*********************************************'
             write(*,*)
             write(*,*)

             cont=cont+1
            endif
       enddo

      else
       cont=1
       do i1=1,plane
        do i2=i1+1,plane
         do i3=i2+1,plane
            occupation2(:)=0

            Sz=spin(i1)+spin(i2)+spin(i3)
            if(Sz.eq.-1)then
             occupation2(:)=0

             phi2(cont,1)=i1
             phi2(cont,2)=i2
             phi2(cont,3)=i3

             occ2(cont,i1)=occ2(cont,i1)+1
             occ2(cont,i2)=occ2(cont,i2)+1
             occ2(cont,i3)=occ2(cont,i3)+1
             occupation2(site(i1))=occupation2(site(i1))+1
             occupation2(site(i2))=occupation2(site(i2))+1
             occupation2(site(i3))=occupation2(site(i3))+1

!          write(*,*)'phi(',cont,') -> ',i1,i2,i3,i4
             do i=1,nsites
              if(occupation2(i).gt.1)then
               flag_u2(cont)=flag_u2(cont)+1
              endif
             enddo

             write(*,*)'*********************************************'
             write(*,*)'determinant ',cont
             write(*,*)'check flag_u '
             write(*,*)'occupation ',occupation3
             write(*,*)'orbitals ',i1,i2,i3
             write(*,*)'sit',site(i1),site(i2),site(i3)
             write(*,*)'sp ',spin(i1),spin(i2),spin(i3)
             write(*,*)'flag_u ',flag_u2(cont)
             write(*,*)'total spin ',Sz
             write(*,*)
             write(*,*)
             write(*,*)'*********************************************'
             write(*,*)
             write(*,*)

             cont=cont+1
            endif
         enddo
        enddo
       enddo
      endif






      write(*,*)'******************************************************'
      write(*,*)'Building hamiltonian matrix for the toy problem ...'
      write(*,*)'Please wait :-) '
      write(*,*)

! Kinetic energy
      do i=1,ord2
       do j=1,ord2
        do m=1,plane
         do t=1,plane
          if(flag_t(m,t).ne.0)then
           do k=1,nptot2
            bra2(k)=phi2(i,k)
            ket2(k)=phi2(j,k)
           enddo
           fattore_bra=1.d0
           fattore_ket=1.d0
           call destroy(ket2,nptot2,t,plane,fattore_ket)
           if(fattore_ket.eq.0)go to 88
           call destroy(bra2,nptot2,m,plane,fattore_bra)
           if(fattore_bra.eq.0)go to 88
           overlap=fattore_bra*fattore_ket
           do k=1,nptot2
            overlap=overlap*delta(bra2(k),ket2(k))
           enddo
           H2(i,j)=H2(i,j)-t_h*overlap
          endif
 88       continue
         enddo
        enddo
       enddo
      enddo

! Potential energy (diagonal in coordinates representation)
      do i=1,ord2
       H2(i,i)=H2(i,i)+U_h*flag_u2(i)
      enddo

      write(*,*)'matrix H ready, see h_ij_n-1.out for details'
      write(*,*)
      call print_matrix('h_ij_n-1.out',H2,ord2)

      write(*,*)'******************************************************'

      allocate(U2(ord2,ord2))
      allocate(epsilon2(ord2))
      lwork=3*ord
      allocate(work(lwork))
      do i=1,ord2
       do j=1,ord2
        U2(i,j)=H2(i,j)
       enddo
      enddo

      call dsyev('V','U',ord2,U2,ord2,epsilon2,work,lwork,info)
      call print_array('epsilon_i_n-1.out',epsilon2,ord2)
      deallocate(work)
      write(*,*)'see epsilon_i_n-1.out for eigenvalues of H'
      write(*,*)'see psi_ij_n-1.out for eigenvectors of H'
      write(*,*)'******************************************************'

      open(14,file='psi_ij_n-1.out')
      write(14,*)'eigenvectors of H'
      write(14,*)'*****************************************************'
      do i=1,ord2
       write(14,*)'psi(',i,') = '
       do j=1,ord2
        if(abs(U2(j,i)).le.eps)then
         write(14,*)'component ',j,0.d0
        else
         write(14,*)'component ',j,U2(j,i)
        endif
       enddo
      enddo
      write(14,*)'*****************************************************'
      close(14)

      deallocate(occupation2)
      deallocate(occ2)
      deallocate(flag_u2)
      deallocate(H2)


! Static Green function G(i,j) = <Psi_0 | a+_i a_j | Psi_0>
      allocate(G(plane,plane))
      G(:,:)=0.d0
      do i=1,plane
       do j=1,plane
        if(spin(i).eq.1.and.spin(j).eq.1)then
         do i1=1,ord
          do i2=1,ord
           do k=1,nptot
            bra(k)=phi(i1,k)
            ket(k)=phi(i2,k)
           enddo
           fattore_bra=1.d0
           fattore_ket=1.d0
           call destroy(ket,nptot,j,plane,fattore_ket)
           call destroy(bra,nptot,i,plane,fattore_bra)
           overlap=fattore_bra*fattore_ket
           do k=1,nptot
            overlap=overlap*delta(bra(k),ket(k))
           enddo
           G(i,j)=G(i,j)+U(i1,1)*overlap*U(i2,1)
          enddo
         enddo
        endif
       enddo
      enddo
      open(36,file='green_ij.out')
      do i=1,plane
       do j=1,plane     
        if(spin(i).eq.1.and.spin(j).eq.1)then 
         if(abs(G(i,j)).le.eps)then
          write(36,*)0.d0,site(i),site(j)
         else
          write(36,*)G(i,j),site(i),site(j)
         endif
        endif
       enddo
      enddo
      close(36)
      deallocate(G)

! Static Green function Gp(pi,j) = <Psi_0 | a_i a+_j | Psi_0>
      allocate(Gp(plane,plane))
      Gp(:,:)=0.d0
      do i=1,plane
       do j=1,plane
        if(spin(i).eq.1.and.spin(j).eq.1)then
         do i1=1,ord
          do i2=1,ord
           do k=1,nptot
            bra(k)=phi(i1,k)
            ket(k)=phi(i2,k)
           enddo
           fattore_bra=1.d0
           fattore_ket=1.d0
           call create(ket,ket3,nptot,j,plane,fattore_ket)
           call create(bra,bra3,nptot,i,plane,fattore_bra)
           overlap=fattore_bra*fattore_ket
           do k=1,nptot3
            overlap=overlap*delta(bra3(k),ket3(k))
           enddo
           Gp(i,j)=Gp(i,j)+U(i1,1)*overlap*U(i2,1)
          enddo
         enddo
        endif
       enddo
      enddo
      open(36,file='greenp_ij.out')
      do i=1,plane
       do j=1,plane
        if(spin(i).eq.1.and.spin(j).eq.1)then
         if(abs(Gp(i,j)).le.eps)then
          write(36,*)0.d0,site(i),site(j)
         else
          write(36,*)Gp(i,j),site(i),site(j)
         endif
        endif
       enddo
      enddo
      close(36)
      deallocate(Gp)
      

! Static Density-Density correlation function 
! Rho(i,j) =  <Psi_0 | a+_i a_i  a+_j a_j | Psi_0>
      allocate(Rho(plane,plane))
      Rho(:,:)=0.d0
      do i=1,plane
       do j=1,plane
        do i1=1,ord
         Rho(i,j)=Rho(i,j)+U(i1,1)*occ(i1,i)*occ(i1,j)*U(i1,1)
        enddo
       enddo
      enddo
      open(37,file='densdens_ij.out')
      do i=1,plane
       do j=1,plane
        if(abs(Rho(i,j)).le.eps)then
         write(37,*)0.d0,site(i),spin(i),site(j),spin(j)
        else
         write(37,*)Rho(i,j),site(i),spin(i),site(j),spin(j)
        endif
       enddo
      enddo
      close(37)
      deallocate(Rho)


! Dynamical  Green function G(i,j) = <Psi_0 | a+_i exp(-t H) a_j | Psi_0>

      allocate(Gd(plane,plane,ntimes))
      Gd(:,:,:)=0.d0

      mu=epsilon(1) !epsilon(1)-epsilon2(1)
      
      do i=1,plane
       do j=1,plane
        if(spin(i).eq.1.and.spin(j).eq.1)then

        do it=1,ntimes
         time=dt*dble(it-1)

          do alpha=1,ord
           do beta=1,ord2

            do k=1,nptot
             bra(k)=phi(alpha,k)
            enddo
            
            do k=1,nptot
             ket2(k)=phi2(beta,k)
            enddo

            fattore_bra=1.d0
            call destroy(bra,nptot,i,plane,fattore_bra)
            ovrlp=fattore_bra
            do k=1,nptot2
             ovrlp=ovrlp*delta(bra(k),ket2(k))
            enddo

            if(ovrlp.ne.0.d0)then
!             write(*,*)'Stati .... ',i,j
!             write(*,*)'bra ov ',bra
!             write(*,*)'ket ov ',ket2
!             write(*,*)'abo',alpha,beta,ovrlp
!             write(*,*)
!             write(*,*)


             do gamma=1,ord2
              do eta=1,ord

               do k=1,nptot2
                bra2(k)=phi2(gamma,k)
               enddo
               do k=1,nptot
                ket(k)=phi(eta,k)
               enddo

               fattore_ket=1.d0
               call destroy(ket,nptot,j,plane,fattore_ket)
               overlap=fattore_ket
               do k=1,nptot2
                overlap=overlap*delta(bra2(k),ket(k))
               enddo

               if(overlap.ne.0.d0)then
!                write(*,*)'FASE 2 '
!                write(*,*)'bra ov ',bra2
!                write(*,*)'ket ov ',ket
!                write(*,*)'gdo',gamma,eta,overlap
!                write(*,*)
!                write(*,*)
    
                do n=1,ord2
 
                 Gd(i,j,it)=Gd(i,j,it)+                                   &
                    exp(-time*(epsilon2(n)-mu))*overlap*ovrlp*    &
                    U(alpha,1)*U(eta,1)*U2(beta,n)*U2(gamma,n)
                enddo
               endif
              enddo
             enddo
            endif
           enddo
          enddo
        enddo
       endif
       enddo
      enddo
      open(41,file='green_t_ij.out')
       do i=1,plane
        do j=1,plane
         if(spin(i).eq.1.and.spin(j).eq.1)then
         do it=1,ntimes
          time=dt*dble(it-1)  
          if(abs(Gd(i,j,it)).le.eps)then
           write(41,*)time,0.d0,site(i),site(j)
          else
           write(41,*)time,Gd(i,j,it),site(i),site(j)
          endif
         enddo
         write(41,*)
         write(41,*)
        endif
        enddo
       enddo
      close(41)
       

      allocate(Q(0:ndim,nq))
! wave vectors
      open(15,file='q.dat',status='old')
      do i=1,nq
       read(15,*)Q(1,i),Q(2,i)
       Q(1,i)=Q(1,i)*pi
       Q(2,i)=Q(2,i)*pi
      enddo
      close(15)
      do i=1,nq
       Q(0,i)=0.d0
       do j=1,ndim
        Q(0,i)=Q(0,i)+Q(j,i)**2
       enddo
       Q(0,i)=sqrt(Q(0,i))
      enddo

      open(42,file='green_t_q.out')
      do i=1,nq
       do it=1,ntimes
        time=dt*dble(it-1)

        temp=dcmplx(0.d0,0.d0)
        do i1=1,plane
         do i2=1,plane
          if(spin(i1).eq.1.and.spin(i2).eq.1)then
! build exp(-i Q (R-S))
           dot=0.d0
           do idim=1,ndim
            dot=dot+Q(idim,i)*(RS(idim,site(i1))-RS(idim,site(i2)))
           enddo
           eiqr=exp(-dcmplx(0.d0,1.d0)*dot)
           temp=temp+eiqr*Gd(i1,i2,it)
          endif
         enddo
        enddo
        write(42,*)time,real(temp)/dble(nsites),Q(0,i)

       enddo
       write(42,*)
       write(42,*)
      enddo
      close(42)
      deallocate(Gd)



! Dynamical  Green function G(i,j) = <Psi_0 | a_i exp(-t H) a+_j | Psi_0>

      allocate(Gpd(plane,plane,ntimes))
      Gpd(:,:,:)=0.d0

      allocate(bra3a(nptot3))
      allocate(ket3e(nptot3))
      bra3a(:)=0
      ket3e(:)=0

      mu=epsilon(1) !epsilon(1)-epsilon2(1)

      do i=1,plane
       do j=1,plane

        if(spin(i).eq.1.and.spin(j).eq.1)then

        do it=1,ntimes
         time=dt*dble(it-1)

          do alpha=1,ord
           do beta=1,ord3

            fattore_bra=1
            do k=1,nptot
             bra(k)=phi(alpha,k)
            enddo
            call create(bra,bra3a,nptot,i,plane,fattore_bra)
            do k=1,nptot3
             ket3(k)=phi3(beta,k)
            enddo
            ovrlp=fattore_bra
            do k=1,nptot3
             ovrlp=ovrlp*delta(bra3a(k),ket3(k))
            enddo

            if(ovrlp.ne.0.d0)then
             if(i.eq.j)then
             write(*,*)'Stati .... ',i,j
             write(*,*)'bra ov ',bra3a
             write(*,*)'ket ov ',ket3
             write(*,*)'abo',alpha,beta,ovrlp
             write(*,*)
             write(*,*)
             endif

             do gamma=1,ord3
              do eta=1,ord
                
               fattore_ket=1
               do k=1,nptot3
                bra3(k)=phi3(gamma,k)
               enddo
               do k=1,nptot
                ket(k)=phi(eta,k)
               enddo
               call create(ket,ket3e,nptot,j,plane,fattore_ket)
               overlap=fattore_ket
               do k=1,nptot3
                overlap=overlap*delta(bra3(k),ket3e(k))
               enddo

               if(overlap.ne.0.d0)then
                if(i.eq.j)then
                write(*,*)'FASE 2 '
                write(*,*)'bra ov ',bra3
                write(*,*)'ket ov ',ket3e
                write(*,*)'gdo',gamma,eta,overlap
                write(*,*)
                write(*,*)
                endif

                do n=1,ord3

                 Gpd(i,j,it)=Gpd(i,j,it)+                         &
                    exp(-time*(epsilon3(n)-mu))*overlap*ovrlp*    &
                    U(alpha,1)*U(eta,1)*U3(beta,n)*U3(gamma,n)
                enddo
               endif
              enddo
             enddo
            endif
           enddo
          enddo
        enddo
       endif
       enddo
      enddo
      open(41,file='greenp_t_ij.out')
       do i=1,plane
        do j=1,plane
         if(spin(i).eq.1.and.spin(j).eq.1)then
         do it=1,ntimes
          time=dt*dble(it-1)
          if(abs(Gpd(i,j,it)).le.eps)then
           write(41,*)time,0.d0,site(i),site(j)
          else
           write(41,*)time,Gpd(i,j,it),site(i),site(j)
          endif
         enddo
         write(41,*)
         write(41,*)
        endif
        enddo
       enddo
      close(41)
      open(42,file='greenp_t_q.out')
      do i=1,nq
       do it=1,ntimes
        time=dt*dble(it-1)

        temp=dcmplx(0.d0,0.d0)
        do i1=1,plane
         do i2=1,plane
          if(spin(i1).eq.1.and.spin(i2).eq.1)then
! build exp(-i Q (R-S))
           dot=0.d0
           do idim=1,ndim
            dot=dot+Q(idim,i)*(RS(idim,site(i1))-RS(idim,site(i2)))
           enddo
           eiqr=exp(-dcmplx(0.d0,1.d0)*dot)
           temp=temp+eiqr*Gpd(i1,i2,it)
          endif
         enddo
        enddo
        write(42,*)time,real(temp)/dble(nsites),Q(0,i)

       enddo
       write(42,*)
       write(42,*)
      enddo
      close(42)
      deallocate(Gpd)







! Dynamical Density-Density correlation function 
! Rho(i,j; t) =  <Psi_0 | a+_i a_i exp(-t (H - E_0)) a+_j a_j | Psi_0>
      allocate(Rhod(plane,plane,ntimes))
      Rhod(:,:,:)=0.d0
      do it=1,ntimes
       time=dt*dble(it-1)

       do i=1,plane
        do j=1,plane

         do i1=1,ord
          do i2=1,ord
           do t=1,ord
            Rhod(i,j,it)=Rhod(i,j,it)+                                 &
             U(i1,1)*U(i1,t)*U(i2,t)*U(i2,1)*occ(i1,i)*occ(i2,j)       &
             *exp(-time*(epsilon(t)-epsilon(1)))    
           enddo
          enddo
         enddo

        enddo
       enddo
      enddo

      open(38,file='dens_t_dens_ij.out')
      do it=1,ntimes
       time=dt*dble(it-1)
       do i=1,plane
        do j=1,plane
         if(abs(Rhod(i,j,it)).le.eps)then
          write(38,*)time,0.d0,i,j
         else
          write(38,*)time,Rhod(i,j,it),i,j
         endif
        enddo
       enddo
       write(38,*)
       write(38,*)
      enddo
      close(38)

      open(39,file='dens_t_dens_q.out')    
      do i=1,nq  
       do it=1,ntimes
        time=dt*dble(it-1)

        temp=dcmplx(0.d0,0.d0)
        do i1=1,plane
         do i2=1,plane
! build exp(-i Q (R-S))
          dot=0.d0
          do idim=1,ndim
           dot=dot+Q(idim,i)*(RS(idim,site(i1))-RS(idim,site(i2)))
          enddo
          eiqr=exp(-dcmplx(0.d0,1.d0)*dot)
          temp=temp+eiqr*Rhod(i1,i2,it)
         enddo
        enddo
        write(39,*)time,real(temp)/dble(nptot),Q(0,i)

       enddo
       write(39,*)
       write(39,*) 
      enddo
      close(39)
 

      end program codicino


      subroutine nn(system,i,j,si,sj,f)

      implicit none
      integer, intent(in) :: system,i,j,si,sj
      integer, intent(out):: f
  
      f=0

      if(system.eq.1)then

       if(si.eq.sj)then
        if(i.eq.1.and.j.eq.2)f=1
        if(i.eq.1.and.j.eq.3)f=1
        if(i.eq.2.and.j.eq.4)f=1
        if(i.eq.3.and.j.eq.4)f=1
       endif

      elseif(system.eq.2)then

       if(si.eq.sj)then
        if(i.eq.1.and.j.eq.2)f=1
        if(i.eq.2.and.j.eq.3)f=1
        if(i.eq.3.and.j.eq.4)f=1
        if(i.eq.4.and.j.eq.4)f=1
       endif      

      elseif(system.eq.3)then

       if(si.eq.sj)then
        if(i.eq.1.and.j.eq.2)f=1
       endif

      else
  
       write(*,*)'Problems with system '
       stop      

      endif

      return
      end
