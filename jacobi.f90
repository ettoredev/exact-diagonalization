      subroutine jacobi(H,order,epsilon,U,eps,quality)

      implicit none

      integer    :: i,j,k,l,m,imax,jmax,order
      real *8    :: eps,temp1,temp2,Smax,quality,theta,pi
      parameter  (pi=dacos(-1.d0))
      real*8,    dimension(order)        :: epsilon
      real*8,    dimension(order,order)  :: H,U,H2
      real*8,    dimension(order)        :: col1,col2

!     trascrizione di H su H2, inizializzazione di U
      do i=1,order
       do j=1,order
        H2(i,j)=H(i,j)
       enddo
      enddo

      do i=1,order
       do j=1,order
        U(i,j)=0.d0
       enddo
       U(i,i)=1.d0
      enddo

      open(14,file='jacobi.out')
      do m=1,5000000
!      ricerca dell'elemento fuori diagonale di massimo modulo
       Smax=abs(H2(1,2))
       imax=1
       jmax=2
       do i=1,order-1
        do j=i+1,order
         if(abs(H2(i,j)).ge.smax)then
          Smax=abs(H2(i,j))
          imax=i
          jmax=j
         endif
        enddo
       enddo
       write(14,*)'diagonalization in progress Smax = ',Smax
       if(Smax.lt.eps)then
        write(14,*)'end of diagonalization' 
        go to 4
       endif
!      costruzione dell'angolo di Givens
       theta=0.d0
       if(H2(imax,imax).eq.H2(jmax,jmax))then
        theta=acos(-1.d0)/4.d0
       else
        theta=0.5d0*atan(2.d0*H2(imax,jmax)/(H2(imax,imax)-H2(jmax,jmax)))
       endif
!      salvataggio delle colonne imax e jmax di U su col1,col2
       do i=1,order
        col1(i)=U(i,imax)
        col2(i)=U(i,jmax)
       enddo
!      update di U
       do i=1,order
        U(i,imax)=cos(theta)*col1(i)+sin(theta)*col2(i)
        U(i,jmax)=cos(theta)*col2(i)-sin(theta)*col1(i)
       enddo
!      salvataggio delle colonne imax e jmax di H2 su col1,col2
       do i=1,order
        col1(i)=H2(i,imax)
        col2(i)=H2(i,jmax)
       enddo
!      update di H2
       do i=1,order
        H2(i,imax)=cos(theta)*col1(i)+sin(theta)*col2(i)
        H2(i,jmax)=cos(theta)*col2(i)-sin(theta)*col1(i)
       enddo
       do i=1,order
        H2(imax,i)=H2(i,imax)
        H2(jmax,i)=H2(i,jmax)
       enddo
       H2(imax,imax)=cos(theta)**2*col1(imax)+sin(theta)**2*col2(jmax)  &
                   +2.d0*sin(theta)*cos(theta)*col1(jmax)
       H2(jmax,jmax)=cos(theta)**2*col2(jmax)+sin(theta)**2*col1(imax)  &
                   -2.d0*sin(theta)*cos(theta)*col1(jmax)
       H2(imax,jmax)=0.d0
       H2(jmax,imax)=0.d0
      enddo

4     continue

      do i=1,order
       epsilon(i)=H2(i,i)
      enddo

      call print_matrix("H2_ij.out",H2,order)

      quality=0.d0
      do i=1,order
       do j=1,order
        temp1=0.d0
        do k=1,order
         temp1=temp1+H(j,k)*U(k,i)
        enddo
        temp1=temp1-epsilon(i)*U(j,i)
        quality=quality+temp1**2
       enddo
      enddo
      quality=dble(order)/quality
      write(14,*)'quality ',quality

      return

      end subroutine jacobi
