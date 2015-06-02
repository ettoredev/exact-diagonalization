      subroutine destroy(stato,nptot,index,plane,fattore)

      integer  :: nptot,index,plane
      integer  :: i,j,temp
      real*8   :: fattore

      integer, dimension(nptot) :: stato

      j=0
      do i=1,nptot
       if(stato(i).eq.index)j=i
      enddo
      if(j.eq.0)then
       fattore=0.d0
       go to 53
      endif
      do i=j,nptot-1
       stato(i)=stato(i+1)
      enddo
      stato(nptot)=plane+100
      fattore=fattore*(-1.d0)**dble(j+1)
53    continue

      return

      end subroutine destroy

      subroutine create(stato,statof,nptot,index,plane,fattore)

      integer  :: nptot,index,plane
      integer  :: i,j,temp
      real*8   :: fattore

      integer, dimension(nptot)   :: stato
      integer, dimension(nptot+1) :: statof

!      write(*,*)'stato iniziale -> ',(stato(i),i=1,nptot)
!      if(stato(nptot).lt.plane+1)then
!       write(*,*)'stato finale -> stato a ',nptot+1,'particelle'
!       stop
!      endif

      j=0
      k=1
      do i=1,nptot-1
       if(stato(i).eq.index)j=i
       if(stato(i).lt.index.and.index.lt.stato(i+1))k=i+1
      enddo
      if(stato(nptot).eq.index)then
       j=nptot
      endif
      if(stato(nptot).lt.index)then
       k=nptot+1
      endif
      if(j.ne.0)then
       fattore=0.d0
       statof(:)=0
      else
       fattore=1.d0
       do i=1,k-1
        statof(i)=stato(i)
       enddo
       statof(k)=index
       do i=k+1,nptot+1
        statof(i)=stato(i-1)
       enddo
       fattore=fattore*(-1.d0)**dble(k+1)
      endif
!      write(*,*)'fattore = ',fattore
!      write(*,*)'stato finale -> ',(stato(i),i=1,nptot)

      return

      end subroutine create
