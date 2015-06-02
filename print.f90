      subroutine print(v,sz)

!      use common

      integer i,j
      integer, dimension(:), intent(in out) :: v

      write(*,*)'called'
      do i=1,sz
       v(i)=i
      enddo

      end subroutine print
