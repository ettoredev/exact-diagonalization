      subroutine print_array(filename,array,order)

      integer                      :: order,i
      character*48                 :: filename
      real*8,    dimension(order)  :: array
      real*8                       :: eps
      parameter ( eps=1.d-10 )

      i=index(filename,'t')
      open(14,file=filename(1:i))
      do i=1,order
       if(abs(array(i)).le.eps)then
        write(14,*)0.d0,i
       else
        write(14,*)array(i),i
       endif
      enddo
     close(14)

end subroutine print_array
