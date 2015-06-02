      function delta(i,j)

      implicit none

      integer i,j
      real*8 delta
      delta=0.d0
      if(i.eq.j)delta=1.d0

      return

      end function
