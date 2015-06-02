      function factorial(n) result(res)
 
      implicit none
      integer, intent (in) :: n
      integer*8            :: res
      integer              :: i
 
      if(n.lt.0)then
       write(*,*)'invalid operand in factorial'
       stop
      endif

      res=1
      if(n.eq.0)then
       go to 53
      endif

      do i=1,n
       res=res*i
      enddo

53    continue
 
      end function factorial

      function binomial(n,k) result(bin)

      implicit none
      integer, intent (in) :: n,k
      integer*8            :: factorial,bin

      if(k.eq.2) then
       bin=(n*(n-1))/2
      else
       bin=factorial(n)/(factorial(k)*factorial(n-k))
      endif

      end function binomial

