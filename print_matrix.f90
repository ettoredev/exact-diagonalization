      subroutine print_matrix(filename,matrix,order)

      implicit none
      integer                            :: order,i,j
      character*20                       :: filename
      real*8,    dimension(order,order)  :: matrix

      write(*,*)'Sono in print_matrix ',filename
      i=index(filename,'t')
      write(*,*)'Sono in print_matrix ',i,filename(1:i)
      write(*,*)
      open(14,file=filename(1:i))
      do i=1,order
       do j=1,order
        write(14,*)matrix(i,j),i,j
       enddo
      enddo
     close(14)

end subroutine print_matrix
