program core_pp
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
 integer :: i
  print *, 'Hello world'
 do i = 1, ao_num
  write(*,'(100(F16.10,X))')ao_val_normed_overlap(i,:)
 enddo
end
