program test_ao_cart_to_sphe_conversion

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none
  call main()

end

! ---

subroutine main()
  implicit none
  double precision :: accu
  integer :: i_ao, j_ao, i, j

  write(*,*) "ao_cart_to_sphe_num: ", ao_cart_to_sphe_num
  write(*,*) "ao_cart_to_sphe_num_old: ", ao_cart_to_sphe_num_old

! write(*,*) "Difference between ao_cart_to_sphe_coef and ao_cart_to_sphe_coef_old"
! write(*,*) sum(ao_cart_to_sphe_coef(1:ao_num,1:ao_cart_to_sphe_num))!-ao_cart_to_sphe_coef_old)
! write(*,*) sum(ao_cart_to_sphe_coef)!-ao_cart_to_sphe_coef_old)
! write(*,*) sum(ao_cart_to_sphe_coef_old)
! write(*,*) sum(ao_cart_to_sphe_coef(1:ao_num,1:ao_cart_to_sphe_num)-ao_cart_to_sphe_coef_old)
 accu = 0.d0
 do j = 1, ao_cart_to_sphe_num
  do i = 1, ao_num
   accu += dabs(ao_cart_to_sphe_coef(i,j) - ao_cart_to_sphe_coef_old(i,j))
  enddo
 enddo
 print*,'accu = ',accu

 provide  overlap_cart

end
