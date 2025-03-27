program test_core_valence_pseudo

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  call main()

end

! ---

subroutine main()
  implicit none
  integer :: i_ao, j_mo, i, j
  double precision :: trace

   print*, "ao_num: ", ao_num
   print*, "n_valence_pseudo_orb: ", n_valence_pseudo_orb
!
!  write(*,*) "TEST ON: mo_val_coef"
!  write(*,*) "SHAPE: ", shape(mo_val_coef)
!  write(*,*) "Coefficients"
!  print*, "MO		AO		Coeff"
!  do j_mo = 1, n_valence_pseudo_orb
!    do i_ao = 1, ao_num
!      print*, j_mo, i_ao, mo_val_coef(i_ao, j_mo)
!    enddo
!  enddo      

!  write(*,*)
!  write(*,*) "TEST ON proj_mo_val"
!  write(*,*) "SHAPE: ", shape(proj_mo_val)
!  do i = 1, ao_num
!    write(*,*) (proj_mo_val(i,j)," ",j=1,n_valence_pseudo_orb)
!  enddo
!
!  write(*,*)
!  write(*,*) "TEST ON proj_mo_val_sphe"
!  write(*,*) "SHAPE: ", shape(proj_mo_val_sphe)
!  do i = 1, ao_sphe_num
!    write(*,*) (proj_mo_val(i,j)," ",j=1,n_valence_pseudo_orb)
!  enddo

  !!write(*,*) "AO_CORE_PSEUDO_COEF" 
  !!do i = 1, ao_sphe_num !ao_num
  !!  write(*,'(100(F16.10,X))') ao_val_coef(i,:)
  !!enddo

   write(*,*)
   write(*,*) "TEST ON ao_val_overlap" 
   write(*,*) "SHAPE: ", shape(ao_val_overlap)
   trace = 0.d0
   do i = 1, ao_num !ao_sphe_num !ao_num
     !write(*,'(100(F16.10,X))') ao_val_overlap(i,:)
     trace += ao_val_overlap(i,i)
   enddo
   write(*,*) "TRACE: ", trace
 
   write(*,*)
   write(*,*) "TEST ON ao_val_overlap_sphe" 
   write(*,*) "SHAPE: ", shape(ao_val_overlap_sphe)
   trace = 0.d0
   do i = 1, ao_sphe_num 
     !write(*,'(100(F16.10,X))') ao_val_overlap_sphe(i,:)
     trace += ao_val_overlap_sphe(i,i)
   enddo
   write(*,*) "TRACE: ", trace
  
  provide mo_overlap_sphe

  provide ao_cart_to_sphe_overlap_inv


  !write(*,*) "ao_pp_overlap_eigval" 
  !do i = 1, ao_num
  ! print*,i,ao_pp_overlap_eigval(i)
  !enddo
  !write(*,*) "ao_pp_overlap_eigvec" 
  !do i = 1, ao_num
  ! write(*,'(100(F16.10,X))')ao_pp_overlap_eigvec(i,:)
  !enddo
  !provide overlap_cart

  !trace = 0.d0
  !do i_ao = 1, ao_num
  !  trace += ao_val_overlap(i_ao,i_ao) 
  !enddo
  !print*, "TRACE OF AO_PP_OVERLAP", trace 
end
