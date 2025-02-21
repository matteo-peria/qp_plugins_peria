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
!  write(*,*) "TEST ON: mo_valence_pseudo_coef"
!  write(*,*) "SHAPE: ", shape(mo_valence_pseudo_coef)
!  write(*,*) "Coefficients"
!  print*, "MO		AO		Coeff"
!  do j_mo = 1, n_valence_pseudo_orb
!    do i_ao = 1, ao_num
!      print*, j_mo, i_ao, mo_valence_pseudo_coef(i_ao, j_mo)
!    enddo
!  enddo      

!  write(*,*)
!  write(*,*) "TEST ON P_mo_pp_valence"
!  write(*,*) "SHAPE: ", shape(P_mo_pp_valence)
!  do i = 1, ao_num
!    write(*,*) (P_mo_pp_valence(i,j)," ",j=1,n_valence_pseudo_orb)
!  enddo
!
!  write(*,*)
!  write(*,*) "TEST ON P_mo_pp_valence_sphe"
!  write(*,*) "SHAPE: ", shape(P_mo_pp_valence_sphe)
!  do i = 1, ao_cart_to_sphe_num
!    write(*,*) (P_mo_pp_valence(i,j)," ",j=1,n_valence_pseudo_orb)
!  enddo

  !!write(*,*) "AO_CORE_PSEUDO_COEF" 
  !!do i = 1, ao_cart_to_sphe_num !ao_num
  !!  write(*,'(100(F16.10,X))') ao_pp_valence_coef(i,:)
  !!enddo

   write(*,*)
   write(*,*) "TEST ON ao_pp_overlap" 
   write(*,*) "SHAPE: ", shape(ao_pp_overlap)
   trace = 0.d0
   do i = 1, ao_num !ao_cart_to_sphe_num !ao_num
     !write(*,'(100(F16.10,X))') ao_pp_overlap(i,:)
     trace += ao_pp_overlap(i,i)
   enddo
   write(*,*) "TRACE: ", trace
 
   write(*,*)
   write(*,*) "TEST ON ao_pp_overlap_sphe" 
   write(*,*) "SHAPE: ", shape(ao_pp_overlap_sphe)
   trace = 0.d0
   do i = 1, ao_cart_to_sphe_num 
     !write(*,'(100(F16.10,X))') ao_pp_overlap_sphe(i,:)
     trace += ao_pp_overlap_sphe(i,i)
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
  !  trace += ao_pp_overlap(i_ao,i_ao) 
  !enddo
  !print*, "TRACE OF AO_PP_OVERLAP", trace 
end
