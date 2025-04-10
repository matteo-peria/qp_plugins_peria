program aos_valence_test

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  call main()

end program

! ---

subroutine main()
  implicit none

  integer :: row,col
  integer :: i,j
  double precision :: accu = 0.d0
  double precision :: trace = 0.d0

!!  print*, '' 
!!  print*, 'VALENCE COEFFICIENTS (SPHERICAL)'
!!  print*, '' 
!!
!!  provide ao_val_coef_sphe
!!  write(*,*) "AO_VAL_COEF, shape = ", shape(ao_val_coef_sphe)
!!  do row = 1, ao_sphe_num
!!    write(*,'(100F12.7)'), (ao_val_coef_sphe(row,col), col = 1, ao_sphe_num)
!!  end do
!!  print*, ''
!!
!!  provide ao_val_sphe_coef_normed
!!  write(*,*) "AO_VAL_SPHE_COEF_NORMED, shape = ", shape(ao_val_sphe_coef_normed)
!!  do row = 1, ao_sphe_num
!!    write(*,'(100F12.7)'), (ao_val_sphe_coef_normed(row,col), col = 1,ao_sphe_num)
!!  end do
!!  print*, ''

  ! START TESTING OVERLAPS 
  print*, '' 
  print*, 'OVERLAP TESTS'
  print*, '' 
  
  print*, '' 
  print*, 'OVERLAP TESTS VALENCE SPACE (SPHERICAL)'
  print*, '' 

  provide ao_val_overlap_sphe
  print*,'ao_val_overlap_sphe'
  do row = 1, ao_sphe_num
    write(*,'(100F12.7)'), (ao_val_overlap_sphe(row,col), col = 1, ao_sphe_num)
    trace += ao_val_overlap_sphe(row,row)
  end do
  print*, "Trace: ", trace 
  trace = 0.d0
  print*, ''

!  provide ao_val_overlap_sphe_evec_overlap
!  print*,'ao_val_overlap_sphe_evec_overlap'
!  do row = 1, ao_sphe_num
!   write(*,'(100F12.7)'), (ao_val_overlap_sphe_evec_overlap(row,col), col = 1, ao_sphe_num)
!   trace += ao_val_overlap_sphe_evec_overlap(row,row)
!  enddo 
!  print*, "Trace: ", trace 
!  trace = 0.d0
!  print*, ''

  provide ao_val_overlap_sphe_as_matprod
  print*,'ao_val_overlap_sphe_as_matprod'
  do row = 1, ao_sphe_num
   write(*,'(100F12.7)'), (ao_val_overlap_sphe_as_matprod(row,col), col = 1, ao_sphe_num)
   trace += ao_val_overlap_sphe_as_matprod(row,row)
  enddo 
  print*, "Trace: ", trace 
  trace = 0.d0
  print*, ''

  provide ao_val_sphe_normed_overlap
  print*,'ao_val_sphe_normed_overlap'
  do row = 1, ao_sphe_num
   write(*,'(100F12.7)'), (ao_val_sphe_normed_overlap(row,col), col = 1, ao_sphe_num)
   trace += ao_val_sphe_normed_overlap(row,row)
  enddo 
  print*, "Trace: ", trace 
  trace = 0.d0
  print*, ''
 
  provide ao_val_overlap_sphe_normed_as_matprod
  print*,'ao_val_overlap_sphe_normed_as_matprod'
  do row = 1, ao_sphe_num
   write(*,'(100F12.7)'), (ao_val_overlap_sphe_normed_as_matprod(row,col), col = 1, ao_sphe_num)
   trace += ao_val_overlap_sphe_normed_as_matprod(row,row)
  enddo 
  print*, "Trace: ", trace 
  trace = 0.d0
  print*, ''
 
  !---

end subroutine
