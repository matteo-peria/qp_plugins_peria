program test_ao_valence_in_r

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call write_int(6, my_n_pt_r_grid, 'radial  external grid over')
  call write_int(6, my_n_pt_a_grid, 'angular external grid over')

  my_extra_grid_becke  = .True.
  PROVIDE tc_grid2_a tc_grid2_r
  my_n_pt_r_extra_grid = tc_grid2_r
  my_n_pt_a_extra_grid = tc_grid2_a
  touch my_extra_grid_becke my_n_pt_r_extra_grid my_n_pt_a_extra_grid

  call write_int(6, my_n_pt_r_extra_grid, 'radial  internal grid over')
  call write_int(6, my_n_pt_a_extra_grid, 'angular internal grid over')

  
  read_ao_integrals_n_e = .false.
  touch read_ao_integrals_n_e

  call main()

end program

! ---

subroutine main()
  implicit none

  integer :: row,col
  integer :: i,j
  double precision :: accu = 0.d0
  double precision :: trace = 0.d0

  read_ao_integrals_n_e = .false.
  write_ao_integrals_n_e = .false.


  ! START TESTING VALENCE INTEGRALS
  provide ao_integrals_ne_to_ao_val_numeric
  print*, "AO_INTEGRALS_NE_TO_AO_VAL_NUMERIC"
  !print*, "Shape of ao_integrals_ne_to_ao_val_numeric: ", shape(ao_integrals_ne_to_ao_val_numeric)
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_integrals_ne_to_ao_val_numeric(row,col), col = 1, ao_num)
    !write(*,*), (ao_integrals_ne_to_ao_val_numeric(row,col), col = 1, ao_num)
  end do
  
  provide ao_integrals_ne_to_ao_val_matprod
  print*, "AO_INTEGRALS_NE_TO_AO_VAL_MATPROD"
  !print*, "Shape of ao_integrals_ne_to_ao_val_matprod: ", shape(ao_integrals_ne_to_ao_val_matprod)
  do row = 1, ao_num
    !write(*,'(100F12.7)'), (ao_integrals_ne_to_ao_val_numeric(row,col), col = 1, ao_num)
    write(*,'(100F12.7)'), (ao_integrals_ne_to_ao_val_matprod(row,col), col = 1, ao_num)
  end do
  print*, ''

!  provide ao_val_coef_normed
!  write(*,*) "AO_VAL_COEF_NORMED, shape = ", shape(ao_val_coef_normed)
!  do row = 1, ao_num
!    write(*,'(100F12.7)'), (ao_val_coef_normed(row,col), col = 1,ao_num)
!  end do
!  print*, ''

  provide ao_val_integrals_ne
  print*, "AO_VAL_INTEGRALS_NE"
  !print*, "Shape of ao_val_integrals_ne: ", shape(ao_val_integrals_ne)
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_val_integrals_ne(row,col), col = 1, ao_num)
  end do


  provide ao_val_integrals_ne_numeric
  print*, "AO_VAL_INTEGRALS_NE_NUMERIC"
  !print*, "Shape of ao_val_integrals_ne_numeric: ", shape(ao_val_integrals_ne_numeric)
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_val_integrals_ne_numeric(row,col), col = 1, ao_num)
  end do


  ! COMPARE ELEMENT BY ELEMENT

  do row = 1, ao_num
    do col = 1, ao_num
      write(*,'(2I, 4F12.7)'), row, col, &
        & ao_integrals_ne_to_ao_val_numeric(row,col), &
        & ao_integrals_ne_to_ao_val_matprod(row,col), &
        & ao_val_integrals_ne(row,col), &
        & ao_val_integrals_ne_numeric(row,col)
    end do
  end do
  ! END TESTING VALENCE INTEGRALS

end subroutine
