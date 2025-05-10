program core_xpot_test

  BEGIN_DOC
  ! Testing the core exchange potential computed with different numerical grids
  END_DOC

  implicit none

  my_grid_becke  = .True.
  my_extra_grid_becke  = .True.
  my_grid_adapt        = .True.

  my_n_pt_r_grid       = my_n_pt_r_grid
  my_n_pt_a_grid       = my_n_pt_a_grid

  my_n_pt_r_extra_grid = my_n_pt_r_extra_grid
  my_n_pt_a_extra_grid = my_n_pt_a_extra_grid
  
  my_n_pt_r_grid_adapt = my_n_pt_r_grid_adapt
  my_n_pt_a_grid_adapt = my_n_pt_a_grid_adapt

  touch my_extra_grid_becke 
  touch my_grid_adapt my_n_pt_r_grid_adapt my_n_pt_a_grid_adapt 

  call main()

end program

! ---

subroutine main()
  !use io_array !io_module
  use io_test_interface
  implicit none

  double precision :: difference


  !print*,'EXACT'
  !call print_db_array(core_xpot_exact, list_valence_pseudo, list_valence_pseudo)
  
  write(*,'(A)') repeat('=', 70)

  write(*,*) 'n_points_rad_extra_grid', n_points_rad_extra_grid

  print*,'CORE_XPOT_NUMERIC'
  call compute_dp_array_diff(core_xpot_numeric, &
                           & core_xpot_exact, &
                           & show = .False.,              &
                           & row_indx=list_valence_pseudo,&
                           & col_indx=list_valence_pseudo)

  write(*,'(A)') repeat('=', 70)
!
!  print*,'CORE_XPOT_NUMERIC_FULL_EXTRA_GRID'
!  call compute_dp_array_diff(core_xpot_numeric_full_extra_grid, &
!                           & core_xpot_exact, &
!                           & show = .False.,              &
!                           & row_indx=list_valence_pseudo,&
!                           & col_indx=list_valence_pseudo)
!
!  write(*,'(A)') repeat('=', 70)
!
!
!  print*,'CORE_XPOT_NUMERIC_ADAPT_OLD'
!  call compute_dp_array_diff(core_xpot_numeric_adapt_old, &
!                           & core_xpot_exact, &
!                           & show = .False.,              &
!                           & row_indx=list_valence_pseudo,       &
!                           & col_indx=list_valence_pseudo)
!
!  write(*,'(A)') repeat('=', 70)
!
!  print*,'CORE_XPOT_NUMERIC_FULL_ADAPT_GRID'
!  call compute_dp_array_diff(core_xpot_numeric_full_adapt_grid, &
!                           & core_xpot_exact, &
!                           & show = .False.,              &
!                           & row_indx=list_valence_pseudo,       &
!                           & col_indx=list_valence_pseudo)
!
!  write(*,'(A)') repeat('=', 70)

end subroutine
