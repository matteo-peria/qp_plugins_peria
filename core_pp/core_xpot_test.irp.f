program core_xpot_test

  BEGIN_DOC
  ! Testing the core exchange potential computed with different numerical grids:
  !  CORE_XPOT_NUMERIC
  !  CORE_XPOT_NUMERIC_FULL_EXTRA_GRID
  !  CORE_XPOT_NUMERIC_ADAPT_OLD
  !  CORE_XPOT_NUMERIC_FULL_ADAPT_GRID
  END_DOC

  implicit none

  character(len=100) :: filename
  logical :: has_arg

  write(*,*) 'filename', filename

  has_arg = command_argument_count() >= 1

  if (has_arg) then
    call get_command_argument(1, filename)
    call test_grids_core_xpot(trim(filename))
  else
    call test_grids_core_xpot()
  end if

end program


subroutine test_grids_core_xpot(filename)
  use iso_fortran_env, only: output_unit
  use io_test_interface
  implicit none

  character(len=*), intent(in), optional :: filename
  character(len=100) :: output

  double precision :: diff_prun1_full2
  double precision :: diff_prun1_prun2
  double precision :: diff_prun1_full3
  double precision :: diff_prun1_full3_old
  double precision :: diff_prun1_prun3
  double precision :: diff_prun1_fuzzy

  integer :: ios
  integer :: unit_out

  !print*,'EXACT'
  !call print_db_array(core_xpot_exact, list_valence_pseudo, list_valence_pseudo)
  
  write(*,'(A)') repeat('=', 70)

  print*,'CORE_XPOT_NUMERIC'
  call compute_dp_array_diff(core_xpot_numeric,            &
                           & core_xpot_exact,              &
                           & show = .False.,               &
                           & row_indx=list_valence_pseudo, &
                           & col_indx=list_valence_pseudo, &
                           & diff=diff_prun1_prun2)

  write(*,'(A)') repeat('=', 70)

  print*,'CORE_XPOT_NUMERIC_FULL_EXTRA_GRID'
  call compute_dp_array_diff(core_xpot_numeric_full_extra_grid, &
                           & core_xpot_exact, &
                           & show = .False.,              &
                           & row_indx=list_valence_pseudo,&
                           & col_indx=list_valence_pseudo, &
                           & diff=diff_prun1_full2)

  write(*,'(A)') repeat('=', 70)

  diff_prun1_prun2=-1.0
  diff_prun1_full2=-1.0

  print*,'CORE_XPOT_NUMERIC_ADAPT_OLD'

  call compute_dp_array_diff(core_xpot_numeric_adapt_old, &
                           & core_xpot_exact, &
                           & show = .False.,              &
                           & row_indx=list_valence_pseudo,       &
                           & col_indx=list_valence_pseudo, &
                           & diff=diff_prun1_full3_old)

  write(*,'(A)') repeat('=', 70)

  print*,'CORE_XPOT_NUMERIC_FULL_ADAPT_GRID'
  call compute_dp_array_diff(core_xpot_numeric_full_adapt_grid, &
                         & core_xpot_exact, &
                         & show = .False.,              &
                         & row_indx=list_valence_pseudo,       &
                         & col_indx=list_valence_pseudo, &
                         & diff=diff_prun1_full3)

  write(*,'(A)') repeat('=', 70)

  print*,'CORE_XPOT_FUZZY_GRID'
  call compute_dp_array_diff(core_xpot_fuzzy_grid, &
                         & core_xpot_exact, &
                         & show = .False.,              &
                         & row_indx=list_valence_pseudo,       &
                         & col_indx=list_valence_pseudo, &
                         & diff=diff_prun1_fuzzy)

  write(*,'(A)') repeat('=', 70)


end subroutine
