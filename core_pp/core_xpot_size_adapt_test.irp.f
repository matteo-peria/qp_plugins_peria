program core_xpot_size_adapt_test

  BEGIN_DOC
  ! Testing the core exchange potential computed with different numerical grids
  END_DOC

  implicit none

  character(len=100) :: filename
  logical :: has_arg


  has_arg = command_argument_count() >= 1

  if (has_arg) then
    call get_command_argument(1, filename)
    call test_grids_core_xpot_adapt(trim(filename))
  else
    call test_grids_core_xpot_adapt()
  end if

end program


subroutine test_grids_core_xpot_adapt(filename)
  use iso_fortran_env, only: output_unit
  use io_test_interface
  implicit none

  character(len=*), intent(in), optional :: filename
  character(len=100) :: output

  double precision :: diff_prun1_full3

  integer :: ios
  integer :: unit_out

  !print*,'EXACT'
  !call print_db_array(core_xpot_exact, list_valence_pseudo, list_valence_pseudo)
  
  write(*,'(A)') repeat('=', 70)

  print*,'CORE_XPOT_NUMERIC_FULL_ADAPT_GRID'
  call compute_dp_array_diff(core_xpot_numeric_full_adapt_grid, &
                           & core_xpot_exact, &
                           & show = .False.,              &
                           & row_indx=list_valence_pseudo,       &
                           & col_indx=list_valence_pseudo, &
                           & diff=diff_prun1_full3)

  write(*,'(A)') repeat('=', 70)

  logical :: file_exists
  if (present(filename)) then
    output = trim(filename)//'/core_xpot_size_adapt.out'
    inquire(file=output, exist=file_exists)
    
    if (file_exists) then
      open(newunit=unit_out, & 
         & file=output,    & 
         & status='old',     & 
         & position='append',&  
         & action='write',   & 
         & iostat=ios        & 
         &)
    else
      open(newunit=unit_out, & 
         & file=output,    & 
         & status='new',     & 
         & action='write',   & 
         & iostat=ios        & 
         &)
    end if

    if (ios /= 0) then
      print *, "Error opening file: ", trim(filename)
      stop 1
    end if

  else
    unit_out=output_unit
  end if

  character(len=100) :: fmt_str
  fmt_str = '(13A8, A12, 2A8, A12)'
  if (.not.file_exists) then
    write(unit_out, fmt_str) &
                       & 'N_rad1', 'N_ang1', 'N_sph1', 'N_ful1', 'N_pru1', &
                       & 'N_rad2', 'N_ang2', 'N_sph2', 'N_ful2',           &
                       & 'N_rad3', 'N_ang3', 'N_sph3', 'N_ful3',           &
                       & 'TOT(p1+f3)', &
                       & 'radius', 'alpha', &
                       & 'RES(p1+f3)'
  end if

  character(len=100) :: fmt_num
  fmt_num = '(13I8, I12, 2F6.3, E12.3)'

  write(unit_out, fmt_num)                            &
      & n_points_rad_grid,                            &
      & n_points_ang_grid,                            &
      & n_points_rad_grid*n_points_ang_grid,          &
      & nucl_num*n_points_rad_grid*n_points_ang_grid, &
      & n_points_final_grid,                          &
      !
      & n_points_rad_extra_grid,                                  &
      & n_points_ang_extra_grid,                                  &
      & n_points_rad_extra_grid*n_points_ang_extra_grid,          &
      & nucl_num*n_points_rad_extra_grid*n_points_ang_extra_grid, &
      !
      & n_points_rad_adapt_grid,                                  &                 
      & n_points_ang_adapt_grid,                                  &
      & n_points_rad_adapt_grid*n_points_ang_adapt_grid,          &
      & nucl_num*n_points_rad_adapt_grid*n_points_ang_adapt_grid, &
      !
      & nucl_num*(n_points_final_grid*(                          &
          & 1 + n_points_rad_adapt_grid*n_points_ang_adapt_grid) &
          & + n_points_rad_extra_grid*n_points_ang_extra_grid),  &
      !
      & my_alpha_knowles, my_radii_ua_av, &
      & diff_prun1_full3

  if (present(filename).and.ios==0) then
    close(unit_out)
  end if


end subroutine
