program core_xpot_test

  BEGIN_DOC
  ! Testing the core exchange potential computed with different numerical grids
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


  !if (command_argument_count() < 1) then
  !  print *, "Usage: ./append_line_main <filename>"
  !  stop 1
  !end if

  !call get_command_argument(1, filename)


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

  integer :: ios
  integer :: unit_out

  !print*,'EXACT'
  !call print_db_array(core_xpot_exact, list_valence_pseudo, list_valence_pseudo)
  
  write(*,'(A)') repeat('=', 70)

  print*,'CORE_XPOT_NUMERIC'
!  call compute_dp_array_diff(core_xpot_numeric,            &
!                           & core_xpot_exact,              &
!                           & show = .False.,               &
!                           & row_indx=list_valence_pseudo, &
!                           & col_indx=list_valence_pseudo, &
!                           & diff=diff_prun1_prun2)
!
  write(*,'(A)') repeat('=', 70)

!  print*,'CORE_XPOT_NUMERIC_FULL_EXTRA_GRID'
!  call compute_dp_array_diff(core_xpot_numeric_full_extra_grid, &
!                           & core_xpot_exact, &
!                           & show = .False.,              &
!                           & row_indx=list_valence_pseudo,&
!                           & col_indx=list_valence_pseudo, &
!                           & diff=diff_prun1_full2)
!
!  write(*,'(A)') repeat('=', 70)


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

  ! not computed yet
  !print*,'CORE_XPOT_NUMERIC_PRUN_ADAPT_GRID'
  !call compute_dp_array_diff(core_xpot_numeric_prun_adapt_grid, &
  !                         & core_xpot_exact, &
  !                         & show = .False.,              &
  !                         & row_indx=list_valence_pseudo,       &
  !                         & col_indx=list_valence_pseudo, &
  !                         & diff=diff_prun1_prun3)

  !write(*,'(A)') repeat('=', 70)
  !diff_prun1_prun3 = -0.123456789

!!!!  logical :: file_exists
!!!!  if (present(filename)) then
!!!!    output = trim(filename)//'/core_xpot.out'
!!!!    inquire(file=output, exist=file_exists)
!!!!    
!!!!    if (file_exists) then
!!!!      open(newunit=unit_out, & 
!!!!         & file=output,    & 
!!!!         & status='old',     & 
!!!!         & position='append',&  
!!!!         & action='write',   & 
!!!!         & iostat=ios        & 
!!!!         &)
!!!!    else
!!!!      open(newunit=unit_out, & 
!!!!         & file=output,    & 
!!!!         & status='new',     & 
!!!!         & action='write',   & 
!!!!         & iostat=ios        & 
!!!!         &)
!!!!    end if
!!!!
!!!!    if (ios /= 0) then
!!!!      print *, "Error opening file: ", trim(filename)
!!!!      stop 1
!!!!    end if
!!!!
!!!!  else
!!!!    unit_out=output_unit
!!!!  end if
!!!!
!!!!  character(len=100) :: fmt_str
!!!!  fmt_str = '(14A8, 6A12)'
!!!!  if (.not.file_exists) then
!!!!    write(unit_out, fmt_str) &
!!!!                       & 'N_rad1', 'N_ang1', 'N_sph1', 'N_ful1', 'N_pru1', &
!!!!                       & 'N_rad2', 'N_ang2', 'N_sph2', 'N_ful2', 'N_pru2', &
!!!!                       & 'N_rad3', 'N_ang3', 'N_sph3', 'N_ful3',           &
!!!!                       & 'TOT(p1+f2)', &
!!!!                       & 'RES(p1+f2)', &
!!!!                       & 'TOT(p1+p2)', &
!!!!                       & 'RES(p1+p2)', &
!!!!                       & 'TOT(p1+f3)', &
!!!!                       & 'RES(p1+f3)'
!!!!  end if
!!!!
!!!!  character(len=100) :: fmt_num
!!!!  fmt_num = '(14I8, I12, E12.3, I12, E12.3, I12, E12.3)'
!!!!
!!!!  write(unit_out, fmt_num)                                &
!!!!      & n_points_rad_grid,                            &
!!!!      & n_points_ang_grid,                            &
!!!!      & n_points_rad_grid*n_points_ang_grid,          &
!!!!      & nucl_num*n_points_rad_grid*n_points_ang_grid, &
!!!!      & n_points_final_grid,                          &
!!!!      !
!!!!      & n_points_rad_extra_grid,                                  &
!!!!      & n_points_ang_extra_grid,                                  &
!!!!      & n_points_rad_extra_grid*n_points_ang_extra_grid,          &
!!!!      & nucl_num*n_points_rad_extra_grid*n_points_ang_extra_grid, &
!!!!      & n_points_extra_final_grid,                                &
!!!!      !
!!!!      & n_points_rad_adapt_grid,                                  &                 
!!!!      & n_points_ang_adapt_grid,                                  &
!!!!      & n_points_rad_adapt_grid*n_points_ang_adapt_grid,          &
!!!!      & nucl_num*n_points_rad_adapt_grid*n_points_ang_adapt_grid, &
!!!!      !
!!!!      & nucl_num*(n_points_final_grid                               &
!!!!             &  + n_points_rad_extra_grid*n_points_ang_extra_grid), &
!!!!      & diff_prun1_full2,                                           &
!!!!      !
!!!!      & nucl_num*(n_points_final_grid + n_points_extra_final_grid), &
!!!!      & diff_prun1_prun2,                                           &
!!!!      !
!!!!      & nucl_num*(n_points_final_grid*(                          &
!!!!          & 1 + n_points_rad_adapt_grid*n_points_ang_adapt_grid) &
!!!!          & + n_points_rad_extra_grid*n_points_ang_extra_grid),  &
!!!!      & diff_prun1_full3
!!!!
!!!!  if (present(filename).and.ios==0) then
!!!!    close(unit_out)
!!!!  end if
!!!!
!!!!!
!!!!!
!!!!!  write(*,'(A14, A14, A14, A14)') &
!!!!!      &'Type', 'Grid 1', 'Grid 2', 'Grid 3'
!!!!!
!!!!!  write(*,'(A14, I14, I14, I14)') &
!!!!!      & 'Radial',                 &
!!!!!      & n_points_rad_grid,        &
!!!!!      & n_points_rad_extra_grid,  &
!!!!!      & n_points_rad_adapt_grid
!!!!!
!!!!!  write(*,'(A14, I14, I14, I14)') &
!!!!!      & 'Angular',                &
!!!!!      & n_points_ang_grid,        &
!!!!!      & n_points_ang_extra_grid,  &
!!!!!      & n_points_ang_adapt_grid
!!!!!
!!!!!  write(*,'(A14, I14, I14, I14)')                        &
!!!!!      & 'Tot sphe',                                      & 
!!!!!      & n_points_rad_grid*n_points_ang_grid,             &
!!!!!      & n_points_rad_extra_grid*n_points_ang_extra_grid, &
!!!!!      & n_points_rad_adapt_grid*n_points_ang_adapt_grid
!!!!!
!!!!!  write(*,'(A14, I14, I14, I14)')                                   &
!!!!!      & 'Tot molec',                                                & 
!!!!!      & nucl_num*n_points_rad_grid*n_points_ang_grid,                 &
!!!!!      & nucl_num*n_points_rad_extra_grid*n_points_ang_extra_grid,     &
!!!!!      & nucl_num*n_points_rad_grid*n_points_ang_grid*(1               &
!!!!!          + nucl_num*n_points_rad_adapt_grid*n_points_ang_adapt_grid) &
!!!!!          + nucl_num*n_points_rad_extra_grid*n_points_ang_extra_grid
!!!!!
!!!!!  write(*,'(A14, I14, I14, A14)')  &
!!!!!      & 'Tot pruned',              &
!!!!!      & n_points_final_grid,       &
!!!!!      & n_points_extra_final_grid, &
!!!!!      & 'no pruning'
!!!!!
!!!!!
!!!!!  write(*,'(A)') repeat('=', 70)
!!!!!
!!!!!  write(*,'(A20, A20, A20, A20, A20)') &
!!!!!      & 'Grid combo',             &
!!!!!      & 'Pruned_1 + Full_2',    &
!!!!!      & 'Pruned_1 + Pruned_2',    &
!!!!!      & 'Pruned_1 + Full_3',    &
!!!!!      & 'Pruned_1 + Pruned_3'
!!!!!  
!!!!!  write(*,'(A20, I20, I20, I20, I20)')                                &
!!!!!      & 'N of points',                                                &
!!!!!      & n_points_final_grid                  &
!!!!!          + nucl_num*n_points_rad_extra_grid*n_points_ang_extra_grid, &
!!!!!      & n_points_final_grid                                           &
!!!!!          + n_points_extra_final_grid, &
!!!!!      & n_points_final_grid*(1                                      &
!!!!!      &   + nucl_num*n_points_rad_adapt_grid*n_points_ang_adapt_grid) &
!!!!!          + nucl_num*n_points_rad_extra_grid*n_points_ang_extra_grid, &
!!!!!      & -9999999
!!!!!
!!!!!  write(*,'(A20, E20.10, E20.10, E20.10, E20.10)') & 
!!!!!      & 'Diff with exact',             &
!!!!!      & diff_prun1_full2,              &
!!!!!      & diff_prun1_prun2,              &
!!!!!      & diff_prun1_full3,              &
!!!!!      & diff_prun1_prun3
!!!!!
!!!!
!!!!!!!!!  write(*,*) 'grid_points_radial', grid_points_radial
!!!!!!!!!
!!!!!!!!!  integer :: j
!!!!!!!!!  double precision :: x, r
!!!!!!!!!  double precision, external :: knowles_function
!!!!!!!!!  do j = 1, n_points_radial_grid-1                                                                
!!!!!!!!!    x = grid_points_radial(j)                                                                     
!!!!!!!!!    r = knowles_function(alpha_knowles(grid_atomic_number(1)), m_knowles, x)                      
!!!!!!!!!    write(*,*) j, x, r
!!!!!!!!!  end do
!!!!
!!!!
end subroutine
