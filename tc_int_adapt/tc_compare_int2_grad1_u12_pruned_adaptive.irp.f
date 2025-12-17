!program tc_compare_int2_grad1_u12
!  BEGIN_DOC
!  END_DOC
!  implicit none
!  call compare_int2_grad1
!end program tc_compare_int2_grad1_u12
!
!
!subroutine compare_int2_grad1
!  use iso_fortran_env, only: output_unit
!  use io_test_interface
!  implicit none
!
!  integer :: i1
!  double precision :: r1(3)
!
!  ! Arrays necessary for adaptive grid
!  double precision, allocatable :: grid_float_points(:,:,:,:)
!  double precision, allocatable :: grid_fixed_weights(:)
!  double precision, allocatable :: grid_float_weights(:,:,:)
!  ! Dummy variables necessary for pruning on floating grid (not implemented yet)
!  integer :: n_fixed_pts_effective(nucl_num)
!  integer :: n_float_pts_effective
!  integer :: n_pts_effective_max
!
!  ! Arrays computed on the adaptive grid
!  double precision, allocatable :: int2_grad1_u12_vect_ao_at_r1(:,:,:)
!  double precision, allocatable :: int2_grad1_u12_sqrd_ao_at_r1(:,:)
!
!  ! Differences
!  double precision :: diff_numeric_adapt_AND_semi_analytic
!  double precision :: diff_numeric_adapt_AND_numeric_fixed
!  double precision :: diff_semi_analytic_AND_numeric_fixed
!
!
!  ! Adaptive grid arrays
!  allocate(grid_float_points(3, n_points_ang_float_grid, n_points_rad_float_grid, 1))
!  allocate(grid_fixed_weights(n_points_extra_final_grid))
!  allocate(grid_float_weights(n_points_ang_float_grid, n_points_rad_float_grid, 1))
!
!  ! Adaptive grid integrals allocation
!  allocate(int2_grad1_u12_vect_ao_at_r1(ao_num, ao_num, 3))
!  allocate(int2_grad1_u12_sqrd_ao_at_r1(ao_num, ao_num))
!
!
!  provide final_grid_points
!  provide final_grid_points_extra
!
!  ! Loop over all points r1 belonging to grid1. 
!  ! Compute and compare integrals for different r1
!  do i1 = 1, size(final_grid_points,2)
!    r1(1:3) = final_grid_points(1:3,i1)
!
!    ! Adaptive grid 
!    call get_adaptive_grid_with_pruned_fixed( r1                                            &
!                          , final_grid_points_extra, grid_float_points &                      
!                          , grid_fixed_weights, grid_float_weights        &         
!                          , n_fixed_pts_effective, n_float_pts_effective  &
!                          , n_pts_effective_max)
!    ! Adaptive grid integrals: initialization
!    int2_grad1_u12_vect_ao_at_r1(:,:,:) = 0.d0 
!    int2_grad1_u12_sqrd_ao_at_r1(:,:) = 0.d0
!    ! Adaptive grid integrals: computation
!    call get_int2_grad1_u12_ao_pruned_adaptive( i1, ao_num     &
!                                              , n_points_extra_final_grid    &
!                                              , n_points_ang_float_grid      &
!                              , n_points_rad_float_grid      &
!                              , 1                            &
!                              , final_grid_points_extra      & 
!                              , grid_float_points            &
!                              , grid_fixed_weights           &
!                              , grid_float_weights           &
!                              , int2_grad1_u12_vect_ao_at_r1 &
!                              , int2_grad1_u12_sqrd_ao_at_r1)
!    ! Output
!    write(*,*) "Iteration: ", i1
!
!    ! Comparison ADAPTIVE and SEMI-ANALYTICAL
!    call compute_dp_array_diff( int2_grad1_u12_sqrd_ao_at_r1              &
!                              , int2_grad1_u12_square_ao(:,:,i1)          &
!                              , show = .False.                            &
!                              , message = "Difference ADAPTIVE and ANALYTITCAL: "  &
!                              , diff=diff_numeric_adapt_AND_semi_analytic )
!    ! Comparison ADAPTIVE and FIXED GRID
!    call compute_dp_array_diff( int2_grad1_u12_sqrd_ao_at_r1           &
!                              , int2_grad1_u12_sqrd_ao_numeric(:,:,i1) &
!                              , show = .False.                         &
!                              , message = "Difference ADAPTIVE and FIXED: "  &
!                              , diff=diff_numeric_adapt_AND_numeric_fixed)
!    ! Comparison SEMI-ANALYTICAL and FIXED GRID
!    call compute_dp_array_diff( int2_grad1_u12_square_ao(:,:,i1)       &
!                              , int2_grad1_u12_sqrd_ao_numeric(:,:,i1) &
!                              , show = .False.                         &
!                              , message = "Difference ANALYTICAL and FIXED: "  &
!                              , diff=diff_semi_analytic_AND_numeric_fixed)
!
!    print*, ''
!  end do
!
!  write(*,'(A)') repeat('=', 70)
!
!end subroutine
