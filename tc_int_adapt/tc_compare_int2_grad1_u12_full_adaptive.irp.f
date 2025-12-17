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
!  double precision, allocatable :: grid_fixed_weights(:,:,:)
!  double precision, allocatable :: grid_float_weights(:,:,:)
!  ! Dummy variables necessary for pruning on floating grid (not implemented yet)
!  integer :: n_fixed_pts_effective(nucl_num)
!  integer :: n_float_pts_effective
!  integer :: n_pts_effective_max
!
!  ! Arrays computed on the adaptive grid
!  double precision, allocatable :: int2_adapt_at_r1(:,:)
!
!  ! Differences
!  double precision :: diff_adapt_exact_sqrd
!  double precision :: diff_adapt_fixed_sqrd
!  double precision :: diff_fixed_exact_sqrd
!
!  double precision :: diff_adapt_exact_vect
!  double precision :: diff_adapt_fixed_vect
!  double precision :: diff_fixed_exact_vect
!
!
!  double precision :: diff_adapt_exact_sqrd_overall
!  double precision :: diff_adapt_fixed_sqrd_overall
!  double precision :: diff_fixed_exact_sqrd_overall
!
!  double precision :: diff_adapt_exact_vect_overall
!  double precision :: diff_adapt_fixed_vect_overall
!  double precision :: diff_fixed_exact_vect_overall
!
!
!  ! Adaptive grid arrays
!  allocate(grid_float_points(3, n_points_ang_float_grid, n_points_rad_float_grid, 1))
!  allocate(grid_fixed_weights(n_points_ang_extra_grid, n_points_rad_extra_grid, nucl_num))
!  allocate(grid_float_weights(n_points_ang_float_grid, n_points_rad_float_grid, 1))
!
!  ! Adaptive grid integrals allocation
!  allocate(int2_adapt_at_r1(ao_num, ao_num))
!
!
!  provide final_grid_points
!
!  diff_adapt_exact_sqrd_overall = 0.d0
!  diff_adapt_fixed_sqrd_overall = 0.d0
!  diff_fixed_exact_sqrd_overall = 0.d0
!
!  diff_adapt_exact_vect_overall = 0.d0
!  diff_adapt_fixed_vect_overall = 0.d0
!  diff_fixed_exact_vect_overall = 0.d0
!
!  integer :: i_rad, i_ang
!  double precision :: dist 
!
!  ! Loop over all points r1 belonging to grid1. 
!  ! Compute and compare integrals for different r1
!  do i1 = 1, size(final_grid_points,2)
!    r1(1:3) = final_grid_points(1:3,i1)
!
!    ! Adaptive grid 
!    call get_adaptive_grid( &
!    !call get_adaptive_grid_nosilence( &
!                           r1                                            &
!                         , grid_points_extra_per_atom, grid_float_points &                      
!                         , grid_fixed_weights, grid_float_weights        &         
!                         , n_fixed_pts_effective, n_float_pts_effective  &
!                         , n_pts_effective_max)
!
!    ! Adaptive grid integrals: initialization
!    int2_adapt_at_r1(:,:) = 0.d0
!    ! Adaptive grid integrals: computation
!
!    call get_int2_adapt_at_r1( i1            &
!                              , ao_num                      &
!                              , n_points_ang_extra_grid      &
!                              , n_points_rad_extra_grid      & 
!                              , nucl_num                     &
!                              , n_points_ang_float_grid      &
!                              , n_points_rad_float_grid      &
!                              , 1                            &
!                              , grid_points_extra_per_atom   & 
!                              , grid_float_points            &
!                              , grid_fixed_weights           &
!                              , grid_float_weights           &
!                              , int2_adapt_at_r1 &
!                              )
!    
!    diff_adapt_exact_sqrd = 0.d0
!    diff_adapt_fixed_sqrd = 0.d0
!    diff_fixed_exact_sqrd = 0.d0
!
!    diff_adapt_exact_vect = 0.d0
!    diff_adapt_fixed_vect = 0.d0
!    diff_fixed_exact_vect = 0.d0
!
!    ! Output
!    write(*,*) "Iteration: ", i1
!    print*, "R1 = ", r1
!    print*, ""
!    print*, "DIFFERENCES IN THE INTEGRALS INVOLVING THE SQUARED GRADIENT"
!
!    ! Comparison ADAPTIVE and SEMI-ANALYTICAL
!    call compute_dp_array_diff( int2_adapt_at_r1              &
!                              , int2_grad1_u12_square_ao(:,:,i1)          &
!                              , show = .False.                            &
!                              , message = "Diff ADAPTIVE and ANALYTICAL: "  &
!                              , diff=diff_adapt_exact_sqrd &
!                              )
!    ! Comparison ADAPTIVE and FIXED GRID
!    call compute_dp_array_diff( int2_adapt_at_r1           &
!                              , int2_grad1_u12_sqrd_ao_numeric(:,:,i1) &
!                              , show = .False.                         &
!                              , message = "Diff ADAPTIVE and NEW FIXED:  "  &
!                              , diff=diff_adapt_fixed_sqrd &
!                              )
!    diff_adapt_fixed_sqrd = 0.d0
!    ! Comparison ADAPTIVE and OLD FIXED GRID (original providers, should be the same as above)
!    call compute_dp_array_diff( int2_adapt_at_r1           &
!                              , int2_grad1_u12_square_ao_num(:,:,i1) &
!                              , show = .False.                         &
!                              , message = "Diff ADAPTIVE and OLD FIXED:  "  &
!                              , diff=diff_adapt_fixed_sqrd &
!                              )
!    ! Comparison SEMI-ANALYTICAL and FIXED GRID
!    call compute_dp_array_diff( int2_grad1_u12_square_ao(:,:,i1)       &
!                              , int2_grad1_u12_sqrd_ao_numeric(:,:,i1) &
!                              , show = .False.                         &
!                              , message = "Diff ANALYTICAL and NEW FIXED: "  &
!                              , diff=diff_fixed_exact_sqrd &
!                              )
!
!
!    print*, ''
!    print*, "DIFFERENCES IN THE INTEGRALS INVOLVING THE VECTORIAL GRADIENT"
!
!    print*, ''
!    ! XXXXXX
!    ! Comparison ADAPTIVE and SEMI-ANALYTICAL
!    call compute_dp_array_diff( int2_grad1_u12_vect_ao_at_r1(:,:,1) &
!                              , int2_grad1_u12_ao(:,:,i1,1) &
!                              , show = .False.                            &
!                              , message = "Diff ADAPTIVE and ANALYTICAL (X):  "  &
!                              , diff=diff_adapt_exact_vect &
!                              )
!    ! Comparison ANALYTICAL and OLD FIXED GRID (original providers, should be the same as above)
!    call compute_dp_array_diff( int2_grad1_u12_ao(:,:,i1,1) &
!                              , int2_grad1_u12_vect_ao_numeric(:,:,i1,1)  &
!                              , show = .False.                             &
!                              , message = "Diff FIXED and    ANALYTICAL (X):  "  &
!                              , diff=diff_fixed_exact_vect &
!                              )
!    print*, ''
!
!    ! YYYYYY
!    ! Comparison ADAPTIVE and SEMI-ANALYTICAL
!    call compute_dp_array_diff( int2_grad1_u12_vect_ao_at_r1(:,:,2) &
!                              , int2_grad1_u12_ao(:,:,i1,2) &
!                              , show = .False.                            &
!                              , message = "Diff ADAPTIVE and ANALYTICAL (Y):  "  &
!                              , diff=diff_adapt_exact_vect &
!                              )
!    ! Comparison ANALYTICAL and OLD FIXED GRID (original providers, should be the same as above)
!    call compute_dp_array_diff( int2_grad1_u12_ao(:,:,i1,2) &
!                              , int2_grad1_u12_vect_ao_numeric(:,:,i1,2)  &
!                              , show = .False.                             &
!                              , message = "Diff FIXED and    ANALYTICAL (Y):  "  &
!                              , diff=diff_fixed_exact_vect &
!                              )
!    print*, ''
!
!    ! ZZZZZZ
!    ! Comparison ADAPTIVE and SEMI-ANALYTICAL
!    call compute_dp_array_diff( int2_grad1_u12_vect_ao_at_r1(:,:,3) &
!                              , int2_grad1_u12_ao(:,:,i1,3) &
!                              , show = .False.                            &
!                              , message = "Diff ADAPTIVE and ANALZTICAL (Z):  "  &
!                              , diff=diff_adapt_exact_vect &
!                              )
!    ! Comparison ANALYTICAL and OLD FIXED GRID (original providers, should be the same as above)
!    call compute_dp_array_diff( int2_grad1_u12_ao(:,:,i1,3) &
!                              , int2_grad1_u12_vect_ao_numeric(:,:,i1,3)  &
!                              , show = .False.                             &
!                              , message = "Diff FIXED and    ANALYTICAL (Z):  "  &
!                              , diff=diff_fixed_exact_vect &
!                              )
!    print*, ''
!
!    ! FULL GRADIENT (X Y and Z)
!    ! Comparison ADAPTIVE and SEMI-ANALYTICAL
!    call compute_dp_array_diff( reshape(int2_grad1_u12_vect_ao_at_r1, [ao_num2, 3]) &
!                              , reshape(int2_grad1_u12_ao(:,:,i1,1:3), [ao_num2,3]) &
!                              , show = .False.                            &
!                              , message = "Diff ADAPTIVE and ANALYTICAL:  "  &
!                              , diff=diff_adapt_exact_vect &
!                              )
!    ! Comparison ADAPTIVE and FIXED GRID
!    call compute_dp_array_diff( reshape(int2_grad1_u12_vect_ao_at_r1, [ao_num2, 3]) &
!                              , reshape(int2_grad1_u12_vect_ao_numeric(:,:,i1,1:3), [ao_num2, 3]) &
!                              , show = .False.                         &
!                              , message = "Diff ADAPTIVE and NEW FIXED:   "  &
!                              , diff=diff_adapt_fixed_vect &
!                              )
!    diff_adapt_fixed_vect = 0.d0
!    ! Comparison ADAPTIVE and FIXED GRID (original providers, should be the same as above)
!    call compute_dp_array_diff( reshape(int2_grad1_u12_vect_ao_at_r1, [ao_num2, 3])      &
!                              , reshape(int2_grad1_u12_ao_num(:,:,i1,1:3), [ao_num2, 3]) &
!                              , show = .False.                                           &
!                              , message = "Diff ADAPTIVE and OLD FIXED:   "              &
!                              , diff=diff_adapt_fixed_vect                               &
!                              )
!    ! Comparison ANALYTICAL and OLD FIXED GRID (original providers, should be the same as above)
!    call compute_dp_array_diff( reshape(int2_grad1_u12_ao(:,:,i1,1:3) , [ao_num2, 3]) &
!                              , reshape(int2_grad1_u12_vect_ao_numeric(:,:,i1,1:3), [ao_num2, 3])  &
!                              , show = .False.                         &
!                              , message = "Diff ANALYTICAL and NEW FIXED: "  &
!                              , diff=diff_fixed_exact_vect &
!                              )
!
!    diff_adapt_exact_sqrd_overall = diff_adapt_exact_sqrd_overall + diff_adapt_exact_sqrd
!    diff_adapt_fixed_sqrd_overall = diff_adapt_fixed_sqrd_overall + diff_adapt_fixed_sqrd
!    diff_fixed_exact_sqrd_overall = diff_fixed_exact_sqrd_overall + diff_fixed_exact_sqrd
!
!    diff_adapt_exact_vect_overall = diff_adapt_exact_vect_overall + diff_adapt_exact_vect
!    diff_adapt_fixed_vect_overall = diff_adapt_fixed_vect_overall + diff_adapt_fixed_vect
!    diff_fixed_exact_vect_overall = diff_fixed_exact_vect_overall + diff_fixed_exact_vect
!
!
!    print*, ''
!    write(*,'(A)') repeat('=', 70)
!  end do
!
!  print*, ''
!  write(*,'(A)') repeat('=', 70)
!  print*, ''
!
!  print*, "Overall differences:"
!
!  print*, "ADAPT - EXACT (sqrd)", diff_adapt_exact_sqrd_overall
!  print*, "ADAPT - FIXED (sqrd)", diff_adapt_fixed_sqrd_overall
!  print*, "FIXED - EXACT (sqrd)", diff_fixed_exact_sqrd_overall
!  print*, ""
!  print*, "ADAPT - EXACT (vect)", diff_adapt_exact_vect_overall
!  print*, "ADAPT - FIXED (vect)", diff_adapt_fixed_vect_overall
!  print*, "FIXED - EXACT (vect)", diff_fixed_exact_vect_overall
!
!
!end subroutine
!
!
!
!subroutine get_int2_adapt_at_r1( i1, &
!                                , n_ang_fixed, n_rad_fixed, n_nuc_fixed  &
!                                , n_ang_float, n_rad_float, n_nuc_float  &
!                                , grid_fixed_points,  grid_float_points  &
!                                , grid_fixed_weights, grid_float_weights &
!                                , int2_adapt_at_r1 &
!  )
!
!  implicit none
!  ! INPUT
!  integer, intent(in) :: i1
!  integer, intent(in) :: n_ang_fixed, n_rad_fixed, n_nuc_fixed
!  integer, intent(in) :: n_ang_float, n_rad_float, n_nuc_float
!  double precision, intent(in) :: grid_float_points(3,n_ang_float,n_rad_float,n_nuc_float)
!  double precision, intent(in) :: grid_fixed_points(3,n_ang_fixed,n_rad_fixed,n_nuc_fixed)
!  double precision, intent(in) :: grid_float_weights(n_ang_float,n_rad_float,n_nuc_float)
!  double precision, intent(in) :: grid_fixed_weights(n_ang_fixed,n_rad_fixed,n_nuc_fixed)
!  ! OUTPUT
!  double precision, intent(out) :: int2_adapt_at_r1(n_aos,n_aos)
!
!  double precision :: r1(3)
!  integer :: j, l, m, i2
! 
!  double precision :: time0, time1
!  
!  double precision :: aos_at_r2(ao_num)
!  double precision, allocatable :: aos_lj_prod_in_r2(:,:,:)
!  double precision, allocatable :: grad1_u12_sqrd_at_r1(:)
!  double precision, allocatable :: grad1_u12_vect_at_r1(:,:)
!  logical, allocatable :: mask(:)
!
!
!  integer :: n_fixed, n_float, n_adaptive
!  integer :: r2_nuc, r2_rad, r2_ang
!  double precision :: ao_l_r2, ao_j_r2
!  double precision :: w2, r2(3)
!
!
!  ! Init the output
!  int2_adapt_at_r1(:,:)   = 0.d0
!
!
!  r1(1:3) = final_grid_points(1:3,i1)
!
!  ! Total number of points belonging to the floating grid
!  n_float = int(size(grid_float_points)/3)
!  ! Total number of points belonging to the fixed grid
!  n_fixed = n_nuc_fixed*(n_rad_fixed-1)*n_ang_fixed
!  ! Total number of points belonging to the fixed (extra) grid +  floating grid
!  n_adaptive = n_fixed + n_float 
!
!
!  call get_grad1_u12_on_full_grid2( i1, 1                       & 
!                             , n_points_rad_float_grid     &
!                             , n_points_ang_float_grid     &
!                             , grid_float_points           &
!                             , grad1_u12_vect_at_r1(:,1:3) &
!                             , grad1_u12_sqrd_at_r1(:)     &
!                             )
!
! call dgemm( "T", "N", ao_num*ao_num, 1, n_float &
!           , -0.5d0, aos_lj_prod_in_r2(:,:,:), n_float             &
!           , grad1_u12_sqrd_at_r1(:), n_float                 &
!           , 0.d0, int2_adapt_at_r1(:,:), ao_num*ao_num) 
!
!  ! This array store the scalars: \phi_i(r2) \phi_j(r2) dr2 for r2 in FIXED GRID
!  allocate(aos_lj_prod_in_r2(n_fixed, ao_num, ao_num))
!  !allocate(aos_lj_prod_in_r2(ao_num, ao_num, n_fixed))
!  aos_lj_prod_in_r2(:,:,:) = 0.d0
!  ! This array store the scalars: |\grad_r1 J(r1,r2)|^2 for r2 in FIXED GRID
!  allocate(grad1_u12_sqrd_at_r1(n_fixed))
!  ! This array store the vectors: \grad_r1 J(r1,r2) for r2 in FIXED GRID
!  allocate(grad1_u12_vect_at_r1(n_fixed,3))
!
!
!
!
!  call get_grad1_u12_on_full_grid2( i1                     &
!                             , n_nuc_fixed    &
!                             , n_rad_fixed    &
!                             , n_ang_fixed    &
!                             , grid_fixed_points           &
!                             , grad1_u12_vect_at_r1(:,1:3) &
!                             , grad1_u12_sqrd_at_r1(:)     &
!                             )
!
!
!
!end subroutine get_int2_adapt_at_r1
!
!
!pure function int2_
