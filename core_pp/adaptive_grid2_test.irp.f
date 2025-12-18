program grid_comparison
  BEGIN_DOC
  ! Program written to compare old and new subroutines about the adaptive grid
  ! Check that when we impose the floating grid to have the same dimension
  ! as the extra grid, both old and new providers give the same results
  END_DOC
  implicit none
  write(*,*) "In this test we are comparing two version of the ADAPTIVE GRID"
  write(*,*) 
  write(*,*) "  - version 1: floating grid + FIXED GRID DEFINED AS EXTRA GRID"
  write(*,*) "  - version 2: floating grid + FIXED GRID DEFINED AS GRID2"
  write(*,*) 
  write(*,*) "********************* WARNING *********************"
  write(*,*) "in order for the test to pass, the 2 grids must have the same dimensions"
  write(*,*) "n_points_ang_extra_grid = n_points_ang_grid2"
  write(*,*) "n_points_rad_float_grid = n_points_rad_grid2"
  write(*,*) 
  write(*,*)  "Here it is a check: "
  write(*,*) 
  write(*,*) "n_points_ang_extra_grid = ", n_points_ang_extra_grid
  write(*,*) "n_points_ang_grid2      = ", n_points_ang_grid2
  write(*,*) 
  write(*,*) "n_points_rad_extra_grid = ", n_points_rad_extra_grid 
  write(*,*) "n_points_rad_grid2      = ", n_points_rad_grid2

  call main()

  write(*,*) "********************* WARNING *********************"
  write(*,*) "in order for the test to pass, the 2 grids must have the same dimensions"
  write(*,*) "n_points_ang_extra_grid = n_points_ang_grid2"
  write(*,*) "n_points_rad_float_grid = n_points_rad_grid2"
  write(*,*) 
  write(*,*)  "Here it is a reminder: "
  write(*,*) 
  write(*,*) "n_points_ang_extra_grid = ", n_points_ang_extra_grid
  write(*,*) "n_points_ang_grid2      = ", n_points_ang_grid2
  write(*,*) 
  write(*,*) "n_points_rad_extra_grid = ", n_points_rad_extra_grid 
  write(*,*) "n_points_rad_grid2      = ", n_points_rad_grid2

end program grid_comparison


subroutine main
  use io_test_interface
  implicit none
  integer :: i_ang, i_rad, i_nuc

  ! NEW ADAPTIVE GRIDS ARRAYS
  double precision, allocatable :: grid_float_points(:,:,:,:)
  double precision, allocatable :: grid_float_weight(:,:,:)
  double precision, allocatable :: grid_fixed_weight(:,:,:)
  ! To be used for pruning (NOT IMPLEMENTED YET BUT NECESSARY FOR THE CALL)
  integer :: n_fixed_pts_effective(nucl_num)
  integer :: n_float_pts_effective
  integer :: n_pts_effective_max

  ! NEW ADAPTIVE GRIDS ARRAYS
  double precision, allocatable :: grid_float_points2(:,:,:,:)
  double precision, allocatable :: grid_float_weight2(:,:,:)
  double precision, allocatable :: grid_fixed_weight2(:,:,:)
  ! To be used for pruning (NOT IMPLEMENTED YET BUT NECESSARY FOR THE CALL)
  integer :: n_fixed_pts_effective2(nucl_num)
  integer :: n_float_pts_effective2
  integer :: n_pts_effective_max2


  integer :: i,j,k,l,ir,d
  double precision :: r(3)
  double precision :: difference
  double precision :: total_difference_float_points
  double precision :: total_difference_float_weights
  double precision :: total_difference_fixed_weights

  ! Allocate adaptive grid defined on FLOAT and EXTRA
  allocate(grid_float_points(3, n_points_ang_float_grid, n_points_rad_float_grid, 1))
  ! Fixed EXTRA
  allocate(grid_float_weight(   n_points_ang_float_grid, n_points_rad_float_grid, 1))
  allocate(grid_fixed_weight(   n_points_ang_extra_grid, n_points_rad_extra_grid, nucl_num))

  ! Allocate adaptive grid defined on FLOAT and GRID2
  allocate(grid_float_points2(3, n_points_ang_float_grid, n_points_rad_float_grid, 1))
  ! Fixed GRID2
  allocate(grid_float_weight2(   n_points_ang_float_grid, n_points_rad_float_grid, 1))
  allocate(grid_fixed_weight2(   n_points_ang_grid2, n_points_rad_grid2, nucl_num))


  ! Initialize adaptive grid defined on FLOAT and EXTRA
  grid_float_points(:,:,:,:) = 0.d0
  grid_float_weight(:,:,:) = 0.d0
  grid_fixed_weight(:,:,:) = 0.d0

  ! Initialize adaptive grid defined on FLOAT and GRID2
  grid_float_points2(:,:,:,:) = 0.d0
  grid_float_weight2(:,:,:) = 0.d0
  grid_fixed_weight2(:,:,:) = 0.d0


  provide n_points_final_grid
  provide grid_points_extra_per_atom
  provide grid_points_per_atom2


  ! Initialize cumulative difference:
  total_difference_float_points  = 0.d0
  total_difference_float_weights = 0.d0
  total_difference_fixed_weights = 0.d0

  ! We do the test for several points in the usual integration grid, named
  ! 'final_grid_points',  by defining an adaptive grid for each point on grid1
  do ir = 1, n_points_final_grid
    r(1:3) = final_grid_points(1:3,ir)
    print*, "r: ", r

    ! Adaptive grid defined on FLOAT + EXTRA
    call get_adaptive_grid( r                          & 
                          , grid_points_extra_per_atom &
                          , grid_float_points      &
                          , grid_fixed_weight      &
                          , grid_float_weight      &
                          , n_fixed_pts_effective  &
                          , n_float_pts_effective  &
                          , n_pts_effective_max)

    ! Adaptive grid defined on FLOAT + GRID2
    call get_adaptive_grid2( r                     & 
                          , grid_points_per_atom2  &
                          , grid_float_points2     &
                          , grid_fixed_weight2     &
                          , grid_float_weight2     &
                          , n_fixed_pts_effective2 &
                          , n_float_pts_effective2 &
                          , n_pts_effective_max2)


    ! Testing function comparing 2D array
    call compute_dp_array_diff( grid_float_weight(:,:,1) &
                              , grid_float_weight2(:,:,1) &
                              , show = .False.        &
                              , diff=difference       &
    )
    print*, "Difference in the floating weights: ", difference
    total_difference_float_weights += difference

    ! At the moment the testing function is restricted to 2D array only, 
    ! so for 3-dimensional fixed weights we need to accumulate the difference
    do i = 1, nucl_num
      call compute_dp_array_diff( grid_fixed_weight(:,:,i) &
                                , grid_fixed_weight2(:,:,i) &
                                , show = .False.        &
                                , diff=difference       &
      )
      total_difference_fixed_weights += difference
    end do

    ! The same goes for comparing 4d array of cartesian coordinates array
    do d = 1,3
      call compute_dp_array_diff( grid_float_points(d,:,:,1) &
                                , grid_float_points2(d,:,:,1) &
                                , show = .False.        &
                                , diff=difference       &
      )
      total_difference_float_points += difference
    end do


  end do 

  write(*,*) "Total difference between the two version of ADAPTIVE GRID"
  write(*,*) "(one defined on EXTRA+FLOAT the other defined on GRID2+FLOAT"
  write(*,*) "total_difference_float_points",  total_difference_float_points
  write(*,*) "total_difference_float_weights", total_difference_float_weights
  write(*,*) "total_difference_fixed_weights",  total_difference_fixed_weights

end subroutine main
