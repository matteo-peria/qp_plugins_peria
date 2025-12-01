program grid_comparison
  BEGIN_DOC
  ! Program written to compare old and new subroutines about the adaptive grid
  ! Check that when we impose the floating grid to have the same dimension
  ! as the extra grid, both old and new providers give the same results
  END_DOC
  implicit none
  write(*,*) "WARNING: for this test you need to set via qp edit"
  write(*,*) "n_points_ang_extra_grid = n_points_ang_float_grid"
  write(*,*) "n_points_rad_float_grid = n_points_rad_extra_grid"
  call main()
  write(*,*) "WARNING: remember that for this test you need to set via qp edit"
  write(*,*) "n_points_ang_extra_grid = n_points_ang_float_grid"
  write(*,*) "n_points_rad_float_grid = n_points_rad_extra_grid"

end program grid_comparison


subroutine main
  use io_test_interface
  implicit none
  integer :: i_ang, i_rad, i_nuc

  ! NEW ADAPTIVE GRIDS ARRAYS
  double precision, allocatable :: new_grid_float_points(:,:,:,:)
  double precision, allocatable :: new_grid_float_weight(:,:,:)
  double precision, allocatable :: new_grid_fixed_weight(:,:,:)
  ! To be used for pruning (NOT IMPLEMENTED YET BUT NECESSARY FOR THE CALL)
  integer :: n_fixed_pts_effective(nucl_num)
  integer :: n_float_pts_effective
  integer :: n_pts_effective_max

  ! OLD ADAPTIVE GRIDS ARRAYS
  double precision, allocatable :: old_grid_adapt_points(:,:,:,:)
  double precision, allocatable :: old_grid_adapt_weight(:,:,:)

  integer :: i,j,k,l,ir,d
  double precision :: r(3)
  double precision :: difference
  double precision :: total_difference_float_points
  double precision :: total_difference_float_weights
  double precision :: total_difference_fixed_weights


  ! Allocate OLD 
  ! (The OLD adaptive grid uses a floating grid with same dimensions as the extra grid)
  allocate(old_grid_adapt_points(3, n_points_ang_extra_grid, n_points_rad_extra_grid, nucl_num+1))
  ! (The OLD adaptive grid store the weights in a single array)
  allocate(old_grid_adapt_weight(   n_points_ang_extra_grid, n_points_rad_extra_grid, nucl_num+1))

  ! Allocate NEW
  ! (The NEW adaptive grid allows for a floating grid with different dimensions,
  !  but since the goal here is to compare the two, they need to have the same size)
  allocate(new_grid_float_points(3, n_points_ang_float_grid, n_points_rad_float_grid, 1))
  ! (The NEW adaptive grid store float and fixed weights on separate arrays)
  allocate(new_grid_float_weight(   n_points_ang_float_grid, n_points_rad_float_grid, 1))
  allocate(new_grid_fixed_weight(   n_points_ang_extra_grid, n_points_rad_extra_grid, nucl_num))

  ! Initialize OLD
  old_grid_adapt_points(:,:,:,:) = 0.d0
  old_grid_adapt_weight(:,:,:) = 0.d0

  ! Initialize NEW
  new_grid_float_points(:,:,:,:) = 0.d0
  new_grid_float_weight(:,:,:) = 0.d0
  new_grid_fixed_weight(:,:,:) = 0.d0

  provide n_points_final_grid
  provide grid_points_extra_per_atom


  ! Initialize cumulative difference:
  total_difference_float_points  = 0.d0
  total_difference_float_weights = 0.d0
  total_difference_fixed_weights = 0.d0

  ! We do the test for several points in the usual grid, 
  ! by defining an adaptive grid for each point on grid1
  do ir = 1, n_points_final_grid
    r(1:3) = final_grid_points(1:3,ir)
    print*, "r: ", r

    ! Old adaptive grid
    call give_adapt_grid_at_r(r, old_grid_adapt_points, old_grid_adapt_weight)
 
    ! New adaptive grid
    call get_adaptive_grid( r                   & 
                          , grid_points_extra_per_atom, new_grid_float_points &
                          , new_grid_fixed_weight, new_grid_float_weight        &
                          , n_fixed_pts_effective, n_float_pts_effective  &
                          , n_pts_effective_max)


    ! Testing function comparing 2D array
    call compute_dp_array_diff( old_grid_adapt_weight(:,:,1) &
                              , new_grid_float_weight(:,:,1) &
                              , show = .False.        &
                              , diff=difference       &
    )
    print*, "Difference in the floating weights: ", difference
    total_difference_float_weights += difference

    ! At the moment the testing function is restricted to 2D array only, 
    ! so for 3-dimensional fixed weights we need to accumulate the difference
    do i = 1, nucl_num
      call compute_dp_array_diff( old_grid_adapt_weight(:,:,i+1) &
                                , new_grid_fixed_weight(:,:,i) &
                                , show = .False.        &
                                , diff=difference       &
      )
      total_difference_fixed_weights += difference
    end do

    ! The same goes for comparing 4d array of cartesian coordinates array
    do d = 1,3
      call compute_dp_array_diff( old_grid_adapt_points(d,:,:,1) &
                                , new_grid_float_points(d,:,:,1) &
                                , show = .False.        &
                                , diff=difference       &
      )
      total_difference_float_points += difference
    end do


  end do 

  write(*,*) "Total difference between old and new definition of the adaptive grid"
  write(*,*) "total_difference_float_points",  total_difference_float_points
  write(*,*) "total_difference_float_weights", total_difference_float_weights
  write(*,*) "total_difference_fixed_weights",  total_difference_fixed_weights




end subroutine main
