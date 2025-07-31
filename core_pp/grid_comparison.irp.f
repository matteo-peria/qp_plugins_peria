! Program written to compare old and new subroutines about the adaptive grid

program grid_comparison
  implicit none
  call main()
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

  integer :: i,j,k,l,ir
  double precision :: r(3)
  double precision :: difference


  ! Allocate OLD 
  allocate(old_grid_adapt_points(3, n_points_ang_extra_grid, n_points_rad_extra_grid, nucl_num+1))
  allocate(old_grid_adapt_weight(   n_points_ang_extra_grid, n_points_rad_extra_grid, nucl_num+1))

  ! Allocate NEW
  allocate(new_grid_float_points(3, n_points_ang_float_grid, n_points_rad_float_grid, 1))
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

  do ir = 1, n_points_final_grid
    !r(1:3) = 1.d0
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

    !do i = 2, n_points_rad_float_grid-1
    !  do j = 1, n_points_ang_float_grid
    !    print*, i, j
    !    print*, old_grid_adapt_points(1:3,j,i,1), new_grid_float_points(1:3,j,i,1)
    !    print*, old_grid_adapt_weight(j,i,1), new_grid_float_weight(j,i,1)
    !  end do
    !end do
    !stop


    call compute_dp_array_diff( old_grid_adapt_weight(:,:,1) &
                              , new_grid_float_weight(:,:,1) &
                              , show = .False.        &
                              , diff=difference       &
    )
    print*, "Difference in the floating weights: ", difference
!    do i = 1, n_points_rad_float_grid
!      do j = 1, n_points_ang_float_grid
!        print*, i, j, old_grid_adapt_weight(j,i,1), new_grid_float_weight(j,i,1)
!      end do
!    end do
!    print*, "Difference in the floating weights ", difference
!     
!    do i = 2, nucl_num
!      call compute_dp_array_diff( old_grid_adapt_weight(:,:,i) &
!                                , new_grid_fixed_weight(:,:,i-1) &
!                                , show = .False.        &
!                                , diff=difference       &
!      )
!      print*, "Iteration ", i, ". Difference in the fixed weights: ", difference
!    end do

  end do 

end subroutine main
