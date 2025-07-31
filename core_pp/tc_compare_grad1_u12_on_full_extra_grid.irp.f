program tc_compare_grad1_u12_on_full_extra_grid
  BEGIN_DOC
  END_DOC
  implicit none
  call compare
end program tc_compare_grad1_u12_on_full_extra_grid


subroutine compare
  use io_test_interface
  implicit none
  integer :: i1
  double precision :: r1(3)

  ! OLD gradients
  double precision, allocatable :: grad1_u12_sqrd_old(:)
  double precision, allocatable :: grad1_u12_vect_old(:,:)
  ! NEW gradients
  double precision, allocatable :: grad1_u12_sqrd_new(:)
  double precision, allocatable :: grad1_u12_vect_new(:,:)

  ! Differences
  double precision :: diff_vect
  double precision :: diff_sqrd

  double precision :: diff_vect_overall
  double precision :: diff_sqrd_overall

  ! Differences with NEW provider
  double precision :: diff_vect_provider
  double precision :: diff_sqrd_provider

  double precision :: diff_vect_overall_provider
  double precision :: diff_sqrd_overall_provider


  integer :: size_extra_grid

  ! Size of extra grid is explictly computed instead of using
  !size_extra_grid = size(grid_points_extra_per_atom)
  ! Because in the current settings of Becke grids the radial dimension has
  ! one useless points (the final one). More over, we have to consider all 
  ! the nuclei in case we are dealing with a multiple atoms system
  size_extra_grid = nucl_num*(n_points_extra_radial_grid-1)*n_points_extra_integration_angular

  ! Allocate OLD gradients
  allocate(grad1_u12_sqrd_old(size_extra_grid))
  allocate(grad1_u12_vect_old(size_extra_grid,3))

  ! Allocate NEW gradients
  allocate(grad1_u12_sqrd_new(size_extra_grid))
  allocate(grad1_u12_vect_new(size_extra_grid,3))

  diff_vect_overall = 0.d0
  diff_sqrd_overall = 0.d0
  diff_vect_overall_provider = 0.d0
  diff_sqrd_overall_provider = 0.d0


  provide final_grid_points
  provide grid_points_extra_per_atom


  print*, "Loop over r1 to compute gradients of u12 for all r2"
  do i1 = 1, size(final_grid_points,2)
    r1(1:3) = final_grid_points(1:3,i1)

    ! NEW grid gradients: initialization
    grad1_u12_sqrd_old(:) = 0.d0
    grad1_u12_vect_old(:,:) = 0.d0

    ! NEW grid gradients: initialization
    grad1_u12_sqrd_new(:) = 0.d0
    grad1_u12_vect_new(:,:) = 0.d0
 
    ! OLD gradients
    print*, "Computing gradients in the old way"
    call get_grad1_u12_on_full_grid2_old( i1                      &
                                        , nucl_num                &
                                        , n_points_rad_extra_grid &
                                        , n_points_ang_extra_grid &
                                        , grad1_u12_vect_old(:,1) &
                                        , grad1_u12_vect_old(:,2) &
                                        , grad1_u12_vect_old(:,3) &
                                        , grad1_u12_sqrd_old(:))
    ! NEW gradients
    print*, "Computing gradients in the new way"
    call get_grad1_u12_on_full_grid2( i1                   &
                              , nucl_num                   &
                              , n_points_rad_extra_grid    &
                              , n_points_ang_extra_grid    &
                              , grid_points_extra_per_atom &
                              , grad1_u12_vect_new(:,1:3)  &
                              , grad1_u12_sqrd_new(:)      &
                              )   

    print*, "Computing difference"

    ! Output
    write(*,*) "Iteration: ", i1

    ! Comparison OLD and NEW gradients VECTORS
    call compute_dp_array_diff( grad1_u12_vect_old(:,1:3) &
                              , grad1_u12_vect_new(:,1:3) &
                              , show = .False.            &
                              , message = "Diff vect: "   &
                              , diff=diff_vect            )
    ! Comparison OLD and NEW gradients SQUARED
    call compute_dp_array_diff( reshape(grad1_u12_sqrd_old(:), [1, size_extra_grid]) &
                              , reshape(grad1_u12_sqrd_new(:), [1, size_extra_grid]) &
                              , show = .False.                                                &
                              , message = "Diff sqrd: "                                       &
                              , diff=diff_sqrd                                                )

    ! Comparison OLD gradients and the NEW PROVIDER computed with the new subroutine (SQUARE)
    call compute_dp_array_diff( reshape(grad1_u12_sqrd_old(:),              [1, size_extra_grid]) &
                              , reshape(grad1_u12_sqrd_on_full_extra(:,i1), [1, size_extra_grid]) &
                              , show = .False.                                                &
                              , message = "Diff sqrd (NEW PROVIDER): "                        &
                              , diff=diff_sqrd_provider        )
    ! Comparison OLD gradients and the NEW PROVIDER computed with the new subroutine (VECTORIAL)
    call compute_dp_array_diff( grad1_u12_vect_old(:,1:3)              &
                              , grad1_u12_vect_on_full_extra(:,i1,1:3) &
                              , show = .False.                         &
                              , message = "Diff vect (NEW PROVIDER): " &
                              , diff=diff_vect_provider        )

    diff_vect_overall = diff_vect_overall + diff_vect
    diff_sqrd_overall = diff_sqrd_overall + diff_sqrd

    diff_vect_overall_provider = diff_vect_overall_provider + diff_vect_provider
    diff_sqrd_overall_provider = diff_sqrd_overall_provider + diff_sqrd_provider

    print*, ''

  end do

  write(*,'(A)') repeat('=', 70)

  write(*,*) "OVERALL Difference in the VECTORIAL gradients: ", diff_vect_overall
  write(*,*) "OVERALL Difference in the SQUARED gradients:   ", diff_sqrd_overall

  write(*,*) "OVERALL Difference in the VECTORIAL gradients: ", diff_vect_overall_provider
  write(*,*) "OVERALL Difference in the SQUARED gradients:   ", diff_sqrd_overall_provider


end subroutine
