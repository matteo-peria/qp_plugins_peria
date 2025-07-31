program tc_compare_grad1_u12_on_pruned_extra_grid
  BEGIN_DOC
  END_DOC
  implicit none
  call compare
end program tc_compare_grad1_u12_on_pruned_extra_grid


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


  ! Allocate OLD gradients
  allocate(grad1_u12_sqrd_old(n_points_extra_final_grid))
  allocate(grad1_u12_vect_old(n_points_extra_final_grid,3))

  ! Allocate NEW gradients
  allocate(grad1_u12_sqrd_new(n_points_extra_final_grid))
  allocate(grad1_u12_vect_new(n_points_extra_final_grid,3))

  diff_vect_overall = 0.d0
  diff_sqrd_overall = 0.d0

  provide final_grid_points

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
    call get_grad1_u12_withsq_r1_seq( i1                        &
                                    , n_points_extra_final_grid &
                                    , grad1_u12_vect_old(:,1)   &
                                    , grad1_u12_vect_old(:,2)   &
                                    , grad1_u12_vect_old(:,3)   &
                                    , grad1_u12_sqrd_old(:))

    ! NEW gradients
    print*, "Computing gradients in the new way"
    call get_grad1_u12_on_pruned_grid2( i1                        &
                                         , n_points_extra_final_grid &
                                         , final_grid_points_extra   &
                                         , grad1_u12_vect_new(:,1:3) &
                                         , grad1_u12_sqrd_new(:)     &
                                         )

    grad1_u12_sqrd_new = sum(grad1_u12_vect_new*grad1_u12_vect_new,dim=2)
    print*, "Computing difference"

    ! Output
    write(*,*) " Iteration:   ", i1
    write(*,*) " r1: ", r1

    ! Comparison OLD and NEW gradients VECTORS
    call compute_dp_array_diff( grad1_u12_vect_old(:,1:3) &
                              , grad1_u12_vect_new(:,1:3) &
                              , show = .False.            &
                              , message = "Diff vect: "   &
                              , diff=diff_vect            )
    ! Comparison OLD and NEW gradients SQUARED
    call compute_dp_array_diff( reshape(grad1_u12_sqrd_old(:), [1, size(grad1_u12_sqrd_new)]) &
                              , reshape(grad1_u12_sqrd_new(:), [1, size(grad1_u12_sqrd_new)]) &
                              , show = .False.            &
                              , message = "Diff sqrd: "   &
                              , diff=diff_sqrd        )

    diff_vect_overall = diff_vect_overall + diff_vect
    diff_sqrd_overall = diff_sqrd_overall + diff_sqrd

    print*, ''
  end do

  write(*,'(A)') repeat('=', 70)

  write(*,*) "OVERALL Difference in the VECTORIAL gradients: ", diff_vect_overall
  write(*,*) "OVERALL Difference in the SQUARED gradients:   ", diff_sqrd_overall



end subroutine
