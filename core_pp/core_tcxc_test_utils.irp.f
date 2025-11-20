! This file contains test-only providers, that is, providers that are just
! used to check that production providers behave well.
! These providers are numerical approximations

BEGIN_PROVIDER [ double precision, ao_overlap_grid1, (ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! TEST-ONLY PROVIDER.
  ! Numerical evaluation of AOs overlap over usual grid
  ! $$ \int d1 \chi_k(1)^* \chi_i(1) 
  END_DOC
  integer :: i, k, i1
  double precision :: w1

  ! Initialization
  ao_overlap_grid1(:,:) = 0.d0

  ! do-loop solution
  do k = 1, ao_num
    do i = 1, ao_num
      do i1 = 1, n_points_final_grid
        w1 = final_weight_at_r_vector(i1)
        ao_overlap_grid1(i,k) += w1 * aos_in_r_array(i,i1) * aos_in_r_array(k,i1)
      end do
    end do
  end do

  !! linear algebra faster solution:
  !! Hadamard product of aos_in_r_array*aos_in_r_array 
  !! and matrix product with grid weights array
  !call get_AB_prod( aos_in_r_array*aos_in_r_array &
  !                , ao_num, n_points_final_grid   &
  !                , final_weight_at_r_vector, 1   &
  !                , ao_overlap_grid1)
END_PROVIDER
