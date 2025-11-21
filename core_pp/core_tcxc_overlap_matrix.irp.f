BEGIN_PROVIDER [ double precision, ao_overlap_mat_grid1, (ao_num, ao_num, n_points_final_grid)]
  implicit none
  BEGIN_DOC
  ! AO overlap matrix weighted at different points on grid1
  ! $$ \chi_k(1)^* \chi_i(1) dr_1 $$
  ! column-wise version, should be faster than the row-wise one
  END_DOC
  integer :: i1, k, i
  double precision :: w1

  ! Initialization
  ao_overlap_mat_grid1(:,:,:) = 0.d0

  do i1 = 1, n_points_final_grid
    do k = 1, ao_num
      do i = 1, ao_num
        w1 = final_weight_at_r_vector(i1) 
        ao_overlap_mat_grid1(i,k,i1) = w1 * aos_in_r_array(i,i1) * aos_in_r_array(k,i1)
      end do
    end do
  end do

END_PROVIDER
