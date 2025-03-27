program test_ao_overlap_

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call write_int(6, my_n_pt_r_grid, 'radial  external grid over')
  call write_int(6, my_n_pt_a_grid, 'angular external grid over')

  my_extra_grid_becke  = .True.
  PROVIDE tc_grid2_a tc_grid2_r
  my_n_pt_r_extra_grid = tc_grid2_r
  my_n_pt_a_extra_grid = tc_grid2_a
  touch my_extra_grid_becke my_n_pt_r_extra_grid my_n_pt_a_extra_grid

  call write_int(6, my_n_pt_r_extra_grid, 'radial  internal grid over')
  call write_int(6, my_n_pt_a_extra_grid, 'angular internal grid over')

  call main()

end

! ---

subroutine main()
  implicit none
  double precision :: weight_r, r(3)
  double precision :: weight_rp, rp(3)
  double precision :: ao_i_r, ao_j_r
  double precision :: numeric
  double precision :: diag_overlap(ao_num)
  double precision :: diag_sum
  integer :: ir 
  integer :: i_ao, j_ao
  
  do i_ao = 1, ao_num
    do j_ao = i_ao, ao_num
    numeric = 0.d0
      do ir = 1, n_points_final_grid
        weight_r = final_weight_at_r_vector(ir)
        ao_i_r = aos_in_r_array_transp(ir,i_ao)
        ao_j_r = aos_in_r_array_transp(ir,j_ao)
        numeric += ao_j_r*ao_i_r*weight_r
      enddo
      print*, i_ao, j_ao
      print*, 'From numerical integration = ', numeric
      print*, 'From ao_overlap subroutine = ', ao_overlap(i_ao,j_ao)
   enddo
  enddo

  diag_overlap = 0.d0
  diag_sum = 0.d0
  do i_ao = 1, ao_num
    print*, "ao_overlap(i_ao,i_ao)", ao_overlap(i_ao,i_ao)
    diag_overlap(i_ao) = ao_overlap(i_ao,i_ao)
    diag_sum += ao_overlap(i_ao,i_ao)
  enddo 
  print*, "Sum over diagonal of overlap matrix ", diag_sum
end
