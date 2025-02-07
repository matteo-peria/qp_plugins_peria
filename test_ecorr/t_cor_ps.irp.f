program test_corr_pseudo_

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
  double precision :: mo_i_r, mo_j_rp, kernel
  double precision :: numeric , exchange_pot_kernel
  integer :: ir,irp
  integer :: i_mo, j_mo
  i_mo = 3
  j_mo = 3
!  ! < phi_i_mo | V_x | phi_j_mo > = -\sum_k in core_pseudo <i_mo k|1/r12|k j_mo>
!                                  = \int dr1 phi_i_mo(r1) \int drp V_x(r1,rp) phi_j_mo(rp)
  numeric = 0.d0
  do ir = 1,  n_points_final_grid
   ! you get x, y and z of the ith grid point
   r(1) = final_grid_points(1,ir)
   r(2) = final_grid_points(2,ir)
   r(3) = final_grid_points(3,ir)
   weight_r = final_weight_at_r_vector(ir)
   mo_i_r = mos_in_r_array_transp(ir,i_mo)
   do irp =1, n_points_extra_final_grid
    rp(1) = final_grid_points_extra(1,irp)
    rp(2) = final_grid_points_extra(2,irp)
    rp(3) = final_grid_points_extra(3,irp)
    weight_rp = final_weight_at_r_vector_extra(irp)
    mo_j_rp = mos_in_r_array_extra_transp(irp,j_mo)
    kernel = exchange_pot_kernel(r,rp)
    numeric+= kernel * mo_j_rp * mo_i_r * weight_rp * weight_r
   enddo
  enddo
  print*,'numeric       = ',numeric
  print*,'x_pot_mo_prov = ',x_pot_mo_prov(i_mo,j_mo)
end
