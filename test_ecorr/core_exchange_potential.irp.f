BEGIN_PROVIDER [ double precision, core_exchange_pot_numeric, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! Numerical evaluation of
 ! < i | V_x^{core} | j > = -\sum_{k\in\text{core}} <ik|1/r12|kj>
 !                        = \int dr1\int drp\phi_i(r1) V_x(r1,rp) phi_j(rp)
 !                        = \int dr1\int drp\phi_i(r1) (\phi_k(r1)\phi_j(rp)/|r-rp|) \phi_j(rp)
 END_DOC
 !
 integer :: ir, irp
 integer :: i, i_core
 integer :: i_mo, j_mo
 double precision :: weight_r, r(3)
 double precision :: weight_rp, rp(3)
 double precision :: mo_i_r, mo_j_rp, kernel
 double precision :: mos_array_r(mo_num), mos_array_rp(mo_num)
 double precision :: distance
 double precision :: v_x_core 
 !
 core_exchange_pot_numeric(:,:) = 0.d0
 !
! do ir =1, n_points_extra_final_grid
!   r(1:3) = final_grid_points_extra(1:3,ir)
!   weight_r = final_weight_at_r_vector_extra(ir) 
 do ir = 1, n_points_final_grid
   r(1:3) = final_grid_points(1:3,ir)
   weight_r = final_weight_at_r_vector(ir)
   call give_all_mos_at_r(r,mos_array_r)
   do irp =1, n_points_extra_final_grid
    rp(1:3) = final_grid_points_extra(1:3,irp)
    weight_rp = final_weight_at_r_vector_extra(irp) 
!   do irp =1, n_points_final_grid
!    rp(1:3) = final_grid_points(1:3,irp)
!    weight_rp = final_weight_at_r_vector(irp) 
    call give_all_mos_at_r(rp,mos_array_rp)
    distance = dsqrt(sum((r-rp)*(r-rp)))
    v_x_core = 0.d0
    do i = 1, n_core_pseudo_orb
      i_core = list_core_pseudo(i)
      v_x_core += mos_array_rp(i_core) * mos_array_r(i_core)
    enddo
    if(distance.gt.1.d-10)then
      v_x_core = -v_x_core/distance 
    else
      v_x_core = 0.d0
    endif
    do j_mo = 1, mo_num
      mo_j_rp = mos_array_rp(j_mo)
      do i_mo = 1, mo_num
        mo_i_r = mos_array_r(i_mo)
        core_exchange_pot_numeric(i_mo,j_mo) += v_x_core * mo_j_rp * mo_i_r * weight_rp * weight_r
      enddo
    enddo
   enddo
 enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, core_exchange_pot_exact, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! Analytical evaluation of
 ! <phi_i | V_x^{core_pseudo} | phi_j >
 END_DOC
 integer :: i_core,j,i,k
 double precision :: get_two_e_integral
 core_exchange_pot_exact(:,:) = 0.d0
 do k = 1, mo_num
  do j = 1, mo_num
   do i = 1, n_core_pseudo_orb
    i_core = list_core_pseudo(i)
    core_exchange_pot_exact(j,k) -= get_two_e_integral(i_core,j,k,i_core,mo_integrals_map) ! <ij|ki>
   enddo
  enddo
 enddo
END_PROVIDER 


