 BEGIN_PROVIDER [ double precision, core_xpot_numeric_full_extra_grid, (mo_num, mo_num)]
  implicit none
  BEGIN_DOC
  ! Numerical evaluation of < i | V_x^{core} | j >
  END_DOC
  !
  integer :: ir 
  !integer :: i_nucl, i_rad, i_ang 
  integer :: j_nucl, j_rad, j_ang 
  integer :: k, k_core
  integer :: i_mo, j_mo
  double precision :: weight_r, r(3)
  double precision :: weight_rp, rp(3)
  double precision :: mo_i_r, mo_j_rp
  double precision :: distance
  double precision :: v_x_core 

  core_xpot_numeric_full_extra_grid(:,:) = 0.d0
  
  ! Loop over first grid (pruned)
  do ir = 1, n_points_final_grid
    r(1:3) = final_grid_points(1:3,ir)
    weight_r = final_weight_at_r_vector(ir)
    ! Loop over full extra grid
    do j_nucl = 1, nucl_num
      do j_rad = 1, n_points_rad_extra_grid
        do j_ang = 1, n_points_ang_extra_grid
          rp(1:3) = grid_points_extra_per_atom(1:3,j_ang,j_rad,j_nucl)
          weight_rp = final_weight_at_r_extra(j_ang,j_rad,j_nucl)
          distance = dsqrt(sum((r-rp)*(r-rp)))
          v_x_core = 0.d0
          ! Loop over all core orbitals
          do k = 1, n_core_pseudo_orb
            k_core = list_core_pseudo(k)
            v_x_core += mos_in_r_full_extra_grid(k_core,j_ang,j_rad,j_nucl) &
                     & *mos_in_r_array(k_core,ir)
          enddo
          if(distance.gt.1.d-10)then
            v_x_core = -v_x_core/distance 
          else
            v_x_core = 0.d0
          endif
          do j_mo = 1, mo_num
            mo_j_rp = mos_in_r_full_extra_grid(j_mo,j_ang,j_rad,j_nucl)              
            do i_mo = 1, mo_num
              mo_i_r = mos_in_r_array(i_mo,ir)
              core_xpot_numeric_full_extra_grid(i_mo,j_mo) += &
                     & v_x_core * mo_j_rp * mo_i_r * weight_rp * weight_r
            enddo
          enddo
        end do
      end do
    end do
  end do
END_PROVIDER


BEGIN_PROVIDER [ double precision, core_xpot_numeric_adapt_old, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! Numerical evaluation of
 ! < i | V_x^{core} | j > = -\sum_{k\in\text{core}} <ik|1/r12|kj>
 !                        = \int dr1\int drp\phi_i(r1) V_x(r1,rp) phi_j(rp)
 !                        = \int dr1\int drp\phi_i(r1) (\phi_k(r1)\phi_j(rp)/|r-rp|) \phi_j(rp)
 END_DOC
 !
 integer :: ir
 integer :: i, i_core,j,k
 integer :: i_mo, j_mo
 integer :: i_nucl
 double precision :: weight_r, r(3)
 double precision :: weight_rp, rp(3)
 double precision :: mo_i_r, mo_j_rp, kernel
 double precision :: mos_array_r(mo_num), mos_array_rp(mo_num)
 double precision :: distance
 double precision :: v_x_core 

 double precision :: adaptive_grid(3,n_points_ang_extra_grid, &
                                  & n_points_rad_extra_grid,  &
                                  & nucl_num+1)
 double precision :: weights_points(n_points_ang_extra_grid, &
                                  & n_points_rad_extra_grid, &
                                  & nucl_num+1)
 adaptive_grid(:,:,:,:) = 0.d0
 weights_points(:,:,:) = 0.d0

 ! Init the provider
 core_xpot_numeric_adapt_old(:,:) = 0.d0

 print*, 'PROVIDING core_xpot_numeric_adapt_old'

 do ir = 1, n_points_final_grid
   r(1:3) = final_grid_points(1:3,ir)
   weight_r = final_weight_at_r_vector(ir)
   call give_adapt_grid_at_r(r, adaptive_grid, weights_points)
   do i_nucl = 1, nucl_num+1
     do j = 1, n_points_rad_extra_grid
       do k = 1, n_points_ang_extra_grid
        rp(1:3) = adaptive_grid(1:3,k,j,i_nucl)
        weight_rp = weights_points(k,j,i_nucl)
        ! call give_all_mos_at_moving_r(rp, mos_array_rp)
        call give_all_mos_at_r(rp,mos_array_rp)
        distance = dsqrt(sum((r-rp)*(r-rp)))
        v_x_core = 0.d0
        do i = 1, n_core_pseudo_orb
          i_core = list_core_pseudo(i)
          v_x_core += mos_array_rp(i_core) * mos_in_r_array(i_core,ir)
        enddo
        if(distance.gt.1.d-10)then
          v_x_core = -v_x_core/distance 
        else
          v_x_core = 0.d0
        endif
        do j_mo = 1, mo_num
          mo_j_rp = mos_array_rp(j_mo)
          do i_mo = 1, mo_num
            mo_i_r = mos_in_r_array(i_mo,ir)
            core_xpot_numeric_adapt_old(i_mo,j_mo) += &
                          & v_x_core * mo_j_rp * mo_i_r * weight_rp * weight_r
          enddo
        enddo
       enddo
     enddo
   enddo
 enddo
END_PROVIDER



BEGIN_PROVIDER [ double precision, core_xpot_numeric_full_adapt_grid, (mo_num, mo_num)]
  implicit none
  BEGIN_DOC
  ! Numerical evaluation of < i | V_x^{core} | j > 
  ! Full adaptive grid (no pruning) 
  ! extra and floating grid have the same resolution
  END_DOC
  !
  integer :: ir
  integer :: i, i_core,j,k
  integer :: i_mo, j_mo
  double precision :: weight_r, r(3)
  double precision :: weight_rp, rp(3)
  double precision :: mo_i_r, mo_j_rp, kernel
  double precision :: mos_array_r(mo_num)
  double precision :: mos_array_rp(mo_num)
  double precision :: distance
  double precision :: v_x_core 

  integer :: i_nucl
  integer :: k_core
  integer :: j_nucl, j_rad, j_ang 

  double precision :: grid_float_points(3, n_points_ang_float_grid, n_points_rad_float_grid, 1)
  double precision :: grid_float_weights(n_points_ang_float_grid, n_points_rad_float_grid, 1)

  double precision :: grid_fixed_weights(n_points_ang_extra_grid,n_points_rad_extra_grid,nucl_num)

  ! To be used for pruning
  integer :: n_fixed_pts_effective(nucl_num)
  integer :: n_float_pts_effective
  integer :: n_pts_effective_max


  print*, 'IM IN THE core_xpot_numeric_full_adapt_grid'

  ! Init the provider
  core_xpot_numeric_full_adapt_grid(:,:) = 0.d0

  grid_float_points(:,:,:,:) = 0.d0
  grid_float_weights(:,:,:) = 0.d0

  do ir = 1, n_points_final_grid
    r(1:3) = final_grid_points(1:3,ir)
    weight_r = final_weight_at_r_vector(ir)
    ! Retrieve the weights on the fixed and the floating grid
    call get_adaptive_grid(r, grid_points_extra_per_atom, grid_float_points, &                      
                         & grid_fixed_weights, grid_float_weights,           &         
                         & n_fixed_pts_effective, n_float_pts_effective,     &
                         & n_pts_effective_max)

    !do k=1, nucl_num
    !  do j_rad = 1, n_points_rad_extra_grid
    !    do j_ang = 1, n_points_ang_extra_grid
    !      print*, k, j_rad, j_ang, grid_fixed_weights(j_ang,j_rad,k)
    !    end do
    !  end do
    !end do

    ! Compute first contributions coming from the floating grid
    do j_rad = 1, n_points_rad_float_grid - 1
      do j_ang = 1, n_points_ang_float_grid
        rp(1:3) = grid_float_points(1:3,j_ang,j_rad,1)
        weight_rp = grid_float_weights(j_ang,j_rad,1)
        distance = dsqrt(sum((r-rp)*(r-rp)))
        call give_all_mos_at_r(rp,mos_array_rp)
        v_x_core = 0.d0
        do k = 1, n_core_pseudo_orb
          k_core = list_core_pseudo(k)
          v_x_core += mos_array_rp(k_core) * mos_in_r_array(k_core,ir)
        enddo
        if (distance.gt.1.d-10) then
          v_x_core = -v_x_core/distance 
        else
          v_x_core = 0.d0
        endif
        do j_mo = 1, mo_num
          mo_j_rp = mos_array_rp(j_mo)
          do i_mo = 1, mo_num
            mo_i_r = mos_in_r_array(i_mo,ir)

              if (isnan(v_x_core * mo_j_rp * mo_i_r * weight_rp * weight_r)) then
                print*, 'floating grid problem'
                print*, 'j_nucl, j_rad, j_ang, k, j_mo, i_mo'
                print*, j_nucl, j_rad, j_ang, k, j_mo, i_mo
                print*, v_x_core, mo_j_rp, mo_i_r, weight_rp, weight_r
                stop
              end if

            core_xpot_numeric_full_adapt_grid(i_mo,j_mo) += &
                         & v_x_core * mo_j_rp * mo_i_r * weight_rp * weight_r
          enddo
        enddo
      enddo
    enddo

    ! Compute contributions coming from the (fixed) extra grid
    do j_nucl = 1, nucl_num
      do j_rad = 1, n_points_rad_extra_grid-1
        do j_ang = 1, n_points_ang_extra_grid
          rp(1:3) = grid_points_extra_per_atom(1:3,j_ang,j_rad,j_nucl)
          weight_rp = grid_fixed_weights(j_ang,j_rad,j_nucl)
          distance = dsqrt(sum((r-rp)*(r-rp)))
          v_x_core = 0.d0
          ! Loop over all core orbitals
          do k = 1, n_core_pseudo_orb
            k_core = list_core_pseudo(k)
            v_x_core += mos_in_r_full_extra_grid(k_core,j_ang,j_rad,j_nucl) &
                     & *mos_in_r_array(k_core,ir)
          enddo
          if(distance.gt.1.d-10)then
            v_x_core = -v_x_core/distance 
          else
            v_x_core = 0.d0
          endif
          do j_mo = 1, mo_num
            mo_j_rp = mos_in_r_full_extra_grid(j_mo,j_ang,j_rad,j_nucl)              
            do i_mo = 1, mo_num
              mo_i_r = mos_in_r_array(i_mo,ir)
              if (isnan(v_x_core * mo_j_rp * mo_i_r * weight_rp * weight_r)) then
                print*, 'fixed grid problem'
                print*, j_nucl, j_rad, j_ang, k, j_mo, i_mo
                print*, v_x_core, mo_j_rp, mo_i_r, weight_rp, weight_r
                stop
              end if

              core_xpot_numeric_full_adapt_grid(i_mo,j_mo) += &
                     & v_x_core * mo_j_rp * mo_i_r * weight_rp * weight_r
            enddo
          enddo
        end do
      end do
    end do
  end do
END_PROVIDER
