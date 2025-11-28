BEGIN_PROVIDER [ double precision, core_tcxc_adapt_grid12aj, (ao_num, ao_num, ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! Numerical evaluation of < kl | V^{\text{TC}}_{x,\text{core}}  | ij >
  !
  ! $$\braket{kl | V^{\text{TC}}_{x,\text{core}} | ij}
  !     \braket{kl | e^{-J} V_{x,\text{core}} e^{+J} | ij}
  !       \int d1 \chi_k(1)^* \chi_i(1) 
  !         \int d2 \chi_l^*(2) e^{-u(1,2)} 
  !           \int d2' v_x^{\text{core}}(2,2') e^{u(1,2')} \chi_j(2') $$
  END_DOC
  !
  integer :: i1, i2
    ! Grids indices
  double precision :: w1, r1(3)
    ! Grid 1 point and weight
  double precision :: w2, r2(3)
    ! Grid 2 point and weight
  double precision :: w2p, r2p(3)
    ! Grid 3 point and weight
  integer :: i, j, k, l
    ! AOs indices
  double precision :: ao_i_r1, ao_k_r1, ao_l_r2, ao_j_r2p
    ! Real values of AOs
  integer :: m, m_core
    ! Core-MOs indices
  double precision :: j_factor, j_r1r2, j_r1r2p
    ! Jastrow factors between r1 and r2 and between r1 and r2p
  double precision, external :: j_mu_env!, core_tcxc_kernel_grid123
    ! Temporary external definition, we should put them in a module an import that
  double precision :: integral, kernel
    ! Temporary variables to stock mid-computation result
  double precision :: distance
    ! Distance for Coulomb integral (exchange)
  integer :: ao_num2

  ! ADAPTIVE GRID VARIABLES
  double precision :: mo_orbitals_in_r2p(mo_num)
    ! MO orbitals on the floating grid
  double precision :: ao_orbitals_in_r2p(ao_num)
    ! AO orbitals on the floating grid
  double precision :: grid_float_points(3, n_points_ang_float_grid, n_points_rad_float_grid, 1)
    ! Floating grid points
  double precision :: grid_float_weights(n_points_ang_float_grid, n_points_rad_float_grid, 1)
    ! Floating grid weights
  !double precision :: grid_fixed_weights(grid3_ang_size, grid3_rad_size, nucl_num)
  ! NOW THE EXTRA GRID IS USED AS THIRD GRID
  double precision :: grid_fixed_weights(n_points_extra_integration_angular,n_points_extra_radial_grid,nucl_num)
    ! Fixed grid weights
  integer :: i2p_nuc, i2p_rad, i2p_ang
    ! Floating grid indices
  integer :: i2p
    ! Grid 3 index
  double precision :: aos_in_r2p(ao_num), mos_in_r2p(mo_num), mos_core_in_r2p(n_core_pseudo_orb)
    ! AOs and MOs computed at points belonging to the floating part of grid3
    ! Depending on how fast we want to go, we either compute all MOs or just the core ones
    
  ! Fixed extra grid weights
  ! To be used for pruning, NOT IMPLEMENTED YET
  integer :: n_fixed_pts_effective(nucl_num)
  integer :: n_float_pts_effective
  integer :: n_pts_effective_max

  !double precision :: sum_w1, sum_ao_i_r1, sum_ao_k_r1, sum_w2, sum_ao_l_r2, sum_integral
  !
  !sum_w1=0.d0
  !sum_ao_i_r1=0.d0
  !sum_ao_k_r1=0.d0
  !sum_w2=0.d0
  !sum_ao_l_r2=0.d0
  !sum_integral=0.d0

  ! Initialize floating grid points and weights
  grid_float_points(:,:,:,:) = 0.d0
  grid_float_weights(:,:,:) = 0.d0

  core_tcxc_adapt_grid12aj(:,:,:,:) = 0.d0

  if (core_tcxc_loops) then
    do i1 = 1, n_points_final_grid
      r1(1:3) = final_grid_points(1:3,i1)
      w1 = final_weight_at_r_vector(i1)
      !write(*,*) "i1 = ", i1, "; w1 = ", w1
      do i2 = 1, n_points_final_grid2
        r2(1:3) = final_grid_points2(1:3,i2)
        w2 = final_weight_at_r_vector2(i2)
        !write(*,*) "i2 = ", i2, "; w2 = ", w2
        ! Compute pair Jastrow factor between r1 and r2
        if (core_tcxc_j0_testing) then
          j_r1r2 = 1.d0
        else
          !j_factor = exp(-j_mu_env(r1,r2,mu_erf))
          j_r1r2 = exp(-j_mu_env(r1,r2,mu_erf))
        end if

        ! Compute adaptive grid centered in r2
        call get_adaptive_grid( r2, grid_points_extra_per_atom, grid_float_points &
                              , grid_fixed_weights, grid_float_weights            &
                              , n_fixed_pts_effective, n_float_pts_effective      &
                              , n_pts_effective_max)

        ! FLOATING grid contributions
        do i2p_rad = 1, n_points_rad_float_grid - 1
          do i2p_ang = 1, n_points_ang_float_grid
            ! Find r2p point from i2p_rad, i2p_ang loop-indices
            r2p(1:3) = grid_float_points(1:3,i2p_ang,i2p_rad,1)
            distance = norm2(r2(1:3) - r2p(1:3))
            ! Compute matrix element only when there is no divergence 
            if (distance.gt.1.d-10) then
              ! Find r'_2 weight from i2p loop-index
              w2p = grid_float_weights(i2p_ang,i2p_rad,1)
              !write(*,*) "i2p_rad = ", i2p_rad, "; i2p_ang = ", i2p_ang, "; w2p = ", w2p
              ! Compute all AOs in r2p 
              call give_all_aos_at_r(r2p, aos_in_r2p)

              ! Get the core-MOs at r2p reciclying the newly computed AOs in r2p
              ! and using only the the sub-matrix of core-MOs coefficients
              call dgemv( 'N', n_core_pseudo_orb, ao_num, 1.d0 &
                        , mo_core_coef_notnorm_transp, n_core_pseudo_orb &
                        , aos_in_r2p, 1, 0.d0, mos_core_in_r2p, 1)

              ! Compute pair Jastrow factor between r1 and r2p
              if (core_tcxc_j0_testing) then
                j_r1r2p = 1.d0
              else
                j_r1r2p = exp(j_mu_env(r1,r2p,mu_erf))
              end if

              ! Initialize kernel
              kernel = 0.d0
              ! Loop over core orbitals to update the kernel
              do m = 1, n_core_pseudo_orb
                m_core = list_core_pseudo(m)
                !! Update kernel using all-MOs(r2p) array
                !kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r2p(m_core)    
                ! Update kernel using core-MOs(r2p) array
                kernel -= mos_in_r_array2_omp(m_core,i2) * mos_core_in_r2p(m)    
              enddo
              kernel = kernel / distance

              ! Loop over all AO (j-th index, variable is r2p)
              do j = 1, ao_num
                ao_j_r2p = aos_in_r2p(j)
                do l = 1, ao_num
                  ao_l_r2 = aos_in_r_array2(l,i2)
                  do k = 1, ao_num
                    ao_k_r1 = aos_in_r_array(k,i1)
                    do i = 1, ao_num
                      ao_i_r1 = aos_in_r_array(i,i1)
        
                      core_tcxc_adapt_grid12aj(i,k,l,j) += w1 * ao_i_r1 * ao_k_r1            &
                                                         * w2 * j_r1r2 * ao_l_r2             &
                                                         * w2p * kernel * ao_j_r2p * j_r1r2p
                    enddo
                  enddo
                enddo
              enddo

            end if

          end do  ! grid 3 floating
        end do    ! grid 3 floating

        ! FIXED GRID3 contributions 
        ! (this grid3 is based on a provider called grid2 because we are recycling)
        do i2p_nuc = 1, nucl_num
          do i2p_rad = 1, n_points_extra_radial_grid - 1
            do i2p_ang = 1, n_points_extra_integration_angular
              r2p(1:3) = grid_points_extra_per_atom(1:3,i2p_ang,i2p_rad,i2p_nuc)
              distance = norm2(r2(1:3) - r2p(1:3))
              ! Compute matrix element only when there is no divergence 
              if (distance.gt.1.d-10) then
                ! NOTICE: THE WEIGHTS ON THE EXTRA GRID ARE UPDATED FOR THE ADAPTIVE GRID 
                w2p = grid_fixed_weights(i2p_ang,i2p_rad,i2p_nuc)

                ! Compute pair Jastrow factor between r1 and r2p
                if (core_tcxc_j0_testing) then
                  !j_factor = 1.d0
                  j_r1r2p = 1.d0
                else
                  !j_factor = exp(j_mu_env(r1,r2p,mu_erf))
                  j_r1r2p = exp(j_mu_env(r1,r2p,mu_erf))
                end if

                ! Initialize kernel
                kernel = 0.d0
                ! Loop over core orbitals to update the kernel
                do m = 1, n_core_pseudo_orb
                  m_core = list_core_pseudo(m)
                  kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r_array_extra_full_omp(m_core,i2p_ang,i2p_rad,i2p_nuc)
                enddo
                kernel = kernel/distance

                ! Loop over all AO (j-th index, variable is r2p)
                do j = 1, ao_num
                  ao_j_r2p = aos_in_r_array_extra_full(j,i2p_ang,i2p_rad,i2p_nuc)
                  do l = 1, ao_num
                    ao_l_r2 = aos_in_r_array2(l,i2)
                    do k = 1, ao_num
                      ao_k_r1 = aos_in_r_array(k,i1)
                      do i = 1, ao_num
                        ao_i_r1 = aos_in_r_array(i,i1)
          
                        core_tcxc_adapt_grid12aj(i,k,l,j) += w1 * ao_i_r1 * ao_k_r1            &
                                                           * w2 * j_r1r2 * ao_l_r2             &
                                                           * w2p * kernel * ao_j_r2p * j_r1r2p
                      enddo   ! i
                    enddo   ! k
                  enddo   ! l
                enddo   ! j

              end if

            end do  ! grid 3 fixed
          end do    ! grid 3 fixed
        end do      ! grid 3 fixed

      enddo ! grid 2
    enddo ! grid 1

  else
    ! dgemm session
    ao_num2 = ao_num*ao_num
    call dgemm( 'N', 'T', ao_num2, ao_num2, n_points_final_grid, 1.d0 &
              , ao_overlap_mat_grid1, ao_num2    &
              , ao_core_xc_mat_grid12a, ao_num2    &
              , 0.d0, core_tcxc_adapt_grid12aj, ao_num2 )
  end if

END_PROVIDER


!BEGIN_PROVIDER [ double precision, ao_core_xc_mat_grid12a, (ao_num, ao_num, n_points_final_grid) ]
!  integer :: i1, i2, i2p
!    ! Grids indices
!  integer :: j, l
!    ! AOs indices
!  double precision :: r1(3)
!    ! Grid 1 point
!  double precision :: w2, r2(3)
!    ! Grid 2 point and weight
!  double precision :: r2p(3), w2p
!    ! Nested grid 2' point (r'_2) and weight (dr'_2)
!  double precision :: ao_l_r2, ao_j_r2p
!    ! Real values of AOs and core-exchange kernel
!  double precision :: j_r1r2, j_r1r2p
!    ! Jastrow factor between r1 and r2
!  double precision, external :: j_mu_env
!    ! Temporary external definition, we should put them in a module an import that
!  double precision :: distance
!    ! Distance for Coulomb integral
!  integer :: m, m_core
!    ! Indices to loop over core orbitals
!  double precision :: kernel
!    ! Kernel of the exchange potential
!
!  ! ADAPTIVE GRID VARIABLES
!  double precision :: mo_orbitals_in_r2p(mo_num)
!    ! MO orbitals on the floating grid
!  double precision :: ao_orbitals_in_r2p(ao_num)
!    ! AO orbitals on the floating grid
!  double precision :: grid_float_points(3, n_points_ang_float_grid, n_points_rad_float_grid, 1)
!    ! Floating grid points
!  double precision :: grid_float_weights(n_points_ang_float_grid, n_points_rad_float_grid, 1)
!    ! Floating grid weights
!  !double precision :: grid_fixed_weights(grid3_ang_size, grid3_rad_size, nucl_num)
!  ! NOW THE EXTRA GRID IS USED AS THIRD GRID
!  double precision :: grid_fixed_weights(n_points_extra_integration_angular,n_points_extra_radial_grid,nucl_num)
!    ! Fixed grid weights
!  integer :: i2p_nuc, i2p_rad, i2p_ang
!    ! Third grid indices
!  double precision :: aos_in_r2p(ao_num), mos_in_r2p(mo_num), mos_core_in_r2p(n_core_pseudo_orb)
!    ! AOs and MOs computed at points belonging to the floating part of grid3
!    ! Depending on how fast we want to go, we either compute all MOs or just the core ones
!    
!  ! Fixed extra grid weights
!  ! To be used for pruning, NOT IMPLEMENTED YET
!  integer :: n_fixed_pts_effective(nucl_num)
!  integer :: n_float_pts_effective
!  integer :: n_pts_effective_max
!
!  ! Initialize floating grid points and weights
!  grid_float_points(:,:,:,:) = 0.d0
!  grid_float_weights(:,:,:) = 0.d0
!  grid_fixed_weights(:,:,:) = 0.d0
!
!  ao_core_xc_mat_grid12a(:,:,:) = 0.d0
!
!  do i1 = 1, n_points_final_grid
!    r1(1:3) = final_grid_points(1:3,i1)
!
!    do i2 = 1, n_points_final_grid2
!      r2(1:3) = final_grid_points2(1:3,i2)
!      w2 = final_weight_at_r_vector2(i2)
!      ! Compute pair Jastrow factor between r1 and r2
!      if (core_tcxc_j0_testing) then
!        j_r1r2 = 1.d0
!      else
!        j_r1r2 = exp(-j_mu_env(r1,r2,mu_erf))
!      end if
!
!      ! Compute adaptive grid centered in r2
!      call get_adaptive_grid( r2, grid_points_extra_per_atom, grid_float_points &
!                            , grid_fixed_weights, grid_float_weights            &
!                            , n_fixed_pts_effective, n_float_pts_effective      &
!                            , n_pts_effective_max)
!
!      ! FLOATING grid contributions
!      do i2p_rad = 1, n_points_rad_float_grid - 1
!        do i2p_ang = 1, n_points_ang_float_grid
!          ! Find r2p point from i2p_rad, i2p_ang loop-indices
!          r2p(1:3) = grid_float_points(1:3,i2p_ang,i2p_rad,1)
!          distance = norm2(r2(1:3) - r2p(1:3))
!          ! Compute matrix element only when there is no divergence 
!          if (distance.gt.1.d-10) then
!            ! Find r'_2 weight from i2p loop-index
!            w2p = grid_float_weights(i2p_ang,i2p_rad,1)
!            !write(*,*) "i2p_rad = ", i2p_rad, "; i2p_ang = ", i2p_ang, "; w2p = ", w2p
!            ! Compute all AOs in r2p 
!            call give_all_aos_at_r(r2p, aos_in_r2p)
!
!            ! Get the core-MOs at r2p reciclying the newly computed AOs in r2p
!            ! and using only the the sub-matrix of core-MOs coefficients
!            call dgemv( 'N', n_core_pseudo_orb, ao_num, 1.d0 &
!                      , mo_core_coef_notnorm_transp, n_core_pseudo_orb &
!                      , aos_in_r2p, 1, 0.d0, mos_core_in_r2p, 1)
!
!            ! Compute pair Jastrow factor between r1 and r2p
!            if (core_tcxc_j0_testing) then
!              j_r1r2p = 1.d0
!            else
!              j_r1r2p = exp(j_mu_env(r1,r2p,mu_erf))
!            end if
!
!            ! Initialize kernel
!            kernel = 0.d0
!            ! Loop over core orbitals to update the kernel
!            do m = 1, n_core_pseudo_orb
!              m_core = list_core_pseudo(m)
!              !! Update kernel using all-MOs(r2p) array
!              !kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r2p(m_core)    
!              ! Update kernel using core-MOs(r2p) array
!              kernel -= mos_in_r_array2_omp(m_core,i2) * mos_core_in_r2p(m)    
!            enddo
!            kernel = kernel / distance
!
!            do l = 1, ao_num
!              ao_l_r2 = aos_in_r_array2(l,i2)
!              do j = 1, ao_num
!                ao_j_r2p = aos_in_r2p(j)
!                ! Notice the index order, because this array will be transposed later!
!                ao_core_xc_mat_grid12a(j,l,i1) += w2 * j_r1r2 * ao_l_r2 &
!                                                * w2p * j_r1r2p * kernel * ao_j_r2p
!              enddo
!            enddo
!
!! OMP ALTERNATIVE: 
!!!$OMP PARALLEL DO
!!!$OMP PARALLEL DO PRIVATE(ao_l_r2, ao_j_r2p) SHARED(aos_in_r_array2, aos_in_r2p, ao_core_xc_mat_grid12a, w2, j_r1r2, w2p, kernel, j_r1r2p, ao_num)
!!do l = 1, ao_num
!!    ao_l_r2 = aos_in_r_array2(l, i2)
!!    do j = 1, ao_num
!!        ao_j_r2p = aos_in_r2p(j)
!!        ! Notice the index order, because this array will be transposed later!
!!        ao_core_xc_mat_grid12a(j, l, i1) = ao_core_xc_mat_grid12a(j, l, i1) + w2 * j_r1r2 * ao_l_r2 * w2p * kernel * ao_j_r2p * j_r1r2p
!!    enddo
!!enddo
!!!$OMP END PARALLEL DO
!
!          endif
!
!        enddo !r2p floating
!      enddo   !r2p floating
!
!      ! FIXED GRID3 contributions 
!      ! (this grid3 is based on a provider called grid2 because we are recycling)
!      do i2p_nuc = 1, nucl_num
!        do i2p_rad = 1, n_points_extra_radial_grid - 1
!          do i2p_ang = 1, n_points_extra_integration_angular
!            r2p(1:3) = grid_points_extra_per_atom(1:3,i2p_ang,i2p_rad,i2p_nuc)
!            distance = norm2(r2(1:3) - r2p(1:3))
!            ! Compute matrix element only when there is no divergence 
!            if (distance.gt.1.d-10) then
!              ! NOTICE: THE WEIGHTS ON THE EXTRA GRID ARE UPDATED FOR THE ADAPTIVE GRID 
!              w2p = grid_fixed_weights(i2p_ang,i2p_rad,i2p_nuc)
!
!              ! Compute pair Jastrow factor between r1 and r2p
!              if (core_tcxc_j0_testing) then
!                j_r1r2p = 1.d0
!              else
!                j_r1r2p = exp(j_mu_env(r1,r2p,mu_erf))
!              end if
!
!              ! Initialize kernel
!              kernel = 0.d0
!              ! Loop over core orbitals to update the kernel
!              do m = 1, n_core_pseudo_orb
!                m_core = list_core_pseudo(m)
!                kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r_array_extra_full_omp(m_core,i2p_ang,i2p_rad,i2p_nuc)
!              enddo
!              kernel = kernel/distance
!
!              do l = 1, ao_num
!                ao_l_r2 = aos_in_r_array2(l,i2)
!                do j = 1, ao_num
!                  ao_j_r2p = aos_in_r_array_extra_full(j,i2p_ang,i2p_rad,i2p_nuc)
!                  ! Notice the index order, because this array will be transposed later!
!                  ao_core_xc_mat_grid12a(j,l,i1) += w2 * j_r1r2 * ao_l_r2 &
!                                                  * w2p * kernel * ao_j_r2p * j_r1r2p
!                enddo  ! j
!              enddo ! l
!
!            endif
!
!          enddo !r2p fixed
!        enddo   !r2p fixed
!      enddo     !r2p fixed
!
!    enddo !r2
!
!  end do !r1
!
!END_PROVIDER



! HERE WE TRY TO PRODUCE TWO HEAVY PROVIDER AT ONCE SO THAT THEY SHARE A
! LARGE PART OF THE CALCULATIONS
 BEGIN_PROVIDER [ double precision, core_tcxc_adapt_grid12aj_loop, (ao_num, ao_num, ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_core_xc_mat_grid12a, (ao_num, ao_num, n_points_final_grid) ]
  integer :: i1, i2, i2p
    ! Grids indices
  integer :: i, j, k, l
    ! AOs indices
  double precision :: ao_i_r1, ao_k_r1, ao_l_r2, ao_j_r2p
    ! Real values of AOs and core-exchange kernel
  double precision :: w1, r1(3)
    ! Grid 1 point
  double precision :: w2, r2(3)
    ! Grid 2 point and weight
  double precision :: r2p(3), w2p
    ! Nested grid 2' point (r'_2) and weight (dr'_2)
  double precision :: j_r1r2, j_r1r2p
    ! Jastrow factor between r1 and r2
  double precision, external :: j_mu_env
    ! Temporary external definition, we should put them in a module an import that
  double precision :: distance
    ! Distance for Coulomb integral
  integer :: m, m_core
    ! Indices to loop over core orbitals
  double precision :: kernel
    ! Kernel of the exchange potential

  ! ADAPTIVE GRID VARIABLES
  double precision :: mo_orbitals_in_r2p(mo_num)
    ! MO orbitals on the floating grid
  double precision :: ao_orbitals_in_r2p(ao_num)
    ! AO orbitals on the floating grid
  double precision :: grid_float_points(3, n_points_ang_float_grid, n_points_rad_float_grid, 1)
    ! Floating grid points
  double precision :: grid_float_weights(n_points_ang_float_grid, n_points_rad_float_grid, 1)
    ! Floating grid weights
  !double precision :: grid_fixed_weights(grid3_ang_size, grid3_rad_size, nucl_num)
  ! NOW THE EXTRA GRID IS USED AS THIRD GRID
  double precision :: grid_fixed_weights(n_points_extra_integration_angular,n_points_extra_radial_grid,nucl_num)
    ! Fixed grid weights
  integer :: i2p_nuc, i2p_rad, i2p_ang
    ! Third grid indices
  double precision :: aos_in_r2p(ao_num), mos_in_r2p(mo_num), mos_core_in_r2p(n_core_pseudo_orb)
    ! AOs and MOs computed at points belonging to the floating part of grid3
    ! Depending on how fast we want to go, we either compute all MOs or just the core ones
    
  ! Fixed extra grid weights
  ! To be used for pruning, NOT IMPLEMENTED YET
  integer :: n_fixed_pts_effective(nucl_num)
  integer :: n_float_pts_effective
  integer :: n_pts_effective_max


  ! Intermediate variables
  double precision :: v_xc_core_ao_j_r2p
  double precision :: ao_core_xc_r2

  ! Initialize floating grid points and weights
  grid_float_points(:,:,:,:) = 0.d0
  grid_float_weights(:,:,:) = 0.d0
  grid_fixed_weights(:,:,:) = 0.d0

  ao_core_xc_mat_grid12a(:,:,:) = 0.d0

  do i1 = 1, n_points_final_grid
    r1(1:3) = final_grid_points(1:3,i1)
    w1 = final_weight_at_r_vector(i1) ! this is only necessary for the testing provider
    do i2 = 1, n_points_final_grid2
      r2(1:3) = final_grid_points2(1:3,i2)
      w2 = final_weight_at_r_vector2(i2)
      ! Compute pair Jastrow factor between r1 and r2
      if (core_tcxc_j0_testing) then
        j_r1r2 = 1.d0
      else
        j_r1r2 = exp(-j_mu_env(r1,r2,mu_erf))
      end if

      ! Compute adaptive grid centered in r2
      call get_adaptive_grid( r2, grid_points_extra_per_atom, grid_float_points &
                            , grid_fixed_weights, grid_float_weights            &
                            , n_fixed_pts_effective, n_float_pts_effective      &
                            , n_pts_effective_max)

      ! FLOATING grid contributions
      do i2p_rad = 1, n_points_rad_float_grid - 1
        do i2p_ang = 1, n_points_ang_float_grid
          ! Find r2p point from i2p_rad, i2p_ang loop-indices
          r2p(1:3) = grid_float_points(1:3,i2p_ang,i2p_rad,1)
          distance = norm2(r2(1:3) - r2p(1:3))
          ! Compute matrix element only when there is no divergence 
          if (distance.gt.1.d-10) then
            ! Find r'_2 weight from i2p loop-index
            w2p = grid_float_weights(i2p_ang,i2p_rad,1)
            !write(*,*) "i2p_rad = ", i2p_rad, "; i2p_ang = ", i2p_ang, "; w2p = ", w2p
            ! Compute all AOs in r2p 
            call give_all_aos_at_r(r2p, aos_in_r2p)

            ! Get the core-MOs at r2p reciclying the newly computed AOs in r2p
            ! and using only the the sub-matrix of core-MOs coefficients
            call dgemv( 'N', n_core_pseudo_orb, ao_num, 1.d0 &
                      , mo_core_coef_notnorm_transp, n_core_pseudo_orb &
                      , aos_in_r2p, 1, 0.d0, mos_core_in_r2p, 1)

            ! Compute pair Jastrow factor between r1 and r2p
            if (core_tcxc_j0_testing) then
              j_r1r2p = 1.d0
            else
              j_r1r2p = exp(j_mu_env(r1,r2p,mu_erf))
            end if

            ! Initialize kernel
            kernel = 0.d0
            ! Loop over core orbitals to update the kernel
            do m = 1, n_core_pseudo_orb
              m_core = list_core_pseudo(m)
              !! Update kernel using all-MOs(r2p) array
              !kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r2p(m_core)    
              ! Update kernel using core-MOs(r2p) array
              kernel -= mos_in_r_array2_omp(m_core,i2) * mos_core_in_r2p(m)    
            enddo
            kernel = kernel / distance

            ! Notice the order of the indices (j,l) in ao_core_xc_mat_grid12a: 
            ! this array will be transposed later in the dgemm

            ! Option 1: OMP parallelism
            if (core_tcxc_omp) then
              ! Since ao_core_xc_mat_grid12a is updated at the second nested loop,
              ! only the first two loops can be OMP parallelised

              !$OMP PARALLEL DEFAULT(NONE) &
              !$OMP& PRIVATE (i, j, k, l, ao_core_xc_r2, v_xc_core_ao_j_r2p) &
              !$OMP& SHARED ( ao_core_xc_mat_grid12a                         &
              !$OMP&        , core_tcxc_adapt_grid12aj_loop                  &
              !$OMP&        , aos_in_r_array, aos_in_r_array2, aos_in_r2p    &
              !$OMP&        , ao_num, i1, i2, w1, w2, w2p, j_r1r2, j_r1r2p, kernel)
              !$OMP DO COLLAPSE(2) SCHEDULE(static)
              do j = 1, ao_num
                do l = 1, ao_num
                  v_xc_core_ao_j_r2p = w2p * j_r1r2p * kernel * aos_in_r2p(j)
                  ao_core_xc_r2 = w2 * j_r1r2 * aos_in_r_array2(l,i2) * v_xc_core_ao_j_r2p
                  ao_core_xc_mat_grid12a(j,l,i1) += ao_core_xc_r2

                  do k = 1, ao_num
                    do i = 1, ao_num
                      core_tcxc_adapt_grid12aj_loop(i,k,l,j) += &
                          w1 * aos_in_r_array(i,i1) * aos_in_r_array(k, i1) * ao_core_xc_r2
                    end do
                  end do
                end do
              end do
              !$OMP END DO
              !$OMP END PARALLEL

            ! Option 2: NO parallelism 
            else
              do l = 1, ao_num
                ao_l_r2 = aos_in_r_array2(l,i2)
                do j = 1, ao_num
                  ao_j_r2p = aos_in_r2p(j)
                  ao_core_xc_mat_grid12a(j,l,i1) += w2 * j_r1r2 * ao_l_r2 &
                                                  * w2p * j_r1r2p * kernel * ao_j_r2p
                  do k = 1, ao_num
                    do i = 1, ao_num
                      core_tcxc_adapt_grid12aj_loop(i,k,l,j) += &
                          &   w1  * aos_in_r_array(i,i1) * aos_in_r_array(k, i1) &
                          & * w2  * j_r1r2               * aos_in_r_array2(l, i2) &
                          & * w2p * j_r1r2p * kernel     * aos_in_r2p(j)
                    end do
                  end do

                enddo
              enddo

            end if

          endif

        enddo !r2p floating
      enddo   !r2p floating

      ! FIXED GRID3 contributions 
      ! (this grid3 is based on a provider called grid2 because we are recycling)
      do i2p_nuc = 1, nucl_num
        do i2p_rad = 1, n_points_extra_radial_grid - 1
          do i2p_ang = 1, n_points_extra_integration_angular
            r2p(1:3) = grid_points_extra_per_atom(1:3,i2p_ang,i2p_rad,i2p_nuc)
            distance = norm2(r2(1:3) - r2p(1:3))
            ! Compute matrix element only when there is no divergence 
            if (distance.gt.1.d-10) then
              ! NOTICE: THE WEIGHTS ON THE EXTRA GRID ARE UPDATED FOR THE ADAPTIVE GRID 
              w2p = grid_fixed_weights(i2p_ang,i2p_rad,i2p_nuc)

              ! Compute pair Jastrow factor between r1 and r2p
              if (core_tcxc_j0_testing) then
                j_r1r2p = 1.d0
              else
                j_r1r2p = exp(j_mu_env(r1,r2p,mu_erf))
              end if

              ! Initialize kernel
              kernel = 0.d0
              ! Loop over core orbitals to update the kernel
              do m = 1, n_core_pseudo_orb
                m_core = list_core_pseudo(m)
                kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r_array_extra_full_omp(m_core,i2p_ang,i2p_rad,i2p_nuc)
              enddo
              kernel = kernel/distance

              ! Notice the order of the indices (j,l): this array will be transposed later!
    
              ! Option 1: OMP parallelism
              if (core_tcxc_omp) then
                ! Since ao_core_xc_mat_grid12a is updated at the second nested loop,
                ! only the first two loops can be OMP parallelised

                !$OMP PARALLEL DEFAULT(NONE) &
                !$OMP& PRIVATE (i, j, k, l, ao_core_xc_r2, v_xc_core_ao_j_r2p) &
                !$OMP& SHARED ( ao_core_xc_mat_grid12a                         &
                !$OMP&        , core_tcxc_adapt_grid12aj_loop                  &
                !$OMP&        , aos_in_r_array, aos_in_r_array2                &
                !$OMP&        , aos_in_r_array_extra_full                      &
                !$OMP&        , ao_num, i1, i2, w1, w2, w2p                    &
                !$OMP&        , i2p_ang,i2p_rad,i2p_nuc, j_r1r2, j_r1r2p, kernel)
                !$OMP DO COLLAPSE(2) SCHEDULE(static)
                do j = 1, ao_num
                  do l = 1, ao_num
                    v_xc_core_ao_j_r2p = w2p * j_r1r2p * kernel * aos_in_r_array_extra_full(j,i2p_ang,i2p_rad,i2p_nuc)
                    ao_core_xc_r2 = w2 * j_r1r2 * aos_in_r_array2(l,i2) * v_xc_core_ao_j_r2p
                    ao_core_xc_mat_grid12a(j,l,i1) += ao_core_xc_r2

                    do k = 1, ao_num
                      do i = 1, ao_num
                        core_tcxc_adapt_grid12aj_loop(i,k,l,j) += &
                            w1 * aos_in_r_array(i,i1) * aos_in_r_array(k,i1) * ao_core_xc_r2
                      end do
                    end do
                  end do
                end do
                !$OMP END DO
                !$OMP END PARALLEL

              ! Option 2: NO parallelism 
              else

                do l = 1, ao_num
                  ao_l_r2 = aos_in_r_array2(l,i2)
                  do j = 1, ao_num
                    ao_j_r2p = aos_in_r_array_extra_full(j,i2p_ang,i2p_rad,i2p_nuc)
                    ao_core_xc_mat_grid12a(j,l,i1) += w2 * j_r1r2 * ao_l_r2 &
                                                    * w2p * j_r1r2p * kernel * ao_j_r2p
                    do k = 1, ao_num
                      do i = 1, ao_num
                        core_tcxc_adapt_grid12aj_loop(i,k,l,j) += &
                            &   w1  * aos_in_r_array(i,i1) * aos_in_r_array(k, i1) &
                            & * w2  * j_r1r2               * aos_in_r_array2(l, i2) &
                            & * w2p * j_r1r2p * kernel     * aos_in_r_array_extra_full(j,i2p_ang,i2p_rad,i2p_nuc)
                      end do
                    end do
  
                  enddo
                enddo

              end if

            endif

          enddo !r2p fixed
        enddo   !r2p fixed
      enddo     !r2p fixed

    enddo !r2

  end do !r1

END_PROVIDER
