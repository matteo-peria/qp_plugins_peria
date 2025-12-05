BEGIN_PROVIDER [ double precision, core_xpot_adapt_grid2aj, (ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! Numerical evaluation of
  ! < l | V_x^{core_pseudo} | j >
  ! = -\sum_{m\in\text{core}} < l | \phi_m(1)\phi_m(2)/r12 | k >
  ! = \int dr_1 \int dr2 \chi_l(r_1) V_x(r_1,r_2) \chi_j(r_2)
  ! = - \int dr_1 \int dr_2 \chi_l(r_1) \sum_{m\in\text{core}} (\phi_m(r_1)\phi_j(r_2)/|r-r_2|) \chi_j(r_2)
  !
  ! 1st integral is over grid2
  ! 2nd integral is over adaptive grid (floating + extra)
  END_DOC
  integer :: i2
    ! Grids 2 index
  double precision :: w2, r2(3)
    ! Grid 2 point and weight
  double precision :: w2p, r2p(3)
    ! Grid 3 point and weight
  integer :: j, l
    ! AOs indices
  double precision :: ao_l_r2, ao_j_r2p
    ! Real values of AOs
  integer :: m, m_core
    ! Core-MOs indices
  double precision :: integral, kernel
    ! Temporary variables to stock mid-computation result
  double precision :: distance
    ! Distance for Coulomb integral (exchange)

  ! ADAPTIVE GRID VARIABLES
  double precision :: mo_orbitals_in_r2p(mo_num)
    ! MO orbitals on the floating grid
  double precision :: ao_orbitals_in_r2p(ao_num)
    ! AO orbitals on the floating grid
  double precision :: grid_float_points(3, n_points_ang_float_grid, n_points_rad_float_grid, 1)
    ! Floating grid points
  double precision :: grid_float_weights(n_points_ang_float_grid, n_points_rad_float_grid, 1)
    ! Floating grid weights
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


  ! Initialize floating grid points and weights
  grid_float_points(:,:,:,:) = 0.d0
  grid_float_weights(:,:,:) = 0.d0

  core_xpot_adapt_grid2aj(:,:) = 0.d0

  do i2 = 1, n_points_final_grid2
    r2(1:3) = final_grid_points2(1:3,i2)
    w2 = final_weight_at_r_vector2(i2)

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

          ! Compute all AOs in r2p 
          call give_all_aos_at_r(r2p, aos_in_r2p)
          ! Get the core-MOs at r2p reciclying the newly computed AOs in r2p
          ! and using only the the sub-matrix of core-MOs coefficients
          call dgemv( 'N', n_core_pseudo_orb, ao_num, 1.d0 &
                    , mo_core_coef_notnorm_transp, n_core_pseudo_orb &
                    , aos_in_r2p, 1, 0.d0, mos_core_in_r2p, 1)
          ! Initialize kernel
          kernel = 0.d0
          do m = 1, n_core_pseudo_orb
            ! Update kernel using core-MOs(r2p) array
            m_core = list_core_pseudo(m)
            kernel -= mos_in_r_array2_omp(m_core,i2) * mos_core_in_r2p(m)    
          enddo

          !! Play safe for the moment and compute all of them
          !call give_all_mos_at_r(r2p,mos_core_in_r2p)
          !! Initialize kernel
          !kernel = 0.d0
          !! Loop over core orbitals to update the kernel
          !do m = 1, n_core_pseudo_orb
          !  m_core = list_core_pseudo(m)
          !  ! Update kernel using all-MOs(r2p) array
          !  kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r2p(m_core)    
          !enddo

          kernel = kernel / distance

          ! Loop over all AO (j-th index, variable is r2p)
          do j = 1, ao_num
            ao_j_r2p = aos_in_r2p(j)
            ! Update of the integral at each new point r2p
            integral = w2p * kernel * ao_j_r2p

            do l = 1, ao_num
              ao_l_r2 = aos_in_r_array2(l,i2)
              core_xpot_adapt_grid2aj(l,j) += w2 * ao_l_r2 * integral
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
              ! Update of the integral at each new point r2p
              integral = w2p * kernel * ao_j_r2p

              do l = 1, ao_num
                ao_l_r2 = aos_in_r_array2(l,i2)
                core_xpot_adapt_grid2aj(l,j) += w2 * ao_l_r2 * integral
              enddo
            enddo

          end if

        end do  ! grid 3 fixed
      end do    ! grid 3 fixed
    end do      ! grid 3 fixed

  enddo ! grid 2

END_PROVIDER


BEGIN_PROVIDER [ double precision, core_tcxc_adapt_j0_grid12aj, (ao_num, ao_num, ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! TEST-ONLY PROVIDER.
  ! Check that the TC exchange integral make sense when the Jastrow factor 
  ! is set to 1 (that is, when the Jastrow exponent is 0).

  ! $$ \braket{kl | e^{-J} V_{x,\text{core}} e^{+J} | ij} |_{J=0}
  !      \int d1 \chi_k(1)^* \chi_i(1) 
  !        \int d2 \chi_l^*(2) \int d2' v_x^{\text{core}}(2,2') \chi_j(2')
  !          \braket{k|i} \braket{l| V_{x,\text{core}} |j} $$
  !
  ! 1st integral is computed along usual grid
  ! 2nd integral is computed along grid2
  ! 3rd integral is computed along adaptive grid (floating + extra grid)
  END_DOC
  integer :: i,j,k,l

  ! Initialization
  core_tcxc_adapt_j0_grid12aj(:,:,:,:) = 0.d0

  ! do-loop version of the tensor product
  do j = 1, ao_num
    do l = 1, ao_num
      do k = 1, ao_num
        do i = 1, ao_num
          core_tcxc_adapt_j0_grid12aj(i,k,l,j) = ao_overlap_grid1(i,k) * core_xpot_adapt_grid2aj(j,l)
        end do
      end do
    end do
  end do

  !! vectorized version of the tensor product
  !core_tcxc_j0_grid123 = reshape( ao_overlap_grid1, shape=[ao_num,ao_num,1,1] ) &
  !                     * reshape( core_xpot_grid22, shape=[1,1,ao_num,ao_num] )
END_PROVIDER



!! BEGIN_PROVIDER [ double precision, core_tcxc_adapt_grid12aj_testing, (ao_num, ao_num, ao_num, ao_num)]
!!   implicit none
!!   BEGIN_DOC
!!   ! Numerical evaluation of < kl | V^{\text{TC}}_{x,\text{core}}  | ij >
!!   !
!!   ! $$\braket{kl | V^{\text{TC}}_{x,\text{core}} | ij}
!!   !     \braket{kl | e^{-J} V_{x,\text{core}} e^{+J} | ij}
!!   !       \int d1 \chi_k(1)^* \chi_i(1) 
!!   !         \int d2 \chi_l^*(2) e^{-u(1,2)} 
!!   !           \int d2' v_x^{\text{core}}(2,2') e^{u(1,2')} \chi_j(2') $$
!!   END_DOC
!!   !
!!   integer :: i1, i2
!!     ! Grids indices
!!   double precision :: w1, r1(3)
!!     ! Grid 1 point and weight
!!   double precision :: w2, r2(3)
!!     ! Grid 2 point and weight
!!   double precision :: w2p, r2p(3)
!!     ! Grid 3 point and weight
!!   integer :: i, j, k, l
!!     ! AOs indices
!!   double precision :: ao_i_r1, ao_k_r1, ao_l_r2, ao_j_r2p
!!     ! Real values of AOs
!!   integer :: m, m_core
!!     ! Core-MOs indices
!!   double precision :: j_factor, j_r1r2, j_r1r2p
!!     ! Jastrow factors between r1 and r2 and between r1 and r2p
!!   double precision, external :: j_mu_env!, core_tcxc_kernel_grid123
!!     ! Temporary external definition, we should put them in a module an import that
!!   double precision :: integral, kernel
!!     ! Temporary variables to stock mid-computation result
!!   double precision :: distance
!!     ! Distance for Coulomb integral (exchange)
!!   integer :: ao_num2
!! 
!!   ! ADAPTIVE GRID VARIABLES
!!   double precision :: mo_orbitals_in_r2p(mo_num)
!!     ! MO orbitals on the floating grid
!!   double precision :: ao_orbitals_in_r2p(ao_num)
!!     ! AO orbitals on the floating grid
!!   double precision :: grid_float_points(3, n_points_ang_float_grid, n_points_rad_float_grid, 1)
!!     ! Floating grid points
!!   double precision :: grid_float_weights(n_points_ang_float_grid, n_points_rad_float_grid, 1)
!!     ! Floating grid weights
!!   !double precision :: grid_fixed_weights(grid3_ang_size, grid3_rad_size, nucl_num)
!!   ! NOW THE EXTRA GRID IS USED AS THIRD GRID
!!   double precision :: grid_fixed_weights(n_points_extra_integration_angular,n_points_extra_radial_grid,nucl_num)
!!     ! Fixed grid weights
!!   integer :: i2p_nuc, i2p_rad, i2p_ang
!!     ! Floating grid indices
!!   integer :: i2p
!!     ! Grid 3 index
!!   double precision :: aos_in_r2p(ao_num), mos_in_r2p(mo_num), mos_core_in_r2p(n_core_pseudo_orb)
!!     ! AOs and MOs computed at points belonging to the floating part of grid3
!!     ! Depending on how fast we want to go, we either compute all MOs or just the core ones
!!     
!!   ! Fixed extra grid weights
!!   ! To be used for pruning, NOT IMPLEMENTED YET
!!   integer :: n_fixed_pts_effective(nucl_num)
!!   integer :: n_float_pts_effective
!!   integer :: n_pts_effective_max
!! 
!!   ! Initialize floating grid points and weights
!!   grid_float_points(:,:,:,:) = 0.d0
!!   grid_float_weights(:,:,:) = 0.d0
!!   grid_fixed_weights(:,:,:) = 0.d0
!! 
!!   core_tcxc_adapt_grid12aj_testing(:,:,:,:) = 0.d0
!! 
!!   do i1 = 1, n_points_final_grid
!!     r1(1:3) = final_grid_points(1:3,i1)
!!     w1 = final_weight_at_r_vector(i1)
!!     !write(*,*) "i1 = ", i1, "; w1 = ", w1
!!     do i2 = 1, n_points_final_grid2
!!       r2(1:3) = final_grid_points2(1:3,i2)
!!       w2 = final_weight_at_r_vector2(i2)
!!       !write(*,*) "i2 = ", i2, "; w2 = ", w2
!!       ! Compute pair Jastrow factor between r1 and r2
!!       if (core_tcxc_j0_testing) then
!!         j_r1r2 = 1.d0
!!       else
!!         j_r1r2 = exp(-j_mu_env(r1,r2,mu_erf))
!!       end if
!! 
!!       ! Compute adaptive grid centered in r2
!!       call get_adaptive_grid( r2, grid_points_extra_per_atom, grid_float_points &
!!                             , grid_fixed_weights, grid_float_weights            &
!!                             , n_fixed_pts_effective, n_float_pts_effective      &
!!                             , n_pts_effective_max)
!! 
!!       ! FLOATING grid contributions
!!       do i2p_rad = 1, n_points_rad_float_grid - 1
!!         do i2p_ang = 1, n_points_ang_float_grid
!!           ! Find r2p point from i2p_rad, i2p_ang loop-indices
!!           r2p(1:3) = grid_float_points(1:3,i2p_ang,i2p_rad,1)
!!           distance = norm2(r2(1:3) - r2p(1:3))
!!           ! Compute matrix element only when there is no divergence 
!!           if (distance.gt.1.d-10) then
!!             ! Find r'_2 weight from i2p loop-index
!!             w2p = grid_float_weights(i2p_ang,i2p_rad,1)
!!             !write(*,*) "i2p_rad = ", i2p_rad, "; i2p_ang = ", i2p_ang, "; w2p = ", w2p
!!             ! Compute all AOs in r2p 
!!             call give_all_aos_at_r(r2p, aos_in_r2p)
!! 
!!             ! Get the core-MOs at r2p reciclying the newly computed AOs in r2p
!!             ! and using only the the sub-matrix of core-MOs coefficients
!!             call dgemv( 'N', n_core_pseudo_orb, ao_num, 1.d0 &
!!                       , mo_core_coef_notnorm_transp, n_core_pseudo_orb &
!!                       , aos_in_r2p, 1, 0.d0, mos_core_in_r2p, 1)
!! 
!!             ! Compute pair Jastrow factor between r1 and r2p
!!             if (core_tcxc_j0_testing) then
!!               j_r1r2p = 1.d0
!!             else
!!               j_r1r2p = exp(j_mu_env(r1,r2p,mu_erf))
!!             end if
!! 
!!             ! Initialize kernel
!!             kernel = 0.d0
!!             ! Loop over core orbitals to update the kernel
!!             do m = 1, n_core_pseudo_orb
!!               m_core = list_core_pseudo(m)
!!               !! Update kernel using all-MOs(r2p) array
!!               !kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r2p(m_core)    
!!               ! Update kernel using core-MOs(r2p) array
!!               kernel -= mos_in_r_array2_omp(m_core,i2) * mos_core_in_r2p(m)    
!!             enddo
!!             kernel = kernel / distance
!! 
!! !            ! Loop over all AO (j-th index, variable is r2p)
!! !            do j = 1, ao_num
!! !              ao_j_r2p = aos_in_r2p(j)
!! !              do l = 1, ao_num
!! !                ao_l_r2 = aos_in_r_array2(l,i2)
!! !                do k = 1, ao_num
!! !                  ao_k_r1 = aos_in_r_array(k,i1)
!! !                  do i = 1, ao_num
!! !                    ao_i_r1 = aos_in_r_array(i,i1)
!! !      
!! !                    core_tcxc_adapt_grid12aj_testing(i,k,l,j) += w1 * ao_i_r1 * ao_k_r1            &
!! !                                                       * w2 * j_r1r2 * ao_l_r2             &
!! !                                                       * w2p * kernel * ao_j_r2p * j_r1r2p
!! !                  enddo
!! !                enddo
!! !              enddo
!! !            enddo
!! 
!! !$OMP PARALLEL DEFAULT(NONE) &
!! !$OMP& PRIVATE(i,j,k,l) &
!! !$OMP& SHARED(ao_num, core_tcxc_adapt_grid12aj_testing,            &
!! !$OMP&        aos_in_r_array, aos_in_r_array2, aos_in_r2p,         &
!! !$OMP&        w1, w2, w2p, j_r1r2, j_r1r2p, kernel, i1, i2)
!! !$OMP DO COLLAPSE(4) SCHEDULE(static)
!! do j = 1, ao_num
!!   do l = 1, ao_num
!!     do k = 1, ao_num
!!       do i = 1, ao_num
!!         core_tcxc_adapt_grid12aj_testing(i,k,l,j) += & !core_tcxc_adapt_grid12aj_testing(i,k,l,j)  &
!!           &   + w1  * aos_in_r_array(i,i1) * aos_in_r_array(k, i1) &
!!           &   * w2  * j_r1r2               * aos_in_r_array2(l, i2) &
!!           &   * w2p * j_r1r2p * kernel     * aos_in_r2p(j)
!!       end do
!!     end do
!!   end do
!! end do
!! !$OMP END DO
!! !$OMP END PARALLEL
!! 
!!           end if
!! 
!!         end do  ! grid 3 floating
!!       end do    ! grid 3 floating
!! 
!!       ! FIXED GRID3 contributions 
!!       ! (this grid3 is based on a provider called grid2 because we are recycling)
!!       do i2p_nuc = 1, nucl_num
!!         do i2p_rad = 1, n_points_extra_radial_grid - 1
!!           do i2p_ang = 1, n_points_extra_integration_angular
!!             r2p(1:3) = grid_points_extra_per_atom(1:3,i2p_ang,i2p_rad,i2p_nuc)
!!             distance = norm2(r2(1:3) - r2p(1:3))
!!             ! Compute matrix element only when there is no divergence 
!!             if (distance.gt.1.d-10) then
!!               ! NOTICE: THE WEIGHTS ON THE EXTRA GRID ARE UPDATED FOR THE ADAPTIVE GRID 
!!               w2p = grid_fixed_weights(i2p_ang,i2p_rad,i2p_nuc)
!! 
!!               ! Compute pair Jastrow factor between r1 and r2p
!!               if (core_tcxc_j0_testing) then
!!                 j_r1r2p = 1.d0
!!               else
!!                 j_r1r2p = exp(j_mu_env(r1,r2p,mu_erf))
!!               end if
!! 
!!               ! Initialize kernel
!!               kernel = 0.d0
!!               ! Loop over core orbitals to update the kernel
!!               do m = 1, n_core_pseudo_orb
!!                 m_core = list_core_pseudo(m)
!!                 kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r_array_extra_full_omp(m_core,i2p_ang,i2p_rad,i2p_nuc)
!!               enddo
!!               kernel = kernel/distance
!! 
!! !              ! Loop over all AO (j-th index, variable is r2p)
!! !              do j = 1, ao_num
!! !                ao_j_r2p = aos_in_r_array_extra_full(j,i2p_ang,i2p_rad,i2p_nuc)
!! !                do l = 1, ao_num
!! !                  ao_l_r2 = aos_in_r_array2(l,i2)
!! !                  do k = 1, ao_num
!! !                    ao_k_r1 = aos_in_r_array(k,i1)
!! !                    do i = 1, ao_num
!! !                      ao_i_r1 = aos_in_r_array(i,i1)
!! !        
!! !                      core_tcxc_adapt_grid12aj_testing(i,k,l,j) += w1 * ao_i_r1 * ao_k_r1            &
!! !                                                         * w2 * j_r1r2 * ao_l_r2             &
!! !                                                         * w2p * kernel * ao_j_r2p * j_r1r2p
!! !                    enddo   ! i
!! !                  enddo   ! k
!! !                enddo   ! l
!! !              enddo   ! j
!! 
!! !$OMP PARALLEL DEFAULT(NONE) &
!! !$OMP& PRIVATE(i,j,k,l) &
!! !$OMP& SHARED(ao_num, core_tcxc_adapt_grid12aj_testing,                   &
!! !$OMP&        aos_in_r_array, aos_in_r_array2, aos_in_r_array_extra_full, &
!! !$OMP&        w1, w2, w2p, j_r1r2, j_r1r2p, kernel, i1, i2,               &
!! !$OMP&        i2p_ang, i2p_rad, i2p_nuc)
!! !$OMP DO COLLAPSE(4) SCHEDULE(static)
!! do j = 1, ao_num
!!   do l = 1, ao_num
!!     do k = 1, ao_num
!!       do i = 1, ao_num
!!         core_tcxc_adapt_grid12aj_testing(i,k,l,j) = core_tcxc_adapt_grid12aj_testing(i,k,l,j)  &
!!      &     + w1  * aos_in_r_array(i,i1) * aos_in_r_array(k,  i1) &
!!      &     * w2  * j_r1r2               * aos_in_r_array2(l, i2) &
!!      &     * w2p * j_r1r2p * kernel     * aos_in_r_array_extra_full(j,i2p_ang,i2p_rad,i2p_nuc)
!!       end do
!!     end do
!!   end do
!! end do
!! !$OMP END DO
!! !$OMP END PARALLEL
!! 
!!             end if
!! 
!!           end do  ! grid 3 fixed
!!         end do    ! grid 3 fixed
!!       end do      ! grid 3 fixed
!! 
!!     enddo ! grid 2
!!   enddo ! grid 1
!! 
!! END_PROVIDER
