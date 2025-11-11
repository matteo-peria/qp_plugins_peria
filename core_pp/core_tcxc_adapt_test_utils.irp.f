BEGIN_PROVIDER [ double precision, core_xpot_adapt_grid23, (ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! Numerical evaluation of
  ! < l | V_x^{core_pseudo} | j >
  ! = -\sum_{m\in\text{core}} <lm|1/r12|mk>
  ! = \int dr_1 \int dr' \phi_l(r_1) V_x(r_1,r') \phi_j(r')
  ! = - \int dr_1 \int dr' \phi_l(r_1) \sum_{m\in\text{core}} (\phi_m(r_1)\phi_j(r')/|r-r'|) \phi_j(r')
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

  core_xpot_adapt_grid23(:,:) = 0.d0

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
        ! Find r'_2 point from i2p_rad, i2p_ang loop-indices
        r2p(1:3) = grid_float_points(1:3,i2p_ang,i2p_rad,1)
        distance = norm2(r2(1:3) - r2p(1:3))
        ! Compute matrix element only when there is no divergence 
        if (distance.gt.1.d-10) then
          ! Find r'_2 weight from i2p loop-index
          w2p = grid_float_weights(i2p_ang,i2p_rad,1)
          ! Compute all AOs in r2p 
          call give_all_aos_at_r(r2p, aos_in_r2p)

          ! For MOs, 3 alternatives:
          ! 1) Compute all-MOs at the new floating grid point r2p
          !    This is the LEAST CONVENIENT.
          !    Cons:
          !    - matrix product involves bigger matrices (all MOs)
          !    - internally it calls give_all_aos_at_r, which was just called 
          !      in the line above
          !    Pro:
          !    + recycle tested software
          call give_all_mos_at_r(r2p, mos_in_r2p)

          ! 2) Compute all-MOs starting from already computed AOs
          !    Cons:
          !    - no tested. Are you sure of what you are doing?
          !    - matrix product involves bigger matrices (all MOs)
          !    Pro:
          !    + recycle aos_in_r2p, just called above
          !call dgemv( 'T'
          !          , mo_num, ao_num, 1.d0, mo_coef
          !          , mo_num, aos_in_r2p, 1
          !          , 0.d0, mos_in_r2p, 1)

          ! 3) Compute core-only-MOs starting from already computed AOs
          !    This should be the MOST CONVENIENT.
          !    Cons:
          !    - no tested. Are you sure of what you are doing?
          !    Pro:
          !    + matrix product involves smaller matrices (just core-MOs)
          !    + recycle aos_in_r2p, just called above
          !call dgemv( 'T'                                             &
          !          , n_core_pseudo_orb, ao_num, 1.d0, mo_core_coef,  &
          !          , n_core_pseudo_orb, aos_in_r2p, 1                &
          !          , 0.d0, mos_core_in_r2p, 1)

          ! Initialize Jastrow-weighted nested core-valence exchange integral
          integral = 0.d0
          ! Loop over all AO (j-th index, variable is r2p)
          do j = 1, ao_num
            ! Get value of the j-th orbital at r2p
            ao_j_r2p = aos_in_r2p(j)

            ! Initialize kernel
            kernel = 0.d0
            ! Loop over core orbitals to update the kernel
            do m = 1, n_core_pseudo_orb
              m_core = list_core_pseudo(m)
              ! Update kernel using all-MOs(r2p) array
              kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r2p(m_core)    
              ! Update kernel using core-MOs(r2p) array
              !kernel -= mos_in_r_array2_omp(m_core,i2) * mos_core_in_r2p(m)    
            enddo
            kernel = kernel/distance
            ! Update of the integral at each new point r2p
            integral += kernel * ao_j_r2p * w2p
            ! Loop over all AO (l-th index, variable is r2)
            do l = 1, ao_num
              ao_l_r2 = aos_in_r_array2(l,i2)
              core_xpot_adapt_grid23(l,j) += w2 * ao_l_r2 * integral
            enddo

          enddo
        endif

      enddo
    enddo

    ! FIXED GRID3 contributions 
    ! (BASED ON A FULL GRID EXTRA)
    do i2p_nuc = 1, nucl_num
      do i2p_rad = 1, n_points_extra_radial_grid - 1
        do i2p_ang = 1, n_points_extra_integration_angular
          r2p(1:3) = grid_points_extra_per_atom(1:3,i2p_ang,i2p_rad,i2p_nuc)
          distance = norm2(r2(1:3) - r2p(1:3))
          ! Compute matrix element only when there is no divergence 
          if (distance.gt.1.d-10) then
            ! NOTICE: THE WEIGHTS ON THE EXTRA GRID ARE UPDATED FOR THE ADAPTIVE GRID 
            w2p = grid_fixed_weights(i2p_ang,i2p_rad,i2p_nuc)

            ! Initialise 3rd nested integral
            integral = 0.d0
            ! Loop over all AO (j-th index, variable is r2p)
            do j = 1, ao_num
              ! Get value of the j-th orbital at r2p
              ao_j_r2p = aos_in_r_array_extra_full(j,i2p_ang,i2p_rad,i2p_nuc)

              ! Initialize kernel
              kernel = 0.d0
              ! Loop over core orbitals to update the kernel
              do m = 1, n_core_pseudo_orb
                m_core = list_core_pseudo(m)
                kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r_array_extra_full_omp(m_core,i2p_ang,i2p_rad,i2p_nuc)
              enddo

              ! Update of the integral at each new point r2p
              integral += kernel * ao_j_r2p * w2p / distance 
              do l = 1, ao_num
                ao_l_r2 = aos_in_r_array2(l,i2)
                core_xpot_adapt_grid23(l,j) += w2 * ao_l_r2 * integral
              enddo

            enddo
          endif

        enddo
      enddo
    enddo

  enddo 

END_PROVIDER




BEGIN_PROVIDER [ double precision, core_tcxc_adapt_j0_grid123, (ao_num, ao_num, ao_num, ao_num)]
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
  ! 2nd integral is computed along extra grid
  ! 3rd integral is computed along extra grid
  END_DOC
  integer :: i,j,k,l

  ! Initialization
  core_tcxc_adapt_j0_grid123(:,:,:,:) = 0.d0

  ! do-loop version of the tensor product
  do l = 1, ao_num
    do k = 1, ao_num
      do j = 1, ao_num
        do i = 1, ao_num
          core_tcxc_adapt_j0_grid123(i,j,k,l) = ao_overlap_grid1(i,k) * core_xpot_adapt_grid23(j,l)
        end do
      end do
    end do
  end do
  !! vectorized version of the tensor product
  !core_tcxc_j0_grid123 = reshape( ao_overlap_grid1, shape=[ao_num,ao_num,1,1] ) &
  !                     * reshape( core_xpot_grid22, shape=[1,1,ao_num,ao_num] )
END_PROVIDER
