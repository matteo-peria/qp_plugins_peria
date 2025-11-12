BEGIN_PROVIDER [ double precision, core_tcxc_adapt_grid123, (ao_num, ao_num, ao_num, ao_num)]
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

  ! Initialize floating grid points and weights
  grid_float_points(:,:,:,:) = 0.d0
  grid_float_weights(:,:,:) = 0.d0

  core_tcxc_adapt_grid123(:,:,:,:) = 0.d0

  do i1 = 1, n_points_final_grid
    write(*,*) "Loop 1: i1 = ", i1
    r1(1:3) = final_grid_points(1:3,i1)
    w1 = final_weight_at_r_vector(i1)
! COMMENTED BECAUSE THE EXTRA GRID NOW MUST BE PART OF THE ADAPTIVE
!    do i2 = 1, n_points_extra_final_grid
!      r2(1:3) = final_grid_points_extra(1:3,i2)
!      w2 = final_weight_at_r_vector_extra(i2) 
    do i2 = 1, n_points_final_grid2
      write(*,*) "Loop 2: i2 = ", i2
      r2(1:3) = final_grid_points2(1:3,i2)
      w2 = final_weight_at_r_vector2(i2)

      ! Compute pair Jastrow factor between r1 and r2
      if (core_tcxc_j0_testing) then
        !write(*,*) "WARNING: for testing reason: j_mu_env = 0)"
        j_factor = 1.d0
        j_r1r2 = 1.d0
      else
        j_factor = exp(-j_mu_env(r1,r2,mu_erf))
        j_r1r2 = exp(-j_mu_env(r1,r2,mu_erf))
      end if

      ! Compute adaptive grid centered in r2
      call get_adaptive_grid( r2, grid_points_extra_per_atom, grid_float_points &
                            , grid_fixed_weights, grid_float_weights            &
                            , n_fixed_pts_effective, n_float_pts_effective      &
                            , n_pts_effective_max)

      !call get_adaptive_grid3( r2 &
      !                       , grid3_molec_r, grid_float_points, &
      !                       , grid_fixed_weights, grid_float_weights &
      !                       )

      ! FLOATING grid contributions
      do i2p_rad = 1, n_points_rad_float_grid - 1
        do i2p_ang = 1, n_points_ang_float_grid
          !write(*,*) "Loops 3 (float): i2p_rad = ", i2p_rad, "  i2p_ang = ", i2p_ang
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
              !write(*,*) "Loop 4.1: j-th AO = ", j
              ! Get value of the j-th orbital at r2p
              ao_j_r2p = aos_in_r2p(j)

              ! Compute pair Jastrow factor between r1 and r2p
              if (core_tcxc_j0_testing) then
                j_factor = 1.d0
                j_r1r2p = 1.d0
              else
                j_factor = exp(j_mu_env(r1,r2p,mu_erf))
                j_r1r2p = exp(j_mu_env(r1,r2p,mu_erf))
              end if

              ! Initialize kernel
              kernel = 0.d0
              ! Loop over core orbitals to update the kernel
              do m = 1, n_core_pseudo_orb
                !write(*,*) "Loop 5.a.1: m-th core-MO = ", m
                m_core = list_core_pseudo(m)
                ! Update kernel using all-MOs(r2p) array
                kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r2p(m_core)    
                ! Update kernel using core-MOs(r2p) array
                !kernel -= mos_in_r_array2_omp(m_core,i2) * mos_core_in_r2p(m)    
              enddo
              ! Update of the integral at each new point r2p
              integral += kernel * j_r1r2p * ao_j_r2p * w2p / distance 
              !
              do l = 1, ao_num
                !write(*,*) "Loop 5.b.1: l-th AO = ", l
                ao_l_r2 = aos_in_r_array2(l,i2)
                do k = 1, ao_num
                  !write(*,*) "Loop 6.1: k-th AO = ", k
                  ao_k_r1 = aos_in_r_array(k,i1)
                  do i = 1, ao_num
                    !write(*,*) "Loop 7.1: i-th AO = ", i
                    ao_i_r1 = aos_in_r_array(i,i1)
                    core_tcxc_adapt_grid123(i,k,l,j) += w1 * ao_i_r1 * ao_k_r1 * w2 * j_r1r2 * ao_l_r2 * integral
                  enddo
                enddo
              enddo
            enddo
          endif
        enddo
      enddo

      ! FIXED GRID3 contributions 
      ! (this grid3 is based on a provider called grid2 because we are recycling)

!      ! OLD OPTION BASED ON PROVIDER GRID2
!      do i2p_nuc = 1, nucl_num
!        do i2p_rad = 1, n_points_rad_grid2 - 1
!          do i2p_ang = 1, n_points_ang_grid2
!            r2p(1:3) = grid_points_per_atom2(1:3,i2p_ang,i2p_rad,i2p_nuc)
!            w2p = final_weight_at_r2(i2p_ang,i2p_rad,i2p_nuc)
!            distance = norm2(r2(1:3) - r2p(1:3))
!            ! Compute matrix element only when there is no divergence 
!            if (distance.gt.1.d-10) then
!
!              ! Compute all AOs in r2p 
!              ! EVEN IF IT'S A FIXED GRID, IT IS NOT PRUNED
!              ! YET MAKING IT A BIT CUMBERSOME TO HAVE A PROVIDER
!              call give_all_aos_at_r(r2p, aos_in_r2p)

! ANOTHER OLD OPTION 
!      do i2p_nuc = 1, nucl_num
!        do i2p_rad = 1, grid3_rad_size-1
!          do i2p_ang = 1, grid3_ang_size
!            ! Find r'_2 point from i2p_rad, i2p_ang loop-indices
!            r2p(1:3) = grid3_molec_r(1:3,i2p_ang,i2p_rad,i2p_nuc)
!            distance = norm2(r2(1:3) - r2p(1:3))
!            ! Compute matrix element only when there is no divergence 
!            if (distance.gt.1.d-10) then
!              ! Find r2p weight from i2p loop-indices
!              w2p = grid3_molec_w(i2p_ang,i2p_rad,i2p_nuc)
!
!              ! Compute all AOs in r2p 
!              ! EVEN IF IT'S A FIXED GRID, IT IS NOT PRUNED
!              ! YET MAKING IT A BIT CUMBERSOME TO HAVE A PROVIDER
!              call give_all_aos_at_r(r2p, aos_in_r2p)
!  
!              ! For MOs, 2 alternatives:
!              ! 1) Compute all-MOs at the new floating grid point r2p
!              !    This is LESS CONVENIENT because
!              !    - matrix product involves bigger matrices (all MOs)
!              !    - internally it calls give_all_aos_at_r, which was just called 
!              !      in the line above
!              !    + recycle tested software
!              !call give_all_mos_at_r(r2p, mos_in_r2p)
!              ! 2) Compute core-only-MOs starting from core-only-AOs
!              !    This is MORE CONVENIENT because
!              !    + matrix product involves smaller matrices (just core-MOs)
!              !    + recycle ao_orbitals_in_r2p, just called above
!              !    - no tested
!              ! Compute core-only MO
!              call dgemv( 'T'                                            &
!                        , n_core_pseudo_orb, ao_num, 1.d0, mo_core_coef, &
!                        , n_core_pseudo_orb, aos_in_r2p, 1               &
!                        , 0.d0, mos_core_in_r2p, 1)
!              !! Recompute all MOs recycling aos_in_r2p
!              !call dgemv( 'T'
!              !          , mo_num, ao_num, 1.d0, mo_coef
!              !          , mo_num, aos_in_r2p, 1
!              !          , 0.d0, mos_in_r2p, 1)
!
!              ! Initialise 3rd nested integral
!              integral = 0.d0
!              ! Loop over all AO (j-th index, variable is r2p)
!              do j = 1, ao_num
!                ! Get value of the j-th orbital at r2p
!                ao_j_r2p = aos_in_r2p(j)
!                !ao_j_r2p = aos_in_r_array_extra(j,i2p) !NOW WE DONT USE EXTRA
!                ! Compute pair Jastrow factor between r1 and r2p
!                if (core_tcxc_j0_testing) then
!                  j_factor = 1.d0
!                else
!                  j_factor = exp(j_mu_env(r1,r2p,mu_erf))
!                end if
!
!                ! Initialize kernel
!                kernel = 0.d0
!                ! Loop over core orbitals to update the kernel
!                do m = 1, n_core_pseudo_orb
!                  m_core = list_core_pseudo(m)
!                  !kernel -= mos_in_r_array_extra_omp(m_core,i2) * mos_in_r2p(m_core)
!                  kernel -= mos_in_r_array_extra_omp(m_core,i2) * mos_core_in_r2p(m)
!                enddo
!
!                ! Update of the integral at each new point r2p
!                integral += kernel * j_factor * ao_j_r2p * w2p / distance 
!                do l = 1, ao_num
!                  ao_l_r2 = aos_in_r_array_extra(l,i2)
!                  do k = 1, ao_num
!                    ao_k_r1 = aos_in_r_array(k,i1)
!                    do i = 1, ao_num
!                      ao_i_r1 = aos_in_r_array(i,i1)
!                      core_tcxc_adapt_grid123(i,k,l,j) += w1*ao_i_r1*ao_k_r1 * w2*ao_l_r2*integral
!                    enddo
!                  enddo
!                enddo
!              enddo
!            endif

      ! NEW OPTION BASED ON A FULL GRID EXTRA
      do i2p_nuc = 1, nucl_num
        do i2p_rad = 1, n_points_extra_radial_grid - 1
          do i2p_ang = 1, n_points_extra_integration_angular
            !write(*,*) "Loops 3 (extra): i2p_nuc = ", i2p_nuc, "  i2p_rad = ", i2p_rad, "  i2p_ang = ", i2p_ang
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
                !write(*,*) "Loop 4.2: j-th AO = ", j
                ! Get value of the j-th orbital at r2p
                ao_j_r2p = aos_in_r_array_extra_full(j,i2p_ang,i2p_rad,i2p_nuc)
                ! Compute pair Jastrow factor between r1 and r2p
                if (core_tcxc_j0_testing) then
                  j_factor = 1.d0
                  j_r1r2p = 1.d0
                else
                  j_factor = exp(j_mu_env(r1,r2p,mu_erf))
                  j_r1r2p = exp(j_mu_env(r1,r2p,mu_erf))
                end if

                ! Initialize kernel
                kernel = 0.d0
                ! Loop over core orbitals to update the kernel
                do m = 1, n_core_pseudo_orb
                  m_core = list_core_pseudo(m)
                  kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r_array_extra_full_omp(m_core,i2p_ang,i2p_rad,i2p_nuc)
                enddo

                ! Update of the integral at each new point r2p
                integral += kernel * j_r1r2p * ao_j_r2p * w2p / distance 
                do l = 1, ao_num
                  !write(*,*) "Loop 5.b.2: l-th AO = ", l
                  !ao_l_r2 = aos_in_r_array_extra(l,i2)
                  ao_l_r2 = aos_in_r_array2(l,i2)
                  do k = 1, ao_num
                    !write(*,*) "Loop 6.2: k-th AO = ", k
                    ao_k_r1 = aos_in_r_array(k,i1)
                    do i = 1, ao_num
                      !write(*,*) "Loop 7.2: i-th AO = ", i
                      ao_i_r1 = aos_in_r_array(i,i1)
                      core_tcxc_adapt_grid123(i,k,l,j) += w1 * ao_i_r1 * ao_k_r1 * w2 * j_r1r2 * ao_l_r2 * integral
                    enddo
                  enddo
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
    enddo 
  enddo

END_PROVIDER

! ALL OF THE FOLLOWING IS COMMENTED BECAUSE I'VE PUT EVERYTHING IN ONE
! SINGLE GIANT HEAVY PROVIDER. LOST ONLY SOME READABILITY, NOT EFFICIENCY
!! remember: pure functions are not allowed in irp.f files
!!pure function core_tcxc_kernel_grid23(r1,i2,j) result(res)
!function core_tcxc_kernel_grid123(r1,i2,j) result(res)
!  BEGIN_DOC 
!  ! Compute the integral
!  ! \int d2' v_x^{\text{core}}(2,2') e^{u(1,2')} \chi_j(2)
!  ! where 2' runs over the adaptive
!  END_DOC
!  implicit none
!  ! INPUT
!  double precision, intent(in) :: r1
!    ! Point of grid 1 at which we are computing the nested integral (parameter)
!  integer, intent(in) :: i2
!    ! Index of the point on grid 2
!  integer, intent(in) :: j
!    ! Index of the MO involved in the integral
!  ! OUTPUT
!  double precision :: res
!  ! INTERNALS
!  integer :: i2p 
!    ! Nested grid 2' point index 
!  double precision :: r2p(3), w2p
!    ! Nested grid 2' point (r'_2) and weight (dr'_2)
!  double precision :: r2(3)
!    ! Grid 2' point (r_2), index i2 retrieved from passed argument i2
!  double precision :: ao_j_r2p
!    ! Real value of MO $\chi_j(r'_2)$, index retrieved from passed argument j
!  double precision :: distance
!    ! Distance for Coulomb integral
!  integer :: m, m_core
!    ! Indices to loop over core orbitals
!  double precision :: kernel
!    ! Kernel of the exchange potential
!  double precision :: j_factor
!    ! Jastrow factor between r1 and r'2
!  double precision, external :: j_mu_env
!    ! 
!  ! ADAPTIVE GRID VARIABLES
!  double precision :: mo_orbitals_in_r2p(mo_num)
!    ! MO orbitals on the floating grid
!  double precision :: ao_orbitals_in_r2p(ao_num)
!    ! AO orbitals on the floating grid
!  double precision :: grid_float_points(3, n_points_ang_float_grid, n_points_rad_float_grid, 1)
!    ! Floating grid points
!  double precision :: grid_float_weights(n_points_ang_float_grid, n_points_rad_float_grid, 1)
!    ! Floating grid weights
!  double precision :: grid_fixed_weights(n_points_ang_extra_grid,n_points_rad_extra_grid,nucl_num)
!  integer :: i2p_nuc, i2p_rad, i2p_ang
!
!    ! Fixed extra grid weights
!  ! To be used for pruning, NOT IMPLEMENTED YET
!  integer :: n_fixed_pts_effective(nucl_num)
!  integer :: n_float_pts_effective
!  integer :: n_pts_effective_max
!
!
!  ! Initialize floating grid points and weights
!  grid_float_points(:,:,:,:) = 0.d0
!  grid_float_weights(:,:,:) = 0.d0
!
!  ! Retrieve r_2 from i2 index passed as argument:
!  r2(1:3) = final_grid_points_extra(1:3,i2)
!
!  ! Compute adaptive grid centered in r2
!  call get_adaptive_grid( r2, grid_points_extra_per_atom, grid_float_points, &
!                        , grid_fixed_weights, grid_float_weights,            &
!                        , n_fixed_pts_effective, n_float_pts_effective,      &
!                        , n_pts_effective_max)
!
!  ! Initialize Jastrow-weighted nested core-valence exchange integral
!  res = 0.d0
!
!  ! Loop over the total, finer grid is splitted into a 
!  ! floating contribution and a fixed one
!
!  ! FLOATING grid contributions
!  do i2p_rad = 1, n_points_rad_float_grid - 1
!    do i2p_ang = 1, n_points_ang_float_grid
!      ! Find r'_2 point from i2p_rad, i2p_ang loop-indices
!      r2p(1:3) = grid_float_points(1:3,i2p_ang,i2p_rad,1)
!      distance = norm2(r2(1:3) - r2p(1:3))
!      ! Compute matrix element only when there is no divergence 
!      if (distance.gt.1.d-10) then
!        ! Find r'_2 weight from i2p loop-index
!        w2p = grid_float_weights(i2p_ang,i2p_rad,1)
!        !! Need to compute all the MOs in the new floating point
!        !! DO WE REALLY NEED TO COMPUTE ALL THE MOs? SHOULDN'T WE CALL JUST THE CURRENT j PASSED AS ARGUMENT?
!
!        ! Need to compute core-MOs in the new floating point
!        call give_all_mos_at_r(r2p, mo_orbitals_in_r2p)
!
!        ! Get value of the j-th orbital at r'_2
!        ao_j_r2p = ao_orbitals_in_r2p(j)
!        ! Compute pair Jastrow factor between r1 and r'_2
!        j_factor = exp(j_mu_env(r1,r2p,mu_erf))
!        ! Initialize kernel
!        kernel = 0.d0
!        ! Loop over core orbitals to update the kernel
!        do m = 1, n_core_pseudo_orb
!          m_core = list_core_pseudo(m)
!          kernel -= mos_in_r_array_extra_omp(m_core,i2) * mo_orbitals_in_r2p(m_core)    
!        enddo
!        ! Update of the integral at each new point r2p
!        res += kernel * j_factor * ao_j_r2p * w2p / distance 
!      endif
!    enddo
!  enddo
!
!
!  ! FULL EXTRA FIXED grid contributions
!  do i2p_nuc = 1, nucl_num
!    do i2p_rad = 1, n_points_rad_extra_grid - 1
!      do i2p_ang = 1, n_points_ang_extra_grid
!        ! Find r'_2 point from i2p_rad, i2p_ang loop-indices
!        r2p(1:3) = grid_points_extra_per_atom(1:3,i2p_ang,i2p_rad,i2p_nuc)
!        distance = norm2(r2(1:3) - r2p(1:3))
!        ! Compute matrix element only when there is no divergence 
!        if (distance.gt.1.d-10) then
!          ! Find r'_2 weight from i2p loop-index
!          w2p = grid_float_weights(i2p_ang,i2p_rad,1)
!          ! Need to compute all the MOs in the new floating point
!          call give_all_mos_at_r(r2p, orbitals_in_r2p)
!          ! Get value of the j-th orbital at r'_2
!          ao_j_r2p = orbitals_in_r2p(j)
!          ! Compute pair Jastrow factor between r1 and r'_2
!          j_factor = exp(j_mu_env(r1,r2p,mu_erf))
!          ! Initialize kernel
!          kernel = 0.d0
!          ! Loop over core orbitals to update the kernel
!          do m = 1, n_core_pseudo_orb
!            m_core = list_core_pseudo(m)
!            kernel -= mos_in_r_array_extra_omp(m_core,i2) * orbitals_in_r2p(m_core)    
!          enddo
!          ! Update of the integral at each new point r2p
!          res += kernel * j_factor * ao_j_r2p * w2p / distance 
!        endif
!      enddo
!    enddo
!  enddo
!end function


