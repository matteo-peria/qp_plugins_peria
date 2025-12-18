 BEGIN_PROVIDER [ double precision, int2b_core_xc_ao_grid2a, (ao_num, ao_num)]
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
  double precision :: wall0,wall1
  print*,'providing int2b_core_xc_ao_grid2a ...'
  call wall_time(wall0)


  ! Initialize floating grid points and weights
  grid_float_points(:,:,:,:) = 0.d0
  grid_float_weights(:,:,:) = 0.d0

  int2b_core_xc_ao_grid2a(:,:) = 0.d0

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
              int2b_core_xc_ao_grid2a(l,j) += w2 * ao_l_r2 * integral
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
                int2b_core_xc_ao_grid2a(l,j) += w2 * ao_l_r2 * integral
              enddo
            enddo

          end if

        end do  ! grid 3 fixed
      end do    ! grid 3 fixed
    end do      ! grid 3 fixed

  enddo ! grid 2
  call wall_time(wall1)
  print*,'time to provide int2b_core_xc_ao_grid2a: ',wall1-wall0

END_PROVIDER


 BEGIN_PROVIDER [ double precision, int3b_ao_overlap_grid1_w_corexc_grid2a, (ao_num, ao_num, ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! TEST-ONLY PROVIDER.
  ! Check that the TC exchange integral make sense when the Jastrow factor 
  ! is set to 1 (that is, when the Jastrow exponent is 0).

  ! $$ \braket{kl | e^{-J} V_{x,\text{core}} e^{+J} | ij} |_{J=0}
  !      \int d1 \chi_k(1)^* \chi_i(1) 
  !        \int d2 \chi_l^*(2) \int d2' v_x^{\text{core}}(2,2') \chi_j(2') $$
  !
  ! which, due to independence of the integrals, can be seen as the tensor 
  ! product between the AO overlap matrix and the AO-based 
  ! core-restrice-exchange-potential:
  !
  !  $$ \braket{k|i} \braket{l| V_{x,\text{core}} |j} $$
  !
  ! 1st integral is computed along usual grid
  ! 2nd integral is computed along grid2
  ! 3rd integral is computed along adaptive grid (floating + extra grid)
  END_DOC
  integer :: i,j,k,l
  double precision :: wall0,wall1
  print*,'providing int3b_ao_overlap_grid1_w_corexc_grid2a ...'

  ! Initialization
  int3b_ao_overlap_grid1_w_corexc_grid2a(:,:,:,:) = 0.d0

  if (core_tcxc_omp) then
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP& PRIVATE (i, j, k, l) &
    !$OMP& SHARED ( ao_num, int3b_ao_overlap_grid1_w_corexc_grid2a          &
    !$OMP&        , ao_overlap_grid1, int2b_core_xc_ao_grid2a )
    !$OMP DO COLLAPSE(4) SCHEDULE(static)
    do j = 1, ao_num
      do l = 1, ao_num
        do k = 1, ao_num
          do i = 1, ao_num
            !int3b_ao_overlap_grid1_w_corexc_grid2a(i,k,l,j) = ao_overlap_grid1(i,k) * int2b_core_xc_ao_grid2a(j,l)
            ! test:
            int3b_ao_overlap_grid1_w_corexc_grid2a(i,k,l,j) = ao_overlap_grid1(i,k) * int2b_core_xc_ao_grid2a(l,j)
          end do
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  else
    do j = 1, ao_num
      do l = 1, ao_num
        do k = 1, ao_num
          do i = 1, ao_num
            !int3b_ao_overlap_grid1_w_corexc_grid2a(i,k,l,j) = ao_overlap_grid1(i,k) * int2b_core_xc_ao_grid2a(j,l)
            ! test:
            int3b_ao_overlap_grid1_w_corexc_grid2a(i,k,l,j) = ao_overlap_grid1(i,k) * int2b_core_xc_ao_grid2a(l,j)
          end do
        end do
      end do
    end do
  end if
  call wall_time(wall1)
  print*,'time to provide int3b_ao_overlap_grid1_w_corexc_grid2a : ',wall1-wall0

END_PROVIDER


! Here we compute two computation-heavy providers at once to recycle the same
! loops and calculations. Since these providers are used only in tests where 
! both are needed, alternative, indepdent versions of these providers exist too 
! in case we need to test only one of them
 BEGIN_PROVIDER [ double precision, int3b_core_tcxc_ao_grid12a_loops_test, (ao_num, ao_num, ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, int2b_core_tcxc_ao_grid2a_at_r1_test, (ao_num, ao_num, n_points_final_grid) ]
  integer :: i1, i2, i2p
    ! Grids indices
  integer :: i, j, k, l
    ! AOs indices
  double precision :: ao_i_r1, ao_k_r1, ao_l_r2, ao_j_r2p
    ! Real values of AOs and core-exchange kernel
  double precision :: r1(3), w1
    ! Grid 1 point and weight
  double precision :: r2(3), w2
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

  ! Adaptive grid variables
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

  int2b_core_tcxc_ao_grid2a_at_r1_test(:,:,:) = 0.d0
  int3b_core_tcxc_ao_grid12a_loops_test(:,:,:,:) = 0.d0

  do i1 = 1, n_points_final_grid
    r1(1:3) = final_grid_points(1:3,i1)
    w1 = final_weight_at_r_vector(i1) ! necessary only to int3b_core_tcxc_ao_grid12a_loops_test
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

            ! Notice the order of the indices: 
            ! - in AO_CORE_XC_MAT_GRID12ADAPT we have (j,l,i1) 
            !   because it will be transposed later!
            ! - in CORE_TCXC_GRID12ADAPT_LOOPS we have (i,k,j,l) 
            !   because it must be already in the 'right' order

            ! Option 1: OMP parallelism
            if (core_tcxc_omp) then
              ! Since int2b_core_tcxc_ao_grid2a_at_r1_test is updated at the second nested loop,
              ! only the first two loops can be OMP-parallelised

              !$OMP PARALLEL DEFAULT(NONE) &
              !$OMP& PRIVATE (i, j, k, l, ao_core_xc_r2, v_xc_core_ao_j_r2p) &
              !$OMP& SHARED ( int2b_core_tcxc_ao_grid2a_at_r1_test                         &
              !$OMP&        , int3b_core_tcxc_ao_grid12a_loops_test                  &
              !$OMP&        , aos_in_r_array, aos_in_r_array2, aos_in_r2p    &
              !$OMP&        , ao_num, i1, i2, w1, w2, w2p, j_r1r2, j_r1r2p, kernel)
              !$OMP DO COLLAPSE(2) SCHEDULE(static)
              do j = 1, ao_num
                do l = 1, ao_num
                  v_xc_core_ao_j_r2p = w2p * j_r1r2p * kernel * aos_in_r2p(j)
                  ao_core_xc_r2 = w2 * j_r1r2 * aos_in_r_array2(l,i2) * v_xc_core_ao_j_r2p
                  int2b_core_tcxc_ao_grid2a_at_r1_test(j,l,i1) += ao_core_xc_r2
                  ! start loops necessary only to int3b_core_tcxc_ao_grid12a_loops_test
                  do k = 1, ao_num
                    do i = 1, ao_num
                      !int3b_core_tcxc_ao_grid12a_loops_test(i,k,l,j) += &
                      int3b_core_tcxc_ao_grid12a_loops_test(i,k,j,l) += &
                          w1 * aos_in_r_array(i,i1) * aos_in_r_array(k,i1) * ao_core_xc_r2
                    end do
                  end do
                  ! end loops necessary only to int3b_core_tcxc_ao_grid12a_loops_test
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
                  int2b_core_tcxc_ao_grid2a_at_r1_test(j,l,i1) += w2 * j_r1r2 * ao_l_r2 &
                                                  * w2p * j_r1r2p * kernel * ao_j_r2p
                  ! start loops necessary only to int3b_core_tcxc_ao_grid12a_loops_test
                  do k = 1, ao_num
                    do i = 1, ao_num
                      !int3b_core_tcxc_ao_grid12a_loops_test(i,k,l,j) += &
                      int3b_core_tcxc_ao_grid12a_loops_test(i,k,j,l) += &
                          &   w1  * aos_in_r_array(i,i1) * aos_in_r_array(k, i1) &
                          & * w2  * j_r1r2               * aos_in_r_array2(l, i2) &
                          & * w2p * j_r1r2p * kernel     * aos_in_r2p(j)
                    end do
                  end do
                  ! end loops necessary only to int3b_core_tcxc_ao_grid12a_loops_test
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

              ! Notice the order of the indices: 
              ! - in AO_CORE_XC_MAT_GRID12ADAPT we have (j,l,i1) 
              !   because it will be transposed later!
              ! - in CORE_TCXC_GRID12ADAPT_LOOPS we have (i,k,j,l) 
              !   because it must be already in the 'right' order
    
              ! Option 1: OMP parallelism
              if (core_tcxc_omp) then
                ! Since int2b_core_tcxc_ao_grid2a_at_r1_test is updated at the second nested loop,
                ! only the first two loops can be OMP-parallelised

                !$OMP PARALLEL DEFAULT(NONE) &
                !$OMP& PRIVATE (i, j, k, l, ao_core_xc_r2, v_xc_core_ao_j_r2p) &
                !$OMP& SHARED ( int2b_core_tcxc_ao_grid2a_at_r1_test           &
                !$OMP&        , int3b_core_tcxc_ao_grid12a_loops_test          &
                !$OMP&        , aos_in_r_array, aos_in_r_array2                &
                !$OMP&        , aos_in_r_array_extra_full                      &
                !$OMP&        , ao_num, i1, i2, w1, w2, w2p                    &
                !$OMP&        , i2p_ang,i2p_rad,i2p_nuc, j_r1r2, j_r1r2p, kernel)
                !$OMP DO COLLAPSE(2) SCHEDULE(static)
                do j = 1, ao_num
                  do l = 1, ao_num
                    v_xc_core_ao_j_r2p = w2p * j_r1r2p * kernel * aos_in_r_array_extra_full(j,i2p_ang,i2p_rad,i2p_nuc)
                    ao_core_xc_r2 = w2 * j_r1r2 * aos_in_r_array2(l,i2) * v_xc_core_ao_j_r2p
                    int2b_core_tcxc_ao_grid2a_at_r1_test(j,l,i1) += ao_core_xc_r2
                    ! start loops necessary only to int3b_core_tcxc_ao_grid12a_loops_test
                    do k = 1, ao_num
                      do i = 1, ao_num
                        !int3b_core_tcxc_ao_grid12a_loops_test(i,k,l,j) += &
                        int3b_core_tcxc_ao_grid12a_loops_test(i,k,j,l) += &
                            w1 * aos_in_r_array(i,i1) * aos_in_r_array(k,i1) * ao_core_xc_r2
                      end do
                    end do
                    ! end loops necessary only to int3b_core_tcxc_ao_grid12a_loops_test
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
                    int2b_core_tcxc_ao_grid2a_at_r1_test(j,l,i1) += w2 * j_r1r2 * ao_l_r2 &
                                                    * w2p * j_r1r2p * kernel * ao_j_r2p

                    ! start loops necessary only to int3b_core_tcxc_ao_grid12a_loops_test
                    do k = 1, ao_num
                      do i = 1, ao_num
                        !int3b_core_tcxc_ao_grid12a_loops_test(i,k,l,j) += &
                        int3b_core_tcxc_ao_grid12a_loops_test(i,k,j,l) += &
                            &   w1  * aos_in_r_array(i,i1) * aos_in_r_array(k, i1) &
                            & * w2  * j_r1r2               * aos_in_r_array2(l, i2) &
                            & * w2p * j_r1r2p * kernel     * aos_in_r_array_extra_full(j,i2p_ang,i2p_rad,i2p_nuc)
                      end do
                    end do
                    ! end loops necessary only to int3b_core_tcxc_ao_grid12a_loops_test
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
