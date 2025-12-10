subroutine get_all_adaptive_grid( all_float_grids_points &
                                , all_float_grids_weight &
                                , all_fixed_grids_weight &
                                , all_float_grids_aos    &
                                , all_float_grids_mos_core)
  integer :: i2
    ! Grids indices
  integer :: j, l
    ! AOs indices
  double precision :: ao_l_r2, ao_j_r2p
    ! Real values of AOs and core-exchange kernel
  double precision :: w2, r2(3)
    ! Grid 2 point and weight
  double precision :: r2p(3), w2p
    ! Nested grid 2' point (r'_2) and weight (dr'_2)
  double precision :: distance
    ! Distance for Coulomb integral
  integer :: m, m_core
    ! Indices to loop over core orbitals
  double precision :: kernel
    ! Kernel of the exchange potential

  ! ADAPTIVE GRID VARIABLES
  integer :: i2p_nuc, i2p_rad, i2p_ang
    ! Floating grid indices
  integer :: i2p
    ! Pruned floating grid indices (to be done)

  ! OUTPUT
  double precision, intent(out) :: all_float_grids_points(   &        
                                     3                       &
                                   , n_points_ang_float_grid &
                                   , n_points_rad_float_grid &
                                   , n_points_final_grid2    &
  )
  double precision, intent(out) :: all_float_grids_weight(   &
                                     n_points_ang_float_grid &
                                   , n_points_rad_float_grid &
                                   , n_points_final_grid2    & 
  )
  double precision, intent(out) :: all_fixed_grids_weight(              &
                                     n_points_extra_integration_angular &
                                   , n_points_extra_radial_grid         &
                                   , nucl_num                           &
                                   , n_points_final_grid2               &
  )
  double precision, intent(out) :: all_float_grids_aos(      &
                                     ao_num                  &
                                   , n_points_ang_float_grid &
                                   , n_points_rad_float_grid &
                                   , n_points_final_grid2    &
  )
  double precision, intent(out) :: all_float_grids_mos_core( &
                                     n_core_pseudo_orb       &
                                   , n_points_ang_float_grid &
                                   , n_points_rad_float_grid &
                                   , n_points_final_grid2    &
  )


  !double precision, intent(out), allocatable :: all_float_grids_points(:,:,:,:)
  !double precision, intent(out), allocatable :: all_float_grids_weight(:,:,:)
  !double precision, intent(out), allocatable :: all_fixed_grids_weight(:,:,:,:)
  !double precision, intent(out), allocatable :: all_float_grids_aos(:,:,:,:)
  !double precision, intent(out), allocatable :: all_float_grids_mos_core(:,:,:,:)

  ! Fixed extra grid weights
  ! To be used for pruning, NOT IMPLEMENTED YET
  integer :: n_fixed_pts_effective(nucl_num)
  integer :: n_float_pts_effective
  integer :: n_pts_effective_max

!  allocate(all_float_grids_points(   3                       &
!                                   , n_points_ang_float_grid &
!                                   , n_points_rad_float_grid &
!                                   , n_points_final_grid2    &
!  ))
!  allocate(all_float_grids_weight(   &
!                                     n_points_ang_float_grid &
!                                   , n_points_rad_float_grid &
!                                   , n_points_final_grid2    & 
!  ))
!  allocate(all_fixed_grids_weight(              &
!                                     n_points_extra_integration_angular &
!                                   , n_points_extra_radial_grid         &
!                                   , nucl_num                           &
!                                   , n_points_final_grid2               &
!  ))
!  allocate(all_float_grids_aos(      &
!                                     ao_num                  &
!                                   , n_points_ang_float_grid &
!                                   , n_points_rad_float_grid &
!                                   , n_points_final_grid2    &
!  ))
!  allocate(all_float_grids_mos_core( &
!                                     n_core_pseudo_orb       &
!                                   , n_points_ang_float_grid &
!                                   , n_points_rad_float_grid &
!                                   , n_points_final_grid2    &
!  ))
!
!
!
  do i2 = 1, n_points_final_grid2
    r2(1:3) = final_grid_points2(1:3,i2)
    w2 = final_weight_at_r_vector2(i2)
    ! Compute adaptive grid centered in r2
    call get_adaptive_grid( r2                                 &
                          , grid_points_extra_per_atom         &
                          , all_float_grids_points(:,:,:,i2) &
                          , all_fixed_grids_weight(:,:,:,i2)   &
                          , all_float_grids_weight(:,:,i2)   &
                          , n_fixed_pts_effective              &
                          , n_float_pts_effective              &
                          , n_pts_effective_max)

    ! FLOATING grid contributions
    do i2p_rad = 1, n_points_rad_float_grid - 1
      do i2p_ang = 1, n_points_ang_float_grid
        ! Find r2p point from i2p_rad, i2p_ang loop-indices
        r2p(1:3) = all_float_grids_points(1:3,i2p_ang,i2p_rad,i2)
        distance = norm2(r2(1:3) - r2p(1:3))
        ! Compute matrix element only when there is no divergence 
        if (distance.gt.1.d-10) then
          ! Find r'_2 weight from i2p loop-index
          w2p = all_float_grids_weight(i2p_ang,i2p_rad,i2)
          ! Compute all AOs in r2p 
          call give_all_aos_at_r(r2p, all_float_grids_aos(:,i2p_ang,i2p_rad,i2))
          ! Get the core-MOs at r2p reciclying the newly computed AOs in r2p
          ! and using only the the sub-matrix of core-MOs coefficients
          call dgemv( 'N', n_core_pseudo_orb, ao_num, 1.d0           &
                    , mo_core_coef_notnorm_transp, n_core_pseudo_orb &
                    , all_float_grids_aos(:,i2p_ang,i2p_rad,i2), 1   &
                    , 0.d0, all_float_grids_mos_core(:,i2p_ang,i2p_rad,i2), 1)
          ! Initialize kernel
          kernel = 0.d0
          ! Loop over core orbitals to update the kernel
          do m = 1, n_core_pseudo_orb
            m_core = list_core_pseudo(m)
            ! Update kernel using core-MOs(r2p) array
            kernel -= mos_in_r_array2_omp(m_core,i2) &
                  & * all_float_grids_mos_core(m,i2p_ang,i2p_rad,i2)
          enddo
          kernel = kernel / distance
        end if
      end do
    end do
  end do
end subroutine get_all_adaptive_grid


 BEGIN_PROVIDER [ double precision, int2b_core_tcxc_ao_grid2a_stored_at_r1, (ao_num, ao_num, n_points_final_grid) ]
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

  !! ALL ADAPTIVE GRID VARIABLES
  !double precision, intent(out) :: all_float_grids_points( 3                       &
  !                                          , n_points_ang_float_grid &
  !                                          , n_points_rad_float_grid &
  !                                          , n_points_final_grid2    &
  !)
  !double precision, intent(out) :: all_float_grids_weight( n_points_final_grid2    & 
  !                                          , n_points_rad_float_grid &
  !                                          , n_points_ang_float_grid &
  !)
  !double precision, intent(out) :: all_fixed_grids_weight( n_points_extra_integration_angular &
  !                                          , n_points_extra_radial_grid         &
  !                                          , nucl_num                           &
  !                                          , n_points_final_grid2               &
  !)
  !double precision, intent(out) :: all_float_grids_aos( ao_num                  &
  !                                       , n_points_ang_float_grid &
  !                                       , n_points_rad_float_grid &
  !                                       , n_points_final_grid2    &
  !)
  !double precision, intent(out) :: all_float_grids_mos_core( n_core_pseudo_orb       &
  !                                            , n_points_ang_float_grid &
  !                                            , n_points_rad_float_grid &
  !                                            , n_points_final_grid2    &
  !)
  double precision, allocatable :: all_float_grids_points(:,:,:,:)
  double precision, allocatable :: all_float_grids_weight(:,:,:)
  double precision, allocatable :: all_fixed_grids_weight(:,:,:,:)
  double precision, allocatable :: all_float_grids_aos(:,:,:,:)
  double precision, allocatable :: all_float_grids_mos_core(:,:,:,:)

  integer :: i2p_nuc, i2p_rad, i2p_ang
    ! Third grid indices

  ! Intermediate variables
  double precision :: v_xc_core_ao_j_r2p
  double precision :: ao_core_xc_r2

  ! Initialize provider
  int2b_core_tcxc_ao_grid2a_stored_at_r1(:,:,:) = 0.d0

  allocate(all_float_grids_points(   3                       &
                                   , n_points_ang_float_grid &
                                   , n_points_rad_float_grid &
                                   , n_points_final_grid2    &
  ))
  allocate(all_float_grids_weight(   &
                                     n_points_ang_float_grid &
                                   , n_points_rad_float_grid &
                                   , n_points_final_grid2    & 
  ))
  allocate(all_fixed_grids_weight(              &
                                     n_points_extra_integration_angular &
                                   , n_points_extra_radial_grid         &
                                   , nucl_num                           &
                                   , n_points_final_grid2               &
  ))
  allocate(all_float_grids_aos(      &
                                     ao_num                  &
                                   , n_points_ang_float_grid &
                                   , n_points_rad_float_grid &
                                   , n_points_final_grid2    &
  ))
  allocate(all_float_grids_mos_core( &
                                     n_core_pseudo_orb       &
                                   , n_points_ang_float_grid &
                                   , n_points_rad_float_grid &
                                   , n_points_final_grid2    &
  ))




  ! Precompute all the adaptive points (float), weights (fixed and float), AOs, core MOs
  call get_all_adaptive_grid( all_float_grids_points &
                            , all_float_grids_weight &
                            , all_fixed_grids_weight &
                            , all_float_grids_aos    &
                            , all_float_grids_mos_core)

  do i1 = 1, n_points_final_grid
    r1(1:3) = final_grid_points(1:3,i1)
    w1 = final_weight_at_r_vector(i1) ! this is only necessary for int3b_core_tcxc_ao_grid12a_loops
    do i2 = 1, n_points_final_grid2
      r2(1:3) = final_grid_points2(1:3,i2)
      w2 = final_weight_at_r_vector2(i2)
      ! Compute pair Jastrow factor between r1 and r2
      if (core_tcxc_j0_testing) then
        j_r1r2 = 1.d0
      else
        j_r1r2 = exp(-j_mu_env(r1,r2,mu_erf))
      end if

      !! Compute adaptive grid centered in r2
      !call get_adaptive_grid( r2, grid_points_extra_per_atom, grid_float_points &
      !                      , grid_fixed_weights, grid_float_weights            &
      !                      , n_fixed_pts_effective, n_float_pts_effective      &
      !                      , n_pts_effective_max)

      ! FLOATING grid contributions
      do i2p_rad = 1, n_points_rad_float_grid - 1
        do i2p_ang = 1, n_points_ang_float_grid
          ! Find r2p point from i2p_rad, i2p_ang loop-indices
          r2p(1:3) = all_float_grids_points(1:3,i2p_ang,i2p_rad,i2)
          distance = norm2(r2(1:3) - r2p(1:3))
          ! Compute matrix element only when there is no divergence 
          if (distance.gt.1.d-10) then
            ! Find r'_2 weight from i2p loop-index
            w2p = all_float_grids_weight(i2p_ang,i2p_rad,i2)
            !! Compute all AOs in r2p 
            !call give_all_aos_at_r(r2p, aos_in_r2p)

            !! Get the core-MOs at r2p reciclying the newly computed AOs in r2p
            !! and using only the the sub-matrix of core-MOs coefficients
            !call dgemv( 'N', n_core_pseudo_orb, ao_num, 1.d0 &
            !          , mo_core_coef_notnorm_transp, n_core_pseudo_orb &
            !          , aos_in_r2p, 1, 0.d0, mos_core_in_r2p, 1)

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
              kernel -= mos_in_r_array2_omp(m_core,i2) &
                     &* all_float_grids_mos_core(m,i2p_ang,i2p_rad,i2)    
            enddo
            kernel = kernel / distance

            ! Notice the order of the indices: 
            ! - in AO_CORE_XC_MAT_GRID12ADAPT we have (j,l,i1) 
            !   because it will be transposed later!
            ! - in CORE_TCXC_GRID12ADAPT_LOOPS we have (i,k,j,l) 
            !   because it must be already in the 'right' order

            ! Option 1: OMP parallelism
            if (core_tcxc_omp) then
              ! Since int2b_core_tcxc_ao_grid2a_stored_at_r1 is updated at the second nested loop,
              ! only the first two loops can be OMP-parallelised

              !$OMP PARALLEL DEFAULT(NONE) &
              !$OMP& PRIVATE (j, l) &
              !$OMP& SHARED ( int2b_core_tcxc_ao_grid2a_stored_at_r1 &
              !$OMP&        , aos_in_r_array2        &
              !$OMP&        , all_float_grids_aos             &
              !$OMP&        , ao_num, i1, i2, i2p_ang, i2p_rad        &
              !$OMP&        , w1, w2, w2p, j_r1r2, j_r1r2p, kernel)
              !$OMP DO COLLAPSE(2) SCHEDULE(static)
              do j = 1, ao_num
                do l = 1, ao_num
                  int2b_core_tcxc_ao_grid2a_stored_at_r1(j,l,i1) += &
                      &    w2 * j_r1r2 * aos_in_r_array2(l,i2) &
                      & * w2p * j_r1r2p * kernel * all_float_grids_aos(j,i2p_ang,i2p_rad,i2)
                end do
              end do
              !$OMP END DO
              !$OMP END PARALLEL

            ! Option 2: NO parallelism 
            else
              do l = 1, ao_num
                ao_l_r2 = aos_in_r_array2(l,i2)
                do j = 1, ao_num
                  ao_j_r2p = all_float_grids_aos(j,i2p_ang,i2p_rad,i2)
                  int2b_core_tcxc_ao_grid2a_stored_at_r1(j,l,i1) += &
                      &    w2 * j_r1r2 * ao_l_r2 &
                      & * w2p * j_r1r2p * kernel * ao_j_r2p
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
              w2p = all_fixed_grids_weight(i2p_ang,i2p_rad,i2p_nuc,i2)
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
                kernel -= mos_in_r_array2_omp(m_core,i2) &
                       & *mos_in_r_array_extra_full_omp(m_core,i2p_ang,i2p_rad,i2p_nuc)
              enddo
              kernel = kernel/distance

              ! Notice the order of the indices: 
              ! - in AO_CORE_XC_MAT_GRID12ADAPT we have (j,l,i1) 
              !   because it will be transposed later!
              ! - in CORE_TCXC_GRID12ADAPT_LOOPS we have (i,k,j,l) 
              !   because it must be already in the 'right' order
    
              ! Option 1: OMP parallelism
              if (core_tcxc_omp) then
                ! Since int2b_core_tcxc_ao_grid2a_stored_at_r1 is updated at the second nested loop,
                ! only the first two loops can be OMP-parallelised

                !$OMP PARALLEL DEFAULT(NONE) &
                !$OMP& PRIVATE (j, l) &
                !$OMP& SHARED ( int2b_core_tcxc_ao_grid2a_stored_at_r1                         &
                !$OMP&        , aos_in_r_array2                &
                !$OMP&        , aos_in_r_array_extra_full                      &
                !$OMP&        , ao_num, i1, i2, w1, w2, w2p                    &
                !$OMP&        , i2p_ang,i2p_rad,i2p_nuc, j_r1r2, j_r1r2p, kernel)
                !$OMP DO COLLAPSE(2) SCHEDULE(static)
                do j = 1, ao_num
                  do l = 1, ao_num
                    int2b_core_tcxc_ao_grid2a_stored_at_r1(j,l,i1) += &
                      &   w2 * j_r1r2 * aos_in_r_array2(l,i2)  &
                      & * w2p * j_r1r2p * kernel * aos_in_r_array_extra_full(j,i2p_ang,i2p_rad,i2p_nuc)
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
                    int2b_core_tcxc_ao_grid2a_stored_at_r1(j,l,i1) += w2 * j_r1r2 * ao_l_r2 &
                                                    * w2p * j_r1r2p * kernel * ao_j_r2p
  
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
