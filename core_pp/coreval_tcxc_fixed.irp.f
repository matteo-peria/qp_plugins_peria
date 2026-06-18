BEGIN_PROVIDER [ double precision, int3b_core_tcxc_ao_grid12e, (ao_num, ao_num, ao_num, ao_num)]
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
  integer :: i1, i2, i2p
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
  double precision, external :: j_mu_env
    ! Temporary external definition, we should put them in a module an import that
  double precision :: integral_j_r2p, kernel
    ! Temporary variables to stock mid-computation result
  double precision :: distance
    ! Distance for Coulomb integral (exchange)
  integer :: ao_num2

  int3b_core_tcxc_ao_grid12e(:,:,:,:) = 0.d0

  ! dgemm session
  ao_num2 = ao_num*ao_num
  call dgemm( 'N', 'T', ao_num2, ao_num2, n_points_final_grid, 1.d0 &
            , ao_overlap_grid1_at_r1, ao_num2    &
            , int2b_core_tcxc_ao_grid2e_at_r1, ao_num2    &
            , 0.d0, int3b_core_tcxc_ao_grid12e, ao_num2 )

END_PROVIDER


BEGIN_PROVIDER [ double precision, int2b_core_tcxc_ao_grid2e_at_r1, (ao_num, ao_num, n_points_final_grid) ]
  integer :: i1, i2, i2p
    ! Grids indices
  integer :: j, l
    ! AOs indices
  double precision :: r1(3)
    ! Grid 1 point
  double precision :: w2, r2(3)
    ! Grid 2 point and weight
  double precision :: r2p(3), w2p
    ! Nested grid 2' point (r'_2) and weight (dr'_2)
  double precision :: ao_l_r2, ao_j_r2p
    ! Real values of AOs and core-exchange kernel
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

  int2b_core_tcxc_ao_grid2e_at_r1(:,:,:) = 0.d0

  do i1 = 1, n_points_final_grid
    r1(1:3) = final_grid_points(1:3,i1)

    do i2 = 1, n_points_final_grid2
      r2(1:3) = final_grid_points2(1:3,i2)
      w2 = final_weight_at_r_vector2(i2)
      ! Compute pair Jastrow factor between r1 and r2
      if (core_tcxc_j0_testing) then
        j_r1r2 = 1.d0
      else
        j_r1r2 = exp(-j_mu_env(r1,r2,mu_erf))
      end if

      ! Loop over the 3rd grid
      do i2p = 1, n_points_extra_final_grid
        r2p(1:3) = final_grid_points_extra(1:3,i2p)
        distance = norm2(r2(1:3) - r2p(1:3))
        ! Compute matrix element only when there is no divergence 
        if (distance.gt.1.d-10) then
          w2p = final_weight_at_r_vector_extra(i2p)
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
            kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r_array_extra_omp(m_core,i2p)
          enddo
          kernel = kernel / distance

          !$OMP PARALLEL DEFAULT(NONE) &
          !$OMP& PRIVATE (j, l) &
          !$OMP& SHARED ( int2b_core_tcxc_ao_grid2e_at_r1         &
          !$OMP&        , aos_in_r_array2, aos_in_r_array_extra   &
          !$OMP&        , ao_num, i1, i2, i2p &
          !$OMP&        , w2, w2p, j_r1r2, j_r1r2p, kernel)
          !$OMP DO COLLAPSE(2) SCHEDULE(static)
          do l = 1, ao_num
            do j = 1, ao_num
              ! Notice the index order, because this array will be transposed later!
              int2b_core_tcxc_ao_grid2e_at_r1(j,l,i1) += w2 * j_r1r2 * aos_in_r_array2(l,i2) * w2p * kernel * aos_in_r_array_extra(j,i2p) * j_r1r2p
            enddo
          enddo
          !$OMP END DO
          !$OMP END PARALLEL

        endif

      enddo !r2p
    enddo !r2

  end do !r1

END_PROVIDER
