BEGIN_PROVIDER [ double precision, core_tcxc_grid12ej_prun, (ao_num, ao_num, ao_num, ao_num)]
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

  core_tcxc_grid12ej_prun(:,:,:,:) = 0.d0

  do i1 = 1, n_points_final_grid
    !write(*,*) "Loop 1: i1 = ", i1
    r1(1:3) = final_grid_points(1:3,i1)
    w1 = final_weight_at_r_vector(i1)

    do i2 = 1, n_points_final_grid2
      !write(*,*) "Loop 2: i2 = ", i2
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
            kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r_array_extra_omp(m_core,i2p)
          enddo
          kernel = kernel / distance

          ! Loop over all AO (j-th index, variable is r2p)
          do j = 1, ao_num
            ! Get value of the j-th orbital at r2p
            ao_j_r2p = aos_in_r_array_extra(j,i2p)
            integral_j_r2p = w2p * kernel * ao_j_r2p * j_r1r2p
            do l = 1, ao_num
              !write(*,*) "Loop 5.b.2: l-th AO = ", l
              ao_l_r2 = aos_in_r_array2(l,i2)
              do k = 1, ao_num
                !write(*,*) "Loop 6.2: k-th AO = ", k
                ao_k_r1 = aos_in_r_array(k,i1)
                do i = 1, ao_num
                  !write(*,*) "Loop 7.2: i-th AO = ", i
                  ao_i_r1 = aos_in_r_array(i,i1)
                  core_tcxc_grid12ej_prun(i,k,l,j) += w1 * ao_i_r1 * ao_k_r1 &
                                                   * w2 * j_r1r2 * ao_l_r2 * integral_j_r2p
                enddo ! loop i
              enddo  ! loop k
            enddo  ! loop l
          enddo  ! loop j


        endif

      end do !loop i2p

    enddo  ! loop i2
  enddo  ! loop i1

END_PROVIDER
