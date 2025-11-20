! This file contains test-only providers, that is, providers that are just
! used to check that production providers behave well.
! These providers are numerical approximations

BEGIN_PROVIDER [ double precision, core_tcxc_j0_grid1eje, (ao_num, ao_num, ao_num, ao_num)]
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
  core_tcxc_j0_grid1eje(:,:,:,:) = 0.d0

  ! do-loop version of the tensor product
  do j = 1, ao_num
    do l = 1, ao_num
      do k = 1, ao_num
        do i = 1, ao_num
          core_tcxc_j0_grid1eje(i,k,l,j) = ao_overlap_grid1(i,k) * core_xpot_eje(j,l)
        end do
      end do
    end do
  end do

  !! vectorized version of the tensor product
  !core_tcxc_j0_grid1eje = reshape( ao_overlap_grid1, shape=[ao_num,ao_num,1,1] ) &
  !                     * reshape( core_xpot_eej, shape=[1,1,ao_num,ao_num] )
END_PROVIDER




BEGIN_PROVIDER [ double precision, core_xpot_eje, (ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! Numerical evaluation of < kl | V^{\text{TC}}_{x,\text{core}}  | ij >
  !
  ! $$\braket{kl | V^{\text{TC}}_{x,\text{core}} | ij}
  !     \braket{kl | e^{-J} V_{x,\text{core}} e^{+J} | ij}
  !       \int d1 \chi_k(1)^* \chi_i(1) 
  !         \int d2 \chi_l^*(2) e^{-u(1,2)} 
  !           \int d2' v_x^{\text{core}}(2,2') e^{u(1,2')} \chi_j(2') $$
  !
  ! 1st integral is over usual grid
  ! 2nd integral is over extra grid
  ! 3rd integral is over extra grid
  END_DOC
  !
  integer :: i1, i2
    ! Grid indices
  integer :: i, j, k, l
    ! AOs indices
  double precision :: w1, r1(3)
    ! Grid 1 point and weight
  double precision :: w2, r2(3)
    ! Grid 2 point and weight
  double precision :: ao_i_r1, ao_k_r1, ao_l_r2, v_x_core_ao_j_tc
    ! Real values of MOs and core-exchange kernel
  double precision :: j_r1r2, j_r1r2p
    ! Jastrow factor between r1 and r2 and between r1 and r2p
  double precision, external :: j_mu_env, core_tcxc_kernel_grid1eje
    ! Temporary external definition, we should put them in a module an import that

  ! INTERNALS
  integer :: i2p 
    ! Nested grid 2' point index 
  double precision :: r2p(3), w2p
    ! Nested grid 2' point (r'_2) and weight (dr'_2)
  double precision :: ao_j_r2p
    ! Real value of MO $\chi_j(r'_2)$, index retrieved from passed argument j
  double precision :: distance
    ! Distance for Coulomb integral
  integer :: m, m_core
    ! Indices to loop over core orbitals
  double precision :: kernel
    ! Kernel of the exchange potential

  core_xpot_eje(:,:) = 0.d0

  do i1 = 1, n_points_final_grid
    r1(1:3) = final_grid_points(1:3,i1)
    w1 = final_weight_at_r_vector(i1)
    do i2 = 1, n_points_extra_final_grid
      r2(1:3) = final_grid_points_extra(1:3,i2)
      w2 = final_weight_at_r_vector_extra(i2) 
      ! Compute pair Jastrow factor between r1 and r2
      if (core_tcxc_j0_testing) then
        j_r1r2 = 1.d0
      else
        j_r1r2 = exp(-j_mu_env(r1,r2,mu_erf))
      end if
      do j = 1, ao_num
        ! Compute nested integral in r2' (r1, i2 and AO index j-th as parameters)
        v_x_core_ao_j_tc = 0.d0

        do i2p = 1, n_points_extra_final_grid
          ! Find r'_2 point from i2p loop-index
          r2p(1:3) = final_grid_points_extra(1:3,i2p)
          distance = norm2(r2(1:3) - r2p(1:3))
          ! Compute matrix element only when there is no divergence 
          if (distance.gt.1.d-10) then
            ! Find r'_2 weight from i2p loop-index
            w2p = final_weight_at_r_vector_extra(i2p) 
            ! Get value of the j orbital at r'_2
            ao_j_r2p = aos_in_r_array_extra(j,i2p)
            ! Compute pair Jastrow factor between r1 and r'_2
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
              kernel -= mos_in_r_array_extra_omp(m_core,i2) * mos_in_r_array_extra_omp(m_core,i2p)
            enddo
            kernel = kernel / distance
            ! Update of the integral at each new point r2p
            v_x_core_ao_j_tc += w2p * kernel * j_r1r2p * ao_j_r2p
          endif
        enddo

        do l = 1, ao_num
          ao_l_r2 = aos_in_r_array_extra(l,i2)
          core_xpot_eje(l,j) += w2 * ao_l_r2 * v_x_core_ao_j_tc
        enddo

      enddo
    enddo
  enddo
END_PROVIDER
