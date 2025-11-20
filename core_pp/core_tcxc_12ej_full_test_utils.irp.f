BEGIN_PROVIDER [ double precision, core_xpot_grid2ej_full, (ao_num, ao_num)]
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
  integer :: i1, i2, i2p_nuc, i2p_rad, i2p_ang
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
  double precision, external :: j_mu_env!, core_tcxc_kernel_grid123
    ! Temporary external definition, we should put them in a module an import that
  double precision :: integral, kernel
    ! Temporary variables to stock mid-computation result
  double precision :: distance
    ! Distance for Coulomb integral (exchange)

  core_xpot_grid2ej_full(:,:) = 0.d0

  do i2 = 1, n_points_final_grid2
    !write(*,*) "Loop 2: i2 = ", i2
    r2(1:3) = final_grid_points2(1:3,i2)
    w2 = final_weight_at_r_vector2(i2)

    ! Loop over the 3rd grid
    do i2p_nuc = 1, nucl_num
      do i2p_rad = 1, n_points_extra_radial_grid - 1
        do i2p_ang = 1, n_points_extra_integration_angular
          !write(*,*) "Loops 3 (extra): i2p_nuc = ", i2p_nuc, "  i2p_rad = ", i2p_rad, "  i2p_ang = ", i2p_ang
          r2p(1:3) = grid_points_extra_per_atom(1:3,i2p_ang,i2p_rad,i2p_nuc)
          distance = norm2(r2(1:3) - r2p(1:3))
          ! Compute matrix element only when there is no divergence 
          if (distance.gt.1.d-10) then
            w2p = final_weight_at_r_extra(i2p_ang,i2p_rad,i2p_nuc)

            ! Initialize kernel
            kernel = 0.d0
            ! Loop over core orbitals to update the kernel
            do m = 1, n_core_pseudo_orb
              m_core = list_core_pseudo(m)
              kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r_array_extra_full_omp(m_core,i2p_ang,i2p_rad,i2p_nuc)
            enddo
            kernel = kernel / distance

            ! Loop over all AO (j-th index, variable is r2p)
            do j = 1, ao_num
              ! Get value of the j-th orbital at r2p
              ao_j_r2p = aos_in_r_array_extra_full(j,i2p_ang,i2p_rad,i2p_nuc)
              ! Update of the integral at each new point r2p
              integral = w2p * kernel * ao_j_r2p
      
              do l = 1, ao_num
                !write(*,*) "Loop 5.b.2: l-th AO = ", l
                ao_l_r2 = aos_in_r_array2(l,i2)
                core_xpot_grid2ej_full(l,j) += w2 * ao_l_r2 * integral
              enddo  ! loop l
      
            enddo  ! loop j

          endif

        end do ! loop 3rd grid
      end do   ! loop 3rd grid
    end do     ! loop 3rd grid

  enddo  ! loop i2

END_PROVIDER


BEGIN_PROVIDER [ double precision, core_tcxc_j0_grid12ej_full, (ao_num, ao_num, ao_num, ao_num)]
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
  core_tcxc_j0_grid12ej_full(:,:,:,:) = 0.d0

  ! do-loop version of the tensor product
  do j = 1, ao_num
    do l = 1, ao_num
      do k = 1, ao_num
        do i = 1, ao_num
          core_tcxc_j0_grid12ej_full(i,k,l,j) = ao_overlap_grid1(i,k) * core_xpot_grid2ej_full(j,l)
        end do
      end do
    end do
  end do

  !! vectorized version of the tensor product
  !core_tcxc_j0_grid12ej_full = reshape( ao_overlap_grid1, shape=[ao_num,ao_num,1,1] ) &
  !                     * reshape( core_xpot_grid22, shape=[1,1,ao_num,ao_num] )
END_PROVIDER
