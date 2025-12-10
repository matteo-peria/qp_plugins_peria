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
  double precision :: j_factor, j_r1r2, j_r1r2p
    ! Jastrow factors between r1 and r2 and between r1 and r2p
  double precision, external :: j_mu_env!, core_tcxc_kernel_grid123
    ! Temporary external definition, we should put them in a module an import that
  double precision :: integral, kernel
    ! Temporary variables to stock mid-computation result
  double precision :: distance
    ! Distance for Coulomb integral (exchange)
  integer :: ao_num2

  int3b_core_tcxc_ao_grid12e(:,:,:,:) = 0.d0

  if (core_tcxc_loops) then
    do i1 = 1, n_points_final_grid
      r1(1:3) = final_grid_points(1:3,i1)
      w1 = final_weight_at_r_vector(i1)

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
        do i2p_nuc = 1, nucl_num
          do i2p_rad = 1, n_points_extra_radial_grid - 1
            do i2p_ang = 1, n_points_extra_integration_angular

              r2p(1:3) = grid_points_extra_per_atom(1:3,i2p_ang,i2p_rad,i2p_nuc)
              distance = norm2(r2(1:3) - r2p(1:3))
              ! Compute matrix element only when there is no divergence 
              if (distance.gt.1.d-10) then
                w2p = final_weight_at_r_extra(i2p_ang,i2p_rad,i2p_nuc)

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
                kernel = kernel / distance

                ! Loop over all AO (j-th index, variable is r2p)
                do j = 1, ao_num
                  ! Get value of the j-th orbital at r2p
                  ao_j_r2p = aos_in_r_array_extra_full(j,i2p_ang,i2p_rad,i2p_nuc)
                  ! Update of the integral at each new point r2p
                  integral = w2p * kernel * j_r1r2p * ao_j_r2p  ! BETTER RENAMING v_xc_core_ao_j_r2p
                  do l = 1, ao_num
                    ao_l_r2 = aos_in_r_array2(l,i2)
                    do k = 1, ao_num
                      ao_k_r1 = aos_in_r_array(k,i1)
                      do i = 1, ao_num
                        ao_i_r1 = aos_in_r_array(i,i1)
                        int3b_core_tcxc_ao_grid12e(i,k,l,j) += w1 * ao_i_r1 * ao_k_r1 &
                                                          * w2 * j_r1r2 * ao_l_r2 * integral
                      enddo ! loop i
                    enddo  ! loop k
                  enddo  ! loop l
                enddo  ! loop j

              endif

            end do ! loop 3rd grid
          end do   ! loop 3rd grid
        end do     ! loop 3rd grid

      enddo  ! loop i2
    enddo  ! loop i1

  else
    ! dgemm session
    ao_num2 = ao_num*ao_num
    call dgemm( 'N', 'T', ao_num2, ao_num2, n_points_final_grid, 1.d0 &
              , ao_overlap_grid1_at_r1, ao_num2    &
              , ao_core_xc_mat_grid12e_full, ao_num2    &
              , 0.d0, int3b_core_tcxc_ao_grid12e, ao_num2 )
  end if

END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_core_xc_mat_grid12e_full, (ao_num, ao_num, n_points_final_grid) ]
 
  integer :: i1, i2, i2p_nuc, i2p_rad, i2p_ang
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

  ao_core_xc_mat_grid12e_full(:,:,:) = 0.d0

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
      do i2p_nuc = 1, nucl_num
        do i2p_rad = 1, n_points_extra_radial_grid - 1
          do i2p_ang = 1, n_points_extra_integration_angular

            r2p(1:3) = grid_points_extra_per_atom(1:3,i2p_ang,i2p_rad,i2p_nuc)
            distance = norm2(r2(1:3) - r2p(1:3))
            ! Compute matrix element only when there is no divergence 
            if (distance.gt.1.d-10) then
              w2p = final_weight_at_r_extra(i2p_ang,i2p_rad,i2p_nuc)

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
              kernel = kernel / distance

              do l = 1, ao_num
                ao_l_r2 = aos_in_r_array2(l,i2)
                do j = 1, ao_num
                  ! Get value of the j-th orbital at r2p
                  ao_j_r2p = aos_in_r_array_extra_full(j,i2p_ang,i2p_rad,i2p_nuc)
                  ! Notice the index order, because this array will be transposed later!
                  ao_core_xc_mat_grid12e_full(j,l,i1) += w2 * j_r1r2 * ao_l_r2 * w2p * kernel * ao_j_r2p * j_r1r2p
                enddo
              enddo
            endif

          enddo !r2p
        enddo   !r2p
      enddo     !r2p

    enddo !r2

  end do !r1

END_PROVIDER


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
  integer :: ao_num2

  core_tcxc_grid12ej_prun(:,:,:,:) = 0.d0

  if (core_tcxc_loops) then
    do i1 = 1, n_points_final_grid
      r1(1:3) = final_grid_points(1:3,i1)
      w1 = final_weight_at_r_vector(i1)

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

  else
    ! dgemm session
    ao_num2 = ao_num*ao_num
    call dgemm( 'N', 'T', ao_num2, ao_num2, n_points_final_grid, 1.d0 &
              , ao_overlap_grid1_at_r1, ao_num2    &
              , ao_core_xc_mat_grid12e_pruned, ao_num2    &
              , 0.d0, core_tcxc_grid12ej_prun, ao_num2 )
  end if

END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_core_xc_mat_grid12e_pruned, (ao_num, ao_num, n_points_final_grid) ]
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

  ao_core_xc_mat_grid12e_pruned(:,:,:) = 0.d0

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

          do l = 1, ao_num
            ! Get value of the l-th orbital at r2
            ao_l_r2 = aos_in_r_array2(l,i2)
            do j = 1, ao_num
              ! Get value of the j-th orbital at r2p
              ao_j_r2p = aos_in_r_array_extra(j,i2p)
              ! Notice the index order, because this array will be transposed later!
              ao_core_xc_mat_grid12e_pruned(j,l,i1) += w2 * j_r1r2 * ao_l_r2 * w2p * kernel * ao_j_r2p * j_r1r2p
            enddo
          enddo
        endif

      enddo !r2p
    enddo !r2

  end do !r1

END_PROVIDER
