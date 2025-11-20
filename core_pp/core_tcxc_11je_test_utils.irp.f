BEGIN_PROVIDER [ double precision, core_tcxc_j0_grid11je, (ao_num, ao_num, ao_num, ao_num)]
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
  ! 2nd integral is computed along usual grid
  ! 3rd integral is computed along extra grid
  END_DOC
  integer :: i,j,k,l

  ! Initialization
  core_tcxc_j0_grid11ej(:,:,:,:) = 0.d0

  ! do-loop version of the tensor product
  do j = 1, ao_num
    do l = 1, ao_num
      do k = 1, ao_num
        do i = 1, ao_num
          core_tcxc_j0_grid11je(i,k,l,j) = ao_overlap_grid1(i,k) * core_xpot_1je(j,l)
        end do
      end do
    end do
  end do

  !! vectorized version of the tensor product
  !core_tcxc_j0_grid11ej = reshape( ao_overlap_grid1, shape=[ao_num,ao_num,1,1] ) &
  !                     * reshape( core_xpot_1ej, shape=[1,1,ao_num,ao_num] )
END_PROVIDER


BEGIN_PROVIDER [ double precision, core_xpot_1je, (ao_num, ao_num)]
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
  ! 2nd integral is over usual grid
  ! 3rd integral is over extra grid
  END_DOC
  integer :: i2, i2p
    ! Indices of the points on grid 2 and 2' respectively
  double precision :: r2(3), w2
    ! Point (r_2) and weight (dr_2) on grid 2
  double precision :: r2p(3), w2p
    ! Point (r'_2) and weight (dr'_2) on nested grid 2'
  integer :: j , l
    ! AO indices
  double precision :: ao_l_r2, ao_j_r2p, v_x_core_ao_j_tc
    ! Real values of AOs $\chi_l(r_2)$ and $\chi_j(r'_2)$
  double precision :: distance
    ! Distance for Coulomb integral
  integer :: m, m_core
    ! Indices to loop over core MOs
  double precision :: kernel
    ! Kernel of the exchange potential
  core_xpot_1je(:,:) = 0.d0

  do i2 = 1, n_points_final_grid
    r2(1:3) = final_grid_points(1:3,i2)
    w2 = final_weight_at_r_vector(i2) 
    do j = 1, ao_num
      ! Compute nested integral in r2' (r1, i2 and AO index j-th as parameters)
      v_x_core_ao_j_tc = 0.d0
      ! Loop over the finer grid
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
          ! Initialize kernel
          kernel = 0.d0
          ! Loop over core orbitals to update the kernel
          do m = 1, n_core_pseudo_orb
            m_core = list_core_pseudo(m)
            kernel -= mos_in_r_array(m_core,i2) * mos_in_r_array_extra_omp(m_core,i2p)    
          enddo
          !kernel = kernel / distance
          ! Update of the integral at each new point r2p
          v_x_core_ao_j_tc += w2p * kernel * ao_j_r2p / distance
        endif
      enddo

      do l = 1, ao_num
        ao_l_r2 = aos_in_r_array(l,i2)
        core_xpot_1je(l,j) += w2 * ao_l_r2 * v_x_core_ao_j_tc
      enddo

    enddo
  enddo
END_PROVIDER
