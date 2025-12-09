BEGIN_PROVIDER [ double precision, core_xpot_grid2ej_prun, (ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! Numerical evaluation of
  ! < l | V_x^{core_pseudo} | j >
  ! = -\sum_{m\in\text{core}} <lm|1/r12|mk>
  ! = \int dr_1 \int dr' \phi_l(r_1) V_x(r_1,r') \phi_j(r')
  ! = - \int dr_1 \int dr' \phi_l(r_1) \sum_{m\in\text{core}} (\phi_m(r_1)\phi_j(r')/|r-r'|) \phi_j(r')
  !
  ! 1st integral is over grid2
  ! 2nd integral is over extra grid (FULL)
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

  integer :: i2p
    ! 3rd grid indices

  core_xpot_grid2ej_prun(:,:) = 0.d0

  do i2 = 1, n_points_final_grid2
    r2(1:3) = final_grid_points2(1:3,i2)
    w2 = final_weight_at_r_vector2(i2)

!    ! Loop over all AO (j-th index, variable is r2p)
!    ! ARE YOU SURE YOU WANT TO KEEP THIS LOOP EXTERNAL?
!    do j = 1, ao_num
!      ! Initialise 3rd nested integral
!      integral = 0.d0

      ! PRUNED FIXED GRID contributions 
      do i2p = 1, n_points_extra_final_grid
        r2p(1:3) = final_grid_points_extra(1:3,i2p)
        distance = norm2(r2(1:3) - r2p(1:3))
        ! Compute matrix element only when there is no divergence 
        if (distance.gt.1.d-10) then
          w2p = final_weight_at_r_vector_extra(i2p)
          !! Get value of the j-th orbital at r2p
          !ao_j_r2p = aos_in_r_array_extra(j,i2p)

          ! Initialize kernel
          kernel = 0.d0
          ! Loop over core orbitals to update the kernel
          do m = 1, n_core_pseudo_orb
            m_core = list_core_pseudo(m)
            kernel -= mos_in_r_array2_omp(m_core,i2) * mos_in_r_array_extra_omp(m_core,i2p)
          enddo

          !! Update of the integral at each new point r2p
          !integral += kernel * ao_j_r2p * w2p / distance 

          ! Loop over all AO (j-th index, variable is r2p)
          ! ARE YOU SURE YOU WANT TO KEEP THIS LOOP EXTERNAL?
          do j = 1, ao_num
            ! Retrieve j-th AO value at position r'
            ao_j_r2p = aos_in_r_array_extra(j,i2p)
            do l = 1, ao_num
              ! Retrieve l-th AO value at position r'
              ao_l_r2 = aos_in_r_array2(l,i2)
              !core_xpot_grid2ej_prun(l,j) += w2 * ao_l_r2 * integral
              core_xpot_grid2ej_prun(l,j) += w2 * ao_l_r2 * w2p * kernel * ao_j_r2p / distance 
            enddo
  
          end do  ! loop over j

        end if

      end do  ! loop oveer r2p

  end do  ! loop oveer r2

END_PROVIDER


BEGIN_PROVIDER [ double precision, core_tcxc_j0_grid12extra_pruned, (ao_num, ao_num, ao_num, ao_num)]
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

  ! Initialization
  core_tcxc_j0_grid12extra_pruned(:,:,:,:) = 0.d0

  ! do-loop version of the tensor product
  do j = 1, ao_num
    do l = 1, ao_num
      do k = 1, ao_num
        do i = 1, ao_num
          core_tcxc_j0_grid12extra_pruned(i,k,l,j) = ao_overlap_grid1(i,k) * core_xpot_grid2ej_prun(j,l)
        end do
      end do
    end do
  end do

  !! vectorized version of the tensor product
  !core_tcxc_j0_grid12extra_pruned = reshape( ao_overlap_grid1, shape=[ao_num,ao_num,1,1] ) &
  !                     * reshape( core_xpot_grid22, shape=[1,1,ao_num,ao_num] )
END_PROVIDER
