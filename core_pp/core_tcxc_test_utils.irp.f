! This file contains test-only providers, that is, providers that are just
! used to check that production providers behave well


!!!!!!!!!!!!!!!!!!!!!!!! START OF NUMERICAL PROVIDERS !!!!!!!!!!!!!!!!!!!!!!!!
BEGIN_PROVIDER [ double precision, core_tcxc_j0_grid122, (ao_num, ao_num, ao_num, ao_num)]
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
  core_tcxc_j0_grid122(:,:,:,:) = 0.d0

  ! do-loop version of the tensor product
  ! WRONG ORDER?
  !do l = 1, ao_num
  !  do k = 1, ao_num
  !    do j = 1, ao_num
  !      do i = 1, ao_num
  !        core_tcxc_j0_grid122(i,j,k,l) = ao_overlap_grid1(i,k) * core_xpot_grid22(j,l)
  !      end do
  !    end do
  !  end do
  !end do

  do j = 1, ao_num
    do l = 1, ao_num
      do k = 1, ao_num
        do i = 1, ao_num
          core_tcxc_j0_grid122(i,k,l,j) = ao_overlap_grid1(i,k) * core_xpot_grid22(j,l)
        end do
      end do
    end do
  end do

  !! vectorized version of the tensor product
  !core_tcxc_j0_grid122 = reshape( ao_overlap_grid1, shape=[ao_num,ao_num,1,1] ) &
  !                     * reshape( core_xpot_grid22, shape=[1,1,ao_num,ao_num] )
END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_overlap_grid1, (ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! TEST-ONLY PROVIDER.
  ! Numerical evaluation of AOs overlap over usual grid
  ! $$ \int d1 \chi_k(1)^* \chi_i(1) 
  END_DOC
  integer :: i, k, i1
  double precision :: w1

  ! Initialization
  ao_overlap_grid1(:,:) = 0.d0

  ! do-loop solution
  do k = 1, ao_num
    do i = 1, ao_num
      do i1 = 1, n_points_final_grid
        w1 = final_weight_at_r_vector(i1)
        ao_overlap_grid1(i,k) += w1 * aos_in_r_array(i,i1) * aos_in_r_array(k,i1)
      end do
    end do
  end do

  !! linear algebra faster solution:
  !! Hadamard product of aos_in_r_array*aos_in_r_array 
  !! and matrix product with grid weights array
  !call get_AB_prod( aos_in_r_array*aos_in_r_array &
  !                , ao_num, n_points_final_grid   &
  !                , final_weight_at_r_vector, 1   &
  !                , ao_overlap_grid1)
END_PROVIDER


BEGIN_PROVIDER [ double precision, core_xpot_grid22, (ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! Numerical evaluation of
  ! < l | V_x^{core_pseudo} | j >
  ! = -\sum_{m\in\text{core}} <lm|1/r12|mk>
  ! = \int dr_1 \int dr' \phi_l(r_1) V_x(r_1,r') \phi_j(r')
  ! = - \int dr_1 \int dr' \phi_l(r_1) \sum_{m\in\text{core}} (\phi_m(r_1)\phi_j(r')/|r-r'|) \phi_j(r')
  !
  ! 1st integral is over extra grid
  ! 2nd integral is over extra grid
  END_DOC
  integer :: i2, i2p
    ! Indices of the points on grid 2 and 2' respectively
  double precision :: r2(3), w2
    ! Point (r_2) and weight (dr2) on grid 2
  double precision :: r2p(3), w2p
    ! Point (r'_2) and weight (dr2') on nested grid 2'
  integer :: j , l
    ! AO indices
  double precision :: ao_l_r2, ao_j_r2p
    ! Real values of AOs $\chi_l(r_2)$ and $\chi_j(r2')$
  double precision :: distance
    ! Distance for Coulomb integral
  integer :: m, m_core
    ! Indices to loop over core MOs
  double precision :: kernel
    ! Kernel of the exchange potential

  ! Initialization
  core_xpot_grid22(:,:) = 0.d0
  
  ! Loop over the finer grid
  do i2 = 1, n_points_extra_final_grid
    ! Find r_2 point from i2 loop-index
    r2(1:3) = final_grid_points_extra(1:3,i2)
    ! Find dr_2 weight from i2 loop-index
    w2 = final_weight_at_r_vector_extra(i2)
    do i2p = 1, n_points_extra_final_grid
      ! Find r'_2 point from i2p loop-index
      r2p(1:3) = final_grid_points_extra(1:3,i2p)
      distance = norm2(r2(1:3)-r2p(1:3))
      ! Compute matrix element only when there is no divergence 
      if(distance.gt.1.d-10) then
        ! Find r'_2 weight from i2p loop-index
        w2p = final_weight_at_r_vector_extra(i2p) 
        ! Initialize kernel
        kernel = 0.d0
        ! Loop over core orbitals to update the kernel
        do m = 1, n_core_pseudo_orb
          m_core = list_core_pseudo(m)
          kernel -= mos_in_r_array_extra_omp(m_core,i2p) * mos_in_r_array_extra_omp(m_core,i2)
        enddo
        do j = 1, ao_num
          ! Retrieve j-th AO value at position r'
          ao_j_r2p = aos_in_r_array_extra(j,i2p)
          do l = 1, ao_num
            ! Retrieve l-th AO value at position r
            ao_l_r2 = aos_in_r_array_extra(l,i2)
            ! Update
            core_xpot_grid22(l,j) += w2 * ao_l_r2 * w2p * kernel * ao_j_r2p / distance
          enddo
        enddo
      endif
    enddo
  enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, core_tcxc_j0_grid112, (ao_num, ao_num, ao_num, ao_num)]
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
  core_tcxc_j0_grid112(:,:,:,:) = 0.d0

  ! do-loop version of the tensor product
  ! WRONG ORDER?
  !do l = 1, ao_num
  !  do k = 1, ao_num
  !    do j = 1, ao_num
  !      do i = 1, ao_num
  !        core_tcxc_j0_grid112(i,j,k,l) = ao_overlap_grid1(i,k) * core_xpot_grid12(j,l)
  !      end do
  !    end do
  !  end do
  !end do

  do j = 1, ao_num
    do l = 1, ao_num
      do k = 1, ao_num
        do i = 1, ao_num
          core_tcxc_j0_grid122(i,k,l,j) = ao_overlap_grid1(i,k) * core_xpot_grid12(j,l)
        end do
      end do
    end do
  end do


  !! vectorized version of the tensor product
  !core_tcxc_j0_grid112 = reshape( ao_overlap_grid1, shape=[ao_num,ao_num,1,1] ) &
  !                     * reshape( core_xpot_grid12, shape=[1,1,ao_num,ao_num] )
END_PROVIDER


BEGIN_PROVIDER [ double precision, core_xpot_grid12, (ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! Numerical evaluation of
  ! < l | V_x^{core_pseudo} | j >
  ! = -\sum_{m\in\text{core}} <lm|1/r12|mk>
  ! = \int dr_1 \int dr' \phi_l(r_1) V_x(r_1,r') \phi_j(r')
  ! = - \int dr_1 \int dr' \phi_l(r_1) \sum_{m\in\text{core}} (\phi_m(r_1)\phi_j(r')/|r-r'|) \phi_j(r')
  !
  ! 1st integral is over usual grid
  ! 2nd integral is over extra grid
  END_DOC
  integer :: i2, i2p
    ! Indices of the points on grid 2 and 2' respectively
  double precision :: r2(3), w2
    ! Point (r_2) and weight (dr_2) on grid 2
  double precision :: r2p(3), w2p
    ! Point (r'_2) and weight (dr'_2) on nested grid 2'
  integer :: j , l
    ! AO indices
  double precision :: ao_l_r2, ao_j_r2p
    ! Real values of AOs $\chi_l(r_2)$ and $\chi_j(r'_2)$
  double precision :: distance
    ! Distance for Coulomb integral
  integer :: m, m_core
    ! Indices to loop over core MOs
  double precision :: kernel
    ! Kernel of the exchange potential

  ! Initialization
  core_xpot_grid12(:,:) = 0.d0
  
  ! Loop over the finer grid
  do i2 = 1, n_points_final_grid
    ! Find r_2 point from i2 loop-index
    r2(1:3) = final_grid_points(1:3,i2)
    ! Find dr_2 weight from i2 loop-index
    w2 = final_weight_at_r_vector(i2)
    do i2p = 1, n_points_extra_final_grid
      ! Find r'_2 point from i2p loop-index
      r2p(1:3) = final_grid_points_extra(1:3,i2p)
      distance = norm2(r2(1:3)-r2p(1:3))
      ! Compute matrix element only when there is no divergence 
      if(distance.gt.1.d-10) then
        ! Find r'_2 weight from i2p loop-index
        w2p = final_weight_at_r_vector_extra(i2p) 
        ! Initialize kernel
        kernel = 0.d0
        ! Loop over core orbitals to update the kernel
        do m = 1, n_core_pseudo_orb
          m_core = list_core_pseudo(m)
          kernel -= mos_in_r_array_extra_omp(m_core,i2p) * mos_in_r_array(m_core,i2)
        enddo
        do j = 1, ao_num
          ! Retrieve j-th AO value at position r'
          ao_j_r2p = aos_in_r_array_extra(j,i2p)
          do l = 1, ao_num
            ! Retrieve l-th AO value at position r'
            ao_l_r2 = aos_in_r_array(l,i2)
            ! Update
            core_xpot_grid12(l,j) += w2 * ao_l_r2 * w2p * kernel * ao_j_r2p / distance
          enddo
        enddo
      endif
    enddo
  enddo
END_PROVIDER



!!!!!!!!!!!!!!!!!!!!!!!! END OF NUMERICAL PROVIDERS !!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!! START OF EXACT PROVIDERS !!!!!!!!!!!!!!!!!!!!!!!!!

BEGIN_PROVIDER [ double precision, core_tcxc_j0_exact, (ao_num, ao_num, ao_num, ao_num)]
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
  ! No grid is used here because this integral can be computed analytically
  END_DOC
  integer :: i,j,k,l

  ! Init
  core_tcxc_j0_exact = 0.d0

  ! do-loop tensor product
  ! WRONG ORDER?
  !do l = 1, ao_num
  !  do k = 1, ao_num
  !    do j = 1, ao_num
  !      do i = 1, ao_num
  !        core_tcxc_j0_exact(i,j,k,l) = ao_overlap(i,k) * core_xpot_exact_ao(j,l)
  !      end do
  !    end do
  !  end do
  !end do

  do j = 1, ao_num
    do l = 1, ao_num
      do k = 1, ao_num
        do i = 1, ao_num
          core_tcxc_j0_exact(i,k,l,j) = ao_overlap(i,k) * core_xpot_exact_ao(j,l)
        end do
      end do
    end do
  end do

  !! vectorized tensor product
  !core_tcxc_jastrow_zero = reshape( ao_overlap, shape=[ao_num,ao_num,1,1] ) &
  !                       * reshape( core_xpot_exact, shape=[1,1,ao_num,ao_num] )
END_PROVIDER


BEGIN_PROVIDER [ double precision, core_xpot_exact_ao, (ao_num, ao_num)]
 implicit none
 BEGIN_DOC
 ! Analytical evaluation of
 ! < l | V_x^{core_pseudo} | j >
 END_DOC
 integer :: j,l
 integer :: m_core,m
 double precision :: get_two_e_integral

 core_xpot_exact_ao(:,:) = 0.d0
 do l = 1, ao_num
   do j = 1, ao_num
     do m = 1, n_core_pseudo_orb
      m_core = list_core_pseudo(m)
      core_xpot_exact_ao(j,l) -= get_two_e_integral(m_core,j,l,m_core,mo_integrals_map) ! <mj|lm>
     enddo
   enddo
 enddo
END_PROVIDER
