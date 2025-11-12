! In other parts of the QP, 2-electrons integrals involving 4 indices
! are stored in 4-dimensional arrays with the following order:
!
! array:      j,l,k,i
! coordinate: 2,2,1,1 
!
! and the do-loops are nested accordingly, that is, in the opposite way
! (from the slowest running index to the fastest: i,k,l,j). 
! In our case, when we are in the loop of the index of particle 2, let's
! say 'j', we need to make a call to another computation-heavy function.
! Therefore, it is better to put that loop as the most external one,
! so that the call is made only once and used several times.
! It may be convenient to access the 4d storing array accordingly, and 
! put the order of the index as
!
! array:      i,k,l,j
! coordinate: 1,1,2,2
!
!  ?


! do a variant grid112

BEGIN_PROVIDER [ double precision, core_tcxc_grid122, (ao_num, ao_num, ao_num, ao_num)]
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
  double precision :: j_factor
    ! Jastrow factor between r1 and r2
  double precision, external :: j_mu_env, core_tcxc_kernel_grid122
    ! Temporary external definition, we should put them in a module an import that

  core_tcxc_grid122(:,:,:,:) = 0.d0

  do i1 = 1, n_points_final_grid
    r1(1:3) = final_grid_points(1:3,i1)
    w1 = final_weight_at_r_vector(i1)
    do i2 = 1, n_points_extra_final_grid
      r2(1:3) = final_grid_points_extra(1:3,i2)
      w2 = final_weight_at_r_vector_extra(i2) 
      ! Compute pair Jastrow factor between r1 and r2
      if (core_tcxc_j0_testing) then
        !write(*,*) "WARNING: for testing reason: j_mu_env = 0)"
        j_factor = 1.d0
      else
        j_factor = exp(-j_mu_env(r1,r2,mu_erf))
      end if
      do j = 1, ao_num
        ! Compute nested integral in r2' (r1, i2 and AO index j-th as parameters)
        v_x_core_ao_j_tc = core_tcxc_kernel_grid122(r1,i2,j) 
        do l = 1, ao_num
          ao_l_r2 = aos_in_r_array_extra(l,i2)
          do k = 1, ao_num
            ao_k_r1 = aos_in_r_array(k,i1)
            do i = 1, ao_num
              ao_i_r1 = aos_in_r_array(i,i1)
              core_tcxc_grid122(i,k,l,j) += w1 * ao_i_r1 * ao_k_r1  *  w2 * j_factor * ao_l_r2 * v_x_core_ao_j_tc
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
END_PROVIDER


! remember: pure functions are not allowed in irp.f files
!pure function core_tcxc_kernel_grid122(r1,i2,j) result(res)
function core_tcxc_kernel_grid122(r1,i2,j) result(res)
  BEGIN_DOC 
  ! Compute the integral
  ! \int d2' v_x^{\text{core}}(2,2') e^{u(1,2')} \chi_j(2)
  ! where 2' runs over the extra_grid
  END_DOC
  implicit none
  ! INPUT
  double precision, intent(in) :: r1
    ! Point of grid 1 at which we are computing the nested integral (parameter)
  integer, intent(in) :: i2
    ! Index of the point on grid 2
  integer, intent(in) :: j
    ! Index of the MO involved in the integral
  ! OUTPUT
  double precision :: res
  ! INTERNALS
  integer :: i2p 
    ! Nested grid 2' point index 
  double precision :: r2p(3), w2p
    ! Nested grid 2' point (r'_2) and weight (dr'_2)
  double precision :: r2(3)
    ! Grid 2' point (r_2), index i2 retrieved from passed argument i2
  double precision :: ao_j_r2p
    ! Real value of MO $\chi_j(r'_2)$, index retrieved from passed argument j
  double precision :: distance
    ! Distance for Coulomb integral
  integer :: m, m_core
    ! Indices to loop over core orbitals
  double precision :: kernel
    ! Kernel of the exchange potential
  double precision :: j_factor
    ! Jastrow factor between r1 and r'2
  double precision, external :: j_mu_env

  ! Retrieve r_2 from i2 index passed as argument:
  ! if r2 belongs to the normal grid
  !r2(1:3) = final_grid_points(1:3,i2)
  ! if r2 belongs to the finer grid
  r2(1:3) = final_grid_points_extra(1:3,i2)

  ! Initialize Jastrow-weighted nested core-valence exchange integral
  res = 0.d0
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
      ! Compute pair Jastrow factor between r1 and r'_2
      if (core_tcxc_j0_testing) then
        !write(*,*) "WARNING: for testing reason: j_mu_env = 0)"
        j_factor = 1.d0
      else
        j_factor = exp(j_mu_env(r1,r2p,mu_erf))
      end if
      ! Initialize kernel
      kernel = 0.d0
      ! Loop over core orbitals to update the kernel
      do m = 1, n_core_pseudo_orb
        m_core = list_core_pseudo(m)
        kernel -= mos_in_r_array_extra_omp(m_core,i2) * mos_in_r_array_extra_omp(m_core,i2p)    
      enddo
      ! Update of the integral at each new point r2p
      res += kernel * j_factor * ao_j_r2p * w2p / distance 
    endif
  enddo
end function


BEGIN_PROVIDER [ double precision, core_tcxc_grid112, (ao_num, ao_num, ao_num, ao_num)]
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
  double precision :: j_factor
    ! Jastrow factor between r1 and r2
  double precision, external :: j_mu_env, core_tcxc_kernel_grid112
    ! Temporary external definition, we should put them in a module an import that

  core_tcxc_grid122(:,:,:,:) = 0.d0

  do i1 = 1, n_points_final_grid
    r1(1:3) = final_grid_points(1:3,i1)
    w1 = final_weight_at_r_vector(i1)
    do i2 = 1, n_points_final_grid
      r2(1:3) = final_grid_points(1:3,i2)
      w2 = final_weight_at_r_vector(i2) 
      ! Compute pair Jastrow factor between r1 and r2
      if (core_tcxc_j0_testing) then
        !write(*,*) "WARNING: for testing reason: j_mu_env = 0)"
        j_factor = 1.d0
      else
        j_factor = exp(-j_mu_env(r1,r2,mu_erf))
      end if
      do j = 1, ao_num
        ! Compute nested integral in r2' (r1, i2 and AO index j-th as parameters)
        v_x_core_ao_j_tc = core_tcxc_kernel_grid112(r1,i2,j) 
        do l = 1, ao_num
          ao_l_r2 = aos_in_r_array(l,i2)
          do k = 1, ao_num
            ao_k_r1 = aos_in_r_array(k,i1)
            do i = 1, ao_num
              ao_i_r1 = aos_in_r_array(i,i1)
              core_tcxc_grid122(i,k,l,j) += w1 * ao_i_r1 * ao_k_r1 * w2 * j_factor * ao_l_r2 * v_x_core_ao_j_tc
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
END_PROVIDER

! remember: pure functions are not allowed in irp.f files
!pure function core_tcxc_kernel_grid112(r1,i2,j) result(res)
function core_tcxc_kernel_grid112(r1,i2,j) result(res)
  BEGIN_DOC 
  ! Compute the integral
  ! \int d2' v_x^{\text{core}}(2,2') e^{u(1,2')} \chi_j(2)
  ! where 2' runs over the extra_grid
  END_DOC
  implicit none
  ! INPUT
  double precision, intent(in) :: r1
    ! Point of grid 1 at which we are computing the nested integral (parameter)
  integer, intent(in) :: i2
    ! Index of the point on grid 1
  integer, intent(in) :: j
    ! Index of the MO involved in the integral
  ! OUTPUT
  double precision :: res
  ! INTERNALS
  integer :: i2p 
    ! Nested grid 2' point index 
  double precision :: r2p(3), w2p
    ! Nested grid 2' point (r'_2) and weight (dr'_2)
  double precision :: r2(3)
    ! Grid 2' point (r_2), index i2 retrieved from passed argument i2
  double precision :: ao_j_r2p
    ! Real value of MO $\chi_j(r'_2)$, index retrieved from passed argument j
  double precision :: distance
    ! Distance for Coulomb integral
  integer :: m, m_core
    ! Indices to loop over core orbitals
  double precision :: kernel
    ! Kernel of the exchange potential
  double precision :: j_factor
    ! Jastrow factor between r1 and r'2
  double precision, external :: j_mu_env

  ! Retrieve r_2 from i2 index passed as argument:
  r2(1:3) = final_grid_points(1:3,i2)
  
  ! Initialize Jastrow-weighted nested core-valence exchange integral
  res = 0.d0
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
      !ao_j_r2p = aos_in_r_array_extra(j,i2p)
      ao_j_r2p = aos_in_r_array_extra(j,i2p)
      ! Compute pair Jastrow factor between r1 and r'_2
      if (core_tcxc_j0_testing) then
        !write(*,*) "WARNING: for testing reason: j_mu_env = 0)"
        j_factor = 1.d0
      else
        j_factor = exp(j_mu_env(r1,r2p,mu_erf))
      end if
      ! Initialize kernel
      kernel = 0.d0
      ! Loop over core orbitals to update the kernel
      do m = 1, n_core_pseudo_orb
        m_core = list_core_pseudo(m)
        !! Update of the exchange kernel if r2 runs on the normal grid
        !!core_tcxc_kernel += aos_in_r_array(m_core,i2) * aos_in_r_array_extra(m_core,i2p)    
        !! Update of the exchange kernel if r2 runs on the finer grid
        !!kernel -= aos_in_r_array_extra(m_core,i2) * aos_in_r_array_extra(m_core,i2p)    
        ! MO variant
        kernel -= mos_in_r_array(m_core,i2) * mos_in_r_array_extra_omp(m_core,i2p)    
      enddo
      ! Update of the integral at each new point r2p
      res += kernel * j_factor * ao_j_r2p * w2p / distance 
    endif
  enddo
end function
