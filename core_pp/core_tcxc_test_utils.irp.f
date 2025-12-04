! This file contains test-only providers, that is, providers that are just
! used to check that production providers behave well.
! These providers are numerical approximations

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
END_PROVIDER


BEGIN_PROVIDER [ double precision, core_density_matrix, (ao_num, ao_num)]
  BEGIN_DOC
  ! Compute the density matrix of the core orbitals
  END_DOC
  implicit none
  integer :: nu, mu, mo_c
  double precision :: c_mu, c_nu

  core_density_matrix(:,:) = 0.d0
  do nu = 1, ao_num
    do mu = 1, ao_num
      do mo_c = 1, n_core_pseudo_orb
        c_mu = mo_coef(mu, list_core_pseudo(mo_c))
        c_nu = mo_coef(nu, list_core_pseudo(mo_c))
        core_density_matrix(mu,nu) += c_mu * c_nu 
      end do
    end do
  end do
END_PROVIDER


BEGIN_PROVIDER [ double precision, core_xpot_exact_ao, (ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! Numerical evaluation of
  ! < l | V_x^{core_pseudo} | j >
  ! = -\sum_{m\in\text{core}} < l | \phi_m(1)\phi_m(2)/r12 | k >
  ! = \int dr_1 \int dr2 \chi_l(r_1) V_x(r_1,r_2) \chi_j(r_2)
  ! = - \int dr_1 \int dr_2 \chi_l(r_1) \sum_{m\in\text{core}} (\phi_m(r_1)\phi_j(r_2)/|r-r_2|) \chi_j(r_2)
  END_DOC
  integer :: j, l, nu, mu
  core_xpot_exact_ao(:,:) = 0.d0

  do j = 1, ao_num
    do l = 1, ao_num
      do nu = 1, n_core_pseudo_orb
        do mu = 1, n_core_pseudo_orb
          core_xpot_exact_ao(l,j) -= ao_two_e_integral(l,nu,mu,j) * core_density_matrix(mu,nu)
        enddo
      enddo
    enddo
  enddo

END_PROVIDER



