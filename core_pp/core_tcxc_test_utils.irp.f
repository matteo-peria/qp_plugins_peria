! This file contains test-only providers, that is, providers that are just
! used to check that production providers behave well.
! These providers are numerical approximations

BEGIN_PROVIDER [ double precision, ao_overlap_grid1, (ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! TEST-ONLY PROVIDER.
  ! Numerical evaluation of AOs overlap over usual grid1
  ! $$ \int d1 \chi_k(1)^* \chi_i(1) 
  END_DOC
  integer :: i, k, i1
  double precision :: w1
  ao_overlap_grid1(:,:) = 0.d0
  do k = 1, ao_num
    do i = 1, ao_num
      do i1 = 1, n_points_final_grid
        w1 = final_weight_at_r_vector(i1)
        ao_overlap_grid1(i,k) += w1 * aos_in_r_array(i,i1) * aos_in_r_array(k,i1)
      end do
    end do
  end do
END_PROVIDER


subroutine test_ao_overlap_EXACT_vs_NUMERIC
  use iso_fortran_env, only: out_unit => output_unit
  implicit none
  BEGIN_DOC
  ! Check that the AO overlap matrix computed numerically is reasonable
  ! The following test is useful to get an idea of the 
  ! discretization error coming from the AO OVERLAP integral 
  ! contribution on the grid 1
  END_DOC
  double precision :: difference

  write(*,*) "ao_overlap_grid1 VS ao_overlap"
  write(*,*) "Depending on ao_num..."
  call write_int(out_unit, ao_num, 'ao_num')                   
  write(*,*) "... Expected sizes:"
  call write_int(out_unit, ao_num*ao_num   , 'ao_num^2')                   
  call write_int(out_unit, size(ao_overlap), 'size(ao_overlap)')                   
  write(*,*) 
  difference = sum(abs(ao_overlap_grid1(:,:) - ao_overlap(:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(ao_overlap)
  write(*,*) 
end subroutine test_ao_overlap_EXACT_vs_NUMERIC



BEGIN_PROVIDER [ double precision, core_density_matrix, (ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! Compute the density matrix of the core orbitals in AO basis.
  ! This will be used to compare the core exchange potential evaluated in the 
  ! AO basis numerically with the exact one (see core_xc_aob_exact)
  END_DOC
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


BEGIN_PROVIDER [ double precision, core_xc_aob_exact, (ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! Numerical evaluation of
  ! < l | V_x^{core_pseudo} | j >
  ! = -\sum_{m\in\text{core}} < l | \phi_m(1)\phi_m(2)/r12 | k >
  ! = \int dr_1 \int dr2 \chi_l(r_1) V_x(r_1,r_2) \chi_j(r_2)
  ! = - \int dr_1 \int dr_2 \chi_l(r_1) \sum_{m\in\text{core}} (\phi_m(r_1)\phi_j(r_2)/|r-r_2|) \chi_j(r_2)
  END_DOC
  integer :: j, l, nu, mu
  double precision, external :: ao_two_e_integral
  core_xc_aob_exact(:,:) = 0.d0
  do j = 1, ao_num
    do mu = 1, n_core_pseudo_orb
      do nu = 1, n_core_pseudo_orb
        do l = 1, ao_num
          core_xc_aob_exact(l,j) -= ao_two_e_integral(l,nu,mu,j) * core_density_matrix(nu,mu)
        enddo
      enddo
    enddo
  enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, int3b_ao_overlap_w_corexc_exact, (ao_num, ao_num, ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! TEST-ONLY PROVIDER.
  ! Check that the TC exchange integral make sense when the Jastrow factor 
  ! is set to 1 (that is, when the Jastrow exponent is 0).

  ! $$ \braket{kl | e^{-J} V_{x,\text{core}} e^{+J} | ij} |_{J=0}
  !      \int d1 \chi_k(1)^* \chi_i(1) 
  !        \int d2 \chi_l^*(2) \int d2' v_x^{\text{core}}(2,2') \chi_j(2')
  !          \braket{k|i} \braket{l| V_{x,\text{core}} |j} $$
  END_DOC
  integer :: i,j,k,l
  ! Initialization
  int3b_ao_overlap_w_corexc_exact(:,:,:,:) = 0.d0
  do j = 1, ao_num
    do l = 1, ao_num
      do k = 1, ao_num
        do i = 1, ao_num
          int3b_ao_overlap_w_corexc_exact(i,k,l,j) = ao_overlap(i,k) * core_xc_aob_exact(j,l)
        end do
      end do
    end do
  end do
END_PROVIDER
