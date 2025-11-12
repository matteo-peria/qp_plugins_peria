! This file contains test-only providers, that is, providers that are just
! used to check that production providers behave well.
! These providers should be exact (analytics)

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


! THE FOLLOWING MIGHT NOT BE CORRECT. FIRST, CHECK INDICES ORDER
! SECOND, get_two_e_integrals uses MOs only, but here we need an exchange
! integral computed on mixed AO and MO!!!!
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
      ! This is wrong, how to mix AO and MO in 2e int?
      core_xpot_exact_ao(j,l) -= get_two_e_integral(m_core,j,l,m_core,mo_integrals_map) ! <mj|lm>
     enddo
   enddo
 enddo
END_PROVIDER
