 BEGIN_PROVIDER [ double precision, int3b_core_tcxc_ao_bch_terms_grid12a, (ao_num, ao_num, ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! Compute the BCH expansion total Nth-order correction of the integral
  ! <kl| e^{-J} V e^{+J} |ij>
  ! where V is the exchange potential with core-only electrons,
  ! and the integral is evaluated numerically.
  !
  ! Since the term can be expanded as
  ! e^{-J}Ve^{+J} = V + [V,J] + [[V,J],J] + [[V,J],J],J] + ...
  !
  ! the correction is computed as
  ! THIS_PROVIDER =  e^{-J}Ve^{+J} - V
  !
  ! Numerical evaluation of the integral is based on the following grids:
  ! int r_1  --> grid1
  ! int r_2  --> grid2
  ! int r_2' --> adaptive_grid
  END_DOC
  int3b_core_tcxc_ao_bch_terms_grid12a = int3b_core_tcxc_ao_grid12a_dgemm - int3b_ao_overlap_grid1_w_corexc_grid2a
END_PROVIDER


 BEGIN_PROVIDER [ double precision, int3b_core_tcxc_ao_bch_terms_grid12e, (ao_num, ao_num, ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! Compute the BCH expansion total Nth-order correction of the integral
  ! <kl| e^{-J} V e^{+J} |ij>
  ! where V is the exchange potential with core-only electrons,
  ! and the integral is evaluated numerically.
  !
  ! Since the term can be expanded as
  ! e^{-J}Ve^{+J} = V + [V,J] + [[V,J],J] + [[V,J],J],J] + ...
  !
  ! the correction is computed as
  ! THIS_PROVIDER =  e^{-J}Ve^{+J} - V
  !
  ! Numerical evaluation of the integral is based on the following grids:
  ! int r_1  --> grid1
  ! int r_2  --> grid2
  ! int r_2' --> extra_grid
  END_DOC

  int3b_core_tcxc_ao_bch_terms_grid12e = int3b_core_tcxc_ao_grid12e - int3b_ao_overlap_grid1_w_corexc_grid2e
END_PROVIDER
