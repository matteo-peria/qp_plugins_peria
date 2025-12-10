 BEGIN_PROVIDER [ double precision, int3b_core_tcxc_ao_bch_terms_grid12a, (ao_num, ao_num, ao_num, ao_num)]
  implicit none
  BEGIN_DOC
  ! Compute 
  ! Used grid:
  ! int r_1 grid1
  ! int r_2 grid2
  ! int r_2' adaptive grid
  END_DOC
  int3b_core_tcxc_ao_bch_terms_grid12a = int3b_core_tcxc_ao_grid12a_dgemm - int3b_ao_overlap_grid1_w_corexc_grid2a
END_PROVIDER


 BEGIN_PROVIDER [ double precision, int3b_core_tcxc_ao_bch_terms_grid12e, (ao_num, ao_num, ao_num, ao_num)]
  implicit none
  int3b_core_tcxc_ao_bch_terms_grid12e = int3b_core_tcxc_ao_grid12e - int3b_ao_overlap_grid1_w_corexc_grid2e
END_PROVIDER
