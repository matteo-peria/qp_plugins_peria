program test_ecorr_tc
  implicit none
  BEGIN_DOC
  ! Test of the plugin on correlation energy 
  END_DOC
  read_wf = .true.      ! Setting by default |true|: read the wave function from the |EZFIO| file 
  touch read_wf         ! Update all the subroutines using read_wf on its new value
  call routine_pouet_tc
end

subroutine routine_pouet_tc
  implicit none
  integer :: i
  double precision :: norm,h00
  call htilde_mu_mat_opt_bi_ortho_tot(ref_bitmask,ref_bitmask,N_int,h00)
  norm = 0.d0
  print*,'n_det_cisd = ',n_det_cisd
  do i = 1,n_det_cisd
!   print*,i
!   call debug_det(psi_det_cisd_proj(1,1,i),N_int)
!   print*,psi_coef_cisd_proj(i)
    norm += psi_coef_cisd_proj_tc_l(i)*psi_coef_cisd_proj_tc_r(i)
  enddo
  print*,'norm = ', norm 
  print*,'E corr from eigval_right  = ' ,  eigval_right_tc_bi_orth - h00 
  print*,'E corr from eigval_left   = ' ,  eigval_left_tc_bi_orth - h00 
  ! can be higher than one
  print*,'e_corr_cisd_proj_tc_r = ', e_corr_cisd_proj_tc_r
  print*,'e_corr_cisd_proj_tc_l = ', e_corr_cisd_proj_tc_l
  print*,'e_corr_cisd_proj_tc (core) = ', e_corr_cisd_proj_tc_r_cc
  print*,'e_corr_cisd_proj_tc (valence) = ', e_corr_cisd_proj_tc_r_vv
  print*,'e_corr_cisd_proj_tc (core and valence) = ', e_corr_cisd_proj_tc_r_cv

  print*, ""
  print*, "Correlation energy from S and D contributions decomposition:"
  write(*,'(A30,A15,A15)')     'Correlation from',                   'Value [LEFT]',             'Value [RIGHT]'
  write(*,'(A30,F15.7,F15.7)') 'Core only (S+D)',            e_corr_cisd_proj_tc_l_cc,   e_corr_cisd_proj_tc_r_cc
  write(*,'(A30,F15.7,F15.7)') 'Valence only (S+D)',         e_corr_cisd_proj_tc_l_vv,   e_corr_cisd_proj_tc_r_vv
  write(*,'(A30,F15.7,F15.7)') 'Core and valence (D)',       e_corr_cisd_proj_tc_l_cv,   e_corr_cisd_proj_tc_r_cv
  write(*,'(A30,F15.7,F15.7)') 'All core and/or val (CISD)', e_corr_cisd_proj_tc_l,      e_corr_cisd_proj_tc_r
  write(*,'(A30,F15.7,F15.7)') 'True total correlation',     eigval_left_tc_bi_orth-h00, eigval_right_tc_bi_orth-h00  

end subroutine
