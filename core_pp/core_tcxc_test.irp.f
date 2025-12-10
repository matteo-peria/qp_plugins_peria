program core_tcxc_test
  implicit none
  BEGIN_DOC
  ! Testing the core exchange potential integral in TC
  ! Using grid1, grid2 and the extra grid in its two versions (full, pruned)
  ! This means that the third integral is computed on a FIXED grid
  END_DOC

  ! Provide this stuff to avoid output later
  PROVIDE final_grid_points final_weight_at_r_vector
  PROVIDE final_grid_points2 final_weight_at_r_vector2
  PROVIDE ao_num mo_num n_core_pseudo_orb n_valence_pseudo_orb
  PROVIDE list_core_pseudo list_valence_pseudo list_core_pseudo_reverse list_valence_pseudo_reverse

  write(*,'(A)') repeat('=', 70)

  write(*,*) "1st grid: grid1,           size = ", n_points_final_grid
  write(*,*) "2nd grid: grid2,           size = ", n_points_final_grid2
  write(*,*) "3rd grid: adaptive"
  write(*,*) "      --> grid_extra_full, size = ", n_points_rad_extra_grid * n_points_ang_extra_grid
  write(*,*) "      --> grid_float,      size = ", n_points_ang_float_grid * n_points_rad_float_grid

  write(*,*) "Jastrow e^{±J}=1 imposed throught the EZFIO interface param"
  write(*,*) "core_tcxc_j0_testing = ", core_tcxc_j0_testing
  write(*,*) "(true is expected when testing)"
  write(*,*) "mu_erf = ", mu_erf
  write(*,*) "(value is ignored when Jastrow is 1 when testing)"

  ! Checking cross providers
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 0"
  call test_tcxc_grid12extra_VS_grid12adapt_j0

  ! Checking cross providers
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 1"
  call test_tcxc_grid12extra_VS_grid12adapt

end program core_tcxc_test


subroutine test_tcxc_grid12extra_VS_grid12adapt_j0
  implicit none
  double precision :: difference
  write(*,*) 
  write(*,*) "int3b_ao_overlap_grid1_w_corexc_grid2a  VS. int3b_ao_overlap_grid1_w_corexc_grid2e"
  difference = sum(abs(int3b_ao_overlap_grid1_w_corexc_grid2a(:,:,:,:) - int3b_ao_overlap_grid1_w_corexc_grid2e(:,:,:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(int3b_ao_overlap_grid1_w_corexc_grid2a)
end subroutine test_tcxc_grid12extra_VS_grid12adapt_j0



subroutine test_tcxc_grid12extra_VS_grid12adapt
  implicit none
  double precision :: difference
  write(*,*) 
  write(*,*) "int3b_core_tcxc_ao_grid12a_dgemm  VS. int3b_core_tcxc_ao_grid12e"
  difference = sum(abs(int3b_core_tcxc_ao_grid12a_dgemm(:,:,:,:) - int3b_core_tcxc_ao_grid12e(:,:,:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(int3b_core_tcxc_ao_grid12e)
end subroutine test_tcxc_grid12extra_VS_grid12adapt

