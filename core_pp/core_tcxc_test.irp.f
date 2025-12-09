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
  write(*,*) "core_tcxc_j0_grid12adapt  VS. core_tcxc_j0_grid12extra"
  difference = sum(abs(core_tcxc_j0_grid12adapt(:,:,:,:) - core_tcxc_j0_grid12extra(:,:,:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid12adapt)
end subroutine test_tcxc_grid12extra_VS_grid12adapt_j0



subroutine test_tcxc_grid12extra_VS_grid12adapt
  implicit none
  double precision :: difference
  write(*,*) 
  write(*,*) "core_tcxc_grid12adapt_dgemm  VS. core_tcxc_grid12extra"
  difference = sum(abs(core_tcxc_grid12adapt_dgemm(:,:,:,:) - core_tcxc_grid12extra(:,:,:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_grid12extra)
end subroutine test_tcxc_grid12extra_VS_grid12adapt

