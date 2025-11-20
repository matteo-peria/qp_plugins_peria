program core_tcxc_12je_prun_test
  implicit none
  BEGIN_DOC
  ! Testing the core exchange potential integral in TC
  ! Using grid1, grid2 and the extra grid in its two version (full, pruned)
 END_DOC

  write(*,*) "1st grid: Grid 1,     size = ", n_points_final_grid
  write(*,*) "2nd grid: Grid 2,     size = ", n_points_final_grid2
  write(*,*) "3rd grid: Grid e prun, size = ", n_points_extra_final_grid

  write(*,*) "Jastrow e^{±J}=1 imposed throught the EZFIO interface param"
  write(*,*) "core_tcxc_j0_testing = ", core_tcxc_j0_testing
  write(*,*) "(true is expected when testing)"
  write(*,*) "mu_erf = ", mu_erf
  write(*,*) "(value is ignored when Jastrow is 1 when testing)"


  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 1"
  call test_tcxc_j0_grid12je_prun

end program core_tcxc_12je_prun_test


subroutine test_tcxc_j0_grid12je_prun
  implicit none
  double precision :: difference

  write(*,*) 
  write(*,*) "... computing the difference between the providers"
  write(*,*) "CORE_TCXC_GRID12je_prun, CORE_TCXC_J0_GRID12je_prun"

  difference = sum(abs(core_tcxc_grid12je_prun(:,:,:,:) - core_tcxc_j0_grid12je_prun(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid12je_prun)

end subroutine test_tcxc_j0_grid12je_prun
