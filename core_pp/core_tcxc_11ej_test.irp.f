program core_tcxc_11ej_test

  implicit none

  write(*,*) "1st grid: Grid 1,        size = ", n_points_final_grid
  write(*,*) "2nd grid: Grid 1,        size = ", n_points_final_grid
  write(*,*) "3rd grid: Grid e pruned, size = ", n_points_extra_final_grid

  write(*,*) "Jastrow e^{±J}=1 imposed throught the EZFIO interface param"
  write(*,*) "core_tcxc_j0_testing = ", core_tcxc_j0_testing
  write(*,*) "(true is expected when testing)"
  write(*,*) "mu_erf = ", mu_erf
  write(*,*) "(value is ignored when Jastrow is 1 when testing)"

  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 1 (grid1 = grid extra)"
  call test_tcxc_j0_grid11ej

end program core_tcxc_11ej_test

subroutine test_tcxc_j0_grid11ej
  implicit none
  double precision :: difference

  write(*,*) 
  write(*,*) "... computing the difference between the providers"
  write(*,*) "CORE_TCXC_GRID11ej, CORE_TCXC_J0_GRID11ej"

  difference = sum(abs(core_tcxc_grid11ej(:,:,:,:) - core_tcxc_j0_grid11ej(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid11ej)

end subroutine test_tcxc_j0_grid11ej

