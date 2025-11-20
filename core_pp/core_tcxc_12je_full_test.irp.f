program core_tcxc_12je_full_test
  implicit none

  write(*,*) "1st grid: Grid 1,  size = ", n_points_final_grid
  write(*,*) "2nd grid: Grid 2,  size = ", n_points_final_grid2
  write(*,*) "3rd grid: Grid ef, size = ", nucl_num * n_points_extra_radial_grid * n_points_extra_integration_angular

  write(*,*) "Jastrow e^{±J}=1 imposed throught the EZFIO interface param"
  write(*,*) "core_tcxc_j0_testing = ", core_tcxc_j0_testing
  write(*,*) "(true is expected when testing)"
  write(*,*) "mu_erf = ", mu_erf
  write(*,*) "(value is ignored when Jastrow is 1 when testing)"


  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 1"
  call test_tcxc_j0_grid12je_full

end program core_tcxc_12je_full_test


subroutine test_tcxc_j0_grid12je_full
  implicit none
  BEGIN_DOC
  ! Check that, when the Jastrow factor is equal to zero,
  ! the numerically-evaluated core exchange potential in a TC context
  ! is the same as a product of standard overlap and exchange integrals (non-TC)
  END_DOC
  double precision :: difference

  write(*,*) "core_tcxc_grid12je_full VS core_tcxc_j0_grid12je_full"

  write(*,*) 
  write(*,*) "... computing the difference between the providers"
  write(*,*) "CORE_TCXC_GRID12je_full, CORE_TCXC_J0_GRID12je_full"

  difference = sum(abs(core_tcxc_grid12je_full(:,:,:,:) - core_tcxc_j0_grid12je_full(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid12je_full)

end subroutine test_tcxc_j0_grid12je_full
