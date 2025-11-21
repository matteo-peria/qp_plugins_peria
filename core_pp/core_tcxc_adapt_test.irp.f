program core_tcxc_adapt_test

  BEGIN_DOC
  ! Testing the core exchange potential integral in TC
  END_DOC

  implicit none

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

  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 1"
  call test_tcxc_j0_12a_12e

  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 2"
  call test_tcxc_12aj_VS_j0

  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 3"
  call test_tcxc_grid12aj_grid12ej_full

  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 4"
  call test_tcxc_adapt_grid12aj_OMP_vs_BLAS

  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 5"
  call test_tcxc_adapt_grid12aj_OMP_vs_grid12ej_full


end program core_tcxc_adapt_test

subroutine test_tcxc_j0_12a_12e
  implicit none
  double precision :: difference

  write(*,*) 
  write(*,*) "CORE_TCXC_ADAPT_J0_GRID12aj VS. CORE_TCXC_ADAPT_J0_GRID12ej "

  difference = sum(abs(core_tcxc_adapt_j0_grid12aj(:,:,:,:) - core_tcxc_j0_grid12ej_full(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_adapt_j0_grid12aj)

end subroutine test_tcxc_j0_12a_12e


subroutine test_tcxc_12aj_VS_j0
  implicit none
  double precision :: difference

  write(*,*) 
  write(*,*) "CORE_TCXC_ADAPT_GRID12aj VS. core_tcxc_adapt_j0_grid12aj"

  difference = sum(abs(core_tcxc_adapt_grid12aj(:,:,:,:) - core_tcxc_adapt_j0_grid12aj(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_adapt_j0_grid12aj)

end subroutine test_tcxc_12aj_VS_j0


subroutine test_tcxc_grid12aj_grid12ej_full
  implicit none
  double precision :: difference

  write(*,*) 
  write(*,*) "CORE_TCXC_ADAPT_GRID12aj VS. CORE_TCXC_GRID12ej_FULL"

  difference = sum(abs(core_tcxc_adapt_grid12aj(:,:,:,:) - core_tcxc_grid12ej_full(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_adapt_j0_grid12aj)

end subroutine test_tcxc_grid12aj_grid12ej_full


subroutine test_tcxc_adapt_grid12aj_OMP_vs_BLAS
  implicit none
  double precision :: difference

  write(*,*) 
  write(*,*) "CORE_TCXC_ADAPT_GRID12aj VS. CORE_TCXC_ADAPT_GRID12aj_testing"

  difference = sum(abs(core_tcxc_adapt_grid12aj(:,:,:,:) - core_tcxc_adapt_grid12aj_testing(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_adapt_grid12aj)

end subroutine test_tcxc_adapt_grid12aj_OMP_vs_BLAS


subroutine test_tcxc_adapt_grid12aj_OMP_vs_grid12ej_full
  implicit none
  double precision :: difference

  write(*,*) 
  write(*,*) "CORE_TCXC_GRID12ej_FULL VS. CORE_TCXC_ADAPT_GRID12aj_testing"

  difference = sum(abs(core_tcxc_grid12ej_full(:,:,:,:) - core_tcxc_adapt_grid12aj_testing(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_grid12ej_full)

end subroutine test_tcxc_adapt_grid12aj_OMP_vs_grid12ej_full


