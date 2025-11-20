program core_tcxc_3grid_cross_provider

  implicit none

  write(*,*) "1st grid: Grid 1,  size = ", n_points_final_grid
  write(*,*) "2nd grid: Grid 2,  size = ", n_points_final_grid2
  write(*,*) "3rd grid: Grid ep, size = ", n_points_extra_final_grid

  write(*,*) "NB: here we use only the pruned version of grid extra"

  write(*,*) "Jastrow e^{±J}=1 imposed throught the EZFIO interface param"
  write(*,*) "core_tcxc_j0_testing = ", core_tcxc_j0_testing
  write(*,*) "(true is expected when testing)"
  write(*,*) "mu_erf = ", mu_erf
  write(*,*) "(value is ignored when Jastrow is 1 when testing)"

  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 1 (grid1 = grid extra)"
  call test_tcxc_cross_provider_12ej_1eej

  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 2 (grid1 = grid extra)"
  call test_tcxc_cross_provider_12je_1eje

  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 3 (grid1 = grid2)"
  call test_tcxc_cross_provider_11je_12je

  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 4 (grid1 = grid2)"
  call test_tcxc_cross_provider_11ej_12ej_pruned

end program core_tcxc_3grid_cross_provider


subroutine test_tcxc_cross_provider_12ej_1eej
  implicit none
  double precision :: difference


  write(*,*) "Expect ZERO DIFFERENCE when usual grid2 = grid_extra"

  write(*,*) 
  write(*,*) "CORE_TCXC_GRID1eej VS. CORE_TCXC_GRID12ej_prun"

  difference = sum(abs(core_tcxc_grid1eej(:,:,:,:) - core_tcxc_grid12ej_prun(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_grid12ej_prun)


  write(*,*) 
  write(*,*) "CORE_TCXC_J0_GRID1eej VS. CORE_TCXC_J0_GRID12ej"

  difference = sum(abs(core_tcxc_j0_grid12ej_prun(:,:,:,:) - core_tcxc_j0_grid1eej(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid12ej_prun)

end subroutine test_tcxc_cross_provider_12ej_1eej


subroutine test_tcxc_cross_provider_12je_1eje
  implicit none
  double precision :: difference

  write(*,*) "Expect ZERO DIFFERENCE when usual grid2 = grid_extra"

  write(*,*) 
  write(*,*) "CORE_TCXC_GRID1eje VS. CORE_TCXC_GRID12je_prun"

  difference = sum(abs(core_tcxc_grid1eje(:,:,:,:) - core_tcxc_grid12je_prun(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_grid12je_prun)


  write(*,*) 
  write(*,*) "CORE_TCXC_J0_GRID1eje VS. CORE_TCXC_J0_GRID12ej_prun"

  difference = sum(abs(core_tcxc_j0_grid12je_prun(:,:,:,:) - core_tcxc_j0_grid1eje(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid12je_prun)

end subroutine test_tcxc_cross_provider_12je_1eje


subroutine test_tcxc_cross_provider_11je_12je
  implicit none
  double precision :: difference

  write(*,*) "Expect ZERO DIFFERENCE when usual grid1 = grid2"

  write(*,*) 
  write(*,*) "CORE_TCXC_GRID12je_prun VS. CORE_TCXC_GRID11je"

  difference = sum(abs(core_tcxc_grid12je_prun(:,:,:,:) - core_tcxc_grid11je(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_grid12je_prun)


  write(*,*) 
  write(*,*) "CORE_TCXC_J0_GRID12je_prun VS. CORE_TCXC_J0_GRID11je"

  difference = sum(abs(core_tcxc_j0_grid12je_prun(:,:,:,:) - core_tcxc_j0_grid11je(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid12je_prun)

end subroutine test_tcxc_cross_provider_11je_12je


subroutine test_tcxc_cross_provider_11ej_12ej_pruned
  implicit none
  double precision :: difference

  write(*,*) "Expect ZERO DIFFERENCE when usual grid2 = grid_extra"

  write(*,*) 
  write(*,*) "CORE_TCXC_GRID12ej_prun VS. CORE_TCXC_GRID11ej"

  difference = sum(abs(core_tcxc_grid12ej_prun(:,:,:,:) - core_tcxc_grid11ej(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_grid12ej_prun)


  write(*,*) 
  write(*,*) "CORE_TCXC_J0_GRID12ej_prun VS. CORE_TCXC_J0_GRID11ej"

  difference = sum(abs(core_tcxc_j0_grid12ej_prun(:,:,:,:) - core_tcxc_j0_grid11ej(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid12ej_prun)

end subroutine test_tcxc_cross_provider_11ej_12ej_pruned
