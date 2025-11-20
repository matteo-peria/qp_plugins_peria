program core_tcxc_12ej_test
  implicit none
  BEGIN_DOC
  ! Testing the core exchange potential integral in TC
  ! Using grid1, grid2 and the extra grid in its two version (full, pruned)
  END_DOC

  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 1"
  call test_tcxc_j0_grid12ej_full

  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 2"
  call test_tcxc_j0_grid12ej_prun

end program core_tcxc_12ej_test


subroutine test_tcxc_j0_grid12ej_full
  implicit none
  BEGIN_DOC
  ! Check that, when the Jastrow factor is equal to zero,
  ! the numerically-evaluated core exchange potential in a TC context
  ! is the same as a product of standard overlap and exchange integrals (non-TC)
  END_DOC
  double precision :: difference

  write(*,*) "core_tcxc_grid12ej_full VS core_tcxc_j0_grid12ej_full"

  write(*,*) "Grid 1: usual,      size = ", n_points_final_grid
  write(*,*) "Grid 2: grid2,      size = ", n_points_final_grid2
  write(*,*) "Grid 3: full extra, size = ", n_points_extra_integration_angular * n_points_extra_radial_grid * nucl_num

  write(*,*) "Expected sizes:"
  write(*,*) "ao_num              = ", ao_num
  write(*,*) "ao_num^4            = ", ao_num*ao_num*ao_num*ao_num
  write(*,*) "size(core_tcxc...)  = ", size(core_tcxc_grid12ej_full)
  write(*,*) 

  write(*,*) "Jastrow e^{±J}=1 imposed throught the EZFIO interface param"
  write(*,*) "core_tcxc_j0_testing = ", core_tcxc_j0_testing
  write(*,*) "(true is expected when testing)"
  write(*,*) "mu_erf = ", mu_erf
  write(*,*) "(value is ignored when Jastrow is 1 when testing)"

  write(*,*) 
  write(*,*) "... computing the difference between the providers"
  write(*,*) "CORE_TCXC_GRID12ej_full, CORE_TCXC_J0_GRID12ej_full"

  difference = sum(abs(core_tcxc_grid12ej_full(:,:,:,:) - core_tcxc_j0_grid12ej_full(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid12ej_full)

end subroutine test_tcxc_j0_grid12ej_full


subroutine test_tcxc_j0_grid12ej_prun
  implicit none
  BEGIN_DOC
  ! Check that, when the Jastrow factor is equal to zero,
  ! the numerically-evaluated core exchange potential in a TC context
  ! is the same as a product of standard overlap and exchange integrals (non-TC)
  END_DOC
  double precision :: difference

  write(*,*) "core_tcxc_grid12ej_prun VS core_tcxc_j0_grid12ej_prun"

  write(*,*) "Grid 1: usual,      size = ", n_points_final_grid
  write(*,*) "Grid 2: grid2,      size = ", n_points_final_grid2
  write(*,*) "Grid 3: prun extra, size = ", n_points_extra_final_grid

  write(*,*) "Expected sizes:"
  write(*,*) "ao_num              = ", ao_num
  write(*,*) "ao_num^4            = ", ao_num*ao_num*ao_num*ao_num
  write(*,*) "size(core_tcxc...)  = ", size(core_tcxc_grid12ej_prun)
  write(*,*) 

  write(*,*) "Jastrow e^{±J}=1 imposed throught the EZFIO interface param"
  write(*,*) "core_tcxc_j0_testing = ", core_tcxc_j0_testing
  write(*,*) "(true is expected when testing)"
  write(*,*) "mu_erf = ", mu_erf
  write(*,*) "(value is ignored when Jastrow is 1 when testing)"

  write(*,*) 
  write(*,*) "... computing the difference between the providers"
  write(*,*) "CORE_TCXC_GRID12ej_prun, CORE_TCXC_J0_GRID12ej_prun"

  difference = sum(abs(core_tcxc_grid12ej_prun(:,:,:,:) - core_tcxc_j0_grid12ej_prun(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid12ej_prun)

end subroutine test_tcxc_j0_grid12ej_prun
