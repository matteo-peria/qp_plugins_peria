program core_tcxc_12ej_test
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
  write(*,*) "3rd grid: 1 of the following"
  write(*,*) "      --> grid_extra_full, size = ", n_points_rad_extra_grid * n_points_ang_extra_grid
  !write(*,*) "      --> grid_extra_prun, size = ", n_points_extra_final_grid

  write(*,*) "Jastrow e^{±J}=1 imposed throught the EZFIO interface param"
  write(*,*) "core_tcxc_j0_testing = ", core_tcxc_j0_testing
  write(*,*) "(true is expected when running a basic test witout jastrow factor)"
  write(*,*) "mu_erf = ", mu_erf
  write(*,*) "(value is ignored when Jastrow is 1 when testing)"
  
  write(*,'(A)') repeat('=', 70)

  !---------------------------------- TESTS ----------------------------------!
  ! The following test is useful to get an idea of the 
  ! discretization error coming from the AO OVERLAP integral 
  ! contribution on the grid 1
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 0"
  call test_ao_overlap_EXACT_vs_NUMERIC

  ! The following test is useful to check that numerically everything is sound
  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 1"
  call test_tcxc_j0_grid12ej_full

  ! The following is just a variation of the full extra grid
  !write(*,*) 
  !write(*,'(A)') repeat('=', 70)
  !write(*,*) "TEST 2"
  !call test_tcxc_j0_grid12ej_prun

  ! This is THE TRUE test
  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 3"
  call test_tcxc_j0_grid12ej_full_VS_EXACT

end program core_tcxc_12ej_test


subroutine test_tcxc_j0_grid12ej_full
  implicit none
  double precision :: difference
  write(*,*) 
  write(*,*) "... computing the difference between the providers"
  write(*,*) "CORE_TCXC_GRID12ej_full, CORE_TCXC_J0_GRID12ej_full"
  difference = sum(abs(core_tcxc_grid12ej_full(:,:,:,:) - core_tcxc_j0_grid12ej_full(:,:,:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid12ej_full)
end subroutine test_tcxc_j0_grid12ej_full


subroutine test_tcxc_j0_grid12ej_prun
  implicit none
  double precision :: difference
  !"TRANSCORRELATED CORE-EXCHANGE TRIPLE INTEGRAL WHEN JASTROW=1"
  write(*,*) "... computing the difference between the providers"
  write(*,*) "CORE_TCXC_GRID12ej_prun, CORE_TCXC_J0_GRID12ej_prun"
  difference = sum(abs(core_tcxc_grid12ej_prun(:,:,:,:) - core_tcxc_j0_grid12ej_prun(:,:,:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid12ej_prun)
end subroutine test_tcxc_j0_grid12ej_prun


subroutine test_tcxc_j0_grid12ej_full_VS_EXACT
  implicit none
  double precision :: difference
  write(*,*) 
  write(*,*) "... computing the difference between the providers"
  write(*,*) "CORE_TCXC_GRID12ej_full, CORE_TCXC_J0_EXACT"
  difference = sum(abs(core_tcxc_grid12ej_full(:,:,:,:) - core_tcxc_j0_exact(:,:,:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_exact)
end subroutine test_tcxc_j0_grid12ej_full_VS_EXACT
