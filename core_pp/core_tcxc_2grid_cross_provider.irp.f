program core_tcxc_2grid_cross_provider

  implicit none

  write(*,*) "Grid 1,     size = ", n_points_final_grid
  write(*,*) "Grid extra, size = ", n_points_extra_final_grid

  write(*,*) "NB: here we use only the pruned version of grid extra"

  write(*,*) "Jastrow e^{±J}=1 imposed throught the EZFIO interface param"
  write(*,*) "core_tcxc_j0_testing = ", core_tcxc_j0_testing
  write(*,*) "(true is expected when testing)"
  write(*,*) "mu_erf = ", mu_erf
  write(*,*) "(value is ignored when Jastrow is 1 when testing)"

  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 1"
  call test_tcxc_cross_provider_11ej_1eej

  !write(*,*) 
  !write(*,'(A)') repeat('=', 70)
  !write(*,*) "TEST 2"
  !call test_tcxc_cross_provider_11je_1eje

end program core_tcxc_2grid_cross_provider


subroutine test_tcxc_cross_provider_11ej_1eej
  implicit none
  BEGIN_DOC
  ! FURTHER TESTING CROSS PROVIDER
  END_DOC
  double precision :: difference


  write(*,*) "Expect ZERO DIFFERENCE when usual grid1 = grid_extra"

  write(*,*) 
  write(*,*) "... computing the difference between the providers"
  write(*,*) "CORE_TCXC_GRID1ee, CORE_TCXC_GRID11e"

  difference = sum(abs(core_tcxc_grid1eej(:,:,:,:) - core_tcxc_grid11ej(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_grid11ej)


  write(*,*) 
  write(*,*) "... computing the difference between the providers"
  write(*,*) "CORE_TCXC_J0_GRID1ee, CORE_TCXC_J0_GRID11e"

  difference = sum(abs(core_tcxc_j0_grid11ej(:,:,:,:) - core_tcxc_j0_grid1eej(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid11ej)

end subroutine test_tcxc_cross_provider_11ej_1eej


!subroutine test_tcxc_cross_provider_11je_1eje
!  implicit none
!  BEGIN_DOC
!  ! FURTHER TESTING CROSS PROVIDER
!  END_DOC
!  double precision :: difference
!
!
!  write(*,*) "Expect ZERO DIFFERENCE when usual grid1 = grid_extra"
!
!  write(*,*) 
!  write(*,*) "... computing the difference between the providers"
!  write(*,*) "CORE_TCXC_GRID1ee, CORE_TCXC_GRID11e"
!
!  difference = sum(abs(core_tcxc_grid1eje(:,:,:,:) - core_tcxc_grid11je(:,:,:,:)))
!
!  write(*,*) "Difference =           ", difference
!  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_grid11je)
!
!
!  write(*,*) 
!  write(*,*) "... computing the difference between the providers"
!  write(*,*) "CORE_TCXC_J0_GRID1ee, CORE_TCXC_J0_GRID11e"
!
!  difference = sum(abs(core_tcxc_j0_grid11je(:,:,:,:) - core_tcxc_j0_grid1eje(:,:,:,:)))
!
!  write(*,*) "Difference =           ", difference
!  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid11je)
!
!end subroutine test_tcxc_cross_provider_11je_1eje

