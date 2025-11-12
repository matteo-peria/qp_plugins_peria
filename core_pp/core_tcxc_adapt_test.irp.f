program core_tcxc_adapt_test

  BEGIN_DOC
  ! Testing the core exchange potential integral in TC
  END_DOC

  implicit none

  call test_tcxc_adapt_j0_grid123
  !call test_tcxc_j0_grid123

end program core_tcxc_adapt_test


subroutine test_tcxc_adapt_j0_grid123
  implicit none
  BEGIN_DOC
  ! Check that, when the Jastrow factor is equal to zero,
  ! the numerically-evaluated core exchange potential in a TC context
  ! is the same as a product of standard overlap and exchange integrals (non-TC)
  ! Expected value is also evaluated numerically
  END_DOC
  double precision :: difference

  write(*,*) 

  provide core_tcxc_adapt_grid123

  write(*,'(A)') repeat('=', 70)

  write(*,*) "TEST 1"
  write(*,*) "core_tcxc_adapt_grid123 VS core_tcxc_adapt_j0_grid123"

  write(*,*) "Grid 1: usual,      size = ", n_points_final_grid
  write(*,*) "Grid 2: grid2,      size = ", n_points_final_grid2
  write(*,*) "Grid 3: full extra, size = ", n_points_rad_extra_grid * n_points_ang_extra_grid
  write(*,*) "Grid 3: floating,   size = ", n_points_ang_float_grid * n_points_rad_float_grid


  write(*,*) "Expected sizes:"
  write(*,*) "ao_num              = ", ao_num
  write(*,*) "ao_num^4            = ", ao_num*ao_num*ao_num*ao_num
  write(*,*) "size(core_tcxc...)  = ", size(core_tcxc_adapt_grid123)
  write(*,*) 

  write(*,*) "Jastrow e^{±J}=1 imposed throught the EZFIO interface param"
  write(*,*) "core_tcxc_j0_testing = ", core_tcxc_j0_testing
  write(*,*) "(true is expected when testing)"
  write(*,*) "mu_erf = ", mu_erf
  write(*,*) "(value is ignored when Jastrow is 1 when testing)"

  write(*,*) 
  write(*,*) "... computing the difference between the providers"
  write(*,*) "CORE_TCXC_ADAPT_GRID123, CORE_TCXC_ADAPT_J0_GRID123"

  difference = sum(abs(core_tcxc_adapt_grid123(:,:,:,:) - core_tcxc_adapt_j0_grid123(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_adapt_j0_grid123)

!  ! THE EXACT INTEGRAL IS NOT READY YET
!  write(*,*) 
!  write(*,*) "... computing the difference between the providers"
!  write(*,*) "CORE_TCXC_ADAPT_GRID123, CORE_TCXC_J0_EXACT"
!
!  difference = sum(abs(core_tcxc_adapt_grid123(:,:,:,:) - core_tcxc_j0_exact(:,:,:,:)))
!
!  write(*,*) "Difference =           ", difference
!  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_adapt_j0_grid123)

end subroutine test_tcxc_adapt_j0_grid123
