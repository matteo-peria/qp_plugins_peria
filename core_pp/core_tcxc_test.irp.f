program core_tcxc_test

  BEGIN_DOC
  ! Testing the core exchange potential integral in TC
  END_DOC

  implicit none

  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 0"
  call test_ao_overlap_exact_vs_numeric

  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 1"
  call test_tcxc_j0_grid122

  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 2"
  call test_tcxc_j0_grid112

end program core_tcxc_test


subroutine test_ao_overlap_exact_vs_numeric
  use iso_fortran_env, only: out_unit => output_unit
  implicit none
  BEGIN_DOC
  ! Check that the AO overlap matrix computed numerically is reasonable
  END_DOC
  double precision :: difference


  write(*,*) "ao_overlap_grid1 VS ao_overlap"
  
  write(*,*) "Depending on ao_num..."
  call write_int(out_unit, ao_num, 'ao_num')                   
  write(*,*) "... Expected sizes:"
  call write_int(out_unit, ao_num*ao_num   , 'ao_num^2')                   
  call write_int(out_unit, size(ao_overlap), 'size(ao_overlap)')                   

  write(*,*) 

  difference = sum(abs(ao_overlap_grid1(:,:) - ao_overlap(:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(ao_overlap)
  
  write(*,*) 

end subroutine test_ao_overlap_exact_vs_numeric


subroutine test_tcxc_j0_grid122
  implicit none
  BEGIN_DOC
  ! Check that, when the Jastrow factor is equal to zero,
  ! the numerically-evaluated core exchange potential in a TC context
  ! is the same as a product of standard overlap and exchange integrals (non-TC)
  ! Expected value is also evaluated numerically
  END_DOC
  double precision :: difference

  write(*,*) "core_tcxc_grid122 VS core_tcxc_j0_grid122"

  write(*,*) "Grid 1: usual, size = ", n_points_final_grid
  write(*,*) "Grid 2: extra, size = ", n_points_extra_final_grid
  write(*,*) "Grid 3: extra, size = ", n_points_extra_final_grid

  write(*,*) "Expected sizes:"
  write(*,*) "ao_num              = ", ao_num
  write(*,*) "ao_num^4            = ", ao_num*ao_num*ao_num*ao_num
  write(*,*) "size(core_tcxc...)  = ", size(core_tcxc_grid122)
  write(*,*) 

  write(*,*) "Jastrow e^{±J}=1 imposed throught the EZFIO interface param"
  write(*,*) "core_tcxc_j0_testing = ", core_tcxc_j0_testing
  write(*,*) "(true is expected when testing)"
  write(*,*) "mu_erf = ", mu_erf
  write(*,*) "(value is ignored when Jastrow is 1 when testing)"

  write(*,*) 
  write(*,*) "... computing the difference between the providers"
  write(*,*) "CORE_TCXC_GRID122, CORE_TCXC_J0_GRID122"

  difference = sum(abs(core_tcxc_grid122(:,:,:,:) - core_tcxc_j0_grid122(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid122)

! THE EXACT PROVIDER DOES NOT MAKE SENSE. NEED TO BE FIXED
!  write(*,*) 
!  write(*,*) "... computing the difference between the providers"
!  write(*,*) "CORE_TCXC_GRID122, CORE_TCXC_J0_EXACT"
!
!  difference = sum(abs(core_tcxc_grid122(:,:,:,:) - core_tcxc_j0_exact(:,:,:,:)))
!
!  write(*,*) "Difference =           ", difference
!  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid122)

end subroutine test_tcxc_j0_grid122



subroutine test_tcxc_j0_grid112
  implicit none
  BEGIN_DOC
  ! Check that, when the Jastrow factor is equal to zero,
  ! the numerically-evaluated core exchange potential in a TC context
  ! is the same as a product of standard overlap and exchange integrals (non-TC)
  ! Expected value is also evaluated numerically
  END_DOC
  double precision :: difference

  write(*,*) "core_tcxc_grid112 VS core_tcxc_j0_grid112"

  write(*,*) "Grid 1: usual, size = ", n_points_final_grid
  write(*,*) "Grid 2: usual, size = ", n_points_final_grid
  write(*,*) "Grid 3: extra, size = ", n_points_extra_final_grid

  write(*,*) "Expected sizes:"
  write(*,*) "ao_num              = ", ao_num
  write(*,*) "ao_num^4            = ", ao_num*ao_num*ao_num*ao_num
  write(*,*) "size(core_tcxc...)  = ", size(core_tcxc_grid112)
  write(*,*) 

  write(*,*) "Jastrow e^{±J}=1 imposed throught the EZFIO interface param"
  write(*,*) "core_tcxc_j0_testing = ", core_tcxc_j0_testing
  write(*,*) "(true is expected when testing)"
  write(*,*) "mu_erf = ", mu_erf
  write(*,*) "(value is ignored when Jastrow is 1 when testing)"

  write(*,*) 
  write(*,*) "... computing the difference between the providers"
  write(*,*) "CORE_TCXC_GRID112, CORE_TCXC_J0_GRID112"

  difference = sum(abs(core_tcxc_grid112(:,:,:,:) - core_tcxc_j0_grid112(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid112)

! THE EXACT PROVIDER DOES NOT MAKE SENSE. NEED TO BE FIXED
!  write(*,*) 
!  write(*,*) "... computing the difference between the providers"
!  write(*,*) "CORE_TCXC_GRID112, CORE_TCXC_J0_EXACT"
!
!  difference = sum(abs(core_tcxc_grid112(:,:,:,:) - core_tcxc_j0_exact(:,:,:,:)))
!
!  write(*,*) "Difference =           ", difference
!  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid122)
!
!  write(*,*) 

end subroutine test_tcxc_j0_grid112
