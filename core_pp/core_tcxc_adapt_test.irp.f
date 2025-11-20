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
  call test_tcxc_adapt_j0_grid12aj

  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 2"
  call test_tcxc_adapt_j0_grid12ja

end program core_tcxc_adapt_test


subroutine test_tcxc_adapt_j0_grid12aj
  implicit none
  BEGIN_DOC
  ! Check that, when the Jastrow factor is equal to zero,
  ! the numerically-evaluated core exchange potential in a TC context
  ! is the same as a product of standard overlap and exchange integrals (non-TC)
  ! Expected value is also evaluated numerically
  END_DOC
  double precision :: difference

  write(*,*) 

  write(*,*) "core_tcxc_adapt_grid12aj VS core_tcxc_adapt_j0_grid12aj"

  write(*,*) 
  write(*,*) "... computing the difference between the providers"
  write(*,*) "CORE_TCXC_ADAPT_GRID12aj, CORE_TCXC_ADAPT_J0_GRID12aj"

  difference = sum(abs(core_tcxc_adapt_grid12aj(:,:,:,:) - core_tcxc_adapt_j0_grid12aj(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_adapt_j0_grid12aj)

  !integer :: i, j, k, l
  !do i = 1, ao_num
  !do j = 1, ao_num
  !do k = 1, ao_num
  !do l = 1, ao_num
  !if ((core_tcxc_adapt_grid12a(l,k,j,i) /= 0.d0).or.( core_tcxc_adapt_j0_grid12aj(l,k,j,i) /= 0.d0)) then
  !write(*,*) l,k,j,i, core_tcxc_adapt_grid12a(l,k,j,i), core_tcxc_adapt_j0_grid12aj(l,k,j,i)
  !end if
  !end do
  !end do
  !end do
  !end do

end subroutine test_tcxc_adapt_j0_grid12aj


subroutine test_tcxc_adapt_j0_grid12ja
  implicit none
  BEGIN_DOC
  ! Check that, when the Jastrow factor is equal to zero,
  ! the numerically-evaluated core exchange potential in a TC context
  ! is the same as a product of standard overlap and exchange integrals (non-TC)
  ! Expected value is also evaluated numerically
  END_DOC
  double precision :: difference

  write(*,*) 

  write(*,*) "core_tcxc_adapt_grid12ja VS core_tcxc_adapt_j0_grid12ja"

  write(*,*) 
  write(*,*) "... computing the difference between the providers"
  write(*,*) "CORE_TCXC_ADAPT_GRID12ja, CORE_TCXC_ADAPT_J0_GRID12ja"

  difference = sum(abs(core_tcxc_adapt_grid12ja(:,:,:,:) - core_tcxc_adapt_j0_grid12ja(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_adapt_j0_grid12ja)

end subroutine test_tcxc_adapt_j0_grid12ja
