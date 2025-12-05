program core_tcxc_adapt_test
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
  write(*,*) "3rd grid: adaptive"
  write(*,*) "      --> grid_extra_full, size = ", n_points_rad_extra_grid * n_points_ang_extra_grid
  write(*,*) "      --> grid_float,      size = ", n_points_ang_float_grid * n_points_rad_float_grid

  write(*,*) "Jastrow e^{±J}=1 imposed throught the EZFIO interface param"
  write(*,*) "core_tcxc_j0_testing = ", core_tcxc_j0_testing
  write(*,*) "(true is expected when testing)"
  write(*,*) "mu_erf = ", mu_erf
  write(*,*) "(value is ignored when Jastrow is 1 when testing)"


  !write(*,'(A)') repeat('=', 70)
  !write(*,*) "TEST 1"
  !call test_tcxc_j0_12a_12e

  ! The following test is useful to check that numerically everything is sound
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 1"
  call test_tcxc_12aj_VS_j0

  ! This is THE TRUE test
  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 3"
  call test_tcxc_j0_grid12aj_VS_EXACT

  !write(*,'(A)') repeat('=', 70)
  !write(*,*) "TEST 3"
  !call test_tcxc_grid12aj_grid12ej_full

  !write(*,'(A)') repeat('=', 70)
  !write(*,*) "TEST 4"
  !call test_tcxc_adapt_grid12aj_OMP_vs_BLAS

  !write(*,'(A)') repeat('=', 70)
  !write(*,*) "TEST 5"
  !call test_tcxc_adapt_grid12aj_OMP_vs_grid12ej_full

  !write(*,'(A)') repeat('=', 70)
  !write(*,*) "TEST 6"
  !call test_tcxc_adapt_grid12aj_OMP_vs_j0

end program core_tcxc_adapt_test


!subroutine test_tcxc_j0_12a_12e
!  implicit none
!  double precision :: difference
!  write(*,*) 
!  write(*,*) "CORE_TCXC_ADAPT_J0_GRID12aj VS. CORE_TCXC_ADAPT_J0_GRID12ej "
!  difference = sum(abs(core_tcxc_adapt_j0_grid12aj(:,:,:,:) - core_tcxc_j0_grid12ej_full(:,:,:,:)))
!  write(*,*) "Difference =           ", difference
!  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_adapt_j0_grid12aj)
!end subroutine test_tcxc_j0_12a_12e


subroutine test_tcxc_12aj_VS_j0
  implicit none
  double precision :: difference
  write(*,*) 
  write(*,*) "CORE_TCXC_ADAPT_GRID12aj VS. core_tcxc_adapt_j0_grid12aj"
  difference = sum(abs(core_tcxc_adapt_grid12aj(:,:,:,:) - core_tcxc_adapt_j0_grid12aj(:,:,:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_adapt_j0_grid12aj)
end subroutine test_tcxc_12aj_VS_j0

subroutine test_tcxc_j0_grid12aj_VS_EXACT
  implicit none
  double precision :: difference
  write(*,*) 
  write(*,*) "... computing the difference between the providers"
  write(*,*) "CORE_TCXC_GRID12aj, CORE_TCXC_J0_EXACT"
  difference = sum(abs(core_tcxc_adapt_grid12aj(:,:,:,:) - core_tcxc_j0_exact(:,:,:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_exact)
end subroutine test_tcxc_j0_grid12ej_full_VS_EXACT

!! THIS IS CROSS CHECK
!subroutine test_tcxc_grid12aj_grid12ej_full
!  implicit none
!  double precision :: difference
!  write(*,*) 
!  write(*,*) "CORE_TCXC_ADAPT_GRID12aj VS. CORE_TCXC_GRID12ej_FULL"
!  difference = sum(abs(core_tcxc_adapt_grid12aj(:,:,:,:) - core_tcxc_grid12ej_full(:,:,:,:)))
!  write(*,*) "Difference =           ", difference
!  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_adapt_j0_grid12aj)
!end subroutine test_tcxc_grid12aj_grid12ej_full

subroutine test_tcxc_adapt_grid12aj_OMP_vs_BLAS
  implicit none
  double precision :: difference

  write(*,*) 
  write(*,*) "CORE_TCXC_ADAPT_GRID12aj VS. CORE_TCXC_ADAPT_GRID12aj_loop"

  difference = sum(abs(core_tcxc_adapt_grid12aj(:,:,:,:) - core_tcxc_adapt_grid12aj_loop(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_adapt_grid12aj)

  !integer :: i, j, k, l
  !do i = 1, ao_num
  !do j = 1, ao_num
  !do k = 1, ao_num
  !do l = 1, ao_num
  !difference = abs(core_tcxc_adapt_grid12aj(i,j,k,l) - core_tcxc_adapt_grid12aj_loop(i,j,k,l))
  !if (difference > 1e-5) then
  !write(*,*) core_tcxc_adapt_grid12aj(i,j,k,l), core_tcxc_adapt_grid12aj_loop(i,j,k,l), difference
  !end if
  !end do
  !end do
  !end do
  !end do

end subroutine test_tcxc_adapt_grid12aj_OMP_vs_BLAS


subroutine test_tcxc_adapt_grid12aj_OMP_vs_grid12ej_full
  implicit none
  double precision :: difference

  write(*,*) 
  write(*,*) "CORE_TCXC_GRID12ej_FULL VS. CORE_TCXC_ADAPT_GRID12aj_loop"

  difference = sum(abs(core_tcxc_grid12ej_full(:,:,:,:) - core_tcxc_adapt_grid12aj_loop(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_grid12ej_full)

end subroutine test_tcxc_adapt_grid12aj_OMP_vs_grid12ej_full


subroutine test_tcxc_adapt_grid12aj_OMP_vs_j0
  implicit none
  double precision :: difference

  write(*,*) 
  write(*,*) "core_tcxc_adapt_j0_grid12aj VS. CORE_TCXC_ADAPT_GRID12aj_loop"

  difference = sum(abs(core_tcxc_adapt_j0_grid12aj (:,:,:,:) - core_tcxc_adapt_grid12aj_loop(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_adapt_grid12aj_loop)

end subroutine test_tcxc_adapt_grid12aj_OMP_vs_j0
