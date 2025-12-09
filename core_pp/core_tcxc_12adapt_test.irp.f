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

  !----------------------------- SIMPLE OUTPUT ------------------------------!
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

  !---------------------------------- TESTS ----------------------------------!

  ! The following test gives an idea of the discretization error
  ! due to the AO OVERLAP integral contribution on the grid 1.
  ! This is useful to understand what should the minimum size of grid1 to
  ! exclude it has a too high influence on the overall error
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 0"
  call test_ao_overlap_EXACT_vs_NUMERIC

  ! The following test is useful to check that numerically everything is sound
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 1"
  call test_tcxc_j0_12Adgemm_NUMERICAL

  ! The following test is useful to check that numerically everything is sound
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 1-bis"
  call test_tcxc_j0_12Aloops

  ! The following test is useful to check that loops and dgemm versions 
  ! are the same independently from the fact that the Jastrow is equal to 1.0
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 2"
  call test_tcxc_12Aloops_VS_12Adgemm

  ! Check the testing providers
  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 3"
  call test_tcxc_j0_grid12adapt_FORCED_VS_EXACT

  ! These are THE TRUE tests
  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 4 dgemm-based"
  call test_tcxc_j0_grid12Adgemm_VS_EXACT

  write(*,*) 
  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 4 loops-based"
  call test_tcxc_j0_grid12Aloops_VS_EXACT


end program core_tcxc_adapt_test


subroutine test_tcxc_j0_12Adgemm_NUMERICAL
  implicit none
  double precision :: difference
  write(*,*) 
  write(*,*) "core_tcxc_grid12adapt_dgemm VS. core_tcxc_j0_grid12adapt"
  difference = sum(abs(core_tcxc_grid12adapt_dgemm(:,:,:,:) - core_tcxc_j0_grid12adapt(:,:,:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid12adapt)
end subroutine test_tcxc_j0_12Adgemm_NUMERICAL


subroutine test_tcxc_j0_12Aloops
  implicit none
  double precision :: difference
  write(*,*) 
  write(*,*) "core_tcxc_grid12adapt_loops VS. core_tcxc_j0_grid12adapt"
  difference = sum(abs(core_tcxc_grid12adapt_loops(:,:,:,:) - core_tcxc_j0_grid12adapt(:,:,:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid12adapt)
end subroutine test_tcxc_j0_12Aloops


subroutine test_tcxc_12Aloops_VS_12Adgemm
  implicit none
  double precision :: difference
  write(*,*) 
  write(*,*) "core_tcxc_grid12adapt_loops VS. core_tcxc_grid12adapt_dgemm"
  difference = sum(abs(core_tcxc_grid12adapt_loops(:,:,:,:) - core_tcxc_grid12adapt_dgemm(:,:,:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_grid12adapt_dgemm)
end subroutine test_tcxc_12Aloops_VS_12Adgemm


subroutine test_tcxc_j0_grid12Aloops_VS_EXACT
  implicit none
  double precision :: difference
  write(*,*) 
  write(*,*) "core_tcxc_grid12adapt_loops VS. core_tcxc_j0_exact"
  difference = sum(abs(core_tcxc_grid12adapt_loops(:,:,:,:) - core_tcxc_j0_exact(:,:,:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_exact)
end subroutine test_tcxc_j0_grid12Aloops_VS_EXACT


subroutine test_tcxc_j0_grid12Adgemm_VS_EXACT
  implicit none
  double precision :: difference
  write(*,*) 
  write(*,*) "core_tcxc_grid12adapt_dgemm VS. core_tcxc_j0_exact"
  difference = sum(abs(core_tcxc_grid12adapt_dgemm(:,:,:,:) - core_tcxc_j0_exact(:,:,:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_exact)
end subroutine test_tcxc_j0_grid12Adgemm_VS_EXACT


subroutine test_tcxc_j0_grid12adapt_FORCED_VS_EXACT
  implicit none
  double precision :: difference
  write(*,*) 
  write(*,*) "core_tcxc_j0_grid12adapt VS. core_tcxc_j0_exact"
  difference = sum(abs(core_tcxc_j0_grid12adapt(:,:,:,:) - core_tcxc_j0_exact(:,:,:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_exact)
end subroutine test_tcxc_j0_grid12adapt_FORCED_VS_EXACT


!! THIS IS CROSS CHECK
!subroutine test_tcxc_grid12aj_grid12ej_full
!  implicit none
!  double precision :: difference
!  write(*,*) 
!  write(*,*) "CORE_TCXC_ADAPT_GRID12aj VS. CORE_TCXC_GRID12ej_FULL"
!  difference = sum(abs(core_tcxc_grid12adapt_dgemm(:,:,:,:) - core_tcxc_grid12extra(:,:,:,:)))
!  write(*,*) "Difference =           ", difference
!  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_j0_grid12adapt)
!end subroutine test_tcxc_grid12aj_grid12ej_full

subroutine test_tcxc_adapt_grid12aj_OMP_vs_BLAS
  implicit none
  double precision :: difference

  write(*,*) 
  write(*,*) "CORE_TCXC_ADAPT_GRID12aj VS. CORE_TCXC_ADAPT_GRID12aj_loop"

  difference = sum(abs(core_tcxc_grid12adapt_dgemm(:,:,:,:) - core_tcxc_grid12adapt_loops(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_grid12adapt_dgemm)

  !integer :: i, j, k, l
  !do i = 1, ao_num
  !do j = 1, ao_num
  !do k = 1, ao_num
  !do l = 1, ao_num
  !difference = abs(core_tcxc_grid12adapt_dgemm(i,j,k,l) - core_tcxc_grid12adapt_loops(i,j,k,l))
  !if (difference > 1e-5) then
  !write(*,*) core_tcxc_grid12adapt_dgemm(i,j,k,l), core_tcxc_grid12adapt_loops(i,j,k,l), difference
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

  difference = sum(abs(core_tcxc_grid12extra(:,:,:,:) - core_tcxc_grid12adapt_loops(:,:,:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_grid12extra)

end subroutine test_tcxc_adapt_grid12aj_OMP_vs_grid12ej_full


subroutine test_tcxc_adapt_grid12aj_OMP_vs_j0
  implicit none
  double precision :: difference
  write(*,*) 
  write(*,*) "core_tcxc_j0_grid12adapt VS. CORE_TCXC_ADAPT_GRID12aj_loop"
  difference = sum(abs(core_tcxc_j0_grid12adapt (:,:,:,:) - core_tcxc_grid12adapt_loops(:,:,:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_tcxc_grid12adapt_loops)
end subroutine test_tcxc_adapt_grid12aj_OMP_vs_j0
