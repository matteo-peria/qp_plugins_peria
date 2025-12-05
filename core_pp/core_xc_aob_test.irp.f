program core_xc_aob_test

  BEGIN_DOC
  ! Testing the core exchange potential in the AO basis instead of the MO one
  END_DOC

  implicit none

  ! Provide this stuff to avoid output later
  PROVIDE final_grid_points final_grid_points2 
  PROVIDE ao_num mo_num n_core_pseudo_orb n_valence_pseudo_orb
  PROVIDE list_core_pseudo list_valence_pseudo list_core_pseudo_reverse list_valence_pseudo_reverse

  write(*,*) "1st grid: grid1,           size = ", n_points_final_grid
  write(*,*) "2nd grid: grid2,           size = ", n_points_final_grid2
  write(*,*) "3rd grid: adaptive"
  write(*,*) "      --> grid_extra_full, size = ", n_points_rad_extra_grid * n_points_ang_extra_grid
  write(*,*) "      --> grid_float,      size = ", n_points_ang_float_grid * n_points_rad_float_grid

  ! TESTS

  ! The following should be the same as the next test when grid2 = grid extra
  !write(*,'(A)') repeat('=', 70)
  !write(*,*) "TEST 1"
  !call test_core_xc_ao_EXACT_vs_FIXED_1ee

  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 1"
  call test_core_xc_ao_EXACT_vs_FIXED

  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 2"
  call test_core_xc_ao_EXACT_vs_ADAPT

end program core_xc_aob_test


subroutine test_core_xc_ao_EXACT_vs_FIXED_1ee
  implicit none
  double precision :: difference
  integer :: i, j

  write(*,*) 
  write(*,*) "core_xpot_eej VS. core_xc_aob_exact "
  difference = sum(abs(core_xpot_eej(:,:) - core_xc_aob_exact(:,:) ))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_xc_aob_exact)

  do i = 1, ao_num
    do j = 1, ao_num
      write(*,*) core_xpot_eej(i,j), core_xc_aob_exact(i,j) 
    end do
  end do
end subroutine test_core_xc_ao_EXACT_vs_FIXED_1ee


subroutine test_core_xc_ao_EXACT_vs_FIXED
  implicit none
  double precision :: difference
  integer :: i, j
  write(*,*) 
  write(*,*) "core_xpot_grid2ej_full VS. core_xc_aob_exact "
  difference = sum(abs(core_xpot_grid2ej_full(:,:) - core_xc_aob_exact(:,:) ))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_xc_aob_exact)
  do i = 1, ao_num
    do j = 1, ao_num
      write(*,*) core_xpot_grid2ej_full(i,j), core_xc_aob_exact(i,j) 
    end do
  end do
end subroutine test_core_xc_ao_EXACT_vs_FIXED


subroutine test_core_xc_ao_EXACT_vs_ADAPT
  implicit none
  double precision :: difference
  integer :: i, j
  write(*,*) 
  write(*,*) "core_xpot_adapt_grid2aj VS. core_xc_aob_exact "
  difference = sum(abs(core_xpot_adapt_grid2aj(:,:) - core_xc_aob_exact(:,:) ))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_xc_aob_exact)
  do i = 1, ao_num
    do j = 1, ao_num
      write(*,*) core_xpot_adapt_grid2aj(i,j), core_xc_aob_exact(i,j) 
    end do
  end do
end subroutine test_core_xc_ao_EXACT_vs_ADAPT


subroutine test_core_xc_ao_FIXED_vs_ADAPT
  implicit none
  double precision :: difference
  write(*,*) 
  write(*,*) "core_xpot_adapt_grid2aj VS. core_xpot_grid2ej_full"
  difference = sum(abs(core_xpot_adapt_grid2aj(:,:) - core_xpot_grid2ej_full(:,:)))
  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(core_xpot_adapt_grid2aj)
end subroutine test_core_xc_ao_FIXED_vs_ADAPT


