program mos_core_in_r_test

  implicit none
  BEGIN_DOC
  ! Testing orbitals function
  END_DOC

  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 1: all MOs comparison"
  call test_all_mos_from_all_aos

  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 2: core MOs comparison"
  call test_core_mos_from_all_aos


end program mos_core_in_r_test

subroutine test_all_mos_from_all_aos
  implicit none
  BEGIN_DOC
  ! Compute all-MOs starting from already computed AOs
  ! Cons:
  ! - no tested. Are you sure of what you are doing?
  ! - matrix product involves bigger matrices (all MOs)
  ! Pro:
  ! + recycle aos_in_r, just called above
  END_DOC
  integer :: i
  double precision :: r(3)
  double precision :: aos_in_r(ao_num)
  double precision :: mos_in_r_new(mo_num)
  double precision :: diff, diff_tot

  diff_tot = 0.d0
  
  do i = 1, n_points_final_grid
    diff = 0.d0
    r(1:3) = final_grid_points(1:3,i)

    call give_all_aos_at_r(r, aos_in_r)
  
    !call dgemv( 'T', ao_num, mo_num, 1.d0, mo_coef, ao_num, &
    !            aos_in_r, 1, 0.d0, mos_in_r_new, 1)

    call dgemv( 'N', mo_num, ao_num, 1.d0, mo_coef_transp, mo_num, &
                aos_in_r, 1, 0.d0, mos_in_r_new, 1)


    diff = sum(abs(mos_in_r_new(:) - mos_in_r_array(:,i)))
    write(*,*) "Iteration: ", i, ". Difference = ", diff
    diff_tot += diff
  end do

  write(*,*) "Total difference: ", diff_tot
 
end subroutine test_all_mos_from_all_aos


subroutine test_core_mos_from_all_aos
  implicit none
  BEGIN_DOC
  ! Compute core-only-MOs starting from already computed AOs
  ! This should be the MOST CONVENIENT.
  ! Cons:
  ! - no tested. Are you sure of what you are doing?
  ! Pro:
  ! + matrix product involves smaller matrices (just core-MOs)
  ! + recycle aos_in_r, just called above
  END_DOC
  integer :: i, m, m_core
  double precision :: r(3)
  double precision :: aos_in_r(ao_num)
  double precision :: mos_core_in_r(mo_num)
  double precision :: diff, diff_tot

  diff_tot = 0.d0
  
  do i = 1, n_points_final_grid
    r(1:3) = final_grid_points(1:3,i)

    call give_all_aos_at_r(r, aos_in_r)

    !call dgemv( 'T', ao_num, n_core_pseudo_orb, 1.d0, mo_core_coef_notnorm, ao_num, &
    !            aos_in_r, 1, 0.d0, mos_core_in_r, 1)
    call dgemv( 'N', n_core_pseudo_orb, ao_num, 1.d0, mo_core_coef_notnorm_transp, n_core_pseudo_orb, &
                aos_in_r, 1, 0.d0, mos_core_in_r, 1)

    diff = 0.d0
    do m = 1, n_core_pseudo_orb
      m_core = list_core_pseudo(m)
      write(*,*) "m = ", m, "m_core = ", m_core
      diff += abs(mos_core_in_r(m) - mos_in_r_array(m_core,i))
    enddo
    !write(*,*) "i: ", i, "Difference = ", diff
    diff_tot += diff

  end do

  write(*,*) "Total difference: ", diff_tot
 
 end subroutine test_core_mos_from_all_aos


 
