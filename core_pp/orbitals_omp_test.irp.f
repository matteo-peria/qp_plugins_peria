program orbitals_test

  implicit none
  BEGIN_DOC
  ! Testing orbitals function
  END_DOC

  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 1"
  call test_give_all_aos_at_r_new

  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 2"
  call test_give_all_aos_at_r_omp

  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 3"
  call speed_test_give_all_aos_at_r_omp

end program orbitals_test


subroutine test_give_all_aos_at_r_new
  use iso_fortran_env, only: out_unit => output_unit
  implicit none
  BEGIN_DOC
  ! Check that 'give_all_aos_at_r_new' gives the same results as the old subroutine
  END_DOC
  integer :: i
  double precision :: diff, diff_tot, r(3)
  double precision :: aos_old(ao_num)
  double precision :: aos_new(ao_num)

  diff = 0.d0
  diff_tot = 0.d0

  aos_old(:) = 0.d0
  aos_new(:) = 0.d0

  do i = 1, n_points_final_grid
    r(1:3) = final_grid_points(1:3,i)

    call give_all_aos_at_r(r, aos_old(:))
    call give_all_aos_at_r_new(r, aos_new(:))

    diff = sum(abs(aos_old - aos_new))
    diff_tot += diff
    !write(*,*) "Diff at point ", i, ": ", diff

  end do

  write(*,*) "Total difference: ", diff_tot

end subroutine test_give_all_aos_at_r_new


subroutine test_give_all_aos_at_r_omp
  use iso_fortran_env, only: out_unit => output_unit
  implicit none
  BEGIN_DOC
  ! Check that 'give_all_aos_at_r_omp' gives the same results as the old subroutine
  END_DOC
  integer :: i
  double precision :: diff, diff_tot, r(3)
  double precision :: aos_old(ao_num)
  double precision :: aos_omp(ao_num)

  integer :: count_start, count_end, count_rate
  real :: elapsed_old, elapsed_new
  
  elapsed_old = 0.d0
  elapsed_new = 0.d0

  diff = 0.d0
  diff_tot = 0.d0

  aos_old(:) = 0.d0
  aos_omp(:) = 0.d0

  do i = 1, n_points_final_grid
    r(1:3) = final_grid_points(1:3,i)

    call system_clock(count_rate = count_rate)
    call system_clock(count_start)
    call give_all_aos_at_r_omp(r, aos_omp(:))
    call system_clock(count_end)
    elapsed_new += real(count_end - count_start) / real(count_rate)


    call system_clock(count_rate = count_rate)
    call system_clock(count_start)
    call give_all_aos_at_r(r, aos_old(:))
    call system_clock(count_end)    
    elapsed_old += real(count_end - count_start) / real(count_rate)



    diff = sum(abs(aos_old - aos_omp))
    diff_tot += diff
    !write(*,*) "Diff at point ", i, ": ", diff

  end do

  write(*,*) "Total difference: ", diff_tot
  write(*,*) "Elapsed old: ", elapsed_old
  write(*,*) "Elapsed new: ", elapsed_new


end subroutine test_give_all_aos_at_r_omp


subroutine speed_test_give_all_aos_at_r_omp
  use iso_fortran_env, only: out_unit => output_unit
  implicit none
  BEGIN_DOC
  ! Check that 'give_all_aos_at_r_omp' gives the same results as the old subroutine
  END_DOC
  integer :: i
  double precision :: diff, diff_tot, r(3)
  double precision :: aos(ao_num)

  integer :: count_start, count_end, count_rate
  real :: elapsed_old, elapsed_new
  
  elapsed_old = 0.d0
  elapsed_new = 0.d0

  aos(:) = 0.d0
  call system_clock(count_rate = count_rate)
  call system_clock(count_start)
  do i = 1, n_points_final_grid
    r(1:3) = final_grid_points(1:3,i)
    call give_all_aos_at_r(r, aos(:))
  end do
  call system_clock(count_end)
  elapsed_old = real(count_end - count_start) / real(count_rate)


  aos(:) = 0.d0
  call system_clock(count_rate = count_rate)
  call system_clock(count_rate = count_rate)
  call system_clock(count_start)
  do i = 1, n_points_final_grid
    r(1:3) = final_grid_points(1:3,i)
    call give_all_aos_at_r_omp(r, aos(:))
  end do
  call system_clock(count_end)
  elapsed_new = real(count_end - count_start) / real(count_rate)

  write(*,*) "Elapsed old: ", elapsed_old
  write(*,*) "Elapsed new: ", elapsed_new


end subroutine speed_test_give_all_aos_at_r_omp
