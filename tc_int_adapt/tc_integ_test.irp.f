program tc_compare_2e_integ
  BEGIN_DOC
  ! Comparing nested TC integrals
  END_DOC
  implicit none
  character(len=100) :: filename
  logical :: has_arg

  write(*,*) 'filename', filename

  has_arg = command_argument_count() >= 1

  if (has_arg) then
    call get_command_argument(1, filename)
    call test_grids_tc_integ(trim(filename))
  else
    call test_grids_tc_integ()
  end if

end program tc_compare_2e_integ


subroutine test_grids_tc_integ(filename)
  use iso_fortran_env, only: output_unit
  use io_test_interface
  implicit none

  character(len=*), intent(in), optional :: filename
  character(len=100) :: output

  double precision :: diff_adapt_exact
  double precision :: diff_adapt_fixed
  double precision :: diff_fixed_exact

  integer :: ios
  integer :: unit_out

  integer :: i, l, j, k

  !print*,'EXACT'
  !call print_db_array(core_xpot_exact, list_valence_pseudo, list_valence_pseudo)
  
  write(*,'(A)') repeat('=', 70)

  provide ao_two_e_tc_tot 
  provide ao_two_e_tc_tot_numeric
  provide ao_two_e_tc_tot_adapt

  !do j = 1, ao_num
  !  do l = 1, ao_num
  !    do i = 1, ao_num
  !      do k = 1, ao_num
  !        if (dabs(ao_two_e_tc_tot(j,l,i,k) - ao_two_e_tc_tot_adapt(j,l,i,k))>1e-10) then 
  !          write(*,'(4I4,2F15.10)'), j, l, i, k, ao_two_e_tc_tot(j,l,i,k), ao_two_e_tc_tot_adapt(j,l,i,k)
  !        end if
  !      end do
  !    end do
  !  end do
  !end do

  diff_adapt_exact = sum(abs(ao_two_e_tc_tot - ao_two_e_tc_tot_adapt))
  diff_adapt_fixed = sum(abs(ao_two_e_tc_tot_numeric - ao_two_e_tc_tot_adapt))
  diff_fixed_exact = sum(abs(ao_two_e_tc_tot_numeric - ao_two_e_tc_tot))

  print*, "Difference semianalytic and adaptive grid", diff_adapt_exact
  print*, "Difference fixed grid and adaptive grid  ", diff_adapt_fixed
  print*, "Difference semianalytic and fixed grid   ", diff_fixed_exact

  write(*,'(A)') repeat('=', 70)

end subroutine
