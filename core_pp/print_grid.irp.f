program print_grid
  implicit none
  character(len=100) :: filepath
  logical :: has_arg


  has_arg = command_argument_count() >= 1

  if (has_arg) then
    call get_command_argument(1, filepath)
    call main(trim(filepath))
  else
    call main()
  end if

end program print_grid


subroutine main(path)
  use iso_fortran_env, only: output_unit
  implicit none
  character(len=*), intent(in), optional :: path
  character(len=19) :: filename
  character(len=100) :: filepath
  !
  integer :: ios
  integer :: unit_out
  logical :: file_exists
  integer :: i_ang, i_rad, i_nuc

  write(filename, '("rad", I4.4, "_ang", I4.4, ".out")') my_n_pt_r_grid, my_n_pt_a_grid

  if (present(path)) then
    filepath = trim(path)//'/'//filename
    open(newunit=unit_out, & 
       & file=filepath,    & 
       & status='new',     & 
       & action='write',   & 
       & iostat=ios        & 
       &)
    if (ios /= 0) then
      print *, "Error opening file: ", trim(filepath)
      stop 1
    end if
  else
    unit_out=output_unit
  end if

  do i_ang = 1, n_points_integration_angular
    do i_rad = 1, n_points_radial_grid
      do i_nuc = 1, nucl_num
        write(unit_out,*) grid_points_per_atom(1:3,i_ang,i_rad,i_nuc)
      end do
    end do
  end do

  if (present(path).and.ios==0) then
    close(unit_out)
  end if

end subroutine main
