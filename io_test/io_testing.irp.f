subroutine compute_dp_array_diff(arr1, arr2, threshold, show, row_indx, col_indx, message, diff)
  use io_color
  use io_const

  implicit none
  double precision, intent(in) :: arr1(:,:)
  double precision, intent(in) :: arr2(:,:)
  double precision, intent(in), optional :: threshold
  logical, intent(in), optional :: show
  integer, intent(in), optional :: row_indx(:)
  integer, intent(in), optional :: col_indx(:)
  character(len=*), intent(in), optional :: message
  ! OUTPUT
  double precision, intent(out) :: diff

  ! Manage optional arguments
  double precision     :: thr
  logical              :: show_array
  integer, allocatable :: r_indx(:)
  integer, allocatable :: c_indx(:)
  character(len=:), allocatable :: msg

  ! Other 
  double precision, allocatable :: arr(:,:)
  logical, allocatable :: mask(:,:)
  integer :: row, col

  ! Check arguments
  if (any(shape(arr1) /= shape(arr2))) then
    print *, "X Arrays have different shapes."
    return
  end if

  ! Manage optional arguments
  if (present(threshold)) then 
    thr = threshold  
  else
    thr = threshold_default
  end if

  if (present(show)) then 
    show_array = show
  else
    show_array = .True.
  end if

  if (present(row_indx)) then 
    allocate(r_indx,source=row_indx)
    r_indx = row_indx
  else
    r_indx = [(row, row=1,size(arr1,1))]
  end if

  if (present(col_indx)) then 
    allocate(c_indx,source=col_indx)
    c_indx = col_indx
  else
    c_indx = [(col, col=1,size(arr1,2))]
  end if

  if (present(message)) then
    allocate(character(len=len_trim(message)) :: msg)
    msg = trim(message)
  else
    msg = "Difference is: "
  end if

  diff = sum(abs(arr1(r_indx,c_indx) - arr2(r_indx,c_indx)))
  print*, msg, diff

  if (show_array) then
    mask = (abs(arr1(r_indx,c_indx) - arr2(r_indx,c_indx)) <= thr)
    arr = arr1(r_indx,c_indx)
    call print_array_color(arr, mask)
  end if
end subroutine compute_dp_array_diff






subroutine print_array_diff(arr1, arr2, threshold)
  use io_color
  use io_const
  implicit none

  double precision, intent(in) :: arr1(:,:), arr2(:,:)
  double precision, intent(in) :: threshold

  integer :: row, col
  character(len=:), allocatable :: esc_green, esc_red, esc_reset
  double precision :: val1, val2

  esc_green = escape_color(color%green)
  esc_red   = escape_color(color%red)
  esc_reset = escape_color(color%reset)

  if (any(shape(arr1) /= shape(arr2))) then
    print *, "❌ Arrays have different shapes."
    return
  end if

  do row = 1, size(arr1, 1)
    do col = 1, size(arr1, 2)
      val1 = arr1(row, col)
      val2 = arr2(row, col)

      if (abs(val1 - val2) <= threshold) then
        write(*,'(A,F12.7,A)', advance="no") esc_green, val1, esc_reset
      else
        write(*,'(A,F12.7,A)', advance="no") esc_red, val1, esc_reset
      end if
    end do
    print*  ! newline
  end do

  print*  ! extra space after matrix 
end subroutine print_array_diff


subroutine print_db_array(array,indx_row,indx_col)
  use io_color
  use io_const
  implicit none
  double precision, intent(in) :: array(:,:)
  integer, intent(in) :: indx_row(:)
  integer, intent(in) :: indx_col(:)
  !
  integer :: row, col
  integer :: val_indx

  do row = 1, size(indx_row)
    write(*,'(100F13.4)') (array(indx_row(row), indx_col(col)), col=1,size(indx_col))
  end do

end subroutine print_db_array
