module io_test_interface
  implicit none

  interface
    subroutine compute_dp_array_diff( arr1, arr2                          &
                                    , threshold, show, row_indx, col_indx &
                                    , message, diff)
      double precision, intent(in) :: arr1(:,:)
      double precision, intent(in) :: arr2(:,:)
      double precision, intent(in), optional :: threshold
      logical,  intent(in), optional :: show
      integer,  intent(in), optional :: row_indx(:)
      integer,  intent(in), optional :: col_indx(:)
      character(len=*), intent(in), optional :: message
      ! OUTPUT
      double precision, intent(out) :: diff
    end subroutine 
  end interface

end module io_test_interface
