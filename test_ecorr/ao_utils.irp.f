subroutine get_nonzero_indices_1d(array,size_array,nonzero,size_nonzero)
  ! Given a dp 1D ARRAY it returns the indices of the non-zero elements
  implicit none
  ! INPUT
  integer, intent(in) :: size_array
  integer, intent(in) :: size_nonzero
  double precision, dimension(size_array), intent(in) :: array
  ! OUPUT
  integer, dimension(size_nonzero), intent(out) :: nonzero
  !
  integer, dimension(:), allocatable :: indices
  integer :: i
  !
  indices = merge(0, [(i, i=1,size(array))], array==0.0d0)
  nonzero = pack(indices, indices /= 0)
  deallocate(indices)
end subroutine


subroutine get_equivalence_partition_1d_dp(array, size_array, tolerance, partition)
  ! Check that for a given set of double precision numbers (provided as an
  ! ARRAY sorted in either increasing or decreasing order) the 
  ! equivalence relationship is transitive within the given TOLERANCE.
  !
  ! Returns the PARTITION of the ARRAY: elements that are equivalent will have
  ! the same partitioning number. Last element of PARTITION will contain the
  ! total number of set-partitions according to the equivalence relation
  implicit none
  ! INPUT
  integer, intent(in) :: size_array
  double precision, dimension(size_array), intent(in) :: array
  double precision, intent(in) :: tolerance
  ! OUTPUT
  integer, dimension(size_array) :: partition
  !
  integer :: i, j, start, n_partition
  integer :: row, col
  integer :: equiv_cart_prod(size_array,size_array)
  double precision :: dist(size_array,size_array)

  ! Pairwise distance
  dist = 0.d0
  dist = abs(spread(array, dim=1, ncopies=size_array) &
          &- spread(array, dim=2, ncopies=size_array) )
   
  ! Check where distance is < tolerance and store it as set-Cartesian product
  equiv_cart_prod(:,:) = 0
  equiv_cart_prod = merge(1, 0, dist < tolerance)

  ! Fill the partition
  partition = 0
  n_partition = 1
  ! Set the first partition to 1
  partition(1) = n_partition

  start = 1
  row = 1
  do while (row < size_array)
    if (equiv_cart_prod(row,row+1)==0) then 
      ! Check that the equivalence is shared among all elements of the block
      if (product(equiv_cart_prod(start:row,start:row))==1) then 
        continue
      else
        print*, "ERROR"
        print*, "Equality relation between double precision numbers with"
        print*, "current tolerance threshold is not transitive"
        print*, "Input array: "
        print*, array
        print*, "Threshold: ", tolerance
        print*, "Distances: "
        do i = 1,size_array
          write(*,'(100F8.4)') (dist(i,j), j=1,size_array)
        end do
        print*, "Cartesian product of the relation: "
        do i = 1,size_array
          write(*,'(100I3)') (equiv_cart_prod(i,j), j=1,size_array)
        end do
        error stop 1
      end if
      ! Set all the equal numbers to belong to the same partition number
      n_partition = n_partition + 1
      ! Move the start of the next block to the next row
      start = row+1
    else if (row+1==size_array) then
      if (product(equiv_cart_prod(start:row+1,start:row+1))==1) then 
        continue
      else
        print*, "ERROR"
        print*, "Equality relation between double precision numbers with"
        print*, "current tolerance threshold is not transitive"
        print*, "Input array: "
        print*, array
        print*, "Threshold: ", tolerance
        do i = 1,size_array
          write(*,'(100F8.4)') (dist(i,j), j=1,size_array)
        end do
        print*, "Cartesian product of the relation: "
        do i = 1,size_array
          write(*,'(100I3)') (equiv_cart_prod(i,j), j=1,size_array)
        end do
        error stop 1
      end if
      ! Set all the equal numbers to belong to the same partition number
      !n_partition = n_partition + 1
      !partition(row+1) = n_partition
    else
      continue
    end if
    row = row+1
    partition(row) = n_partition
  end do
end subroutine


function are_dp_equivalent(a, b, tolerance)
  implicit none
  double precision, intent(in) :: a, b, tolerance
  logical :: are_dp_equivalent
  are_dp_equivalent = abs(a-b) < tolerance
end function 
