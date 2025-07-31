! This is just a testing program to check 
! - if gradient wrt to r1 on the whole grid2 can be stored
! - if a more general provider (GRAD1_U12_SQRD_ON_PRUNE_GRID2),
!   which computed the that gradient on the whole grid2 for 
!   different type of Jastrow function, gives the same result as
!   an old one (GRAD1_U12_SQUARED_NUM), which does the same but
!   only for numerical Jastrows


program tc_compare_grad1_on_fixed_grid
  BEGIN_DOC
  ! Comparing nested TC integrals
  END_DOC
  implicit none
  call main

end program tc_compare_grad1_on_fixed_grid


subroutine main
  use io_test_interface
  implicit none

  double precision :: difference

  !provide grad1_u12_vect_on_full_extra
  !provide grad1_u12_sqrd_on_full_extra

  print*, "PROVIDING grad1_u12_vect_on_prune_grid2"
  print*, "PROVIDING grad1_u12_sqrd_on_prune_grid2"


  provide grad1_u12_vect_on_prune_grid2
  provide grad1_u12_sqrd_on_prune_grid2

  print*, "grad1_u12_vect_on_prune_grid2 PROVIDED"
  print*, "grad1_u12_sqrd_on_prune_grid2 PROVIDED" 


  provide grad1_u12_num
  provide grad1_u12_squared_num

  !call compute_dp_array_diff( grad1_u12_sqrd_on_full_extra &
  call compute_dp_array_diff( grad1_u12_sqrd_on_prune_grid2 &
                            , grad1_u12_squared_num        &
                            , show = .False.            &
                            , diff = difference            )


  print*, "Difference : ", difference

end subroutine main

