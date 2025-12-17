!! This is just a testing program to check 
!! - if gradient wrt to r1 on the whole grid2 can be stored (small atoms)
!! - if more general providers contained in tc_jast_deriv_on_full_extra.irp.f,
!!
!!     grad1_u12_vect_on_full_extra_grid
!!     grad1_u12_sqrd_on_full_extra
!!     grad1_u12_vect_on_pruned_extra_grid
!!     grad1_u12_sqrd_on_prune_grid2,
!!
!!   which store gradients on the whole grid2 for different of Jastrow functions, 
!!   give the same result as other old providers, such as those contained in non_h_ints_mu/jast_deriv.irp.f, 
!!
!!     grad1_u12_num
!!     grad1_u12_squared_num,
!!
!!   which do the same, but for numerical Jastrows only
!
!program tc_compare_grad1_on_fixed_grid
!  BEGIN_DOC
!  ! Comparing PROVIDERS of Jastrow gradients
!  END_DOC
!  implicit none
!  call main
!
!end program tc_compare_grad1_on_fixed_grid
!
!
!subroutine main
!  use io_test_interface
!  implicit none
!
!  double precision :: difference
!
!  !provide grad1_u12_vect_on_full_extra_grid
!  !provide grad1_u12_sqrd_on_full_extra
!
!  print*, "PROVIDING grad1_u12_vect_on_pruned_extra_grid"
!  print*, "PROVIDING grad1_u12_sqrd_on_prune_grid2"
!
!
!  provide grad1_u12_vect_on_pruned_extra_grid
!  provide grad1_u12_sqrd_on_prune_grid2
!
!  print*, "grad1_u12_vect_on_pruned_extra_grid PROVIDED"
!  print*, "grad1_u12_sqrd_on_prune_grid2 PROVIDED" 
!
!
!  provide grad1_u12_num
!  provide grad1_u12_squared_num
!
!  !call compute_dp_array_diff( grad1_u12_sqrd_on_full_extra &
!  call compute_dp_array_diff( grad1_u12_sqrd_on_prune_grid2 &
!                            , grad1_u12_squared_num        &
!                            , show = .False.            &
!                            , message = "Difference : ", &
!                            , diff = difference            )
!
!end subroutine main
!
