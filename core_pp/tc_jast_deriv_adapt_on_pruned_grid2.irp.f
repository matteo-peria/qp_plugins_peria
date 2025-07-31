! TEST ONLY SUBROUTINES!
subroutine get_grad1_u12_on_pruned_grid2(i1, ngrid2, grid2, grad_vect, grad_sqrd)
  BEGIN_DOC
  ! Temporary test of GET_GRAD1_U12_ON_GRID2 on the pruned grid2
  ! (the original subroutine does not act on the 2d-pruned-grid2)
  END_DOC

  implicit none
  integer, intent(in) :: i1
  integer, intent(in) :: ngrid2
  double precision, intent(in) :: grid2(3,ngrid2)
  double precision, intent(out) :: grad_vect(ngrid2,3)
  double precision, intent(out) :: grad_sqrd(ngrid2)

  double precision :: r1(3)
  

  PROVIDE j1e_type
  PROVIDE j2e_type
  PROVIDE env_type
  
  r1(1:3) = final_grid_points(1:3,i1)

  if (     (j2e_type .eq. "Mu")       &
      .or. (j2e_type .eq. "Mur")      &
      .or. (j2e_type .eq. "Jpsi")     &
      .or. (j2e_type .eq. "Mugauss")  &
      .or. (j2e_type .eq. "Murgauss") &
      .or. (j2e_type .eq. "Bump")     &
      .or. (j2e_type .eq. "Boys")) then
    if(env_type .eq. "None") then
      call get_grad1_j12_on_pruned_grid2(r1, ngrid2, grid2, grad_vect)
    else
      error stop "ERROR: env_type treatment with adaptive grid not implemented yet"
    endif ! env_type
  else
      error stop "ERROR: j2e_type treatment with adaptive grid not implemented yet"
  end if

  grad_sqrd = sum(grad_vect*grad_vect,dim=2)

end subroutine get_grad1_u12_on_pruned_grid2


subroutine get_grad1_j12_on_pruned_grid2(r1, ngrid2, grid2, grad_vect)
  use jastrow_module, only : grad1_j12
  BEGIN_DOC
  ! Temporary test of GET_GRAD1_J12_ON_GRID2 on the normal grid2
  END_DOC
  include 'constants.include.F'
  implicit none
  double precision, intent(in)  :: r1(3)
  integer, intent(in) :: ngrid2
  double precision, intent(in) :: grid2(3,ngrid2)
  double precision, intent(out) :: grad_vect(ngrid2,3)
  !
  integer :: i2
  double precision :: r2(3)

  
  do i2 = 1, ngrid2
    r2(1:3) = grid2(1:3,i2)
    grad_vect(i2,1:3) = grad1_j12(r1,r2,mu_erf)
  end do

end subroutine get_grad1_j12_on_pruned_grid2
