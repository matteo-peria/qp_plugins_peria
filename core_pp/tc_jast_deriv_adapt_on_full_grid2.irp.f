! MAKE THIS PURE TO AVOID SIDE EFFECTS IN PARALLEL EXECUTION?
subroutine get_grad1_u12_on_full_grid2(i1, n_nuc, n_rad, n_ang, grid2, grad_vect, grad_sqrd)

  BEGIN_DOC
  ! For specific Jastrow factors defined as product of pair
  ! correlation functions of the type 
  !
  ! :math: u12 = u(r1-r2) = u(r1,r2), 
  !
  ! where r1 and r2 are the positions of particle 1 and 2 respectively,
  ! this subroutine computes the gradient of u with respect to r1
  ! 
  ! :math: \grad_1 u(r1,r2)
  !
  ! For several provided values of r2.
  !
  ! From a numerical point of view, the gradient computation is based on 
  ! two grids, grid1 and grid2.
  !
  ! i1 is the index of the specific point on grid1 we are interested in; 
  ! grid2 is the array of r2.
  ! 
  ! From a qp-implementation point of view:
  ! * i1:        integer between 1 and size(final_grid_points)
  ! * grid2:     dp array (grid_points_extra_per_atom, or a floating grid)
  ! * grad_vect: is the resulting array of gradients
  ! * grad_sqrd: is the resulting array of gradients' norms
  END_DOC

  implicit none
  integer, intent(in) :: i1
  integer, intent(in) :: n_nuc, n_rad, n_ang
  double precision, intent(in) :: grid2(3,n_ang,n_rad,n_nuc)
  double precision, intent(out) :: grad_vect(n_nuc*(n_rad-1)*n_ang,3)
  double precision, intent(out) :: grad_sqrd(n_nuc*(n_rad-1)*n_ang)

  double precision :: r1(3)

  PROVIDE j1e_type
  PROVIDE j2e_type
  PROVIDE env_type

  !!grad1_u12_vect_adapt_grid(n_pts_float_grid:n_grid2,:,:) = grad1_u12_num(:,:,:)
  !!grad1_u12_norm_adapt_grid(n_pts_float_grid:n_grid2,:,:) = grad1_u12_squared_num(:,:)

  ! Probably we do not need the following providers in the current simple version of the subroutine
  !PROVIDE final_grid_points
  !PROVIDE final_grid_points_extra

  r1(1:3) = final_grid_points(1:3,i1)

  if (     (j2e_type .eq. "Mu")       &
      .or. (j2e_type .eq. "Mur")      &
      .or. (j2e_type .eq. "Jpsi")     &
      .or. (j2e_type .eq. "Mugauss")  &
      .or. (j2e_type .eq. "Murgauss") &
      .or. (j2e_type .eq. "Bump")     &
      .or. (j2e_type .eq. "Boys")) then
    if(env_type .eq. "None") then
      call get_grad1_j12_on_full_grid2(r1, n_nuc, n_rad, n_ang, grid2, grad_vect)
    else
      error stop "ERROR: env_type treatment with adaptive grid not implemented yet"
    endif ! env_type
  else
      error stop "ERROR: j2e_type treatment with adaptive grid not implemented yet"
  end if

  !! THIS LOOP HAS BEEN REPLACED BY THE NORM2 INTRINSTIC
  !do i2 = 1, n_points_extra_final_grid
  !  res(i2) = resx(i2) * resx(i2) + resy(i2) * resy(i2) + resz(i2) * resz(i2)
  !enddo
  ! Substitute by 
  grad_sqrd = sum(grad_vect*grad_vect,dim=2)

end subroutine get_grad1_u12_on_full_grid2


subroutine get_grad1_j12_on_full_grid2(r1, n_nuc, n_rad, n_ang, grid2, grad_vect)
  use jastrow_module, only : grad1_j12
  BEGIN_DOC
  !
  !  d/dx1 j_2e(1,2)
  !  d/dy1 j_2e(1,2)
  !  d/dz1 j_2e(1,2)
  !
  END_DOC
  include 'constants.include.F'
  implicit none
  double precision, intent(in)  :: r1(3)
  integer, intent(in) :: n_nuc, n_rad, n_ang
  double precision, intent(in) :: grid2(3,n_ang,n_rad,n_nuc)
  ! OUTPUT
  double precision, intent(out) :: grad_vect(n_nuc*(n_rad-1)*n_ang,3)
  !
  integer :: i2, i_ang2, i_rad2, i_nuc2
  double precision :: r2(3)
  !double precision, external :: grad1_j12 !now it is imported from the module

  !print*, "IM IN get_grad1_j12_on_grid2"
  i2 = 1
  do i_nuc2 = 1, n_nuc
    do i_rad2 = 1, n_rad-1
      do i_ang2 = 1, n_ang
        r2(1:3) = grid2(1:3,i_ang2,i_rad2,i_nuc2)
        grad_vect(i2,1:3) = grad1_j12(r1,r2,mu_erf)
        i2 = i2 + 1
      end do
    end do
  end do

end subroutine get_grad1_j12_on_full_grid2
