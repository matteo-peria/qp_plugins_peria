! These providers are inspired by those in 
! ${QPROPT}/non_h_ints_mu/jast_deriv.irp.f
!
!   grad1_u12_num
!   grad1_u12_squared_num
!
! with the difference that now the gradients wrt to r1 are computed for all 
! values of r2 of the un-pruned extra grid. 
! This might be a too big array to be stored in practical calculations, 
! but here we are defining provisional providers for testing


 BEGIN_PROVIDER [ double precision, grad1_u12_vect_on_pruned_extra, (n_points_extra_final_grid, n_points_final_grid, 3) ]
&BEGIN_PROVIDER [ double precision, grad1_u12_sqrd_on_pruned_extra, (n_points_extra_final_grid, n_points_final_grid) ]  
  ! Provisory providers following a naming convention
  grad1_u12_vect_on_pruned_extra = grad1_u12_num
  grad1_u12_sqrd_on_pruned_extra = grad1_u12_squared_num
END_PROVIDER




 BEGIN_PROVIDER [ double precision, grad1_u12_vect_on_full_extra_grid, (n_points_tot_extra_grid, n_points_final_grid, 3) ]
&BEGIN_PROVIDER [ double precision, grad1_u12_sqrd_on_full_extra, (n_points_tot_extra_grid, n_points_final_grid) ]  
  implicit none
  integer :: i1
  ! Buffer to fill the provider for different values of r1
  double precision, allocatable :: grad1_u12_vect_at_r1(:,:)
  double precision, allocatable :: grad1_u12_sqrd_at_r1(:)

  ! Initializide providers
  grad1_u12_vect_on_full_extra_grid(:,:,:) = 0.d0
  grad1_u12_sqrd_on_full_extra(:,:) = 0.d0

  allocate(grad1_u12_vect_at_r1(n_points_tot_extra_grid, 3))
  allocate(grad1_u12_sqrd_at_r1(n_points_tot_extra_grid))

  do i1 = 1, n_points_final_grid 
    call get_grad1_u12_on_full_grid2( i1                     &
                               , nucl_num                    & 
                               , n_points_rad_extra_grid     &
                               , n_points_ang_extra_grid     &
                               , grid_points_extra_per_atom  &
                               , grad1_u12_vect_at_r1(:,1:3) &
                               , grad1_u12_sqrd_at_r1(:)     &
                               )
    ! Compute the square of the gradient at r1
    !grad1_u12_sqrd_at_r1(:) = sum(grad1_u12_vect_at_r1(:,1:3)*grad1_u12_vect_at_r1(:,1:3),dim=2)

    grad1_u12_vect_on_full_extra_grid(:,i1,1:3) = grad1_u12_vect_at_r1(:,1:3)
    grad1_u12_sqrd_on_full_extra(:,i1) = grad1_u12_sqrd_at_r1(:)  
  end do

END_PROVIDER


! BEGIN_PROVIDER [ double precision, grad1_u12_vect_on_pruned_extra_grid, (n_points_extra_final_grid, n_points_final_grid, 3) ]
!&BEGIN_PROVIDER [ double precision, grad1_u12_sqrd_on_prune_grid2, (n_points_extra_final_grid, n_points_final_grid) ]  
!  implicit none
!  integer :: i1
!  ! Buffer to fill the provider for different values of r1
!  double precision, allocatable :: grad1_u12_vect_at_r1(:,:)
!  double precision, allocatable :: grad1_u12_sqrd_at_r1(:)
!
!  ! Initialize providers
!  grad1_u12_vect_on_pruned_extra_grid(:,:,:) = 0.d0
!  grad1_u12_sqrd_on_prune_grid2(:,:) = 0.d0
!
!  allocate(grad1_u12_vect_at_r1(n_points_extra_final_grid, 3))
!  allocate(grad1_u12_sqrd_at_r1(n_points_extra_final_grid))
!
!  do i1 = 1, n_points_final_grid 
!    call get_grad1_u12_on_pruned_grid2( i1                        &
!                                      , n_points_extra_final_grid &
!                                      , final_grid_points_extra   &
!                                      , grad1_u12_vect_at_r1(:,1:3) &
!                                      , grad1_u12_sqrd_at_r1(:)     &
!                                      )
!    ! Compute the square of the gradient at r1
!    !grad1_u12_sqrd_at_r1(:) = sum(grad1_u12_vect_at_r1(:,1:3)*grad1_u12_vect_at_r1(:,1:3),dim=2)
!
!    grad1_u12_vect_on_pruned_extra_grid(:,i1,1:3) = grad1_u12_vect_at_r1(:,1:3)
!    grad1_u12_sqrd_on_prune_grid2(:,i1) = grad1_u12_sqrd_at_r1(:)  
!  end do
!
!  !print*, "FINISH WITH  grad1_u12_sqrd_on_prune_grid2"
!
!END_PROVIDER
