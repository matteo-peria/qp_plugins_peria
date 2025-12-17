!! This is an OLD version of the subroutine, where the floating grid is 
!! computed internally. 
!! The NEW versions of this subroutine are
!!
!!   get_int2_grad1_u12_ao_full_adaptive    in tc_integ_k_r2_full_adaptive.irp.f
!!   get_int2_grad1_u12_ao_pruned_adaptive  in tc_integ_k_r2_pruned_adaptive.irp.f 
!!
!! where the full and pruned actually refers to the fixed part of grid2


!subroutine get_int2_grad1_u12_ao_adapt(i1, int2_grad1_u12_ao_vect, int2_grad1_u12_ao_sqrd)
!  BEGIN_DOC
!  ! Compute the integrals in r2 involving the gradients of the Jastrow at 
!  ! a given r1 (indexed by i1 in grid1)
!  !
!  ! int2_grad1_u12_ao_vect_at_r1(j,l,:) = \int dr2 [\grad_r1 J(r1,r2)] \phi_i(r2) \phi_j(r2) 
!  ! int2_grad1_u12_ao_sqrd_at_r1(j,l)   = -0.5 \int dr2 chi_l(r2) chi_j(r2) [grad_1 u(r1,r2)]^2
!  ! 
!  END_DOC
!
!  implicit none
!  ! INPUT
!  integer, intent(in) :: i1
!  ! OUTPUT
!  double precision, intent(out) :: int2_grad1_u12_ao_vect(:,:,:)
!  double precision, intent(out) :: int2_grad1_u12_ao_sqrd(:,:)
!
!  double precision :: r1(3)
!  integer :: j, l, m, i2
!  !integer :: i1, i2, i1_start, i1_final
!  !integer :: i1, i2, i1_start, i1_final
!  
!  !integer :: chunk_size, chunk_rest, num_chunks
!  !integer :: chunk_col, rest_col, chunk_idx
!  
!  double precision :: time0, time1
!  !double precision :: mem, num_double
!  
!  double precision :: aos_at_r2(ao_num)
!  double precision, allocatable :: aos_lj_prod_in_r2(:,:,:)
!!  double precision, allocatable :: grad1_u12_vect(:,:,:)
!!  double precision, allocatable :: grad1_u12_sqrd(:,:)
!!  ! Buffer-array to fill the above ones
!  double precision, allocatable :: grad1_u12_sqrd_at_r1(:)
!  double precision, allocatable :: grad1_u12_vect_at_r1(:,:)
!
!  !! Dimensions of the following arrays are not known at compile time
!  !double precision :: grid_float_points(3, n_points_ang_float_grid, n_points_rad_float_grid, 1)
!  !double precision :: grid_fixed_weights(n_points_ang_extra_grid, n_points_rad_extra_grid, nucl_num)
!  !double precision :: grid_float_weights(n_points_ang_float_grid, n_points_rad_float_grid, 1)
!  !
!  double precision, allocatable :: grid_float_points(:,:,:,:)
!  double precision, allocatable :: grid_fixed_weights(:,:,:)
!  double precision, allocatable :: grid_float_weights(:,:,:)
!
!  ! Dummy variables necessary for pruning on floating grid (not coded yet)
!  integer :: n_fixed_pts_effective(nucl_num)
!  integer :: n_float_pts_effective
!  integer :: n_pts_effective_max
!
!  integer :: n_float, n_adaptive
!  integer :: r2_nuc, r2_rad, r2_ang, ao_l_r2, ao_j_r2
!  double precision :: w2, r2(3)
!
!
!  print*, i1
!  print*, ' Providing int2_grad1_u12_ao_vect & int2_grad1_u12_ao_sqrd ...'
!
!  call wall_time(time0)
!
!  provide n_points_ang_float_grid n_points_rad_float_grid
!  provide n_points_ang_extra_grid n_points_rad_extra_grid nucl_num
!
!
!  allocate(grid_float_points(3, n_points_ang_float_grid, n_points_rad_float_grid, 1))
!  allocate(grid_fixed_weights(n_points_ang_extra_grid, n_points_rad_extra_grid, nucl_num))
!  allocate(grid_float_weights(n_points_ang_float_grid, n_points_rad_float_grid, 1))
!
!
!  !!!! SUBSTITUTE THE FOLLOWING
!  !!!PROVIDE final_weight_at_r_vector_extra 
!  !!!provide aos_in_r_array_extra
!
!  print*, "final_grid_points(1:3,", i1, ") = ", final_grid_points(1:3,i1)
!  
!  r1(1:3) = final_grid_points(1:3,i1)
!
!  ! Retrieve the weights grid 2 composed by the fixed and the floating grid
!  call get_adaptive_grid(r1                                            &
!                       , grid_points_extra_per_atom, grid_float_points &                      
!                       , grid_fixed_weights, grid_float_weights        &         
!                       , n_fixed_pts_effective, n_float_pts_effective  &
!                       , n_pts_effective_max)
!
!  ! Total number of points belonging to the floating grid
!  n_float = size(grid_float_points)
!  ! Total number of points belonging to the floating grid + extra grid
!  n_adaptive = size(grid_points_extra_per_atom) + n_float 
!
!  allocate(aos_lj_prod_in_r2(n_adaptive, ao_num, ao_num))
!  allocate(grad1_u12_sqrd_at_r1(n_adaptive))
!  allocate(grad1_u12_vect_at_r1(n_adaptive,3))
!
!
!  !$OMP PARALLEL                             &
!  !$OMP DEFAULT (NONE)                       &
!  !$OMP PRIVATE ( r2_rad, r2_ang, l, j, i2     &
!  !$OMP         , r2, w2, ao_l_r2, ao_j_r2)  &
!  !$OMP SHARED ( n_points_rad_float_grid, n_points_ang_float_grid, ao_num  &
!  !$OMP        , aos_lj_prod_in_r2, grid_float_points, grid_float_weights           & 
!  !$OMP        , aos_at_r2                                              )
!  !$OMP DO COLLAPSE(2) SCHEDULE (static)
!  ! MOs product at point r2 belonging to the floating grid
!  do r2_rad = 1, n_points_rad_float_grid - 1
!    do r2_ang = 1, n_points_ang_float_grid
!      ! Index running on the total adaptive grid in r2
!      i2 = (r2_rad - 1) * n_points_ang_float_grid + r2_ang
!      ! Point on grid2 (floating part)
!      r2(1:3) = grid_float_points(1:3,r2_ang,r2_rad,1)
!      ! Weight on grid2 (floating part)
!      w2 = grid_float_weights(r2_ang,r2_rad,1)
!      ! Compute all the MOs at the new r2 point of the floating grid
!      call give_all_aos_at_r(r2, aos_at_r2)
!      do l = 1, ao_num
!        ao_l_r2 = aos_at_r2(l)
!        do j = 1, ao_num
!          ao_j_r2 = aos_at_r2(j)
!          ! Update
!          aos_lj_prod_in_r2(i2,j,l) = w2 * ao_l_r2 * ao_j_r2
!        end do
!      end do
!      i2 = i2+1
!    end do
!  end do
!  !$OMP END DO
!  !$OMP END PARALLEL
!
!  !$OMP PARALLEL                                    &
!  !$OMP DEFAULT (NONE)                              &
!  !$OMP PRIVATE ( r2_nuc, r2_rad, r2_ang, l, j, i2     &
!  !$OMP         , r2, w2, ao_l_r2, ao_j_r2       )  &
!  !$OMP SHARED ( nucl_num, n_points_rad_extra_grid, n_points_ang_extra_grid &
!  !$OMP        , ao_num, grid_points_extra_per_atom, grid_fixed_weights     &
!  !$OMP        , aos_in_r_full_extra_grid, aos_lj_prod_in_r2                )
!  !$OMP DO COLLAPSE(3) SCHEDULE (static)
!  ! MOs product at point r2 belonging to the fixed extra grid
!  do r2_nuc = 1, nucl_num
!    do r2_rad = 1, n_points_rad_extra_grid-1
!      do r2_ang = 1, n_points_ang_extra_grid
!        ! Index running on the total adaptive grid in r2
!        i2 = (r2_nuc-1)*(n_points_rad_extra_grid-1)*n_points_ang_extra_grid + &
!             (r2_rad-1)*n_points_ang_extra_grid + r2_ang
!        ! Point on grid2 (fixed part)
!        r2(1:3) = grid_points_extra_per_atom(1:3,r2_ang,r2_rad,r2_nuc)
!        ! Weight on grid2 (fixed part)
!        w2 = grid_fixed_weights(r2_ang,r2_rad,r2_nuc)
!        do l = 1, ao_num
!          ao_l_r2 = aos_in_r_full_extra_grid(l,r2_ang,r2_rad,r2_nuc)              
!          do j = 1, ao_num
!            ao_j_r2 = aos_in_r_full_extra_grid(j,r2_ang,r2_rad,r2_nuc)              
!            ! Update
!            aos_lj_prod_in_r2(i2,j,l) = w2 * ao_l_r2 * ao_j_r2
!          end do
!        end do
!        i2 = i2+1
!      enddo
!    enddo
!  enddo
!  !$OMP END DO
!  !$OMP END PARALLEL
!
!  print*, "IM CALLING get_grad1_u12_on_grid2"
!  call get_grad1_u12_on_grid2( i1, 1                       & 
!                             , n_points_rad_float_grid     &
!                             , n_points_ang_float_grid     &
!                             , grid_float_points           &
!                             , grad1_u12_vect_at_r1(:,1:3) &
!                             , grad1_u12_sqrd_at_r1(:)     &
!                             )
! do m = 1, 3
!   call dgemm( "T", "N", ao_num*ao_num, 1, n_adaptive &
!             , 1.d0, aos_lj_prod_in_r2(1,1,1), n_adaptive               & 
!             , grad1_u12_vect_at_r1(1,m), n_adaptive           &
!             , 0.d0, int2_grad1_u12_ao_vect(1,1,m), ao_num*ao_num) 
! enddo
! call dgemm( "T", "N", ao_num*ao_num, 1, n_adaptive &
!           , -0.5d0, aos_lj_prod_in_r2(1,1,1), n_adaptive             &
!           , grad1_u12_sqrd_at_r1(1), n_adaptive                 &
!           , 0.d0, int2_grad1_u12_ao_sqrd(1,1), ao_num*ao_num) 
! 
!
!  deallocate(aos_lj_prod_in_r2)
!  deallocate(grad1_u12_vect_at_r1)
!  deallocate(grad1_u12_sqrd_at_r1)
!  deallocate(grid_float_points)
!  deallocate(grid_fixed_weights)
!  deallocate(grid_float_weights)
!
!
!
!  call wall_time(time1)
!  print*, ' wall time for int2_grad1_u12_ao_vect & int2_grad1_u12_ao_sqrd = (min)', (time1-time0) / 60.d0
!  call print_memory_usage()
!
!
!end subroutine get_int2_grad1_u12_ao_adapt
