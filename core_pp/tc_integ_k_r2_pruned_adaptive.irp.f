subroutine get_int2_grad1_u12_ao_pruned_adaptive( i1, n_aos         &
                                , n_grid2 &
                                , n_ang_float, n_rad_float, n_nuc_float  &
                                , grid_fixed_points,  grid_float_points  &
                                , grid_fixed_weights, grid_float_weights &
                                , int2_grad1_u12_ao_vect                 &
                                , int2_grad1_u12_ao_sqrd)
  BEGIN_DOC
  !
  ! Compute the integrals in r2 involving the gradients of the Jastrow at 
  ! a given r1, specified with the index i1 running on the grid1
  !
  ! int2_grad1_u12_ao_vect_at_r1(j,l,:) = \int dr2 [\grad_r1 J(r1,r2)] \phi_i(r2) \phi_j(r2) 
  ! int2_grad1_u12_ao_sqrd_at_r1(j,l)   = -0.5 \int dr2 chi_l(r2) chi_j(r2) [grad_1 u(r1,r2)]^2
  !
  END_DOC

  implicit none
  ! INPUT
  integer, intent(in) :: i1
  integer, intent(in) :: n_aos
  integer, intent(in) :: n_grid2
  integer, intent(in) :: n_ang_float, n_rad_float, n_nuc_float
  double precision, intent(in) :: grid_fixed_points(3,n_grid2)
  double precision, intent(in) :: grid_float_points(3,n_ang_float,n_rad_float,n_nuc_float)
  double precision, intent(in) :: grid_fixed_weights(n_grid2)
  double precision, intent(in) :: grid_float_weights(n_ang_float,n_rad_float,n_nuc_float)
  ! OUTPUT
  double precision, intent(out) :: int2_grad1_u12_ao_vect(n_aos,n_aos,3)
  double precision, intent(out) :: int2_grad1_u12_ao_sqrd(n_aos,n_aos)

  double precision :: r1(3)
  integer :: j, l, m, i2
 
  double precision :: time0, time1
  
  double precision :: aos_at_r2(ao_num)
  double precision, allocatable :: aos_lj_prod_in_r2(:,:,:)
  double precision, allocatable :: grad1_u12_sqrd_at_r1(:)
  double precision, allocatable :: grad1_u12_vect_at_r1(:,:)
  logical, allocatable :: mask(:)


  integer :: n_fixed, n_float, n_adaptive
  integer :: r2_nuc, r2_rad, r2_ang
  double precision :: ao_l_r2, ao_j_r2
  double precision :: w2, r2(3)


  ! Init the output
  int2_grad1_u12_ao_vect(:,:,:) = 0.d0
  int2_grad1_u12_ao_sqrd(:,:)   = 0.d0


  r1(1:3) = final_grid_points(1:3,i1)

  ! Total number of points belonging to the floating grid
  n_float = int(size(grid_float_points)/3)
  ! Total number of points belonging to the fixed grid
  n_fixed = n_grid2 !int(size(grid_fixed_points)/3)
  ! Total number of points belonging to the fixed (extra) grid +  floating grid
  n_adaptive = n_fixed + n_float 

  ! This array store the scalars: \phi_i(r2) \phi_j(r2) dr2 for r2 in FLOATING GRID
  allocate(aos_lj_prod_in_r2(n_float, ao_num, ao_num))
  ! This array store the scalars: |\grad_r1 J(r1,r2)|^2 for r2 in FLOATING GRID
  allocate(grad1_u12_sqrd_at_r1(n_float))
  ! This array store the vectors: \grad_r1 J(r1,r2) for r2 in FLOATING GRID
  allocate(grad1_u12_vect_at_r1(n_float,3))

  grad1_u12_sqrd_at_r1(:) = 0.d0
  grad1_u12_vect_at_r1(:,:) = 0.d0


! hard silenced !  !!$OMP PARALLEL                             &
! hard silenced !  !!$OMP DEFAULT (NONE)                       &
! hard silenced !  !!$OMP PRIVATE ( r2_rad, r2_ang, l, j, i2     &
! hard silenced !  !!$OMP         , r2, w2, ao_l_r2, ao_j_r2)  &
! hard silenced !  !!$OMP SHARED ( n_points_rad_float_grid, n_points_ang_float_grid, ao_num  &
! hard silenced !  !!$OMP        , aos_lj_prod_in_r2, grid_float_points, grid_float_weights           & 
! hard silenced !  !!$OMP        , aos_at_r2                                              )
! hard silenced !  !!$OMP DO COLLAPSE(2) SCHEDULE (static)
! hard silenced !  ! MOs product at point r2 belonging to the floating grid
! hard silenced !
! hard silenced !  print*, "FLOATING GRID dr2\phi_i(r2)\phi_j(r2) (SHOULD BE ALL ZERO)"
! hard silenced !  do r2_rad = 1, n_points_rad_float_grid - 1
! hard silenced !    do r2_ang = 1, n_points_ang_float_grid
! hard silenced !      ! Index running on the total adaptive grid in r2
! hard silenced !      i2 = (r2_rad - 1) * n_points_ang_float_grid + r2_ang
! hard silenced !      ! Point on grid2 (floating part)
! hard silenced !      r2(1:3) = grid_float_points(1:3,r2_ang,r2_rad,1)
! hard silenced !      ! Weight on grid2 (floating part)
! hard silenced !      w2 = grid_float_weights(r2_ang,r2_rad,1)
! hard silenced !      ! Compute all the MOs at the new r2 point of the floating grid
! hard silenced !      !print*, "CALLING GIVE_ALL_AOS_AT_R at indices: "
! hard silenced !      !print*, "r2_rad, r2_ang, i2: ", r2_rad, r2_ang, i2
! hard silenced !      call give_all_aos_at_r(r2, aos_at_r2)
! hard silenced !      !!print*, "aos_at_r2, when r2 = ", r2!, aos_at_r2
! hard silenced !      !!do l = 1, ao_num
! hard silenced !      !!  print*, l, aos_at_r2(l)
! hard silenced !      !!end do
! hard silenced !      !!stop
! hard silenced !
! hard silenced !      !print*, "aos_at_r2 (r2 = ", r2, "): "
! hard silenced !      !print*, aos_at_r2
! hard silenced !      do l = 1, ao_num
! hard silenced !        ao_l_r2 = aos_at_r2(l)
! hard silenced !        do j = 1, ao_num
! hard silenced !          ao_j_r2 = aos_at_r2(j)
! hard silenced !          ! Update
! hard silenced !          aos_lj_prod_in_r2(i2,j,l) = w2 * ao_l_r2 * ao_j_r2
! hard silenced !          !write(*,'(3A5,3A14,A25)') 'i2','j','l','w2','ao_l_r2','ao_j_r2','aos_lj_prod_in_r2(i2,j,l)'
! hard silenced !          !write(*,'(3I5, 4F14.8)') i2, j, l, w2, ao_l_r2, ao_j_r2, aos_lj_prod_in_r2(i2,j,l)
! hard silenced !        end do
! hard silenced !      end do
! hard silenced !      !!!!i2 = i2+1
! hard silenced !    end do
! hard silenced !  end do
! hard silenced !  !!$OMP END DO
! hard silenced !  !!$OMP END PARALLEL
! hard silenced !
! hard silenced !  print*, "DEALLOCATE(AOS_LJ_PROD_IN_R2)"
! hard silenced !  print*, "ALLOCATE(AOS_LJ_PROD_IN_R2(N_FIXED, AO_NUM, AO_NUM))"
! hard silenced !
! hard silenced !
! hard silenced !  print*, "DONE WITH THE ALLOCATION"
! hard silenced !
! hard silenced !
! hard silenced !  print*, "IM CALLING get_grad1_u12_on_grid2"
! hard silenced !  !call get_grad1_u12_on_grid2_w_pruning( i1, 1                       & 
! hard silenced !  call get_grad1_u12_on_grid2( i1, 1                       & 
! hard silenced !                             , n_points_rad_float_grid     &
! hard silenced !                             , n_points_ang_float_grid     &
! hard silenced !                             , grid_float_points           &
! hard silenced !                             , grad1_u12_vect_at_r1(:,1:3) &
! hard silenced !                             , grad1_u12_sqrd_at_r1(:)     &
! hard silenced !                             )
! hard silenced !
! hard silenced ! print*, "grad1_u12_vect_at_r1 COMPUTED"
! hard silenced ! !allocate(mask(size(grad1_u12_sqrd_at_r1)))
! hard silenced !
! hard silenced ! !if (any(isnan(grad1_u12_sqrd_at_r1))) then
! hard silenced ! !  print*, "There are NaN in grad1_u12_sqrd_at_r1 at iteration: ", i1
! hard silenced ! !  do m = 1, size(mask)
! hard silenced ! !    if (mask(m)) then
! hard silenced ! !      print '(A,I0,A,ES24.16)', 'NaN at index ', m, ' with value ', grad1_u12_sqrd_at_r1(m)       
! hard silenced ! !      !print*, m, grad1_u12_sqrd_at_r1(m)
! hard silenced ! !    end if
! hard silenced ! !  end do
! hard silenced ! !  !print*, grad1_u12_sqrd_at_r1
! hard silenced ! !end if
! hard silenced !
! hard silenced ! print*, "DGEMMING THE FLOAT GRADIENT"
! hard silenced !
! hard silenced ! do m = 1, 3
! hard silenced !   call dgemm( "T", "N", ao_num*ao_num, 1, n_float &
! hard silenced !             , 1.d0, aos_lj_prod_in_r2(:,:,1), n_float      & 
! hard silenced !             , grad1_u12_vect_at_r1(1,m), n_float           &
! hard silenced !             , 0.d0, int2_grad1_u12_ao_vect(1,1,m), ao_num*ao_num) 
! hard silenced ! enddo
! hard silenced !
! hard silenced ! print*, "DGEMMING THE FLOAT GRADIENT SQUARE"
! hard silenced !
! hard silenced ! call dgemm( "T", "N", ao_num*ao_num, 1, n_float &
! hard silenced !           , -0.5d0, aos_lj_prod_in_r2(1,1,1), n_float             &
! hard silenced !           , grad1_u12_sqrd_at_r1(1), n_float                 &
! hard silenced !           , 0.d0, int2_grad1_u12_ao_sqrd(1,1), ao_num*ao_num) 
! hard silenced !
! hard silenced !  deallocate(aos_lj_prod_in_r2)
  allocate(aos_lj_prod_in_r2(n_fixed, ao_num, ao_num))

  !print*, "NOW EVERYTHING SHOULD BE ZERO:"
  !print*, sum(dabs(int2_grad1_u12_ao_vect))
  !print*, sum(dabs(int2_grad1_u12_ao_sqrd))
  !stop


  !!$OMP PARALLEL                                    &
  !!$OMP DEFAULT (NONE)                              &
  !!$OMP PRIVATE ( r2_nuc, r2_rad, r2_ang, l, j, i2     &
  !!$OMP         , r2, w2, ao_l_r2, ao_j_r2       )  &
  !!$OMP SHARED ( nucl_num, n_points_rad_extra_grid, n_points_ang_extra_grid &
  !!$OMP        , ao_num, grid_fixed_points, grid_fixed_weights     &
  !!$OMP        , aos_in_r_full_extra_grid, aos_lj_prod_in_r2                )
  !!$OMP DO COLLAPSE(3) SCHEDULE (static)

  !! MOs product at point r2 belonging to the fixed extra grid

  !print*, "Limiti dei loop: "
  !print*, nucl_num,    n_points_rad_extra_grid, n_points_ang_extra_grid
  !print*, n_nuc_fixed, n_rad_fixed,             n_ang_fixed
  print*, "FIXED GRID dr2\phi_i(r2)\phi_j(r2) (SHOULD BE GOOD)"
  do i2 = 1, n_grid2
    ! Point on grid2 (fixed part)
    r2(1:3) = grid_fixed_points(1:3,i2)
    ! Weight on grid2 (fixed part)
    w2 = grid_fixed_weights(i2)
    do l = 1, ao_num
      ao_l_r2 = aos_in_r_full_extra_grid(l,r2_ang,r2_rad,r2_nuc)              
      do j = 1, ao_num
        ao_j_r2 = aos_in_r_full_extra_grid(j,r2_ang,r2_rad,r2_nuc)              
        ! Update
        aos_lj_prod_in_r2(i2,j,l) = w2 * ao_l_r2 * ao_j_r2
        !write(*,'(3A5,3A14,A25)') 'i2','j','l','w2','ao_l_r2','ao_j_r2','aos_lj_prod_in_r2(i2,j,l)'
        !write(*,'(3I5, 4F14.8)') i2, j, l, w2, ao_l_r2, ao_j_r2, aos_lj_prod_in_r2(i2,j,l)
      end do
    end do
  enddo

  !!$OMP END DO
  !!$OMP END PARALLEL


  do m = 1, 3
    call dgemm( "T", "N", ao_num*ao_num, 1, n_fixed                &
              , 1.d0, aos_lj_prod_in_r2(1,1,1), n_fixed            & 
              , grad1_u12_vect_on_full_extra(:,i1,m), n_fixed        &
              , 0.d0, int2_grad1_u12_ao_vect(1,1,m), ao_num*ao_num &
              ) 
  enddo
  call dgemm( "T", "N", ao_num*ao_num, 1, n_fixed              &
            , -0.5d0, aos_lj_prod_in_r2(1,1,1), n_fixed        &
            , grad1_u12_sqrd_on_full_extra(:,i1), n_fixed        &
            , 0.d0, int2_grad1_u12_ao_sqrd(1,1), ao_num*ao_num &
           ) 


  deallocate(aos_lj_prod_in_r2)
  deallocate(grad1_u12_vect_at_r1)
  deallocate(grad1_u12_sqrd_at_r1)
  !! OLD DEALLOCATES FROM PREVIOUS VERSION OF THE SUBROUTINE
  !deallocate(grid_float_points)
  !deallocate(grid_fixed_weights)
  !deallocate(grid_float_weights)



!  call wall_time(time1)
!  print*, ' wall time for int2_grad1_u12_ao_vect & int2_grad1_u12_ao_sqrd = (min)', (time1-time0) / 60.d0
!  call print_memory_usage()


end subroutine get_int2_grad1_u12_ao
