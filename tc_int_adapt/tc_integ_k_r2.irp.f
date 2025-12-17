subroutine get_int2_grad1_u12_ao( i1, n_aos                              &
                                , n_ang_fixed, n_rad_fixed, n_nuc_fixed  &
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
  integer, intent(in) :: n_ang_fixed, n_rad_fixed, n_nuc_fixed
  integer, intent(in) :: n_ang_float, n_rad_float, n_nuc_float
  double precision, intent(in) :: grid_float_points(3,n_ang_float,n_rad_float,n_nuc_float)
  double precision, intent(in) :: grid_fixed_points(3,n_ang_fixed,n_rad_fixed,n_nuc_fixed)
  double precision, intent(in) :: grid_float_weights(n_ang_float,n_rad_float,n_nuc_float)
  double precision, intent(in) :: grid_fixed_weights(n_ang_fixed,n_rad_fixed,n_nuc_fixed)
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


  integer :: n_float, n_adaptive
  integer :: r2_nuc, r2_rad, r2_ang
  double precision :: ao_l_r2, ao_j_r2
  double precision :: w2, r2(3)

  !print*, i1
  !print*, ' Providing int2_grad1_u12_ao_vect & int2_grad1_u12_ao_sqrd ...'

!  call wall_time(time0)

  ! OLD STUFF FROM PREVIOUS VERSION OF THIS SUBROUTINE
  !TRY AVOID PROVIDE n_points_ang_float_grid n_points_rad_float_grid
  !TRY AVOID PROVIDE n_points_ang_extra_grid n_points_rad_extra_grid nucl_num

  !TRY AVOID allocate(grid_float_points(3, n_points_ang_float_grid, n_points_rad_float_grid, 1))
  !TRY AVOID allocate(grid_fixed_weights(n_points_ang_extra_grid, n_points_rad_extra_grid, nucl_num))
  !TRY AVOID allocate(grid_float_weights(n_points_ang_float_grid, n_points_rad_float_grid, 1))


  r1(1:3) = final_grid_points(1:3,i1)

  !! Retrieve the weights grid 2 composed by the fixed and the floating grid
  !call get_adaptive_grid(r1                                            &
  !                     , grid_fixed_points, grid_float_points &                      
  !                     , grid_fixed_weights, grid_float_weights        &         
  !                     , n_fixed_pts_effective, n_float_pts_effective  &
  !                     , n_pts_effective_max)

  ! Total number of points belonging to the floating grid
  n_float = int(size(grid_float_points)/3)
  ! Total number of points belonging to the floating grid + extra grid
  n_adaptive = int(size(grid_fixed_points))/3 + n_float 

  !print*, "size(grid_float_points)", int(size(grid_float_points)/3)
  !print*, "size(grid_fixed_points)", int(size(grid_fixed_points)/3)
  !print*, "N_FLOAT: ", N_FLOAT, "N_ADAPTIVE: ", N_ADAPTIVE

  ! This array store the scalars: \phi_i(r2) \phi_j(r2) dr2
  allocate(aos_lj_prod_in_r2(n_adaptive, ao_num, ao_num))
  ! This array store the scalars: |\grad_r1 J(r1,r2)|^2
  allocate(grad1_u12_sqrd_at_r1(n_adaptive))
  ! This array store the vectors: \grad_r1 J(r1,r2)
  allocate(grad1_u12_vect_at_r1(n_adaptive,3))

  grad1_u12_sqrd_at_r1(:) = 0.d0
  grad1_u12_vect_at_r1(:,:) = 0.d0


  !!$OMP PARALLEL                             &
  !!$OMP DEFAULT (NONE)                       &
  !!$OMP PRIVATE ( r2_rad, r2_ang, l, j, i2     &
  !!$OMP         , r2, w2, ao_l_r2, ao_j_r2)  &
  !!$OMP SHARED ( n_points_rad_float_grid, n_points_ang_float_grid, ao_num  &
  !!$OMP        , aos_lj_prod_in_r2, grid_float_points, grid_float_weights           & 
  !!$OMP        , aos_at_r2                                              )
  !!$OMP DO COLLAPSE(2) SCHEDULE (static)
  ! MOs product at point r2 belonging to the floating grid

  print*, "FLOATING GRID dr2\phi_i(r2)\phi_j(r2) (SHOULD BE ALL ZERO)"
  do r2_rad = 1, n_points_rad_float_grid - 1
    do r2_ang = 1, n_points_ang_float_grid
      ! Index running on the total adaptive grid in r2
      i2 = (r2_rad - 1) * n_points_ang_float_grid + r2_ang
      ! Point on grid2 (floating part)
      r2(1:3) = grid_float_points(1:3,r2_ang,r2_rad,1)
      ! Weight on grid2 (floating part)
      w2 = grid_float_weights(r2_ang,r2_rad,1)
      ! Compute all the MOs at the new r2 point of the floating grid
      !print*, "CALLING GIVE_ALL_AOS_AT_R at indices: "
      !print*, "r2_rad, r2_ang, i2: ", r2_rad, r2_ang, i2
      call give_all_aos_at_r(r2, aos_at_r2)
      !!print*, "aos_at_r2, when r2 = ", r2!, aos_at_r2
      !!do l = 1, ao_num
      !!  print*, l, aos_at_r2(l)
      !!end do
      !!stop

      !print*, "aos_at_r2 (r2 = ", r2, "): "
      !print*, aos_at_r2
      do l = 1, ao_num
        ao_l_r2 = aos_at_r2(l)
        do j = 1, ao_num
          ao_j_r2 = aos_at_r2(j)
          ! Update
          aos_lj_prod_in_r2(i2,j,l) = w2 * ao_l_r2 * ao_j_r2
          print*, " i2,  j,  l,        w2,  ao_l_r2,  ao_j_r2, aos_lj_prod_in_r2(i2,j,l): "
          write(*,'(3I4, 4F10.7)') i2, j, l, w2, ao_l_r2, ao_j_r2, aos_lj_prod_in_r2(i2,j,l)
        end do
      end do
      !!!!i2 = i2+1
    end do
  end do
  !!$OMP END DO
  !!$OMP END PARALLEL

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
  do r2_nuc = 1, nucl_num
    do r2_rad = 1, n_points_rad_extra_grid-1
      do r2_ang = 1, n_points_ang_extra_grid
  !do r2_nuc = 1, nucl_num
  !  do r2_rad = 1, n_points_rad_extra_grid-1
  !    do r2_ang = 1, n_points_ang_extra_grid
        ! Index running on the total adaptive grid in r2
        i2 = (r2_nuc-1)*(n_points_rad_extra_grid-1)*n_points_ang_extra_grid + &
             (r2_rad-1)*n_points_ang_extra_grid + r2_ang + n_float
        ! Point on grid2 (fixed part)
        r2(1:3) = grid_fixed_points(1:3,r2_ang,r2_rad,r2_nuc)
        ! Weight on grid2 (fixed part)
        w2 = grid_fixed_weights(r2_ang,r2_rad,r2_nuc)
        do l = 1, ao_num
          ao_l_r2 = aos_in_r_full_extra_grid(l,r2_ang,r2_rad,r2_nuc)              
          do j = 1, ao_num
            ao_j_r2 = aos_in_r_full_extra_grid(j,r2_ang,r2_rad,r2_nuc)              
            ! Update
            aos_lj_prod_in_r2(i2,j,l) = w2 * ao_l_r2 * ao_j_r2
            print*, " i2,  j,  l,        w2,  ao_l_r2,  ao_j_r2, aos_lj_prod_in_r2(i2,j,l): "
            write(*,'(3I4, 4F10.7)') i2, j, l, w2, ao_l_r2, ao_j_r2, aos_lj_prod_in_r2(i2,j,l)
          end do
        end do
        !i2 = i2+1
      enddo
    enddo
  enddo
  !!$OMP END DO
  !!$OMP END PARALLEL

  stop

  !print*, "IM CALLING get_grad1_u12_on_grid2"
  !call get_grad1_u12_on_grid2_w_pruning( i1, 1                       & 
  call get_grad1_u12_on_grid2( i1, 1                       & 
                             , n_points_rad_float_grid     &
                             , n_points_ang_float_grid     &
                             , grid_float_points           &
                             , grad1_u12_vect_at_r1(:,1:3) &
                             , grad1_u12_sqrd_at_r1(:)     &
                             )

 !print*, "grad1_u12_vect_at_r1"
 !if (any(isnan(grad1_u12_vect_at_r1))) then
 !  print*, "There are NaN in grad1_u12_vect_at_r1"
 !  print*, grad1_u12_vect_at_r1
 !end if

 !print*, "grad1_u12_sqrd_at_r1"
 allocate(mask(size(grad1_u12_sqrd_at_r1)))

 if (any(isnan(grad1_u12_sqrd_at_r1))) then
   print*, "There are NaN in grad1_u12_sqrd_at_r1 at iteration: ", i1
   do m = 1, size(mask)
     if (mask(m)) then
       print '(A,I0,A,ES24.16)', 'NaN at index ', m, ' with value ', grad1_u12_sqrd_at_r1(m)       
       !print*, m, grad1_u12_sqrd_at_r1(m)
     end if
   end do
   !print*, grad1_u12_sqrd_at_r1
 end if
 !stop

 do m = 1, 3
   call dgemm( "T", "N", ao_num*ao_num, 1, n_adaptive &
             , 1.d0, aos_lj_prod_in_r2(1,1,1), n_adaptive               & 
             , grad1_u12_vect_at_r1(1,m), n_adaptive           &
             , 0.d0, int2_grad1_u12_ao_vect(1,1,m), ao_num*ao_num) 
 enddo
 call dgemm( "T", "N", ao_num*ao_num, 1, n_adaptive &
           , -0.5d0, aos_lj_prod_in_r2(1,1,1), n_adaptive             &
           , grad1_u12_sqrd_at_r1(1), n_adaptive                 &
           , 0.d0, int2_grad1_u12_ao_sqrd(1,1), ao_num*ao_num) 
 

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
