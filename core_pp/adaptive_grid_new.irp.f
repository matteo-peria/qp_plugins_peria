 subroutine prune_adaptive_grid(fixed_grid, float_grid,                     &
                            & fixed_grid_weights, float_grid_weights, &
                            & n_fixed_pts_effective, n_float_pts_effective, &
                            & n_pts_effective_max, pruned_grid)
  implicit none
  BEGIN_DOC
  ! Prunes an adaptive grid
  END_DOC
  ! INPUT
  double precision, intent(in) :: fixed_grid(3,n_points_ang_extra_grid,n_points_rad_extra_grid,nucl_num)
  double precision, intent(in) :: float_grid(3,n_points_ang_float_grid,n_points_rad_float_grid,1)
  double precision, intent(in) :: fixed_grid_weights(n_points_ang_extra_grid,n_points_rad_extra_grid,nucl_num)
  double precision, intent(in) :: float_grid_weights(n_points_ang_float_grid,n_points_rad_float_grid,1)
  integer, intent(in)          :: n_fixed_pts_effective(nucl_num)
  integer, intent(in)          :: n_float_pts_effective
  integer, intent(in)          :: n_pts_effective_max
  ! OUTPUT
  double precision, intent(out) :: pruned_grid(3, n_pts_effective_max, nucl_num+1)
  integer :: i,k,i_nucl

  pruned_grid(:,:,:) = 0.d0
  i_nucl = 1
  do i = 1, n_points_rad_float_grid-1
    do k = 1, n_points_ang_float_grid
    end do
  end do
end subroutine


subroutine get_pruned_adaptive_grid()
  implicit none
  !call get_adaptive_grid()
  !call prune_adaptive_grid()
end subroutine

subroutine get_adaptive_grid(r, fixed_grid, float_grid, &
                            & fixed_grid_weights, float_grid_weights, &
                            & n_fixed_pts_effective, n_float_pts_effective, &
                            & n_pts_effective_max)
  implicit none
  BEGIN_DOC
  !
  ! Computes the points of the floating part of an adaptive grid,
  ! the updated weights on the total adaptive grid (floating + fixed part)
  ! and the number of effective points for pruning
  !
  END_DOC
  ! INPUT
  double precision, intent(in)  :: r(3)
  double precision, intent(in)  :: fixed_grid(3,n_points_ang_extra_grid,n_points_rad_extra_grid,nucl_num)
  ! OUTPUT
  double precision, intent(out) :: float_grid(3,n_points_ang_float_grid,n_points_rad_float_grid,1)
  double precision, intent(out) :: fixed_grid_weights(n_points_ang_extra_grid,n_points_rad_extra_grid,nucl_num)
  double precision, intent(out) :: float_grid_weights(n_points_ang_float_grid,n_points_rad_float_grid,1)
  integer,          intent(out) :: n_fixed_pts_effective(nucl_num)
  integer,          intent(out) :: n_float_pts_effective
  integer,          intent(out) :: n_pts_effective_max
  !
  double precision :: a_knowles,radii_ua_av,slater_inter_per_atom(nucl_num)
  double precision :: becke_weights_at_fixed_grid(n_points_ang_extra_grid,n_points_rad_extra_grid,nucl_num)
  double precision :: becke_weights_at_float_grid(n_points_ang_float_grid,n_points_rad_float_grid,1)
  double precision :: weights_per_atom(nucl_num+1)
  double precision :: norm,x,derivative_knowles_function,knowles_function
  double precision :: weight_knowles, weight_total
  integer :: i,k,i_nucl

  ! Compute parameters from the new floating grid
  call get_floating_grid_param(r, a_knowles, radii_ua_av, slater_inter_per_atom)

  ! Compute cartesian coordinates of floating grid only 
  call get_floating_grid_at_r(r, a_knowles, angular_adapt_quadrature_points, float_grid)
  
  ! Compute Becke weights of the whole adaptive grid (floating + extra grid)
  call get_adaptive_grid_becke_weights(r, float_grid, fixed_grid,   &
                                     & slater_inter_per_atom,       &
                                     & becke_weights_at_float_grid, &
                                     & becke_weights_at_fixed_grid)

  ! Floating grid weights (grid points of extra ghost-atom)
  n_float_pts_effective = 0
  i_nucl = 1
  do i = 1, n_points_rad_float_grid-1
    x = grid_adapt_points_radial(i)
    do k = 1, n_points_ang_float_grid
      ! Contribution due to Knowles change of variable
      weight_knowles = derivative_knowles_function(a_knowles, m_knowles, x) &
                    & *knowles_function(a_knowles, m_knowles, x)**2
      ! Quadrature weight * Becke partitioning weight
      float_grid_weights(k,i,i_nucl) = becke_weights_at_float_grid(k,i,i_nucl) &
                                    & *weights_angular_adapt_points(k)  &
                                    & *weight_knowles*dr_adapt_radial_integral
      if (isnan(float_grid_weights(k,i,i_nucl))) then
        print*, 'errore nella float_grid_weights'
        print*, 'cycle: ', i_nucl, i, k
        print*, becke_weights_at_float_grid(k,i,i_nucl), weights_angular_points_extra(k), weight_knowles, dr_radial_extra_integral
        stop
      end if

      ! Pruning on the moving grid
      if (dabs(float_grid_weights(k,i,i_nucl)) < thresh_grid) then                    
        cycle                                                                 
      endif                                                                   
      n_float_pts_effective += 1
    enddo                                                                       
  enddo                                                                         
  
  ! Fixed grid weights (grid points associated to real atoms)
  n_fixed_pts_effective(:) = 0
  do i_nucl = 1, nucl_num
    do i = 1, n_points_rad_extra_grid -1
      x = grid_points_extra_radial(i)
      do k = 1, n_points_ang_extra_grid
        weight_knowles = &
            &  derivative_knowles_function(alpha_knowles(grid_atomic_number(i_nucl)), m_knowles, x) &
            & *knowles_function(alpha_knowles(grid_atomic_number(i_nucl)), m_knowles, x)**2
        ! Quadrature weight * Becke partitioning weight
        fixed_grid_weights(k,i,i_nucl) = becke_weights_at_fixed_grid(k,i,i_nucl)   &
                                      & *weights_angular_points_extra(k)           &
                                      & *weight_knowles * dr_radial_extra_integral
        if (isnan(fixed_grid_weights(k,i,i_nucl))) then
          print*, 'errore nella fixed_grid_weights'
          print*, 'cycle: ', i_nucl, i, k
          print*, becke_weights_at_fixed_grid(k,i,i_nucl), weights_angular_points_extra(k), weight_knowles, dr_radial_extra_integral
          stop
        end if
        ! Pruning on the fixed grid
        if (dabs(fixed_grid_weights(k,i,i_nucl)) < thresh_grid) then                    
          cycle                                                                 
        endif                                                                   
        n_fixed_pts_effective(i_nucl) += 1
      enddo
    enddo
  enddo
  n_pts_effective_max = maxval([n_float_pts_effective, n_fixed_pts_effective])
end subroutine
