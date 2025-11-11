 subroutine get_adaptive_grid3( r                                            &
                              , fixed_grid, float_grid                       &
                              , fixed_grid_weights, float_grid_weights)!       &
  implicit none
  BEGIN_DOC
  ! Computes the points of the floating part of an adaptive grid,
  ! the updated weights on the total adaptive grid (floating + fixed part)
  ! and the number of effective points for pruning
  END_DOC
  ! INPUT
  double precision, intent(in)  :: r(3)
  double precision, intent(in)  :: fixed_grid(3,grid3_ang_size,grid3_rad_size,nucl_num)
  ! OUTPUT
  double precision, intent(out) :: float_grid(3,n_points_ang_float_grid,n_points_rad_float_grid,1)
  double precision, intent(out) :: fixed_grid_weights(grid3_ang_size, grid3_rad_size, nucl_num)
  double precision, intent(out) :: float_grid_weights(n_points_ang_float_grid,n_points_rad_float_grid,1)

! These are no longer output
!  integer,          intent(out) :: n_fixed_pts_effective(nucl_num)
!  integer,          intent(out) :: n_float_pts_effective
!  integer,          intent(out) :: n_pts_effective_max
  integer :: n_fixed_pts_effective(nucl_num)
  integer :: n_float_pts_effective
  integer :: n_pts_effective_max

  !
  double precision :: a_knowles
  double precision :: radii_ua_av
  double precision :: slater_inter_per_atom(nucl_num)

  double precision :: becke_weights_at_fixed_grid(grid3_ang_size,grid3_rad_size,nucl_num)
  double precision :: becke_weights_at_float_grid(n_points_ang_float_grid,n_points_rad_float_grid,1)
  double precision :: norm,x,derivative_knowles_function,knowles_function
  double precision :: weight_knowles, weight_total
  integer :: i,k,i_nucl


  !print *, 'Lower bound of 1st dimension: ', lbound(becke_weights_at_fixed_grid, 1)
  !print *, 'Upper bound of 1st dimension: ', ubound(becke_weights_at_fixed_grid, 1)
  !print *, 'Lower bound of 2nd dimension: ', lbound(becke_weights_at_fixed_grid, 2)
  !print *, 'Upper bound of 2nd dimension: ', ubound(becke_weights_at_fixed_grid, 2)
  !print *, 'Lower bound of 3rd dimension: ', lbound(becke_weights_at_fixed_grid, 3)
  !print *, 'Upper bound of 3rd dimension: ', ubound(becke_weights_at_fixed_grid, 3)

  if (my_grid_adapt_param) then 
    ! Adaptive parameters are hard coded into the EZFIO file
    a_knowles = alph_knowles !my_alpha_knowles
    radii_ua_av = radius_ua_av !my_radii_ua_av
    slater_inter_per_atom = slater_rad_ratio_new
  else 
    ! Compute parameters of the new floating grid starting from the one of extra
    print*, "call get_floating_grid_param"
    call get_floating_grid_param(r, a_knowles, radii_ua_av, slater_inter_per_atom)
  end if
  print*, "Parameters obtained"
  print*, "a_knowles                 ", a_knowles
  print*, "radii_ua_av               ", radii_ua_av
  print*, "slater_inter_per_atom     ", slater_inter_per_atom

  ! Compute cartesian coordinates of floating grid only 
  print*, "call get_floating_grid_at_r"
  call get_floating_grid_at_r(r, a_knowles, angular_adapt_quadrature_points, float_grid)

  ! Compute Becke weights of the whole adaptive grid (floating + extra grid)
  print*, "call get_adaptive_grid_becke_weights3"
  call get_adaptive_grid_becke_weights3(r, float_grid, fixed_grid,   &
                                     & slater_inter_per_atom,       &
                                     & becke_weights_at_float_grid, &
                                     & becke_weights_at_fixed_grid)
  print*, becke_weights_at_fixed_grid

  print*, "managing floating grid weights"
  ! Floating grid weights (grid points of extra ghost-atom)
  if (silence_floating_grid.eq..true.) then
    ! no contribution from floating grid
    float_grid_weights(:,:,:) = 0.d0
  else
    n_float_pts_effective = 0
    i_nucl = 1
    do i = 1, n_points_rad_float_grid-1
      print*, "i_rad = ", i
      x = grid_adapt_points_radial(i)
      do k = 1, n_points_ang_float_grid
        print*, "k_rad = ", k
        ! Contribution due to Knowles change of variable
        weight_knowles = derivative_knowles_function(a_knowles, m_knowles, x) &
                      & *knowles_function(a_knowles, m_knowles, x)**2
        ! Quadrature weight * Becke partitioning weight
        !print*,'becke_weights',becke_weights_at_float_grid(k,i,i_nucl),x
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
  endif                                                                   
  
  print*, "managing fixed grid weights"
  ! Fixed grid weights (grid points associated to real atoms)
  n_fixed_pts_effective(:) = 0
  do i_nucl = 1, nucl_num
    print*, "i_nucl = ", i_nucl
    do i = 1, grid3_rad_size -1
      print*, "i = ", i
      !x = grid_points_extra_radial(i)
      x = grid3_rad(i)
      do k = 1, grid3_ang_size
        print*, "k = ", k
        weight_knowles = &
            &  derivative_knowles_function(alpha_knowles(grid_atomic_number(i_nucl)), m_knowles, x) &
            & *knowles_function(alpha_knowles(grid_atomic_number(i_nucl)), m_knowles, x)**2
        print*, "weight_knowles = ", weight_knowles
        ! Quadrature weight * Becke partitioning weight
        print*, "becke_weights_at_fixed_grid(k,i,i_nucl)", becke_weights_at_fixed_grid(k,i,i_nucl), "kind = ", kind(becke_weights_at_fixed_grid(k,i,i_nucl))
        print*, "grid3_lebedev_weights(k)               ", grid3_lebedev_weights(k)               , "kind = ", kind(grid3_lebedev_weights(k)               )
        print*, "weight_knowles                         ", weight_knowles                         , "kind = ", kind(weight_knowles                         )
        print*, "dr3                                    ", dr3                                    , "kind = ", kind(dr3                                    )

        print*, "FUNZIONA?"
        print*, "shape(fixed_grid_weights)          = ", shape(fixed_grid_weights)
        print*, "shape(becke_weights_at_fixed_grid) = ", shape(becke_weights_at_fixed_grid)
print *, "Bounds: ", k, i, i_nucl
if (k < 1 .or. k > size(fixed_grid_weights, 1)) then
   print *, "Index k is out of bounds"
end if
if (i < 1 .or. i > size(fixed_grid_weights, 2)) then
   print *, "Index i is out of bounds"
end if
if (i_nucl < 1 .or. i_nucl > size(fixed_grid_weights, 3)) then
   print *, "Index i_nucl is out of bounds"
end if

        fixed_grid_weights(k,i,i_nucl) = becke_weights_at_fixed_grid(k,i,i_nucl) * grid3_lebedev_weights(k) * weight_knowles * dr3
        print*, "NON FUNZIONA"
        !fixed_grid_weights(k,i,i_nucl) = becke_weights_at_fixed_grid(k,i,i_nucl)   &
        !                              & *weights_angular_points_extra(k)           &
        !                              & *weight_knowles * dr_radial_extra_integral


        print*, "fixed_grid_weights(k,i,i_nucl) = ", fixed_grid_weights(k,i,i_nucl)

        if (isnan(fixed_grid_weights(k,i,i_nucl))) then
          print*, 'errore nella fixed_grid_weights'
          print*, 'cycle: ', i_nucl, i, k
          print*, becke_weights_at_fixed_grid(k,i,i_nucl), grid3_lebedev_weights(k), weight_knowles, dr3
          stop
        end if
        print*, "Pruning, not for real"
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
