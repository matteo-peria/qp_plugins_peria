 BEGIN_PROVIDER [integer, n_points_rad_grid]
&BEGIN_PROVIDER [integer, n_points_ang_grid]  
  implicit none
  n_points_rad_grid = n_points_radial_grid
  n_points_ang_grid = n_points_integration_angular
END_PROVIDER


 BEGIN_PROVIDER [integer, n_points_rad_extra_grid]
&BEGIN_PROVIDER [integer, n_points_ang_extra_grid]
&BEGIN_PROVIDER [integer, n_points_tot_extra_grid]
  implicit none
  n_points_rad_extra_grid = n_points_extra_radial_grid
  n_points_ang_extra_grid = n_points_extra_integration_angular
  ! Tot num of points somehow doesn't consider the last point in the radial direction
  n_points_tot_extra_grid = nucl_num*(n_points_rad_extra_grid-1)*n_points_ang_extra_grid
END_PROVIDER


 BEGIN_PROVIDER [integer, n_points_rad_float_grid]
&BEGIN_PROVIDER [integer, n_points_ang_float_grid]
  implicit none
  n_points_rad_float_grid = n_points_radial_grid_adapt
  n_points_ang_float_grid = n_points_integration_angular_adapt
END_PROVIDER


 BEGIN_PROVIDER [integer, n_points_rad_adapt_grid]
&BEGIN_PROVIDER [integer, n_points_ang_adapt_grid]
  implicit none
  n_points_rad_adapt_grid = n_points_radial_grid_adapt
  n_points_ang_adapt_grid = n_points_integration_angular_adapt
END_PROVIDER
