! grid_adapt_points_radial dr_adapt_radial_integral


 BEGIN_PROVIDER [double precision, alph_knowles]
  alph_knowles = my_alpha_knowles
END_PROVIDER

 BEGIN_PROVIDER [double precision, radius_ua_av]
&BEGIN_PROVIDER [double precision, slater_rad_ratio_new, (nucl_num)]

  radius_ua_av = my_radii_ua_av 

  double precision :: chi, u_ij
  integer :: j
  do j = 1, nucl_num
    chi = radius_ua_av / slater_bragg_radii_per_atom_ua(j)
    !chi =  slater_bragg_radii_per_atom_ua(j) /radius_ua_av 
    u_ij = (chi - 1.d0 ) / (chi + 1.d0)
    slater_rad_ratio_new(j) = u_ij  / (u_ij * u_ij - 1.d0)
    if(slater_rad_ratio_new(j).gt.0.5d0)then
      slater_rad_ratio_new(j) = 0.5d0
    else if( slater_rad_ratio_new(j) .le.-0.5d0)then
      slater_rad_ratio_new(j) = -0.5d0
    endif
  enddo

  print*, "SLATER_RAD_RATIO_NEW, ", SLATER_RAD_RATIO_NEW
END_PROVIDER


 BEGIN_PROVIDER [integer, n_points_radial_grid_adapt]
&BEGIN_PROVIDER [integer, n_points_integration_angular_adapt]
  BEGIN_DOC
  ! n_points_radial_grid_adapt = number of radial grid points per atom
  !
  ! n_points_integration_angular_adapt = number of angular grid points per atom
  !
  ! These numbers are automatically set by setting the grid_adapt_type_sgn parameter
  END_DOC
 
  implicit none

  if(.not. my_grid_adapt) then
    print*, "MY_GRID_ADAPT be careful in modifying these values"
    select case (adapt_grid_type_sgn)
      case(0)
        n_points_radial_grid_adapt = 23
        n_points_integration_angular_adapt = 170  
      case(1)
        n_points_radial_grid_adapt = 50
        n_points_integration_angular_adapt = 194
      case(2)
        n_points_radial_grid_adapt = 30
        n_points_integration_angular_adapt = 194
      case(3)
        n_points_radial_grid_adapt = 99
        n_points_integration_angular_adapt = 590
      case default
        write(*,*) '!!! Quadrature grid not available !!!'
        stop
    end select
  else
    n_points_radial_grid_adapt = my_n_pt_r_grid_adapt
    n_points_integration_angular_adapt = my_n_pt_a_grid_adapt
  endif

  !if (verbose) then
    write(6,*) 'Adaptive grid parameters:'
    call write_int(6, n_points_radial_grid_adapt, 'N radial points ')
    call write_int(6, n_points_integration_angular_adapt, 'N angular points')
    write(6,*) ''
  !end if
END_PROVIDER

! ---

BEGIN_PROVIDER [integer, n_total_adapt_grid]
  BEGIN_DOC
  ! Number of points in the adaptive grid
  END_DOC
  implicit none
  n_total_adapt_grid = n_points_grid_adapt_per_atom*(nucl_num+1) 
END_PROVIDER 


BEGIN_PROVIDER [integer, n_points_grid_adapt_per_atom]
  BEGIN_DOC
  ! Number of grid points per atom
  END_DOC
  implicit none
  n_points_grid_adapt_per_atom = n_points_ang_extra_grid * n_points_rad_extra_grid
END_PROVIDER


 BEGIN_PROVIDER [double precision, grid_adapt_points_radial, (n_points_rad_adapt_grid)]
&BEGIN_PROVIDER [double precision, dr_adapt_radial_integral]
  BEGIN_DOC
  ! points in [0,1] to map the radial integral [0,\infty]
  END_DOC
  implicit none
  integer :: i
  dr_adapt_radial_integral = 1.d0 / dble(n_points_rad_adapt_grid-1)
  !do i = 1, n_points_rad_extra_grid
  do i = 1, n_points_rad_adapt_grid
    grid_adapt_points_radial(i) = dble(i-1) * dr_adapt_radial_integral
  enddo
END_PROVIDER
