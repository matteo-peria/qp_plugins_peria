subroutine give_adapt_grid_at_r(r_input, grid_points, final_adapt_weights)
 implicit none
 BEGIN_DOC
 !
 ! grid_points(1:3,i,j,k) where the indices are
 ! = x/y/z for the ith angular, jth radial points attached to the kth atom
 ! 
 ! but we add an extra atom corresponding, LOCATED FIRST in the list, located at r_input
 !
 ! NB: this is the old adaptive grid where the number of points on the moving
 ! grid is the same as the extra grid. In order to make a subroutine that allow
 ! for real different treatment of floting (moving) and fixed extra grid,
 ! in this subroutine I adopt GRID_POINTS_EXTRA_RADIAL and 
 ! DR_RADIAL_EXTRA_INTEGRAL as points and integral-bin-width even for the
 ! adaptive grid
 !
 !
 END_DOC
 double precision, intent(in) :: r_input(3)
 !
 double precision, intent(out):: grid_points(3,n_points_ang_extra_grid, &
                                           & n_points_rad_extra_grid,   &
                                           & nucl_num+1)
 double precision, intent(out):: final_adapt_weights(n_points_ang_extra_grid, &
                                                   & n_points_rad_extra_grid, &
                                                   & nucl_num+1)
 !
 double precision :: alpha_av,radii_ua_av,slater_inter_per_atom(nucl_num)
 double precision :: weight_tmp(n_points_ang_extra_grid, &
                              & n_points_rad_extra_grid, &
                              & nucl_num+1)
 double precision :: weights_per_atom(nucl_num+1)
 double precision :: norm, x
 double precision :: derivative_knowles_function, knowles_function, contrib_integration

 integer :: i,k,i_nucl

 call get_info_at_r_adapt_grid(r_input,alpha_av,radii_ua_av,slater_inter_per_atom)
 call get_all_points_at_r_adapt_grid(r_input,alpha_av,grid_points)
 call get_voronoi_partition(r_input,grid_points,slater_inter_per_atom,weight_tmp)

 ! Set the extra atom portion of the grid (corresponding to 1st entry)
 i_nucl = 1
 do i = 1, n_points_rad_extra_grid-1
   !x = grid_adapt_points_radial(i)  ! See NB
   x = grid_points_extra_radial(i)
   do k = 1, n_points_ang_extra_grid
     contrib_integration = derivative_knowles_function(alpha_av, m_knowles, x) &
                         * knowles_function(alpha_av, m_knowles, x)**2
     !final_adapt_weights(k,i,i_nucl) = weights_angular_adapt_points(k) & ! See NB
     final_adapt_weights(k,i,i_nucl) = weights_angular_points_extra(k) &
                                     & *weight_tmp(k,i,i_nucl)         &
                                     & *contrib_integration            &
                                     !& *dr_adapt_radial_integral ! See NB
                                     & *dr_radial_extra_integral
   enddo
 enddo
 
 ! Set all the remaining real atoms portions of the grid (starting from the 2nd)
 do i_nucl = 2, nucl_num+1
   do i = 1, n_points_rad_extra_grid-1
     !x = grid_adapt_points_radial(i) ! See NB
     x = grid_points_extra_radial(i)
     do k = 1, n_points_ang_extra_grid
       contrib_integration = &
           & derivative_knowles_function(alpha_knowles(grid_atomic_number(i_nucl-1)), m_knowles, x) &
           & * knowles_function(alpha_knowles(grid_atomic_number(i_nucl-1)), m_knowles, x)**2
       
       !final_adapt_weights(k,i,i_nucl) = weights_angular_adapt_points(k) & ! See NB
       final_adapt_weights(k,i,i_nucl) = weights_angular_points_extra(k) &
                                      & *weight_tmp(k,i,i_nucl)          &
                                      & *contrib_integration             &
                                      !& *dr_adapt_radial_integral
                                      & *dr_radial_extra_integral
     enddo
   enddo
 enddo

end subroutine



 subroutine get_info_at_r_adapt_grid(r_input,alpha_av,radii_ua_av,slater_inter_per_atom)
 implicit none
 double precision, intent(in) :: r_input(3)
 double precision, intent(out):: alpha_av, radii_ua_av, slater_inter_per_atom(nucl_num)
 double precision :: dist,norm,w_i,j
 integer :: i_nucl, k
 alpha_av = 0.d0
 radii_ua_av = 0.d0
 norm = 0.d0
 do i_nucl = 1, nucl_num
  dist = 0.d0
  do k = 1, 3
   dist += (nucl_coord(i_nucl,k)-r_input(k)) * (nucl_coord(i_nucl,k)-r_input(k))
  enddo
  dist = dsqrt(dist)
  w_i = dexp(-dist)  
  norm += w_i
  alpha_av += w_i * alpha_knowles(grid_atomic_number(i_nucl))
  radii_ua_av += w_i * slater_bragg_radii_per_atom_ua(i_nucl)
 enddo
 alpha_av = alpha_av/norm
 radii_ua_av = radii_ua_av/norm
 
 double precision :: xhi_tmp,u_ij
 !! determines the slater_inter distance radii for each atoms with respect to r_input
 do j = 1, nucl_num
  xhi_tmp = radii_ua_av / slater_bragg_radii_per_atom_ua(j)
  u_ij = (xhi_tmp - 1.d0 ) / (xhi_tmp +1.d0)
  slater_inter_per_atom(j) = u_ij  / (u_ij * u_ij - 1.d0)
  if(slater_inter_per_atom(j).gt.0.5d0)then
   slater_inter_per_atom(j) = 0.5d0
  else if( slater_inter_per_atom(j) .le.-0.5d0)then
   slater_inter_per_atom(j) = -0.5d0
  endif
 enddo

 end


subroutine get_all_points_at_r_adapt_grid(r_input,alpha_av,grid_points)
  implicit none
  BEGIN_DOC
  ! subroutine that returns all grid points with respect to r_input, 
  ! with the first spherical grid centered around r_input
  END_DOC
  double precision, intent(in)  :: r_input(3),alpha_av
  !double precision, intent(out) :: grid_points(3,n_points_integration_angular_adapt,n_points_radial_grid_adapt,nucl_num+1)
  double precision, intent(out) :: grid_points(3,n_points_extra_integration_angular,& 
                                             & n_points_extra_radial_grid, &
                                             & nucl_num+1)
  integer :: i_nucl,j,k
  double precision :: x_ref,y_ref,z_ref,x,r,knowles_function
  !
  ! First fill the array of points with first the points centered on r_input
  i_nucl = 1
  !
  x_ref = r_input(1)
  y_ref = r_input(2)
  z_ref = r_input(3)

  ! Loop over all the radial points
  !do j = 1, n_points_radial_grid_adapt-1
  do j = 1, n_points_extra_radial_grid-1
    !x = grid_adapt_points_radial(j)
    x = grid_points_extra_radial(j)
    r = knowles_function(alpha_av, m_knowles, x)
    ! Loop over all the angular points at fixed R
    !do k = 1, n_points_integration_angular_adapt
    do k = 1, n_points_extra_integration_angular
      !grid_points(1,k,j,i_nucl) = x_ref + angular_adapt_quadrature_points(k,1) * r
      !grid_points(2,k,j,i_nucl) = y_ref + angular_adapt_quadrature_points(k,2) * r
      !grid_points(3,k,j,i_nucl) = z_ref + angular_adapt_quadrature_points(k,3) * r
      grid_points(1,k,j,i_nucl) = x_ref + angular_quadrature_points_extra(k,1) * r
      grid_points(2,k,j,i_nucl) = y_ref + angular_quadrature_points_extra(k,2) * r
      grid_points(3,k,j,i_nucl) = z_ref + angular_quadrature_points_extra(k,3) * r
    enddo
  enddo
  !
  ! Then fills the rest of the grid 
  do i_nucl = 1, nucl_num
    x_ref = nucl_coord(i_nucl,1)
    y_ref = nucl_coord(i_nucl,2)
    z_ref = nucl_coord(i_nucl,3)
    !do j = 1, n_points_radial_grid_adapt-1
    do j = 1, n_points_extra_radial_grid-1
      ! x value for the mapping of the [0, +\infty] to [0,1]
      !x = grid_adapt_points_radial(j)
      x = grid_points_extra_radial(j)
      ! value of the radial coordinate for the integration
      r = knowles_function(alpha_knowles(grid_atomic_number(i_nucl)), m_knowles, x)

      ! explicit values of the grid points centered around each atom
      !do k = 1, n_points_integration_angular_adapt
      do k = 1, n_points_extra_integration_angular
        !grid_points(1,k,j,i_nucl+1) = x_ref + angular_adapt_quadrature_points(k,1) * r
        !grid_points(2,k,j,i_nucl+1) = y_ref + angular_adapt_quadrature_points(k,2) * r
        !grid_points(3,k,j,i_nucl+1) = z_ref + angular_adapt_quadrature_points(k,3) * r
        grid_points(1,k,j,i_nucl+1) = x_ref + angular_quadrature_points_extra(k,1) * r
        grid_points(2,k,j,i_nucl+1) = y_ref + angular_quadrature_points_extra(k,2) * r
        grid_points(3,k,j,i_nucl+1) = z_ref + angular_quadrature_points_extra(k,3) * r

      enddo
    enddo
  enddo
end subroutine


double precision function cell_function_becke_general_adapt(r, atom_number, r_input,slater_inter_per_input)
  BEGIN_DOC
  ! 
  ! General Becke function for a list of atom supplied with an extra atom located at "r_input"
  ! and whose slater ratio slater_inter_per_input(i_nucl) is known for each nucl "i_nucl"
  ! 
  ! See Becke (1988, JCP,88(4))
  ! r(1:3)                       :: x,y,z coordinantes of the current point
  END_DOC

  implicit none
  double precision, intent(in) :: r(3)
  double precision, intent(in) :: r_input(3),slater_inter_per_input(nucl_num)
  integer         , intent(in) :: atom_number
  integer                      :: j
  double precision             :: mu_ij, nu_ij
  double precision             :: distance_i, distance_j, step_function_becke

  ! Distance of the point r from the nucleus belonging to the same Voronoi cell
  distance_i  = (r(1) - nucl_coord_transp(1,atom_number) ) * (r(1) - nucl_coord_transp(1,atom_number))
  distance_i += (r(2) - nucl_coord_transp(2,atom_number) ) * (r(2) - nucl_coord_transp(2,atom_number))
  distance_i += (r(3) - nucl_coord_transp(3,atom_number) ) * (r(3) - nucl_coord_transp(3,atom_number))
  distance_i  = dsqrt(distance_i)

  cell_function_becke_general_adapt = 1.d0
  ! Loop for the contributes coming from all atoms' Voronoi cell
  do j = 1, nucl_num
    if(j==atom_number) cycle

    distance_j  = (r(1) - nucl_coord_transp(1,j) ) * (r(1) - nucl_coord_transp(1,j))
    distance_j += (r(2) - nucl_coord_transp(2,j) ) * (r(2) - nucl_coord_transp(2,j))
    distance_j += (r(3) - nucl_coord_transp(3,j) ) * (r(3) - nucl_coord_transp(3,j))
    distance_j  = dsqrt(distance_j)

    mu_ij = (distance_i - distance_j) * nucl_dist_inv(atom_number,j)
    nu_ij = mu_ij + slater_bragg_type_inter_distance_ua(atom_number,j) * (1.d0 - mu_ij*mu_ij)

    cell_function_becke_general_adapt *= step_function_becke(nu_ij)

    write(*,'(I4,100E13.4)'), j, distance_j, mu_ij, nu_ij, cell_function_becke_general_adapt
  enddo

  ! Contribute coming from the extra atom
  distance_j  = (r(1) - r_input(1) ) * (r(1) - r_input(1))
  distance_j += (r(2) - r_input(2) ) * (r(2) - r_input(2))
  distance_j += (r(3) - r_input(3) ) * (r(3) - r_input(3))
  distance_j  = dsqrt(distance_j)

  double precision :: rij
  rij  = (r_input(1) - nucl_coord_transp(1,atom_number) ) * (r_input(1) - nucl_coord_transp(1,atom_number))
  rij += (r_input(2) - nucl_coord_transp(2,atom_number) ) * (r_input(2) - nucl_coord_transp(2,atom_number))
  rij += (r_input(3) - nucl_coord_transp(3,atom_number) ) * (r_input(3) - nucl_coord_transp(3,atom_number))
  rij = dsqrt(rij)

  mu_ij = (distance_i - distance_j) / rij
  nu_ij = mu_ij + slater_inter_per_input(atom_number) * (1.d0 - mu_ij*mu_ij)

  cell_function_becke_general_adapt *= step_function_becke(nu_ij)
 
  !write(*,'(100E13.4)'), distance_j, mu_ij, nu_ij, cell_function_becke_general_adapt

  return
end


double precision function cell_function_becke_extra_atom(r, r_input, slater_inter_per_input)
  BEGIN_DOC
  ! 
  ! Return the Becke-weight of a given point "r" with respect to a point characterized by "r_input" and slater_inter_per_input
  !
  ! 
  ! See Becke (1988, JCP,88(4))
  ! r(1:3)                       :: x,y,z coordinantes of the current point
  END_DOC

  implicit none
  double precision, intent(in) :: r(3)
  double precision, intent(in) :: r_input(3)
  double precision, intent(in) :: slater_inter_per_input(nucl_num)
  !
  integer           :: j
  double precision  :: mu_ij, nu_ij
  double precision  :: distance_i, distance_j, inv_dist
  double precision :: step_function_becke

  distance_i  = (r(1) - r_input(1) ) * (r(1) - r_input(1))
  distance_i += (r(2) - r_input(2) ) * (r(2) - r_input(2))
  distance_i += (r(3) - r_input(3) ) * (r(3) - r_input(3))
  distance_i  = dsqrt(distance_i)

  !write(*,'(100E13.5)') r, r_input, distance_i

  cell_function_becke_extra_atom = 1.d0
  do j = 1, nucl_num

    distance_j  = (r(1) - nucl_coord_transp(1,j) ) * (r(1) - nucl_coord_transp(1,j))
    distance_j += (r(2) - nucl_coord_transp(2,j) ) * (r(2) - nucl_coord_transp(2,j))
    distance_j += (r(3) - nucl_coord_transp(3,j) ) * (r(3) - nucl_coord_transp(3,j))
    distance_j  = dsqrt(distance_j)

    inv_dist  = (r_input(1) - nucl_coord_transp(1,j) ) * (r_input(1) - nucl_coord_transp(1,j))
    inv_dist += (r_input(2) - nucl_coord_transp(2,j) ) * (r_input(2) - nucl_coord_transp(2,j))
    inv_dist += (r_input(3) - nucl_coord_transp(3,j) ) * (r_input(3) - nucl_coord_transp(3,j))
    inv_dist = 1.d0/dsqrt(inv_dist)

    ! Elliptical coordinate (Becke Eq. 11)
    mu_ij = (distance_i - distance_j) * inv_dist
    ! Heteronuclear correction (Becke Eq. A2)
    nu_ij = mu_ij + slater_inter_per_input(j) * (1.d0 - mu_ij*mu_ij)

    cell_function_becke_extra_atom *= step_function_becke(nu_ij)
    !write(*,'(I4,100E13.4)'), j, distance_j, 1.d0/inv_dist, mu_ij, nu_ij, cell_function_becke_extra_atom
  enddo

  return
end function


subroutine get_becke_functions_general(r,r_input,slater_inter_per_input,weights_per_atom)
 implicit none
 BEGIN_DOC
 !
 ! Computes the Becke-weight associated to each atom, 'weights_per_atom',
 ! at the point in real space 'r'.
 ! The 3D space has been already Voronoi-partitioned in atom-centered Becke-cells
 ! plus a new one associated to a ghost atom centered in 'r_input'.
 ! The new Becke-cell is characterized by its center, r_input, and its
 ! a_ij parameter defining offset Voronoi boundaries, 'slater_inter_per_input'
 !
 !
 !
 ! in terms of nucl_num atoms + 1 additional atom characterized by "r_input" and "slater_inter_per_input"
 !
 ! returns the following partition :: weights_per_atom(i) == w(r,i) =  P(r,i) / \sum_{k in partition } P(r,k)
 !
 ! WARNING the first entry corresponds to the fake atom "r_input"
 END_DOC
 double precision, intent(in) :: r(3), r_input(3), slater_inter_per_input(nucl_num)
 double precision, intent(out):: weights_per_atom(nucl_num+1)
 double precision :: norm,cell_function_becke_extra_atom,cell_function_becke_general_adapt
 integer :: i
 !
 norm = 0.d0
 !
 ! weights_per_atom(i_nucl) = P(r,i_nucl)/sum_k P(r,k) 
 ! BUT NOW sum_k P(r,k) = P(r,new_atom) + \sum_k = 1, nucl_num  P(r,k)

 ! you do the first part of the sum corresponding to the extra atom "r_input"
 weights_per_atom(1) = cell_function_becke_extra_atom(r, r_input, slater_inter_per_input)
 norm += weights_per_atom(1)
 ! do the rest of the usual atoms
 ! For each of these points in space, ou need to evaluate the P_n(r)
 do i = 1, nucl_num
   ! function defined for each atom "i" by equation (13) and (21) with k == 3
   weights_per_atom(i+1) = cell_function_becke_general_adapt(r, i, r_input,slater_inter_per_input) ! P_n(r)
   ! Then you compute the summ the P_n(r) function for each of the "r" points
   norm += weights_per_atom(i+1)
 enddo
 do i = 1, nucl_num+1
  weights_per_atom(i) = weights_per_atom(i) / norm
 enddo
end


subroutine get_voronoi_partition(r_input,grid_points,slater_inter_per_input,weight_tmp)
  implicit none
  BEGIN_DOC
  ! Computes the Voronoi Becke-weights for a set of 'grid_points'.
  ! The Voronoi partition is made of nucl_num atoms + 1 cells, where the +1 
  ! stands for an additional spherical cell centered around 'r_input', treated
  ! as an extra atom.
  !
  !extra atom at r_input and characterized by slater_inter_per_input
  END_DOC
  double precision, intent(in) :: r_input(3)
  !double precision, intent(in) :: grid_points(3,n_points_integration_angular_adapt,n_points_radial_grid_adapt,nucl_num+1)
  double precision, intent(in) :: grid_points(3,n_points_extra_integration_angular, &
                                            & n_points_extra_radial_grid,           &
                                            & nucl_num+1)
  double precision, intent(in) :: slater_inter_per_input(nucl_num)
  !double precision, intent(out):: weight_tmp(n_points_integration_angular_adapt,n_points_radial_grid_adapt,nucl_num+1)
  double precision, intent(out):: weight_tmp(n_points_extra_integration_angular, &
                                           & n_points_extra_radial_grid,         &
                                           & nucl_num+1)
  integer :: i_nucl,k,l
  double precision :: r(3)
  double precision :: weights_per_atom(nucl_num+1)
  !
  ! Loop over all grid points fixed in space (those points of the grid that do 
  ! not change because they are referred to each real atom, not the ghost-one)
  do i_nucl = 1, nucl_num
    !do k = 1, n_points_radial_grid_adapt -1
    do k = 1, n_points_extra_radial_grid -1
      !do l = 1, n_points_integration_angular_adapt
      do l = 1, n_points_extra_integration_angular
        r(1:3) = grid_points(1:3,l,k,i_nucl+1)
        call get_becke_functions_general(r,r_input,slater_inter_per_input,weights_per_atom)
        weight_tmp(l,k,i_nucl+1) = weights_per_atom(i_nucl+1)
        !write(*,'(3I4,E13.6)') i_nucl, k, l, weight_tmp(l,k,i_nucl+1)
      enddo
    enddo
  enddo

 ! Loop over all new grid points that belong to the moving grid of the 
 ! ghost-atom centered in 'r_input'
 i_nucl = 1
 !do k = 1, n_points_radial_grid_adapt -1
 do k = 1, n_points_extra_radial_grid -1
   !do l = 1, n_points_integration_angular_adapt
   do l = 1, n_points_extra_integration_angular
     r(1:3) = grid_points(1:3,l,k,i_nucl)
     call get_becke_functions_general(r,r_input,slater_inter_per_input,weights_per_atom)
     weight_tmp(l,k,i_nucl) = weights_per_atom(i_nucl)
     !write(*,'(3I4,E13.6)') i_nucl, k, l, weight_tmp(l,k,i_nucl)
   enddo
 enddo

end subroutine
