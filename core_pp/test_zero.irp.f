program test_zero
 implicit none
 integer :: i_nucl, j,k,n_atoms_max
 double precision, allocatable :: grid_points(:,:,:,:),weights_points(:,:,:)
 double precision :: r_input(3),r_gauss(3),r(3),weights,alpha
 allocate(grid_points(3,n_points_integration_angular_adapt,n_points_radial_grid_adapt,nucl_num+1))
 allocate(weights_points(n_points_integration_angular_adapt,n_points_radial_grid_adapt,nucl_num+1))
 r_gauss = 0.d0
 r_gauss(1) = 5.5d0
 r_gauss(2) = 5.75d0
 alpha = 1.d0 
 r_input = r_gauss
 call give_adapt_grid_at_r(r_input, grid_points, weights_points,n_atoms_max)
 double precision :: integral,dist
 integral = 0.d0
 do i_nucl = 1, nucl_num+1
  do j = 1, n_points_radial_grid_adapt
   do k = 1, n_points_integration_angular_adapt
    r(1:3) = grid_points(1:3,k,j,i_nucl)
    weights = weights_points(k,j,i_nucl)
    dist = dsqrt((r(1)-r_gauss(1))**2 + (r(2)-r_gauss(2))**2 + (r(3)-r_gauss(3))**2)
    integral += weights * dexp(-alpha*dist*dist)
   enddo
  enddo
 enddo
 print*,'integral     = ',integral
 print*,'pi/alpha^3/2 = ',(dacos(-1.d0)/alpha)**(1.5)
 integer :: i
 print*,'Usual grid'
 integral = 0.D0
 do i = 1,  n_points_final_grid
  ! you get x, y and z of the ith grid point
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  weights = final_weight_at_r_vector(i)
  dist = dsqrt((r(1)-r_gauss(1))**2 + (r(2)-r_gauss(2))**2 + (r(3)-r_gauss(3))**2)
  integral += weights * dexp(-alpha*dist*dist)
 enddo
 print*,'integral     = ',integral


end
