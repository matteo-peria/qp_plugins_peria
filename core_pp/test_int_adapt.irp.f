program test
 implicit none
 BEGIN_DOC
 ! int dr 1/|r-C| AO_i(r) AO_j(r) = \int dr 1/r AO_i(r+C) AO_j(r+C)
 END_DOC
 include 'constants.include.F'
 integer :: i,j,k,l,ipoint,i_ao, j_ao
 double precision :: alpha,center(3)
 double precision :: weight, r(3)
 double precision :: distance, f_r, mu
 double precision, allocatable :: aos_tmp(:)
 double precision, allocatable :: integrals_ao(:,:),integral_1(:,:), integral_2(:,:)
 allocate(integrals_ao(ao_num, ao_num),integral_1(ao_num,ao_num), integral_2(ao_num,ao_num))
 allocate(aos_tmp(ao_num))
 ! center 
 center(1:3)=0.d0
 center(1) = 1.75d0
 center(2) = 1.5d0
 integral_1 = 0.d0
 ! you browse all the grid points as a one-dimensional array
 do i = 1,  n_points_final_grid
  ! you get x, y and z of the ith grid point
  r(1) = final_grid_points(1,i)
  r(2) = final_grid_points(2,i)
  r(3) = final_grid_points(3,i)
  weight = final_weight_at_r_vector(i)
  ! you compute the function to be integrated
  distance = dsqrt( (r(1) - center(1))**2 +  (r(2) - center(2))**2 + (r(3) - center(3))**2 )
  call give_all_aos_at_r(r, aos_tmp)
  do i_ao = 1, ao_num
   do j_ao = 1, ao_num
    if(distance.ne.0.d0)then
     f_r = aos_tmp(i_ao) * aos_tmp(j_ao) / distance
    endif
    ! you add the contribution of the grid point to the integral
    integral_1(j_ao,i_ao) += f_r * weight
   enddo
  enddo
 enddo
 print*,'integral_1 = '
 do i_ao = 1, ao_num
  write(*,'(100(F16.10,X))')integral_1(i_ao,:)
 enddo

 integral_2 = 0.d0
 integer :: i_nucl, n_atoms_max
 double precision, allocatable :: grid_points(:,:,:,:),weights_points(:,:,:)
 allocate(grid_points(3,n_points_integration_angular_adapt,n_points_radial_grid_adapt,nucl_num+1))
 allocate(weights_points(n_points_integration_angular_adapt,n_points_radial_grid_adapt,nucl_num+1))
 call give_adapt_grid_at_r(center, grid_points, weights_points,n_atoms_max)
 do i_nucl = 1, nucl_num+1
  do j = 1, n_points_radial_grid_adapt
   do k = 1, n_points_integration_angular_adapt
    r(1:3) = grid_points(1:3,k,j,i_nucl)
    weight = weights_points(k,j,i_nucl)
    distance = dsqrt( (r(1) - center(1))**2 +  (r(2) - center(2))**2 + (r(3) - center(3))**2 )
    call give_all_aos_at_r(r, aos_tmp)
    do i_ao = 1, ao_num
     do j_ao = 1, ao_num
      if(distance.ne.0.d0)then
       f_r = aos_tmp(i_ao) * aos_tmp(j_ao) / distance
      endif
      ! you add the contribution of the grid point to the integral
      integral_2(j_ao,i_ao) += f_r * weight
     enddo
    enddo
   enddo
  enddo
 enddo

 print*,'integral_2 = '
 do i_ao = 1, ao_num
  write(*,'(100(F16.10,X))')integral_2(i_ao,:)
 enddo
 mu = 100000000.d0
 call give_all_erf_kl_ao(integrals_ao,mu,center)
 print*,'integrals_ao = '
 do i_ao = 1, ao_num
  write(*,'(100(F16.10,X))')integrals_ao(i_ao,:)
 enddo
 double precision :: accu_1, accu_2
 accu_1 = 0.d0
 accu_2 = 0.d0
 do i_ao = 1, ao_num
  do j_ao = 1, ao_num
   accu_1 += dabs(integral_1(j_ao, i_ao) - integrals_ao(j_ao,i_ao))
   accu_2 += dabs(integral_2(j_ao, i_ao) - integrals_ao(j_ao,i_ao))
  enddo
 enddo

 print*,'accu_1 = ',accu_1
 print*,'accu_2 = ',accu_2


end
