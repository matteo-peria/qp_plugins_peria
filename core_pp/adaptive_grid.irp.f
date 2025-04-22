subroutine give_adapt_grid_at_r(r_input, grid_points, final_adapt_weights,n_atoms_max)
 implicit none
 BEGIN_DOC
 ! grid_points(1:3,i,j,k) = x/y/z for the ith angular, jth radial points attached to the kth atom
 ! 
 ! but we add an extra atom corresponding, LOCATED FIRST in the list, located at r_input
 !
 ! n_atoms_max is the maximum number of atoms to be considered in the list
 END_DOC
 double precision, intent(in) :: r_input(3)
 double precision, intent(out):: grid_points(3,n_points_angular_grid_adapt,n_points_radial_grid_adapt,nucl_num+1)
 double precision, intent(out):: final_adapt_weights(n_points_angular_grid_adapt,n_points_radial_grid_adapt,nucl_num+1)
 integer, intent(out)         :: n_atoms_max

 double precision ::  alpha_av,radii_ua_av,slater_inter_per_atom(nucl_num)
 double precision :: weight_tmp(n_points_angular_grid_adapt,n_points_radial_grid_adapt,nucl_num+1)
 double precision :: weights_per_atom(nucl_num+1),norm,x,derivative_knowles_function,knowles_function,contrib_integration
 integer :: i,k,i_nucl
!! TO BE CHANGED FOR PRUNING !!
 n_atoms_max = nucl_num+1
!! determines the alpha averaged based on the distance
 call get_info_at_r_adapt_grid(r_input,alpha_av,radii_ua_av,slater_inter_per_atom)
 call get_all_points_at_r_adapt_grid(r_input,alpha_av,grid_points)
 call get_voronoi_partition(r_input,grid_points,slater_inter_per_atom,weight_tmp)

 i_nucl = 1
   do i = 1, n_points_radial_grid_adapt -1 !for each radial grid attached to the "jth" atom
     x = grid_adapt_points_radial(i) ! x value for the mapping of the [0, +\infty] to [0,1]

     do k = 1, n_points_angular_grid_adapt  ! for each angular point attached to the "i_nucl-th" atom
       contrib_integration = derivative_knowles_function(alpha_av, m_knowles, x) &
                           * knowles_function(alpha_av, m_knowles, x)**2

       final_adapt_weights(k,i,i_nucl) = weights_angular_points(k)  * weight_tmp(k,i,i_nucl) * contrib_integration * dr_adapt_radial_integral

!       if(isnan(final_adapt_weights(final_adapt_pointsk,i,i_nucl))) then
!        print*,'isnan(final_weight_at_r_adapt(k,i,i_nucl))' 
!        print*,k,i,i_nucl
!        write(*,'(100(F16.10,X))') weights_angular_points(k), weight_tmp(k,i,i_nucl), contrib_integration
!        stop 
!       endif
     enddo
   enddo
 
 ! run over all points in space
 do i_nucl = 2, nucl_num+1 ! that are referred to each atom
   do i = 1, n_points_radial_grid_adapt -1 !for each radial grid attached to the "jth" atom
     x = grid_adapt_points_radial(i) ! x value for the mapping of the [0, +\infty] to [0,1]

     do k = 1, n_points_angular_grid_adapt  ! for each angular point attached to the "i_nucl-th" atom
       contrib_integration = derivative_knowles_function(alpha_knowles(grid_atomic_number(i_nucl-1)), m_knowles, x) &
                           * knowles_function(alpha_knowles(grid_atomic_number(i_nucl-1)), m_knowles, x)**2

       final_adapt_weights(k,i,i_nucl) = weights_angular_points(k)  * weight_tmp(k,i,i_nucl) * contrib_integration * dr_adapt_radial_integral

!       if(isnan(final_adapt_points(final_adapt_pointsk,i,i_nucl))) then
!        print*,'isnan(final_weight_at_r_adapt(k,i,i_nucl))' 
!        print*,k,i,i_nucl
!        write(*,'(100(F16.10,X))') weights_angular_points(k), weight_tmp(k,i,i_nucl), contrib_integration
!        stop 
!       endif
     enddo
   enddo
 enddo

end

! ---

 BEGIN_PROVIDER [integer, n_points_radial_grid_adapt]
&BEGIN_PROVIDER [integer, n_points_angular_grid_adapt]

  BEGIN_DOC
  ! n_points_radial_grid_adapt = number of radial grid points per atom
  !
  ! n_points_angular_grid_adapt = number of angular grid points per atom
  !
  ! These numbers are automatically set by setting the grid_adapt_type_sgn parameter
  END_DOC
 
  implicit none

!  if(.not. my_grid_adapt_adapt_becke) then
        n_points_radial_grid_adapt = 23
        n_points_angular_grid_adapt = 170
!      case(1)
!        n_points_radial_grid_adapt = 50
!        n_points_angular_grid_adapt = 194
!      case(2)
!        n_points_radial_grid_adapt = 75
!        n_points_angular_grid_adapt = 302
!      case(3)
!        n_points_radial_grid_adapt = 99
!        n_points_angular_grid_adapt = 590
!      case default
!        write(*,*) '!!! Quadrature grid not available !!!'
!        stop
!    end select
!  else
!    n_points_radial_grid_adapt         = my_n_pt_r_grid_adapt
!    n_points_angular_grid_adapt = my_n_pt_a_grid_adapt
!  endif

  print*, " n_points_radial_grid_adapt         = ", n_points_radial_grid_adapt
  print*, " n_points_angular_grid_adapt = ", n_points_angular_grid_adapt

END_PROVIDER

! ---
BEGIN_PROVIDER [integer, n_total_adapt_grid]
 implicit none
 n_total_adapt_grid = n_points_grid_adapt_adapt_per_atom + n_points_angular_grid_adapt * n_points_radial_grid_adapt
END_PROVIDER 

BEGIN_PROVIDER [integer, n_points_grid_adapt_adapt_per_atom]

  BEGIN_DOC
  ! Number of grid points per atom
  END_DOC

  implicit none

  n_points_grid_adapt_adapt_per_atom = n_points_angular_grid_adapt * n_points_radial_grid_adapt

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, grid_adapt_points_radial, (n_points_radial_grid_adapt)]
&BEGIN_PROVIDER [double precision, dr_adapt_radial_integral]

  BEGIN_DOC
  ! points in [0,1] to map the radial integral [0,\infty]
  END_DOC

  implicit none
  integer :: i

  dr_adapt_radial_integral = 1.d0 / dble(n_points_radial_grid_adapt-1)

  do i = 1, n_points_radial_grid_adapt
    grid_adapt_points_radial(i) = dble(i-1) * dr_adapt_radial_integral
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, grid_adapt_points_per_atom, (3,n_points_angular_grid_adapt,n_points_radial_grid_adapt,nucl_num)]

  BEGIN_DOC
  ! x,y,z coordinates of grid points used for integration in 3d space
  END_DOC

  implicit none
  integer                    :: i, j, k
  double precision           :: dr, x_ref, y_ref, z_ref
  double precision           :: x, r, tmp
  double precision, external :: knowles_function

  grid_adapt_points_per_atom = 0.d0


    do i = 1, nucl_num
      x_ref = nucl_coord(i,1)
      y_ref = nucl_coord(i,2)
      z_ref = nucl_coord(i,3)
      do j = 1, n_points_radial_grid_adapt-1

        ! x value for the mapping of the [0, +\infty] to [0,1]
        x = grid_adapt_points_radial(j)
        ! value of the radial coordinate for the integration
        r = knowles_function(alpha_knowles(grid_atomic_number(i)), m_knowles, x)

        ! explicit values of the grid points centered around each atom
        do k = 1, n_points_angular_grid_adapt
          grid_adapt_points_per_atom(1,k,j,i) = x_ref + angular_quadrature_points(k,1) * r
          grid_adapt_points_per_atom(2,k,j,i) = y_ref + angular_quadrature_points(k,2) * r
          grid_adapt_points_per_atom(3,k,j,i) = z_ref + angular_quadrature_points(k,3) * r
        enddo
      enddo
    enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, weight_at_r_adapt, (n_points_angular_grid_adapt,n_points_radial_grid_adapt,nucl_num)]
&BEGIN_PROVIDER [double precision, interm_weight_at_r_adapt, (n_points_angular_grid_adapt,n_points_radial_grid_adapt,nucl_num)]

  BEGIN_DOC
  ! Weight function at grid points : w_n(r) according to the equation (22)
  ! of Becke original paper (JCP, 88, 1988)
  !
  ! The "n" discrete variable represents the nucleis which in this array is
  ! represented by the last dimension and the points are labelled by the
  ! other dimensions.
  END_DOC

  implicit none
  integer                    :: i, j, k, l, m
  double precision           :: r(3), accu
  double precision           :: tmp_array(nucl_num)
  double precision, external :: cell_function_becke

  ! run over all points in space
  ! that are referred to each atom
  do j = 1, nucl_num
    !for each radial grid attached to the "jth" atom
    do k = 1, n_points_radial_grid_adapt -1
      ! for each angular point attached to the "jth" atom
      do l = 1, n_points_angular_grid_adapt
        r(1) = grid_adapt_points_per_atom(1,l,k,j)
        r(2) = grid_adapt_points_per_atom(2,l,k,j)
        r(3) = grid_adapt_points_per_atom(3,l,k,j)

        accu = 0.d0
        ! For each of these points in space, ou need to evaluate the P_n(r)
        do i = 1, nucl_num
          ! function defined for each atom "i" by equation (13) and (21) with k == 3
          tmp_array(i) = cell_function_becke(r, i) ! P_n(r)
          ! Then you compute the summ the P_n(r) function for each of the "r" points
          accu += tmp_array(i)
        enddo
        interm_weight_at_r_adapt(l,k,j) = accu
        accu = 1.d0/accu
        weight_at_r_adapt(l,k,j) = tmp_array(j) * accu

        if(isnan(weight_at_r_adapt(l,k,j))) then
          print*,'isnan(weight_at_r_adapt(l,k,j))'
          print*,l,k,j
          accu = 0.d0
          do i = 1, nucl_num
            ! function defined for each atom "i" by equation (13) and (21) with k == 3
            tmp_array(i) = cell_function_becke(r,i) ! P_n(r)
            print*,i,tmp_array(i)
            ! Then you compute the summ the P_n(r) function for each of the "r" points
            accu += tmp_array(i)
          enddo
          write(*,'(100(F16.10,X))')tmp_array(j) , accu
          stop
        endif
      enddo
    enddo
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, final_weight_at_r_adapt, (n_points_angular_grid_adapt,n_points_radial_grid_adapt,nucl_num)]

  BEGIN_DOC
  ! Total weight on each grid point which takes into account all Lebedev, Voronoi and radial weights.
  END_DOC

  implicit none
  integer                    :: i, j, k, l, m
  double precision           :: r(3)
  double precision           :: tmp_array(nucl_num)
  double precision           :: contrib_integration, x, tmp
  double precision, external :: derivative_knowles_function, knowles_function

  final_weight_at_r_adapt = 0.d0

    ! run over all points in space
    do j = 1, nucl_num  ! that are referred to each atom
      do i = 1, n_points_radial_grid_adapt -1 !for each radial grid attached to the "jth" atom
        x = grid_adapt_points_radial(i) ! x value for the mapping of the [0, +\infty] to [0,1]

        do k = 1, n_points_angular_grid_adapt  ! for each angular point attached to the "jth" atom
          contrib_integration = derivative_knowles_function(alpha_knowles(grid_atomic_number(j)), m_knowles, x) &
                              * knowles_function(alpha_knowles(grid_atomic_number(j)), m_knowles, x)**2

          final_weight_at_r_adapt(k,i,j) = weights_angular_points(k)  * weight_at_r_adapt(k,i,j) * contrib_integration * dr_adapt_radial_integral

          if(isnan(final_weight_at_r_adapt(k,i,j))) then
           print*,'isnan(final_weight_at_r_adapt(k,i,j))' 
           print*,k,i,j
           write(*,'(100(F16.10,X))') weights_angular_points(k), weight_at_r_adapt(k,i,j), contrib_integration
           stop 
          endif
        enddo
      enddo
    enddo

END_PROVIDER


