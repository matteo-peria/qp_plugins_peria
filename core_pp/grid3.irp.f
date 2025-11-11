! Temporary grid needed in 3-nested integrals

 BEGIN_PROVIDER [integer, grid3_rad_size]
&BEGIN_PROVIDER [integer, grid3_ang_size]
  implicit none
  BEGIN_DOC
  ! grid3_rad_size = number of radial grid3 points per atom
  ! grid3_ang_size = number of angular grid3 points per atom
  !
  ! These numbers are automatically set by setting the grid_type_sgn parameter
  END_DOC
 
  if(.not. my_grid3) then
    ! If we do not impose the radial/angular resolution
    ! these are retrieved from the extra grid
    grid3_rad_size = n_points_extra_radial_grid
    grid3_ang_size = n_points_extra_integration_angular
  else
    ! Otherwise, read from EZFIO file
    grid3_rad_size = my_grid3_rad_size
    grid3_ang_size = my_grid3_ang_size
  endif
  
  write(6,*) 'Grid parameters:'
  call write_int(6, grid3_rad_size, 'Grid 3 radial resolution ')
  call write_int(6, grid3_ang_size, 'Grid 3 angular resolution') 
  write(6,*) ''      
END_PROVIDER


BEGIN_PROVIDER [integer, grid3_nuc_size]
  implicit none
  BEGIN_DOC
  ! Number of grid3 points per atom
  END_DOC
  grid3_nuc_size = grid3_ang_size * grid3_rad_size
END_PROVIDER


 BEGIN_PROVIDER [double precision, grid3_rad, (grid3_rad_size)]
&BEGIN_PROVIDER [double precision, dr3]
  BEGIN_DOC
  ! points in [0,1] to map the radial integral [0,\infty]
  END_DOC
  implicit none
  integer :: i

  dr3 = 1.d0 / dble(grid3_rad_size-1)

  do i = 1, grid3_rad_size
    grid3_rad(i) = dble(i-1) * dr3
  enddo

END_PROVIDER


BEGIN_PROVIDER [double precision, grid3_molec_r, (3,grid3_ang_size,grid3_rad_size,nucl_num)]
  implicit none
  BEGIN_DOC
  ! x,y,z coordinates of grid points used for integration in 3d space
  END_DOC
  integer                    :: i, j, k
  double precision           :: pos(3), dr!, x_ref, y_ref, z_ref
  double precision           :: x, r, tmp
  double precision, external :: knowles_function

  grid3_molec_r = 0.d0

  do i = 1, nucl_num
    pos(1:3) = nucl_coord(i,1:3)
    !x_ref = nucl_coord(i,1)
    !y_ref = nucl_coord(i,2)
    !z_ref = nucl_coord(i,3)
    do j = 1, grid3_rad_size-1
      ! x value for the mapping of the [0, +\infty] to [0,1]
      x = grid3_rad(j)
      ! value of the radial coordinate for the integration
      r = knowles_function(alpha_knowles(grid_atomic_number(i)), m_knowles, x)
      ! explicit values of the grid points centered around each atom
      do k = 1, grid3_ang_size
        grid3_molec_r(1:3,k,j,i) = pos(1:3) + grid3_lebedev_points(k,1:3) * r
        !grid3_molec_r(1,k,j,i) = x_ref + grid3_lebedev_points(k,1) * r
        !grid3_molec_r(2,k,j,i) = y_ref + grid3_lebedev_points(k,2) * r
        !grid3_molec_r(3,k,j,i) = z_ref + grid3_lebedev_points(k,3) * r
      enddo
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER [double precision, grid3_molec_qw, (grid3_ang_size,grid3_rad_size,nucl_num)]
  implicit none
  BEGIN_DOC
  ! Quadrature weights
  END_DOC

  integer                    :: i, j, k, l, m
  double precision           :: r(3), accu
  double precision           :: tmp_array(nucl_num)
  double precision, external :: cell_function_becke

  ! run over all points in space
  ! that are referred to each atom
  do j = 1, nucl_num
    !for each radial grid attached to the "jth" atom
    do k = 1, grid3_rad_size -1
      ! for each angular point attached to the "jth" atom
      do l = 1, grid3_ang_size
        r(1:3) = grid3_molec_r(1:3,l,k,j)
        accu = 0.d0
        ! For each of these points in space, ou need to evaluate the P_n(r)
        do i = 1, nucl_num
          ! function defined for each atom "i" by equation (13) and (21) with k == 3
          tmp_array(i) = cell_function_becke(r, i) ! P_n(r)
          ! Then you compute the summ the P_n(r) function for each of the "r" points
          accu += tmp_array(i)
        enddo
        accu = 1.d0/accu
        grid3_molec_qw(l,k,j) = tmp_array(j) * accu

        if(isnan(grid3_molec_qw(l,k,j))) then
          print*,'isnan(grid3_molec_qw(l,k,j))'
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


BEGIN_PROVIDER [double precision, grid3_molec_w, (grid3_ang_size,grid3_rad_size,nucl_num)]

  BEGIN_DOC
  ! Total weight on each grid point which takes into account all Lebedev, Voronoi and radial weights.
  END_DOC

  implicit none
  integer                    :: i, j, k, l, m
  double precision           :: r(3)
  double precision           :: tmp_array(nucl_num)
  double precision           :: contrib_integration, x, tmp
  double precision, external :: derivative_knowles_function, knowles_function

  grid3_molec_w = 0.d0

  ! run over all points in space
  do j = 1, nucl_num  ! that are referred to each atom
    do i = 1, grid3_rad_size -1 !for each radial grid attached to the "jth" atom
      x = grid3_rad(i) ! x value for the mapping of the [0, +\infty] to [0,1]
      do k = 1, grid3_ang_size  ! for each angular point attached to the "jth" atom
        contrib_integration = derivative_knowles_function(alpha_knowles(grid_atomic_number(j)), m_knowles, x) &
                            * knowles_function(alpha_knowles(grid_atomic_number(j)), m_knowles, x)**2

        grid3_molec_w(k,i,j) = grid3_lebedev_weights(k)  * grid3_molec_qw(k,i,j) * contrib_integration * dr3

        if(isnan(grid3_molec_w(k,i,j))) then
         print*,'isnan(grid3_molec_w(k,i,j))' 
         print*,k,i,j
         write(*,'(100(F16.10,X))') grid3_lebedev_weights(k), grid3_molec_qw(k,i,j), contrib_integration
         stop 
        endif
      enddo
    enddo
  enddo
END_PROVIDER


! PRUNED GRID PROVIDERS


BEGIN_PROVIDER [integer, grid3_pruned_size]
  implicit none
  BEGIN_DOC
  ! Number of points which are non zero
  END_DOC

  integer :: i, j, k, l

  write(6,*) 'PROVIDING grid3_pruned_size ...'
  write(6,*) '... through pruning of PROVIDER grid3_molec_w'

  grid3_pruned_size = 0
  do j = 1, nucl_num
    do i = 1, grid3_rad_size -1
      do k = 1, grid3_ang_size
        if(dabs(grid3_molec_w(k,i,j)) < thresh_grid3) then
          cycle
        endif
        grid3_pruned_size += 1
      enddo
    enddo
  enddo
  call write_int(6, grid3_ang_size*(grid3_rad_size*nucl_num-1), 'N points before pruning')
  call write_int(6, grid3_pruned_size, 'N points after pruning')

!!!  print*,' grid3_pruned_size = ', grid3_pruned_size
!!!  print*,' n max point         = ', grid3_ang_size*(grid3_rad_size*nucl_num - 1)
!!! ! no reason to write in the EZFIO file the number of grid points ?
!!!!  call ezfio_set_becke_numerical_grid_grid3_pruned_size(grid3_pruned_size)
!!!
END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, grid3_prune_r,     (3,grid3_pruned_size)]
&BEGIN_PROVIDER [double precision, grid3_prune_w,     (grid3_pruned_size)]
&BEGIN_PROVIDER [integer,          grid3_prune_index, (3,grid3_pruned_size)]
&BEGIN_PROVIDER [integer,  grid3_prune_index_reverse, (grid3_ang_size,grid3_rad_size,nucl_num)]

  BEGIN_DOC
  !  grid3_prune_r(1:3,j) = (/ x, y, z /) of the jth grid point
  !
  ! grid3_prune_w(i) = Total weight function of the ith grid point which contains the Lebedev, Voronoi and radial weights contributions
  !
  ! grid3_prune_index(1:3,i) = gives the angular, radial and atomic indices associated to the ith grid point
  !
  ! grid3_prune_index_reverse(i,j,k) = index of the grid point having i as angular, j as radial and l as atomic indices
  END_DOC

  implicit none
  integer          :: i, j, k, l, i_count
  double precision :: r(3)
  double precision :: wall0, wall1

  call wall_time(wall0)
  print *, ' Providing grid3_prune_r ...'

  i_count = 0
  do j = 1, nucl_num
    do i = 1, grid3_rad_size -1
      do k = 1, grid3_ang_size
        if(dabs(grid3_molec_w(k,i,j)) < thresh_grid3) then
          cycle
        endif
        i_count += 1
        grid3_prune_r(1:3,i_count) = grid3_molec_r(1:3,k,i,j)
        grid3_prune_w(i_count) = grid3_molec_w(k,i,j)
        grid3_prune_index(1,i_count) = k
        grid3_prune_index(2,i_count) = i
        grid3_prune_index(3,i_count) = j
        grid3_prune_index_reverse(k,i,j) = i_count
      enddo
    enddo
  enddo

!  FREE grid3_prune_r
!  FREE grid3_prune_w

  call wall_time(wall1)
  print *, ' wall time for grid3_prune_r,', wall1 - wall0
  call print_memory_usage()

END_PROVIDER


BEGIN_PROVIDER [double precision, grid3_prune_r_transp, (grid3_pruned_size,3)]
  implicit none
  BEGIN_DOC
  ! Transposed grid3_prune_r
  END_DOC

  integer :: i,j

  do j = 1, 3
    do i = 1, grid3_pruned_size
      grid3_prune_r_transp(i,j) = grid3_prune_r(j,i)
    enddo
  enddo

END_PROVIDER
