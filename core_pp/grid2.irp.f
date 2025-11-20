 BEGIN_PROVIDER [integer, n_points_rad_grid2]
&BEGIN_PROVIDER [integer, n_points_ang_grid2]

  BEGIN_DOC
  ! n_points_rad_grid2 = number of radial grid points per atom
  ! n_points_ang_grid2 = number of angular grid points per atom
  !
  ! These numbers are automatically set by setting the grid_type_sgn parameter
  END_DOC
 
  implicit none

  if(.not. my_grid_becke2) then
    select case (grid_type_sgn)
      case(0)
        n_points_rad_grid2 = 23
        n_points_ang_grid2 = 170
      case(1)
        n_points_rad_grid2 = 50
        n_points_ang_grid2 = 194
      case(2)
        n_points_rad_grid2 = 75
        n_points_ang_grid2 = 302
      case(3)
        n_points_rad_grid2 = 99
        n_points_ang_grid2 = 590
      case default
        write(*,*) '!!! Quadrature grid not available !!!'
        stop
    end select
  else
    n_points_rad_grid2 = my_n_pt_r_grid2
    n_points_ang_grid2 = my_n_pt_a_grid2
  endif
  
  !if (verbose) then 
  !  print*, " n_points_rad_grid2         = ", n_points_rad_grid2
  !  print*, " n_points_ang_grid2 = ", n_points_ang_grid2
  write(6,*) 'Grid2 parameters:'
  call write_int(6, n_points_rad_grid2, 'N radial points ')
  call write_int(6, n_points_ang_grid2, 'N angular points') 
  write(6,*) ''      
  !end if

END_PROVIDER


BEGIN_PROVIDER [integer, n_points_grid_per_atom2]
  BEGIN_DOC
  ! Number of grid points per atom
  END_DOC
  implicit none
  n_points_grid_per_atom2 = n_points_ang_grid2 * n_points_rad_grid2
END_PROVIDER


 BEGIN_PROVIDER [double precision, grid_points_radial2, (n_points_rad_grid2)]
&BEGIN_PROVIDER [double precision, dr_radial_integral2]
  BEGIN_DOC
  ! points in [0,1] to map the radial integral [0,\infty]
  END_DOC
  implicit none
  integer :: i
  dr_radial_integral2 = 1.d0 / dble(n_points_rad_grid2-1)
  do i = 1, n_points_rad_grid2
    grid_points_radial2(i) = dble(i-1) * dr_radial_integral2
  enddo
END_PROVIDER


BEGIN_PROVIDER [double precision, grid_points_per_atom2, (3,n_points_ang_grid2,n_points_rad_grid2,nucl_num)]

  BEGIN_DOC
  ! x,y,z coordinates of grid points used for integration in 3d space
  END_DOC

  implicit none
  integer                    :: i, j, k
  double precision           :: dr, x_ref, y_ref, z_ref
  double precision           :: x, r, tmp
  double precision, external :: knowles_function

  grid_points_per_atom2 = 0.d0

  PROVIDE rad_grid_type
  if(rad_grid_type .eq. "KNOWLES") then

    do i = 1, nucl_num
      x_ref = nucl_coord(i,1)
      y_ref = nucl_coord(i,2)
      z_ref = nucl_coord(i,3)
      do j = 1, n_points_rad_grid2-1

        ! x value for the mapping of the [0, +\infty] to [0,1]
        x = grid_points_radial2(j)
        ! value of the radial coordinate for the integration
        r = knowles_function(alpha_knowles(grid_atomic_number(i)), m_knowles, x)

        ! explicit values of the grid points centered around each atom
        do k = 1, n_points_ang_grid2
          grid_points_per_atom2(1,k,j,i) = x_ref + grid2_lebedev_points(k,1) * r
          grid_points_per_atom2(2,k,j,i) = y_ref + grid2_lebedev_points(k,2) * r
          grid_points_per_atom2(3,k,j,i) = z_ref + grid2_lebedev_points(k,3) * r
        enddo
      enddo
    enddo

  elseif(rad_grid_type .eq. "GILL") then
    ! GILL & CHIEN, 2002

    do i = 1, nucl_num
      x_ref = nucl_coord(i,1)
      y_ref = nucl_coord(i,2)
      z_ref = nucl_coord(i,3)
      do j = 1, n_points_rad_grid2-1

        r = R_gill * dble(j-1)**2 / dble(n_points_rad_grid2-j+1)**2 

        ! explicit values of the grid points centered around each atom
        do k = 1, n_points_ang_grid2
          grid_points_per_atom2(1,k,j,i) = x_ref + grid2_lebedev_points(k,1) * r
          grid_points_per_atom2(2,k,j,i) = y_ref + grid2_lebedev_points(k,2) * r
          grid_points_per_atom2(3,k,j,i) = z_ref + grid2_lebedev_points(k,3) * r
        enddo
      enddo
    enddo
 
  else

    print*, " rad_grid_type = ", rad_grid_type, ' is not implemented'
    stop

  endif

END_PROVIDER


BEGIN_PROVIDER [double precision, weight_at_r2, (n_points_ang_grid2,n_points_rad_grid2,nucl_num)]

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
    do k = 1, n_points_rad_grid2 -1
      ! for each angular point attached to the "jth" atom
      do l = 1, n_points_ang_grid2
        r(1:3) = grid_points_per_atom2(1:3,l,k,j)

        accu = 0.d0
        ! For each of these points in space, ou need to evaluate the P_n(r)
        do i = 1, nucl_num
          ! function defined for each atom "i" by equation (13) and (21) with k == 3
          tmp_array(i) = cell_function_becke(r, i) ! P_n(r)
          ! Then you compute the summ the P_n(r) function for each of the "r" points
          accu += tmp_array(i)
        enddo
        accu = 1.d0/accu
        weight_at_r2(l,k,j) = tmp_array(j) * accu

        if(isnan(weight_at_r2(l,k,j))) then
          print*,'isnan(weight_at_r2(l,k,j))'
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

BEGIN_PROVIDER [double precision, final_weight_at_r2, (n_points_ang_grid2,n_points_rad_grid2,nucl_num)]

  BEGIN_DOC
  ! Total weight on each grid point which takes into account all Lebedev, Voronoi and radial weights.
  END_DOC

  implicit none
  integer                    :: i, j, k, l, m
  double precision           :: r(3)
  double precision           :: tmp_array(nucl_num)
  double precision           :: contrib_integration, x, tmp
  double precision, external :: derivative_knowles_function, knowles_function

  final_weight_at_r2 = 0.d0

  PROVIDE rad_grid_type
  if(rad_grid_type .eq. "KNOWLES") then

    ! run over all points in space
    do j = 1, nucl_num  ! that are referred to each atom
      do i = 1, n_points_rad_grid2 -1 !for each radial grid attached to the "jth" atom
        x = grid_points_radial2(i) ! x value for the mapping of the [0, +\infty] to [0,1]

        do k = 1, n_points_ang_grid2  ! for each angular point attached to the "jth" atom
          contrib_integration = derivative_knowles_function(alpha_knowles(grid_atomic_number(j)), m_knowles, x) &
                              * knowles_function(alpha_knowles(grid_atomic_number(j)), m_knowles, x)**2

          final_weight_at_r2(k,i,j) = grid2_lebedev_weights(k)  * weight_at_r2(k,i,j) * contrib_integration * dr_radial_integral2

          if(isnan(final_weight_at_r2(k,i,j))) then
           print*,'isnan(final_weight_at_r2(k,i,j))' 
           print*,k,i,j
           write(*,'(100(F16.10,X))') grid2_lebedev_weights(k), weight_at_r2(k,i,j), contrib_integration
           stop 
          endif
        enddo
      enddo
    enddo

  elseif(rad_grid_type .eq. "GILL") then
    ! GILL & CHIEN, 2002

    tmp = 2.d0 * R_gill * R_gill * R_gill * dble(n_points_rad_grid2)

    ! run over all points in space
    do j = 1, nucl_num  ! that are referred to each atom
      do i = 1, n_points_rad_grid2 - 1 !for each radial grid attached to the "jth" atom
        contrib_integration = tmp * dble(i-1)**5 / dble(n_points_rad_grid2-i+1)**7
        do k = 1, n_points_ang_grid2  ! for each angular point attached to the "jth" atom
          final_weight_at_r2(k,i,j) = grid2_lebedev_weights(k) * weight_at_r2(k,i,j) * contrib_integration

          if(isnan(final_weight_at_r2(k,i,j))) then
           print*,'isnan(final_weight_at_r2(k,i,j))' 
           print*,k,i,j
           write(*,'(100(F16.10,X))') grid2_lebedev_weights(k), weight_at_r2(k,i,j), contrib_integration, dr_radial_integral2
           stop 
          endif
        enddo
      enddo
    enddo

  else

    print*, " rad_grid_type = ", rad_grid_type, ' is not implemented'
    stop

  endif

END_PROVIDER


! this is the vectorial version
BEGIN_PROVIDER [integer, n_points_final_grid2]

  BEGIN_DOC
  ! Number of points which are non zero
  END_DOC

  implicit none
  integer :: i, j, k, l


  write(6,*) 'PROVIDING n_points_final_grid2 ...'
  write(6,*) '... through pruning of PROVIDER final_weight_at_r2'

  n_points_final_grid2 = 0
  do j = 1, nucl_num
    do i = 1, n_points_rad_grid2 -1
      do k = 1, n_points_ang_grid2
        !if(dabs(final_weight_at_r(k,i,j)) < thresh_grid)then
        if(dabs(final_weight_at_r2(k,i,j)) < thresh_grid)then
          cycle
        endif
        n_points_final_grid2 += 1
      enddo
    enddo
  enddo
  call write_int(6, n_points_ang_grid2*(n_points_rad_grid2*nucl_num-1), 'N points before pruning')
  call write_int(6, n_points_final_grid2, 'N points after pruning')

!!!  print*,' n_points_final_grid2 = ', n_points_final_grid2
!!!  print*,' n max point         = ', n_points_ang_grid2*(n_points_rad_grid2*nucl_num - 1)
!!! ! no reason to write in the EZFIO file the number of grid points ?
!!!!  call ezfio_set_becke_numerical_grid_n_points_final_grid2(n_points_final_grid2)
!!!
END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, final_grid_points2,      (3,n_points_final_grid2)]
&BEGIN_PROVIDER [double precision, final_weight_at_r_vector2, (n_points_final_grid2)]
&BEGIN_PROVIDER [integer,          index_final_points2,     (3,n_points_final_grid2)]
&BEGIN_PROVIDER [integer, index_final_points_reverse2, (n_points_ang_grid2,n_points_rad_grid2,nucl_num)]

  BEGIN_DOC
  !  final_grid_points2(1:3,j) = (/ x, y, z /) of the jth grid point
  !
  ! final_weight_at_r_vector2(i) = Total weight function of the ith grid point which contains the Lebedev, Voronoi and radial weights contributions
  !
  ! index_final_points2(1:3,i) = gives the angular, radial and atomic indices associated to the ith grid point
  !
  ! index_final_points_reverse2(i,j,k) = index of the grid point having i as angular, j as radial and l as atomic indices
  END_DOC

  implicit none
  integer          :: i, j, k, l, i_count
  double precision :: r(3)
  double precision :: wall0, wall1

  call wall_time(wall0)
  print *, ' Providing final_grid_points2 ...'

  i_count = 0
  do j = 1, nucl_num
    do i = 1, n_points_rad_grid2 -1
      do k = 1, n_points_ang_grid2
        if(dabs(final_weight_at_r2(k,i,j)) < thresh_grid) then
          cycle
        endif
        i_count += 1
        final_grid_points2(1:3,i_count) = grid_points_per_atom2(1:3,k,i,j)

        final_weight_at_r_vector2(i_count) = final_weight_at_r2(k,i,j)
        index_final_points2(1,i_count) = k
        index_final_points2(2,i_count) = i
        index_final_points2(3,i_count) = j
        index_final_points_reverse2(k,i,j) = i_count
      enddo
    enddo
  enddo

!  FREE grid_points_per_atom
!  FREE final_weight_at_r2

  call wall_time(wall1)
  print *, ' wall time for final_grid_points2,', wall1 - wall0
  call print_memory_usage()

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, final_grid_points_transp2, (n_points_final_grid2,3)]

  BEGIN_DOC
  ! Transposed final_grid_points2
  END_DOC

  implicit none
  integer :: i,j

  do j = 1, 3
    do i = 1, n_points_final_grid2
      final_grid_points_transp2(i,j) = final_grid_points2(j,i)
    enddo
  enddo

END_PROVIDER

