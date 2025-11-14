 BEGIN_PROVIDER[double precision, aos_in_r_array_extra_full, (ao_num, n_points_extra_integration_angular,  n_points_extra_radial_grid, nucl_num)]
  implicit none
  BEGIN_DOC
  ! aos_in_r_array_extra_full(i,j) 
  ! = value of the i-th AO on the jth grid point of the EXTRA grid
  END_DOC
  integer :: i_nuc, i_rad, i_ang
  double precision :: r(3)
 
  !$OMP  PARALLEL DEFAULT(NONE) &
  !$OMP& PRIVATE(i_nuc,i_rad,i_ang,r) &
  !$OMP& SHARED(nucl_num, n_points_extra_radial_grid, n_points_extra_integration_angular, &
  !$OMP&        grid_points_extra_per_atom, aos_in_r_array_extra_full)
  !$OMP  DO COLLAPSE(3) SCHEDULE(static)
 
  do i_nuc = 1, nucl_num
    do i_rad = 1, n_points_extra_radial_grid - 1
      do i_ang = 1, n_points_extra_integration_angular
        r(1:3) = grid_points_extra_per_atom(1:3,i_ang,i_rad,i_nuc)
        call give_all_aos_at_r(r, aos_in_r_array_extra_full(1,i_ang,i_rad,i_nuc))
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL 
 END_PROVIDER


 BEGIN_PROVIDER[double precision, mos_in_r_array_extra_full_omp, (mo_num, n_points_extra_integration_angular,  n_points_extra_radial_grid, nucl_num)]
  implicit none
  BEGIN_DOC
  ! mos_in_r_array_extra(i,j)        = value of the ith mo on the jth grid point on the EXTRA GRID 
  END_DOC
  integer :: i_nuc, i_rad, i_ang
  integer :: j
  double precision :: mos_array_extra(mo_num), r(3)

  !$OMP  PARALLEL DEFAULT(NONE) &
  !$OMP& PRIVATE(i_nuc,i_rad,i_ang,j,r,mos_array_extra) &
  !$OMP& SHARED(nucl_num, n_points_extra_radial_grid, n_points_extra_integration_angular, &
  !$OMP&        grid_points_extra_per_atom, mo_num, mos_in_r_array_extra_full_omp)
  !$OMP DO COLLAPSE(3) SCHEDULE(static)
  do i_nuc = 1, nucl_num
    do i_rad = 1, n_points_extra_radial_grid - 1
      do i_ang = 1, n_points_extra_integration_angular
        r(1:3) = grid_points_extra_per_atom(1:3, i_ang, i_rad, i_nuc)
        ! can't we just
        call give_all_mos_at_r(r, mos_in_r_array_extra_full_omp(1, i_ang, i_rad, i_nuc))
        ! instead of doing the loop 
        !call give_all_mos_at_r(r, mos_array_extra)
        !do j = 1, mo_num
        !  mos_in_r_array_extra_full_omp(j, i_ang, i_rad, i_nuc) = mos_array_extra(j)
        !end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

 END_PROVIDER

