 BEGIN_PROVIDER[double precision, aos_in_r_array_extra_full, (ao_num, n_points_extra_integration_angular,  n_points_extra_radial_grid, nucl_num)]
  implicit none
  BEGIN_DOC
  ! aos_in_r_array_extra(i,j)        = value of the ith ao on the jth grid point of the EXTRA grid
  END_DOC
  integer :: i_nuc, i_rad, i_ang
  double precision :: r(3)
 
  ! CHANGE OMP DIRECTIVES!!!
 
  !$no OMP PARALLEL DO &
  !$no OMP DEFAULT (NONE)  PRIVATE (i) &
  !$no OMP SHARED(aos_in_r_array_extra,n_points_extra_final_grid,final_grid_points_extra)
 
  do i_nuc = 1, nucl_num
    do i_rad = 1, n_points_extra_radial_grid - 1
      do i_ang = 1, n_points_extra_integration_angular
        r(1:3) = grid_points_extra_per_atom(1:3,i_ang,i_rad,i_nuc)
        call give_all_aos_at_r(r, aos_in_r_array_extra_full(1,i_ang,i_rad,i_nuc))
      enddo
    enddo
  enddo
  !$no OMP END PARALLEL DO
 
 END_PROVIDER


 BEGIN_PROVIDER[double precision, mos_in_r_array_extra_full_omp, (mo_num, n_points_extra_integration_angular,  n_points_extra_radial_grid, nucl_num)]
 implicit none
 BEGIN_DOC
 ! mos_in_r_array_extra(i,j)        = value of the ith mo on the jth grid point on the EXTRA GRID 
 END_DOC
 integer :: i_nuc, i_rad, i_ang
 integer :: j
 double precision :: mos_array_extra(mo_num), r(3)

 ! CHANGE THE OMP DIRECTIVES
 !$no OMP PARALLEL DO &
 !$no OMP DEFAULT (NONE)  &
 !$no OMP PRIVATE (i,r,mos_array_extra,j) & 
 !$no OMP SHARED(mos_in_r_array_extra_full_omp,n_points_extra_final_grid,mo_num,final_grid_points_extra)
 !!!!!do i = 1, n_points_extra_final_grid

 do i_nuc = 1, nucl_num
   do i_rad = 1, n_points_extra_radial_grid - 1
     do i_ang = 1, n_points_extra_integration_angular
       r(1:3) = grid_points_extra_per_atom(1:3,i_ang,i_rad,i_nuc)
       call give_all_mos_at_r(r,mos_array_extra)
       do j = 1, mo_num
         mos_in_r_array_extra_full_omp(j,i_ang,i_rad,i_nuc) = mos_array_extra(j)
       enddo
     enddo
   enddo
 enddo
 !$no OMP END PARALLEL DO
 END_PROVIDER

