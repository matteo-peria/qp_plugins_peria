! AOs and MOs IN R2 (grid2)

BEGIN_PROVIDER[double precision, aos_in_r_array2, (ao_num,n_points_final_grid2)]

  BEGIN_DOC
  ! aos_in_r_array2(i,j) = value of the ith ao on the jth grid point
  END_DOC

  implicit none
  integer          :: i

  !$OMP PARALLEL DO               &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i) &
  !$OMP SHARED(aos_in_r_array2,n_points_final_grid2,final_grid_points2)
  do i = 1, n_points_final_grid2
    call give_all_aos_at_r(final_grid_points2(1,i), aos_in_r_array2(1,i))
  enddo
  !$OMP END PARALLEL DO

END_PROVIDER


BEGIN_PROVIDER[double precision, aos_in_r_array2_transp, (n_points_final_grid2,ao_num)]

  BEGIN_DOC
  ! aos_in_r_array2_transp(i,j) = value of the jth ao on the ith grid point
  END_DOC

  implicit none
  integer          :: i, j
  double precision :: aos_array(ao_num), r(3)

  do i = 1, n_points_final_grid2
    do j = 1, ao_num
      aos_in_r_array2_transp(i,j) = aos_in_r_array2(j,i)
    enddo
  enddo

END_PROVIDER


BEGIN_PROVIDER[double precision, mos_in_r_array2, (mo_num,n_points_final_grid2)]
 implicit none
 BEGIN_DOC
 ! mos_in_r_array2(i,j) = value of the ith mo on the jth grid point
 END_DOC
 integer :: i,j
 double precision :: mos_array(mo_num), r(3)
 do i = 1, n_points_final_grid2
  r(1:3) = final_grid_points2(1:3,i)
  call give_all_mos_at_r(r,mos_array)
  do j = 1, mo_num
   mos_in_r_array2(j,i) = mos_array(j)
  enddo
 enddo
 END_PROVIDER


 BEGIN_PROVIDER[double precision, mos_in_r_array2_omp, (mo_num,n_points_final_grid2)]
 implicit none
 BEGIN_DOC
 ! mos_in_r_array2(i,j)        = value of the ith mo on the jth grid point
 END_DOC
 integer :: i
 !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE (i) & 
 !$OMP SHARED(mos_in_r_array2_omp,n_points_final_grid2,mo_num,final_grid_points2)
 do i = 1, n_points_final_grid2
  call give_all_mos_at_r(final_grid_points2(1,i),mos_in_r_array2_omp(1,i))
 enddo
 !$OMP END PARALLEL DO
 END_PROVIDER


 BEGIN_PROVIDER[double precision, mos_in_r_array2_transp,(n_points_final_grid2,mo_num)]
 implicit none
 BEGIN_DOC
 ! mos_in_r_array2_transp(i,j) = value of the jth mo on the ith grid point
 END_DOC
 integer :: i,j
 do i = 1, n_points_final_grid2
  do j = 1, mo_num
   mos_in_r_array2_transp(i,j) = mos_in_r_array2_omp(j,i) 
  enddo
 enddo
 END_PROVIDER


