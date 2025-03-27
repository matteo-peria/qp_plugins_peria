 BEGIN_PROVIDER[double precision, aos_val_in_r_from_matprod, (ao_num,n_points_final_grid)]
&BEGIN_PROVIDER[double precision, aos_val_in_r_from_matprod_transp, (n_points_final_grid,ao_num)]
  implicit none
  BEGIN_DOC
  ! aos_val_in_r_from_matprod(i,j) = value of the i-th valence AO on the j-th grid point
  END_DOC
  !
  !call get_AB_prod(                       &
  !    ao_val_coef_normed, ao_num, ao_num, &
  !    aos_in_r_array, n_points_final_grid,&
  !    aos_val_in_r_from_matprod               &
  !)
  call get_AB_prod(                       &
      ao_val_coef_normed_transp, ao_num, ao_num, &
      aos_in_r_array, n_points_final_grid,&
      aos_val_in_r_from_matprod               &
  )
  ! Compute transpose
  integer :: i, j
  do i = 1,ao_num
    do j = 1,n_points_final_grid
      aos_val_in_r_from_matprod_transp(j,i) = aos_val_in_r_from_matprod(i,j)
    enddo
  enddo
END_PROVIDER
