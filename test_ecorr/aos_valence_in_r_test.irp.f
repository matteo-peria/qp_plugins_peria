BEGIN_PROVIDER[double precision, aos_val_in_r_array, (ao_num,n_points_final_grid)]
  BEGIN_DOC
  ! aos_in_r_array(i,j) = value of the ith ao on the jth grid point
  END_DOC
  implicit none
  integer          :: i, j
  double precision :: tmp_array(ao_num), r(3)
  !$OMP PARALLEL DO               &
  !$OMP DEFAULT (NONE)            &
  !$OMP PRIVATE (i,r,tmp_array,j) & 
  !$OMP SHARED(aos_val_in_r_array,n_points_final_grid,ao_num,final_grid_points)
  do i = 1, n_points_final_grid
    r(1) = final_grid_points(1,i)
    r(2) = final_grid_points(2,i)
    r(3) = final_grid_points(3,i)
    call give_all_aos_val_at_r(r, tmp_array)
    do j = 1, ao_num
      aos_val_in_r_array(j,i) = tmp_array(j)
    enddo
  enddo
  !$OMP END PARALLEL DO
END_PROVIDER


BEGIN_PROVIDER[double precision, aos_val_in_r_array_transp, (n_points_final_grid,ao_num)]
  implicit none
  integer :: row, col
  do row=1,n_points_final_grid 
    do col=1,ao_num
      aos_val_in_r_array_transp(row,col) = aos_val_in_r_array(col,row)
    end do
  end do
END_PROVIDER


subroutine give_all_aos_val_at_r(r, tmp_array)
  BEGIN_dOC
  ! input  : r == r(1) = x and so on
  ! output : tmp_array(i) = aos(i) evaluated in $\textbf{r}$
  END_DOC
  implicit none
  double precision, intent(in)  :: r(3)
  double precision, intent(out) :: tmp_array(ao_num)
  integer                       :: p_ao(3)
  integer                       :: i, j, k, l, m
  double precision              :: dx, dy, dz, r2
  double precision              :: dx2, dy2, dz2
  double precision              :: c_ao(3)
  double precision              :: beta
  !
  do i = 1, nucl_num
    c_ao(1:3) = nucl_coord(i,1:3)
    dx = r(1) - c_ao(1)
    dy = r(2) - c_ao(2)
    dz = r(3) - c_ao(3)
    r2 = dx*dx + dy*dy + dz*dz
    do j = 1, Nucl_N_Aos(i)
      k = Nucl_Aos_transposed(j,i) ! index of the ao in the ordered format
      p_ao(1:3) = ao_power_ordered_transp_per_nucl(1:3,j,i)
      dx2 = dx**p_ao(1)
      dy2 = dy**p_ao(2)
      dz2 = dz**p_ao(3)
      tmp_array(k) = 0.d0
      do l = 1, ao_val_prim_num(k)
        beta = ao_val_prim_expo_transp_per_nucl(l,j,i)
        if(beta*r2.gt.50.d0) cycle
        tmp_array(k) += ao_val_prim_coef_normed_transp_per_nucl(l,j,i) * dexp(-beta*r2)
      enddo
      tmp_array(k) = tmp_array(k) * dx2 * dy2 * dz2
    enddo
  enddo
  return
end subroutine
