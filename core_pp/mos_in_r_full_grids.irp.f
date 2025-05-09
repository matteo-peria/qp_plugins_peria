! These providers are inspired by src/mo_basis/mos_in_r.irp.f
! But instead of providing the MOs on the pruned grids
! we consider the full grids indexed by nuclei, radial points and angular points

BEGIN_PROVIDER[double precision, mos_in_r_full_grid, (mo_num,n_points_ang_grid,n_points_rad_grid,nucl_num)]
  implicit none
  BEGIN_DOC
  ! MOs in r on the FULL EXTRA grid
  END_DOC
  integer :: i_nucl, i_rad, i_ang
  double precision :: r(3)

  do i_nucl = 1, nucl_num
    do i_rad = 1, n_points_rad_grid
      do i_ang = 1, n_points_ang_grid
        r(1:3) =  grid_points_per_atom(1:3,i_ang,i_rad,i_nucl)
        call give_all_mos_at_r(r,mos_in_r_full_grid(:,i_ang,i_rad,i_nucl))
      end do
    end do
  end do
END_PROVIDER


 BEGIN_PROVIDER[ double precision, mos_in_r_full_extra_grid, (mo_num,n_points_ang_extra_grid,n_points_rad_extra_grid,nucl_num) ]
  implicit none
  BEGIN_DOC
  ! MOs in r on the FULL EXTRA grid
  END_DOC
  integer :: i_nucl, i_rad, i_ang
  double precision :: r(3)
  integer :: j

  write(*,*) 'PROVIDING mos_in_r_full_extra_grid'
  write(*,*) '... shape(mos_in_r_full_extra_grid)', shape(mos_in_r_full_extra_grid)

  mos_in_r_full_extra_grid(:,:,:,:) = 0.d0
  do i_nucl = 1, nucl_num
    do i_rad = 1, n_points_rad_extra_grid-1
      do i_ang = 1, n_points_ang_extra_grid
        r(1:3) =  grid_points_extra_per_atom(1:3,i_ang,i_rad,i_nucl)
        call give_all_mos_at_r(r,mos_in_r_full_extra_grid(:,i_ang,i_rad,i_nucl))
      end do
    end do
  end do

END_PROVIDER

