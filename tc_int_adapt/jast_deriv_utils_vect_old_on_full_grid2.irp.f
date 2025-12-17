! TEST ONLY SUBROUTINE OF THE ORIGINAL NON_H_INTS CODE IN ORDER TO SEE
! THEM IN ACTION ON THE FULL GRID2 AND NOT ONLY ON THE 2D-PRUNED ONE

!subroutine get_grad1_u12_on_full_grid2_old(ipoint, n_grid2, resx, resy, resz, res)
subroutine get_grad1_u12_on_full_grid2_old(ipoint, n_nuc, n_rad, n_ang, resx, resy, resz, res)

  implicit none
  integer, intent(in)  :: ipoint
  !integer, intent(in)  :: n_grid2
  integer, intent(in)  :: n_nuc
  integer, intent(in)  :: n_rad
  integer, intent(in)  :: n_ang

  !double precision, intent(out) :: resx(n_grid2), resy(n_grid2), resz(n_grid2), res(n_grid2)
  ! above line substituted by the whole array sizes:
  double precision, intent(out) :: resx(n_nuc*(n_rad-1)*n_ang)
  double precision, intent(out) :: resy(n_nuc*(n_rad-1)*n_ang)
  double precision, intent(out) :: resz(n_nuc*(n_rad-1)*n_ang)
  double precision, intent(out) :: res( n_nuc*(n_rad-1)*n_ang)

  integer                       :: jpoint
  double precision              :: env_r1, tmp
  double precision              :: grad1_env(3), r1(3)
  double precision, allocatable :: env_r2(:)
  double precision, allocatable :: u2b_r12(:), gradx1_u2b(:), grady1_u2b(:), gradz1_u2b(:)
  double precision, allocatable :: u2b_mu(:), gradx1_mu(:), grady1_mu(:), gradz1_mu(:)
  double precision, allocatable :: u2b_nu(:), gradx1_nu(:), grady1_nu(:), gradz1_nu(:)
  double precision, external    :: env_nucl

  PROVIDE j1e_type j2e_type env_type
  PROVIDE mu_erf nu_erf a_boys
  PROVIDE final_grid_points
  PROVIDE final_grid_points_extra

  r1(1) = final_grid_points(1,ipoint)
  r1(2) = final_grid_points(2,ipoint)
  r1(3) = final_grid_points(3,ipoint)

  if( (j2e_type .eq. "Mu")        .or. &
      (j2e_type .eq. "Mur")       .or. &
      (j2e_type .eq. "Jpsi")   .or. &
      (j2e_type .eq. "Mugauss")   .or. &
      (j2e_type .eq. "Murgauss")  .or. &
      (j2e_type .eq. "Bump")      .or. &
      (j2e_type .eq. "Boys") ) then

    if(env_type .eq. "None") then

      !call get_grad1_j12_on_full_grid2_old(r1, n_grid2, resx, resy, resz)
      call get_grad1_j12_on_full_grid2_old(r1, n_nuc, n_rad, n_ang, resx, resy, resz)

    else
      print*, "Problem in get_grad1_u12_on_full_grid2_old: env_type is not NONE"
      stop
    end if
  else
    print*, "Problem in get_grad1_u12_on_full_grid2_old: j2e_type is not MU"
    stop
  end if

  !do jpoint = 1, n_points_extra_final_grid
  !  res(jpoint) = resx(jpoint) * resx(jpoint) + resy(jpoint) * resy(jpoint) + resz(jpoint) * resz(jpoint)
  !enddo


  integer :: irad, iang, inuc
  jpoint = 1
  do inuc = 1, nucl_num
    do irad = 1, n_points_extra_radial_grid-1
      do iang = 1, n_points_extra_integration_angular
        res(jpoint) = resx(jpoint)*resx(jpoint) &
                    + resy(jpoint)*resy(jpoint) &
                    + resz(jpoint)*resz(jpoint)
        jpoint = jpoint + 1
      end do
    end do
  end do

  return
end subroutine


!subroutine get_grad1_j12_on_full_grid2_old(r1, n_grid2, gradx, grady, gradz)
subroutine get_grad1_j12_on_full_grid2_old(r1, n_nuc, n_rad, n_ang, gradx, grady, gradz)
  include 'constants.include.F'

  implicit none
  double precision, intent(in)  :: r1(3)
  !integer, intent(in)  :: n_grid2
  integer, intent(in)  :: n_nuc
  integer, intent(in)  :: n_rad
  integer, intent(in)  :: n_ang
  !double precision, intent(out) :: gradx(n_grid2)
  !double precision, intent(out) :: grady(n_grid2)
  !double precision, intent(out) :: gradz(n_grid2)
  double precision, intent(out) :: gradx(n_nuc*(n_rad-1)*n_ang)
  double precision, intent(out) :: grady(n_nuc*(n_rad-1)*n_ang)
  double precision, intent(out) :: gradz(n_nuc*(n_rad-1)*n_ang)


  integer                       :: jpoint
  integer                       :: irad, iang, inuc
  integer                       :: i_nucl, p, mpA, npA, opA
  double precision              :: r2(3)
  double precision              :: dx, dy, dz, r12, tmp
  double precision              :: mu_val, mu_tmp, mu_der(3)
  double precision              :: rn(3), f1A, grad1_f1A(3), f2A, grad2_f2A(3), g12, grad1_g12(3)
  double precision              :: tmp1, tmp2


  PROVIDE j2e_type

  jpoint = 1
  if(j2e_type .eq. "Mu") then
    !do jpoint = 1, n_points_extra_final_grid
    do inuc = 1, nucl_num
      do irad = 1, n_points_extra_radial_grid-1
        do iang = 1, n_points_extra_integration_angular
          r2(1) = grid_points_extra_per_atom(1,iang,irad,inuc)
          r2(2) = grid_points_extra_per_atom(2,iang,irad,inuc)
          r2(3) = grid_points_extra_per_atom(3,iang,irad,inuc)

          dx = r1(1) - r2(1)
          dy = r1(2) - r2(2)
          dz = r1(3) - r2(3)

          r12 = dsqrt(dx * dx + dy * dy + dz * dz)
          if(r12 .lt. 1d-10) then
            gradx(jpoint) = 0.d0 
            grady(jpoint) = 0.d0 
            gradz(jpoint) = 0.d0 
            cycle
          endif

          tmp = 0.5d0 * (1.d0 - derf(mu_erf * r12)) / r12

          gradx(jpoint) = tmp * dx
          grady(jpoint) = tmp * dy
          gradz(jpoint) = tmp * dz

          jpoint = jpoint + 1
        end do
      end do
    enddo

  else
    print*, "Problem in get_grad1_j12_on_full_grid2_old"
  end if

  return
end subroutine
