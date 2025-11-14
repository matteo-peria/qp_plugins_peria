subroutine give_all_aos_at_r_new(r, tmp_array)
  implicit none
  BEGIN_DOC
  ! input  : r == r(1) = x and so on
  ! output : tmp_array(i) = aos(i) evaluated in r
  END_DOC

  double precision, intent(in)  :: r(3)
  double precision, intent(out) :: tmp_array(ao_num)

  integer          :: p_ao(3)
  integer          :: i, j, l, k
  double precision :: dr(3), r2
  double precision :: beta
  double precision :: angular_factor, radial_factor, radial_exponent

  tmp_array(:) = 0.d0

  ! Loop over all nuclei
  do i = 1, nucl_num

    ! vectorial distance to nucleus i
    dr(1:3) = r(1:3) - nucl_coord(i,1:3)
    ! r^2 = |r - R_i|^2
    r2 = dot_product(dr, dr)   ! slightly nicer than sum(dr*dr)

    ! Loop over all AOs on nucleus i
    do j = 1, Nucl_N_Aos(i)

      ! angular exponents of this AO
      p_ao(1:3) = ao_power_ordered_transp_per_nucl(1:3, j, i)

      ! angular polynomial factor (Cartesian)
      angular_factor = dr(1)**p_ao(1) * dr(2)**p_ao(2) * dr(3)**p_ao(3)

      ! index of this AO in the global ordering
      k = Nucl_Aos_transposed(j,i)

      radial_factor = 0.d0

      ! Loop over all primitives in this contracted AO
      do l = 1, ao_prim_num(k)
        beta = ao_expo_ordered_transp_per_nucl(l, j, i)
        radial_exponent = beta*r2

        if (radial_exponent .gt. 50.d0) cycle

        radial_factor += ao_coef_normalized_ordered_transp_per_nucl(l, j, i) * dexp(-radial_exponent)
      end do

      ! AO value
      tmp_array(k) = radial_factor * angular_factor
    end do
  end do

  return
end subroutine

subroutine give_all_aos_at_r_omp(r, tmp_array)
  implicit none
  double precision, intent(in)  :: r(3)
  double precision, intent(out) :: tmp_array(ao_num)

  integer          :: p_ao(3)
  integer          :: i, j, l, k
  double precision :: dr(3), r2
  double precision :: beta
  double precision :: angular_factor, radial_factor, radial_exponent

  tmp_array(:) = 0.d0

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP SHARED(r, tmp_array, nucl_num, nucl_coord, &
  !$OMP        Nucl_N_Aos, Nucl_Aos_transposed, &
  !$OMP        ao_power_ordered_transp_per_nucl, &
  !$OMP        ao_prim_num, ao_expo_ordered_transp_per_nucl, &
  !$OMP        ao_coef_normalized_ordered_transp_per_nucl) &
  !$OMP PRIVATE(i, j, l, k, dr, r2, p_ao, beta, &
  !$OMP         angular_factor, radial_factor, radial_exponent)
  
  !$OMP DO
  do i = 1, nucl_num

    dr(1:3) = r(1:3) - nucl_coord(i,1:3)
    r2 = dot_product(dr, dr)

    do j = 1, Nucl_N_Aos(i)
      p_ao(1:3) = ao_power_ordered_transp_per_nucl(1:3, j, i)
      angular_factor = product(dr(:)**p_ao(:))

      k = Nucl_Aos_transposed(j,i)
      radial_factor = 0.d0

      do l = 1, ao_prim_num(k)
        beta            = ao_expo_ordered_transp_per_nucl(l, j, i)
        radial_exponent = beta*r2
        if (radial_exponent .gt. 50.d0) cycle

        radial_factor = radial_factor + &
          ao_coef_normalized_ordered_transp_per_nucl(l, j, i) * dexp(-radial_exponent)
      end do

      tmp_array(k) = radial_factor * angular_factor
    end do
  end do
!$OMP END DO
!$OMP END PARALLEL

end subroutine
