! Adapted from ao_one_e_ints/ao_overlap.irp.f to the case of valence orbitals

 BEGIN_PROVIDER [double precision, ao_val_overlap_from_prim, (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_val_overlap_x, (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_val_overlap_y, (ao_num, ao_num)]
&BEGIN_PROVIDER [double precision, ao_val_overlap_z, (ao_num, ao_num)]

  BEGIN_DOC
  ! Overlap between atomic basis functions of the valence subspace:
  !
  ! :math:`\int \chi_i(r) \chi_j(r) dr`
  END_DOC

  implicit none
  integer          :: i, j, n, l, dim1, power_A(3), power_B(3)
  double precision :: overlap, overlap_x, overlap_y, overlap_z
  double precision :: alpha, beta, c
  double precision :: A_center(3), B_center(3)

  ao_val_overlap_from_prim   = 0.d0
  ao_val_overlap_x = 0.d0
  ao_val_overlap_y = 0.d0
  ao_val_overlap_z = 0.d0

  dim1=100
  !$OMP PARALLEL DO SCHEDULE(GUIDED) &
  !$OMP DEFAULT(NONE) &
  !$OMP PRIVATE(A_center,B_center,power_A,power_B,&
  !$OMP  overlap_x,overlap_y, overlap_z, overlap, &
  !$OMP  alpha, beta,i,j,n,l,c) &
  !$OMP SHARED(nucl_coord,ao_power,ao_val_prim_num, &
  !$OMP  ao_val_overlap_x,ao_val_overlap_y,ao_val_overlap_z,ao_val_overlap_from_prim,ao_num,ao_val_prim_coef_normed_transp,ao_nucl, &
  !$OMP  ao_val_prim_expo_transp,dim1)
  do j=1,ao_num
    A_center(1) = nucl_coord( ao_nucl(j), 1 )
    A_center(2) = nucl_coord( ao_nucl(j), 2 )
    A_center(3) = nucl_coord( ao_nucl(j), 3 )
    power_A(1)  = ao_power( j, 1 )
    power_A(2)  = ao_power( j, 2 )
    power_A(3)  = ao_power( j, 3 )
    do i= 1,ao_num
      B_center(1) = nucl_coord( ao_nucl(i), 1 )
      B_center(2) = nucl_coord( ao_nucl(i), 2 )
      B_center(3) = nucl_coord( ao_nucl(i), 3 )
      power_B(1)  = ao_power( i, 1 )
      power_B(2)  = ao_power( i, 2 )
      power_B(3)  = ao_power( i, 3 )
      do n = 1,ao_val_prim_num(j)
        alpha = ao_val_prim_expo_transp(n,j)
        do l = 1, ao_val_prim_num(i)
          beta = ao_val_prim_expo_transp(l,i)
          call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)
          c = ao_val_prim_coef_normed_transp(n,j) * ao_val_prim_coef_normed_transp(l,i)
          ao_val_overlap_from_prim(i,j) += c * overlap
          if (isnan(ao_val_overlap_from_prim(i,j))) then
            print*,'i,j',i,j
            print*,'l,n',l,n
            print*,'c,overlap',c,overlap
            print*,overlap_x,overlap_y,overlap_z
            stop
          endif
          ao_val_overlap_x(i,j) += c * overlap_x
          ao_val_overlap_y(i,j) += c * overlap_y
          ao_val_overlap_z(i,j) += c * overlap_z
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

END_PROVIDER
