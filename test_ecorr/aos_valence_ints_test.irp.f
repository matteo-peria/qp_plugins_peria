BEGIN_PROVIDER [ double precision, ao_integrals_ne_to_ao_val_matprod, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Matrix of the nuclear-electron potential for valence orbitals
  ! obtained as matrix-base change from the cartesian AOs and the valence AOs
  END_DOC
  !
  call ao_cart_to_ao_val_normed(ao_integrals_n_e,           &
                         ao_num,                            &
                         ao_integrals_ne_to_ao_val_matprod, &
                         ao_num                             &
  )
  !call ao_cart_to_ao_val(ao_integrals_n_e,                         &
  !                       ao_num,ao_integrals_ne_to_ao_val_matprod, &
  !                       ao_num                                    &
  !)
END_PROVIDER


 BEGIN_PROVIDER [ double precision, ao_integrals_ne_to_ao_val_numeric, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Numerical matrix of the nuclear-electron potential for valence orbitals
  END_DOC

  integer :: ir 
  integer :: i_ao, j_ao, k
  double precision :: weight_r, r(3)
  double precision :: ao_i_r, ao_j_r

  ao_integrals_ne_to_ao_val_numeric(:,:) = 0.d0
  
  ! Notice on the second do-loop we run over all |AO|s
  ! but it is possible to compute only the upper-right triangle of the matrix
  ! and then symmetrize it by imposibine V(i,j) = V(j,i)
  do i_ao = 1, ao_num
    do j_ao = 1, ao_num !i_ao, ao_num 
      do ir = 1, n_points_final_grid
        r(1:3) = final_grid_points(1:3,ir)                                               
        weight_r = final_weight_at_r_vector(ir)
        ao_i_r = aos_val_in_r_from_matprod_transp(ir,i_ao)
        ao_j_r = aos_val_in_r_from_matprod_transp(ir,j_ao)
        ! Sum over all nuclei 
        do k = 1, nucl_num
          double precision :: Z_nucl, R_nucl(3)
          Z_nucl = nucl_charge(k)
          R_nucl(1:3) = nucl_coord(k,1:3)
          ao_integrals_ne_to_ao_val_numeric(i_ao,j_ao) -= Z_nucl/norm2(r-R_nucl)*ao_j_r*ao_i_r*weight_r
          !ao_integrals_ne_to_ao_val_numeric(j_ao,i_ao) -= Z_nucl/norm2(r-R_nucl)*ao_j_r*ao_i_r*weight_r
        end do
      enddo
    enddo
  enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_val_integrals_ne_numeric, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Numerical matrix of the nuclear-electron potential for valence orbitals
  END_DOC

  integer :: ir 
  integer :: i_ao, j_ao, k
  double precision :: weight_r, r(3)
  double precision :: ao_i_r, ao_j_r

  ao_val_integrals_ne_numeric(:,:) = 0.d0
  
  ! Notice on the second do-loop we run over all |AO|s
  ! but it is possible to compute only the upper-right triangle of the matrix
  ! and then symmetrize it by imposibine V(i,j) = V(j,i)
  do i_ao = 1, ao_num
    do j_ao = 1, ao_num !i_ao, ao_num 
      do ir = 1, n_points_final_grid
        r(1:3) = final_grid_points(1:3,ir)                                               
        weight_r = final_weight_at_r_vector(ir)
        ao_i_r = aos_val_in_r_array_transp(ir,i_ao)
        ao_j_r = aos_val_in_r_array_transp(ir,j_ao)
        ! Sum over all nuclei 
        do k = 1, nucl_num
          double precision :: Z_nucl, R_nucl(3)
          Z_nucl = nucl_charge(k)
          R_nucl(1:3) = nucl_coord(k,1:3)
          ao_val_integrals_ne_numeric(i_ao,j_ao) -= Z_nucl/norm2(r-R_nucl)*ao_j_r*ao_i_r*weight_r
          !ao_integrals_ne_to_ao_val_numeric(j_ao,i_ao) -= Z_nucl/norm2(r-R_nucl)*ao_j_r*ao_i_r*weight_r
        end do
      enddo
    enddo
  enddo
END_PROVIDER
