! Adapted from ao_integrals_n_e
BEGIN_PROVIDER [ double precision, ao_val_integrals_ne, (ao_num,ao_num)]
  BEGIN_DOC
  !  Nucleus-electron interaction, in the valence-projected |AO| basis set.
  !
  !  :math:`\langle \tilde{chi}_i | -\sum_A \frac{1}{|r-R_A|} | \tilde{chi}_j \rangle`
  !
  !  These integrals also contain the pseudopotential integrals.
  END_DOC

  implicit none
  integer          :: num_A, num_B, power_A(3), power_B(3)
  integer          :: i, j, k, l, n_pt_in, m
  double precision :: alpha, beta
  double precision :: A_center(3),B_center(3),C_center(3)
  double precision :: overlap_x,overlap_y,overlap_z,overlap,dx,NAI_pol_mult

  ao_val_integrals_ne = 0.d0


  !if (read_ao_val_integrals_ne) then
  !  call ezfio_get_ao_one_e_ints_ao_val_integrals_ne(ao_val_integrals_ne)
  !  print *,  'AO N-e integrals read from disk'
  !else

    !$OMP PARALLEL                                                   &
        !$OMP DEFAULT (NONE)                                         &
        !$OMP PRIVATE (i,j,k,l,m,alpha,beta,A_center,B_center,C_center,power_A,power_B,&
        !$OMP          num_A,num_B,Z,c,c1,n_pt_in)                      &
        !$OMP SHARED (ao_num,ao_pp_prim_num,ao_pp_expo_transp,ao_power,ao_nucl,nucl_coord,ao_pp_coef_normalized_transp,&
        !$OMP         n_pt_max_integrals,ao_val_integrals_ne,nucl_num,nucl_charge)

    n_pt_in = n_pt_max_integrals

    !$OMP DO SCHEDULE (dynamic)
    do j = 1, ao_num
      num_A = ao_nucl(j)                      ! Nucleus associated to j-th AO
      power_A(1:3)= ao_power(j,1:3)           ! Power over the cartesian polynomial of j-th AO
      A_center(1:3) = nucl_coord(num_A,1:3)   ! Nuclear coordinate
      

      do i = 1, ao_num
        num_B = ao_nucl(i)
        power_B(1:3)= ao_power(i,1:3)
        B_center(1:3) = nucl_coord(num_B,1:3)

        do l=1,ao_pp_prim_num(j)
          alpha = ao_pp_expo_transp(l,j)

          do m=1,ao_pp_prim_num(i)
            beta = ao_pp_expo_transp(m,i)

            double precision :: c, c1
            c = 0.d0
            ! Loop over all postive charges' centers (nuclei)
            do  k = 1, nucl_num
              double precision :: Z
              Z = nucl_charge(k)
              C_center(1:3) = nucl_coord(k,1:3)
              c1 = NAI_pol_mult(A_center, B_center, power_A, power_B &
                               , alpha, beta, C_center, n_pt_in )
              c = c - Z * c1
            enddo
            ao_val_integrals_ne(i,j) = ao_val_integrals_ne(i,j) &
                + ao_pp_coef_normalized_transp(l,j)             &
                * ao_pp_coef_normalized_transp(m,i) * c
          enddo
        enddo
      enddo
    enddo

    !$OMP END DO
    !$OMP END PARALLEL

    IF(do_pseudo) THEN
       ao_val_integrals_ne += ao_pseudo_integrals
    ENDIF
    IF(point_charges) THEN
       ao_val_integrals_ne += ao_integrals_pt_chrg
    ENDIF

  !endif


  !if (write_ao_val_integrals_ne) then
  !  call ezfio_set_ao_one_e_ints_ao_val_integrals_ne(ao_val_integrals_ne)
  !  print *,  'AO N-e integrals written to disk'
  !endif

END_PROVIDER
