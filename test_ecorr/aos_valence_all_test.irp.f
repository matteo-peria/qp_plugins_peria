program test_ao_valence_in_r

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call write_int(6, my_n_pt_r_grid, 'radial  external grid over')
  call write_int(6, my_n_pt_a_grid, 'angular external grid over')

  my_extra_grid_becke  = .True.
  PROVIDE tc_grid2_a tc_grid2_r
  my_n_pt_r_extra_grid = tc_grid2_r
  my_n_pt_a_extra_grid = tc_grid2_a
  touch my_extra_grid_becke my_n_pt_r_extra_grid my_n_pt_a_extra_grid

  call write_int(6, my_n_pt_r_extra_grid, 'radial  internal grid over')
  call write_int(6, my_n_pt_a_extra_grid, 'angular internal grid over')

  
  read_ao_integrals_n_e = .false.
  touch read_ao_integrals_n_e

  call main()

end program

! ---

subroutine main()
  implicit none

  integer :: row,col
  integer :: i,j
  double precision :: accu = 0.d0
  double precision :: trace = 0.d0

  read_ao_integrals_n_e = .false.
  write_ao_integrals_n_e = .false.

!!  ! START TESTING VALENCE-SUBSPACE PROJECTOR
!!  print*, "norm of the valence orbitals"
!!  do col=1,n_valence_pseudo_orb
!!    print*, col, sqrt(dot_product(mo_val_coef(:,col),mo_val_coef(:,col)))
!!  end do
!!
!!  provide mo_val_coef
!!  print*, "mo_val_coef"
!!  do row=1,ao_num
!!    write(*,'(100F12.7)'), (mo_val_coef(row,col), col=1,n_valence_pseudo_orb)
!!  end do
!!  
!!  provide proj_mo_val
!!  print*, "proj_mo_val"
!!  do row=1,ao_num
!!    write(*,'(100F12.7)'), (proj_mo_val(row,col), col=1,ao_num)
!!    trace += proj_mo_val(row,row)
!!  end do
!!  print*, "Trace: ", trace 
!!  ! END TESTING VALENCE-SUBSPACE PROJECTOR


!!  ! START TESTING VALENCE COEFFICIENTS 
!!  print*, '' 
!!  print*, 'VALENCE COEFFICIENTS (CARTESIAN)'
!!  print*, '' 
!!
!!  provide ao_val_coef
!!  write(*,*) "AO_VAL_COEF, shape = ", shape(ao_val_coef)
!!  do row = 1, ao_num
!!    write(*,'(100F12.7)'), (ao_val_coef(row,col), col = 1, ao_num)
!!  end do
!!  print*, ''
!!
!!  provide ao_val_coef_normed
!!  write(*,*) "AO_VAL_COEF_NORMED, shape = ", shape(ao_val_coef_normed)
!!  do row = 1, ao_num
!!    write(*,'(100F12.7)'), (ao_val_coef_normed(row,col), col = 1,ao_num)
!!  end do
!!  print*, ''
!!
!!  print*, '' 
!!  print*, 'VALENCE COEFFICIENTS (SPHERICAL)'
!!  print*, '' 
!!
!!  provide ao_val_coef_sphe
!!  write(*,*) "AO_VAL_COEF, shape = ", shape(ao_val_coef_sphe)
!!  do row = 1, ao_sphe_num
!!    write(*,'(100F12.7)'), (ao_val_coef_sphe(row,col), col = 1, ao_sphe_num)
!!  end do
!!  print*, ''
!!
!!  provide ao_val_sphe_coef_normed
!!  write(*,*) "AO_VAL_SPHE_COEF_NORMED, shape = ", shape(ao_val_sphe_coef_normed)
!!  do row = 1, ao_sphe_num
!!    write(*,'(100F12.7)'), (ao_val_sphe_coef_normed(row,col), col = 1,ao_sphe_num)
!!  end do
!!  print*, ''


!!  ! START TESTING OVERLAPS 
!!  print*, '' 
!!  print*, 'OVERLAP TESTS'
!!  print*, '' 
!!
!!  provide ao_overlap
!!  print*,'ao_overlap'
!!  do row = 1, ao_num
!!    write(*,'(100F12.7)'), (ao_overlap(row,col), col = 1, ao_num)
!!    trace += ao_overlap(row,row)
!!  end do
!!  print*, "Trace: ", trace 
!!  trace = 0.d0
!!  print*, ''
!!
!!  print*, '' 
!!  print*, 'OVERLAP TESTS VALENCE SPACE (CARTESIAN)'
!!  print*, '' 
!!
!!  provide ao_val_overlap
!!  print*,'ao_val_overlap'
!!  do row = 1, ao_num
!!    write(*,'(100F12.7)'), (ao_val_overlap(row,col), col = 1, ao_num)
!!    trace += ao_val_overlap(row,row)
!!  end do
!!  print*, "Trace: ", trace 
!!  trace = 0.d0
!!  print*, ''
!!
!!  provide ao_val_overlap_evec_overlap
!!  print*,'ao_val_overlap_evec_overlap'
!!  do row = 1, ao_num
!!    write(*,'(100F12.7)'), (ao_val_overlap_evec_overlap(row,col), col = 1, ao_num)
!!    trace += ao_val_overlap_evec_overlap(row,row)
!!  end do
!!  print*, "Trace: ", trace 
!!  trace = 0.d0
!!  print*, ''
!!
!!  provide ao_val_overlap_as_matprod
!!  print*,'ao_val_overlap_as_matprod'
!!  do row = 1, ao_num
!!    write(*,'(100F12.7)'), (ao_val_overlap_as_matprod(row,col), col = 1, ao_num)
!!    trace += ao_val_overlap_as_matprod(row,row)
!!  end do
!!  print*, "Trace: ", trace 
!!  trace = 0.d0
!!  print*, ''
!!
!!  provide ao_val_normed_overlap
!!  print*,'ao_val_normed_overlap'
!!  do row = 1, ao_num
!!    write(*,'(100F12.7)'), (ao_val_normed_overlap(row,col), col = 1, ao_num)
!!    trace += ao_val_normed_overlap(row,row)
!!  end do
!!  print*, "Trace: ", trace 
!!  trace = 0.d0
!!  print*, ''
!!
  print*,'overlap'
  provide ao_val_overlap
  do row = 1, ao_num
    write(*,'(100F16.13)'), (ao_val_overlap(row,col), col = 1, ao_num)
  end do
  print*,'eigvec'
  provide ao_val_overlap_eigvec
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_val_overlap_eigvec(row,col), col = 1, ao_num)
  end do
  pause
  
  provide ao_val_overlap_normed_as_matprod
  print*,'ao_val_overlap_normed_as_matprod'
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_val_overlap_normed_as_matprod(row,col), col = 1, ao_num)
    trace += ao_val_overlap_normed_as_matprod(row,row)
  end do
  print*, "Trace: ", trace 
  trace = 0.d0
  print*, ''
!!
!!  !---
!!  
!!  print*, '' 
!!  print*, 'OVERLAP TESTS VALENCE SPACE (SPHERICAL)'
!!  print*, '' 
!!
!!  provide ao_val_overlap_sphe
!!  print*,'ao_val_overlap_sphe'
!!  do row = 1, ao_sphe_num
!!    write(*,'(100F12.7)'), (ao_val_overlap_sphe(row,col), col = 1, ao_sphe_num)
!!    trace += ao_val_overlap_sphe(row,row)
!!  end do
!!  print*, "Trace: ", trace 
!!  trace = 0.d0
!!  print*, ''
!!
!!  provide ao_val_overlap_sphe_evec_overlap
!!  print*,'ao_val_overlap_sphe_evec_overlap'
!!  do row = 1, ao_sphe_num
!!   write(*,'(100F12.7)'), (ao_val_overlap_sphe_evec_overlap(row,col), col = 1, ao_sphe_num)
!!   trace += ao_val_overlap_sphe_evec_overlap(row,row)
!!  enddo 
!!  print*, "Trace: ", trace 
!!  trace = 0.d0
!!  print*, ''
!!
!!  provide ao_val_overlap_sphe_as_matprod
!!  print*,'ao_val_overlap_sphe_as_matprod'
!!  do row = 1, ao_sphe_num
!!   write(*,'(100F12.7)'), (ao_val_overlap_sphe_as_matprod(row,col), col = 1, ao_sphe_num)
!!   trace += ao_val_overlap_sphe_as_matprod(row,row)
!!  enddo 
!!  print*, "Trace: ", trace 
!!  trace = 0.d0
!!  print*, ''
!!
!!  provide ao_val_sphe_normed_overlap
!!  print*,'ao_val_sphe_normed_overlap'
!!  do row = 1, ao_sphe_num
!!   write(*,'(100F12.7)'), (ao_val_sphe_normed_overlap(row,col), col = 1, ao_sphe_num)
!!   trace += ao_val_sphe_normed_overlap(row,row)
!!  enddo 
!!  print*, "Trace: ", trace 
!!  trace = 0.d0
!!  print*, ''
!! 
!!  provide ao_val_overlap_sphe_normed_as_matprod
!!  print*,'ao_val_overlap_sphe_normed_as_matprod'
!!  do row = 1, ao_sphe_num
!!   write(*,'(100F12.7)'), (ao_val_overlap_sphe_normed_as_matprod(row,col), col = 1, ao_sphe_num)
!!   trace += ao_val_overlap_sphe_normed_as_matprod(row,row)
!!  enddo 
!!  print*, "Trace: ", trace 
!!  trace = 0.d0
!!  print*, ''
 
  !---

!!  print*, '' 
!!  print*, 'OVERLAP TESTS CORE SPACE (CARTESIAN)'
!!  print*, '' 
!!
!!  provide ao_core_overlap
!!  print*,'ao_core_overlap'
!!  do row = 1, ao_num
!!    write(*,'(100F12.7)'), (ao_core_overlap(row,col), col = 1, ao_num)
!!    trace += ao_core_overlap(row,row)
!!  end do
!!  print*, "Trace: ", trace 
!!  trace = 0.d0
!!  print*, ''
!!
!!  provide ao_core_overlap_evec_overlap
!!  print*,'ao_core_overlap_evec_overlap'
!!  do row = 1, ao_num
!!    write(*,'(100F12.7)'), (ao_core_overlap_evec_overlap(row,col), col = 1, ao_num)
!!    trace += ao_core_overlap_evec_overlap(row,row)
!!  end do
!!  print*, "Trace: ", trace 
!!  trace = 0.d0
!!  print*, ''
!!
!!  provide ao_core_overlap_as_matprod
!!  print*,'ao_core_overlap_as_matprod'
!!  do row = 1, ao_num
!!    write(*,'(100F12.7)'), (ao_core_overlap_as_matprod(row,col), col = 1, ao_num)
!!    trace += ao_core_overlap_as_matprod(row,row)
!!  end do
!!  print*, "Trace: ", trace 
!!  trace = 0.d0
!!  print*, ''
!!
!!  provide ao_core_normed_overlap
!!  print*,'ao_core_normed_overlap'
!!  do row = 1, ao_num
!!    write(*,'(100F12.7)'), (ao_core_normed_overlap(row,col), col = 1, ao_num)
!!    trace += ao_core_normed_overlap(row,row)
!!  end do
!!  print*, "Trace: ", trace 
!!  trace = 0.d0
!!  print*, ''
!!
!!  provide ao_core_overlap_normed_as_matprod
!!  print*,'ao_core_overlap_normed_as_matprod'
!!  do row = 1, ao_num
!!    write(*,'(100F12.7)'), (ao_core_overlap_normed_as_matprod(row,col), col = 1, ao_num)
!!    trace += ao_core_overlap_normed_as_matprod(row,row)
!!  end do
!!  print*, "Trace: ", trace 
!!  trace = 0.d0
!!  print*, ''

  !---



  ! START TESTING VALENCE INTEGRALS
  provide ao_integrals_ne_to_ao_val_numeric
  print*, "Shape of ao_integrals_ne_to_ao_val_numeric: ", shape(ao_integrals_ne_to_ao_val_numeric)
  print*, "ao_integrals_ne_to_ao_val_numeric"
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_integrals_ne_to_ao_val_numeric(row,col), col = 1, ao_num)
    !write(*,*), (ao_integrals_ne_to_ao_val_numeric(row,col), col = 1, ao_num)
  end do
  
  provide ao_integrals_ne_to_ao_val_matprod
  print*, "Shape of ao_integrals_ne_to_ao_val_matprod: ", shape(ao_integrals_ne_to_ao_val_matprod)
  print*, "ao_integrals_ne_to_ao_val_matprod"
  do row = 1, ao_num
    !write(*,'(100F12.7)'), (ao_integrals_ne_to_ao_val_numeric(row,col), col = 1, ao_num)
    write(*,'(100F12.7)'), (ao_integrals_ne_to_ao_val_matprod(row,col), col = 1, ao_num)
  end do
  print*, ''

  provide ao_val_coef_normed
  write(*,*) "AO_VAL_COEF_NORMED, shape = ", shape(ao_val_coef_normed)
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_val_coef_normed(row,col), col = 1,ao_num)
  end do
  print*, ''


!!  do i=1,size(ao_integrals_ne_to_ao_val_numeric,1)
!!    do j=1,size(ao_integrals_ne_to_ao_val_numeric,2)
!!      accu += ao_integrals_ne_to_ao_val_numeric(i,j) - ao_integrals_ne_to_ao_val_matprod(i,j)
!!    end do
!!  end do
!!  print*, "Difference between numerical and matrix product V_ne: ", accu

  provide ao_val_integrals_ne
  print*, "ao_val_integrals_ne"
  print*, "Shape of ao_val_integrals_ne: ", shape(ao_val_integrals_ne)
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_val_integrals_ne(row,col), col = 1, ao_num)
  end do



  provide ao_val_integrals_ne_numeric
  print*, "Shape of ao_val_integrals_ne_numeric: ", shape(ao_val_integrals_ne_numeric)
  print*, "ao_val_integrals_ne_numeric"
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_val_integrals_ne_numeric(row,col), col = 1, ao_num)
  end do
  ! END TESTING VALENCE INTEGRALS
!!
!!  ! START TESTING VALENCE PRIMITIVES
!!  provide ao_pp_prim_num_temp
!!  write(*,*) "AO_PP_PRIM_NUM_TEMP"
!!  write(*,*) ao_pp_prim_num_temp
!!  print*, ''
!!
!!  provide ao_pp_prim_num_max_temp
!!  write(*,*) "AO_PP_PRIM_NUM_MAX_TEMP"
!!  write(*,*) ao_pp_prim_num_max_temp
!!  print*, ''
!!
!!  provide ao_pp_prim_num
!!  write(*,*) "AO_PP_PRIM_NUM"
!!  write(*,*) ao_pp_prim_num
!!  print*, ''
!!
!!  provide ao_pp_prim_num_max
!!  write(*,*) "AO_PP_PRIM_NUM_MAX"
!!  write(*,*) ao_pp_prim_num_max
!!  print*, ''
!!
!!  provide ao_expo_transp
!!  write(*,*) "AO_EXPO_TRANSP, shape = ", shape(ao_expo_transp)
!!  do row = 1, ao_prim_num_max
!!    write(*,'(100F12.7)'), (ao_expo_transp(row,col), col = 1, ao_num)
!!  end do
!!  print*, ''
!!
!!  provide ao_pp_expo_transp_temp
!!  write(*,*) "AO_PP_EXPO_TRANSP_TEMP, shape = ", shape(ao_pp_expo_transp_temp)
!!  do row = 1, ao_pp_prim_num_max_temp
!!    write(*,'(100F12.7)'), (ao_pp_expo_transp_temp(row,col), col = 1, ao_num)
!!  end do
!!  print*, ''
!!
!!  provide ao_pp_expo_transp
!!  write(*,*) "AO_PP_EXPO_TRANSP, shape = ", shape(ao_pp_expo_transp)
!!  do row = 1, ao_pp_prim_num_max
!!    write(*,'(100F12.7)'), (ao_pp_expo_transp(row,col), col = 1, ao_num)
!!  end do
!!  print*, ''
!!
!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  
!!  provide ao_coef_transp
!!  write(*,*) "AO_COEF_TRANSP, shape = ", shape(ao_coef_transp)
!!  do row = 1, ao_prim_num_max
!!    write(*,'(100F12.7)'), (ao_coef_transp(row,col), col = 1, ao_num)
!!  end do
!!  print*, ''
!!
!!  provide ao_pp_coef_transp_temp
!!  write(*,*) "AO_PP_COEF_TRANSP_TEMP, shape = ", shape(ao_pp_coef_transp_temp)
!!  do row = 1, ao_pp_prim_num_max_temp
!!    write(*,'(100F12.7)'), (ao_pp_coef_transp_temp(row,col), col = 1, ao_num)
!!  end do
!!  print*, ''
!!
!!  provide ao_pp_coef_transp
!!  write(*,*) "AO_PP_COEF_TRANSP, shape = ", shape(ao_pp_coef_transp)
!!  do row = 1, ao_pp_prim_num_max
!!    write(*,'(100F12.7)'), (ao_pp_coef_transp(row,col), col = 1, ao_num)
!!  end do
!!  print*, ''
!!  ! END TESTING VALENCE PRIMITIVES
  

end subroutine


! Test providers for aos_valence.irp.f

BEGIN_PROVIDER [ double precision, ao_val_overlap_as_matprod, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Matrix of the nuclear-electron potential for valence orbitals
  ! obtained as matrix-base change from the cartesian AOs and the valence AOs
  END_DOC
  !
  call ao_cart_to_ao_val(ao_overlap,                &
                         ao_num,                    &
                         ao_val_overlap_as_matprod, &
                         ao_num                     &
  )
END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_val_overlap_normed_as_matprod, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Matrix of the nuclear-electron potential for valence orbitals
  ! obtained as matrix-base change from the cartesian AOs and the valence AOs
  END_DOC
  !
  call ao_cart_to_ao_val_normed(ao_overlap,                &
                                ao_num,                    &
                                ao_val_overlap_normed_as_matprod, &
                                ao_num                     &
  )
END_PROVIDER
