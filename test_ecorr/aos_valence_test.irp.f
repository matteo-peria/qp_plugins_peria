program aos_valence_test

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

  call main()

end program

! ---

subroutine main()
  implicit none

  integer :: row,col
  integer :: i,j
  double precision :: accu = 0.d0
  double precision :: trace = 0.d0

  print*, '' 
  print*, 'VALENCE COEFFICIENTS (CARTESIAN)'
  print*, '' 

  provide ao_val_coef
  write(*,*) "AO_VAL_COEF, shape = ", shape(ao_val_coef)
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_val_coef(row,col), col = 1, ao_num)
  end do
  print*, ''

  provide ao_val_coef_normed
  write(*,*) "AO_VAL_COEF_NORMED, shape = ", shape(ao_val_coef_normed)
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_val_coef_normed(row,col), col = 1,ao_num)
  end do
  print*, ''


  ! START TESTING OVERLAPS 
  print*, '' 
  print*, 'OVERLAP TESTS'
  print*, '' 

  print*, '' 
  print*, 'OVERLAP ON THE FULL SPACE'
  print*, '' 

  provide ao_overlap
  print*,'AO_OVERLAP'
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_overlap(row,col), col = 1, ao_num)
    trace += ao_overlap(row,row)
  end do
  print*, "Trace: ", trace 
  trace = 0.d0
  print*, ''

  print*, '' 
  print*, 'OVERLAP TESTS VALENCE SPACE (CARTESIAN)'
  print*, '' 

  provide ao_val_overlap
  print*,'AO_VAL_OVERLAP computed as matrix product (C^v)^TSC^v'
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_val_overlap(row,col), col = 1, ao_num)
    trace += ao_val_overlap(row,row)
  end do
  print*, "Trace: ", trace 
  trace = 0.d0
  print*, ''
 
  print*,'Now we diagonalize the val-subspace overlap matrix by using the ao_val_overlap eigenvectors:'
  print*, ''
  print*, 'AO_VAL_OVERLAP_EIGVAL'
  provide ao_val_overlap_eigval
  print*, ao_val_overlap_eigval
  print*, ''
  print*, 'AO_VAL_OVERLAP_EIGVEC'
  provide ao_val_overlap_eigvec
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_val_overlap_eigvec(row,col), col = 1, ao_num)
  end do
  print*, ''

  provide ao_val_overlap_diag
  print*,'AO_VAL_OVERLAP_DIAG'
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_val_overlap_diag(row,col), col = 1, ao_num)
    trace += ao_val_overlap_diag(row,row)
  end do
  print*, "Trace: ", trace 
  trace = 0.d0
  print*, ''

  print*, '' 
  print*, 'OVERLAP USING NORMED EIGENVECTORS'
  print*, '' 

  provide ao_val_normed_overlap
  print*,'AO_VAL_NORMED_OVERLAP'
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_val_normed_overlap(row,col), col = 1, ao_num)
    trace += ao_val_normed_overlap(row,row)
  end do
  print*, "Trace: ", trace 
  trace = 0.d0
  print*, ''

  provide ao_val_overlap_from_prim
  print*,'AO_VAL_OVERLAP_FROM_PRIM'
  do row = 1, ao_num
    write(*,'(100F12.7)'), (ao_val_overlap_from_prim(row,col), col = 1, ao_num)
    trace += ao_val_overlap_from_prim(row,row)
  end do
  print*, "Trace: ", trace 
  trace = 0.d0
  print*, ''

!mo_val_coef_transp*ao_val_overlap*mo_val_coef

!  double precision :: temp(ao_num, n_valence_pseudo_orb) 
!  double precision :: temp2(n_valence_pseudo_orb, n_valence_pseudo_orb)

  double precision, allocatable :: temp(:,:) 
  double precision, allocatable :: temp2(:,:)

  allocate(temp(ao_num,n_valence_pseudo_orb))
  allocate(temp2(n_valence_pseudo_orb,n_valence_pseudo_orb))

  call dgemm("N","N", &
        & ao_num, n_valence_pseudo_orb, ao_num, &
        & 1.d0, ao_val_overlap, size(ao_val_overlap,1), &
        & mo_val_coef, size(mo_val_coef,1), 0.d0, &
        & temp, size(temp,1)&
  )

  temp2 = 0.d0
  !do row=1,n_valence_pseudo_orb
  !  temp2(row,row) = -1.d0
  !end do

  call dgemm("T","N", &
        & n_valence_pseudo_orb, n_valence_pseudo_orb, ao_num, &
        & 1.d0, mo_val_coef, size(mo_val_coef,1), &
        & temp, size(temp,1), 0.d0, &
        & temp2, size(temp2,1)&
  )


  print*, n_valence_pseudo_orb
  write(*,*) minval(temp2), maxval(temp2)
  do row = 1, n_valence_pseudo_orb
    write(*,'(100F12.7)'), (temp2(row,col), col = 1, n_valence_pseudo_orb)
  end do

  print*, ''
  print*, list_valence_pseudo
  print*, ''

  do row = 1,ao_num 
    write(*,'(100F12.7)'), (mo_val_coef(row,col), col = 1, n_valence_pseudo_orb)
  end do

  print*, ''
  print*, 'C1.S0.C1T'
  call dgemm("N","N", &
        & ao_num, n_valence_pseudo_orb, ao_num, &
        & 1.d0, ao_overlap, size(ao_overlap,1), &
        & mo_val_coef, size(mo_val_coef,1), 0.d0, &
        & temp, size(temp,1)&
  )
  temp2 = 0.d0
  call dgemm("T","N", &
        & n_valence_pseudo_orb, n_valence_pseudo_orb, ao_num, &
        & 1.d0, mo_val_coef, size(mo_val_coef,1), &
        & temp, size(temp,1), 0.d0, &
        & temp2, size(temp2,1)&
  )
  do row = 1,size(temp2,1)
    write(*,'(100F12.7)'), (temp2(row,col), col = 1,size(temp2,2))
  end do
  
  print*, ''
  print*, 'C0.S0.C0T'
  call dgemm("N","N", &
        & ao_num, mo_num, ao_num, &
        & 1.d0, ao_overlap, size(ao_overlap,1), &
        & mo_val_coef, size(mo_val_coef,1), 0.d0, &
        & temp, size(temp,1)&
  )
  temp2 = 0.d0
  call dgemm("T","N", &
        & mo_num, mo_num, ao_num, &
        & 1.d0, mo_coef, size(mo_val_coef,1), &
        & temp, size(temp,1), 0.d0, &
        & temp2, size(temp2,1)&
  )
  print*, ''
  do row = 1,size(temp2,1)
    write(*,'(100F12.7)'), (temp2(row,col), col = 1,size(temp2,2))
  end do













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

end subroutine
