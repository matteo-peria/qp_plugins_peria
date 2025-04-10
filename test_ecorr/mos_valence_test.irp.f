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
  double precision :: accu = 0.d0
  double precision :: trace = 0.d0

  ! START TESTING VALENCE-SUBSPACE PROJECTOR
  print*, "norm of the valence orbitals"
  do col=1,n_valence_pseudo_orb
    print*, col, sqrt(dot_product(mo_val_coef(:,col),mo_val_coef(:,col)))
  end do

  provide mo_val_coef
  print*, "mo_val_coef"
  do row=1,ao_num
    write(*,'(100F12.7)'), (mo_val_coef(row,col), col=1,n_valence_pseudo_orb)
  end do
  
  provide proj_mo_val
  print*, "proj_mo_val"
  do row=1,ao_num
    write(*,'(100F12.7)'), (proj_mo_val(row,col), col=1,ao_num)
    trace += proj_mo_val(row,row)
  end do
  print*, "Trace: ", trace 
  ! END TESTING VALENCE-SUBSPACE PROJECTOR

end subroutine
