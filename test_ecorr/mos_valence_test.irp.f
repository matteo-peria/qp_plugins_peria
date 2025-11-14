program mos_valence_test
  implicit none
  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  !call old_test()
  call mos_loop_vs_array_slicing

end program mos_valence_test

subroutine mos_loop_vs_array_slicing
  implicit none
  double precision :: diff

  double precision :: mo_core_coef_normed(ao_num, n_core_pseudo_orb)
  double precision :: mo_val_coef_normed(ao_num, n_valence_pseudo_orb)

  double precision :: norm_core_mos(n_core_pseudo_orb)
  double precision :: norm_valence_mos(n_valence_pseudo_orb)
  norm_core_mos(:)    = norm2(mo_core_coef_notnorm(:,:), dim=1)
  norm_valence_mos(:) = norm2(mo_val_coef_notnorm(:,:), dim=1)

  integer :: row
  do row = 1, ao_num 
    mo_core_coef_normed(row,:) = mo_core_coef_notnorm(row,:) / norm_core_mos(:)    
    mo_val_coef_normed(row,:)  = mo_val_coef_notnorm(row,:)  / norm_valence_mos(:) 
  end do

  diff = sum(abs(mo_core_coef_normed - mo_core_coef))
  write(*,*) "Difference MOs core from array slicing and loops: ", diff 
  diff = sum(abs(mo_val_coef_normed - mo_val_coef))
  write(*,*) "Difference MOs valence from array slicing and loops: ", diff 

end subroutine mos_loop_vs_array_slicing


subroutine old_test
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
