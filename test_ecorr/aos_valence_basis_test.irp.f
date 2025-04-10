program test_aos_valence_basis
  implicit none

  call main()
end program test_aos_valence_basis


subroutine main()
  implicit none


!!
!!  ! START TESTING VALENCE PRIMITIVES
!!  provide ao_val_prim_num_temp
!!  write(*,*) "AO_PP_PRIM_NUM_TEMP"
!!  write(*,*) ao_val_prim_num_temp
!!  print*, ''
!!
!!  provide ao_val_prim_num_max_temp
!!  write(*,*) "AO_PP_PRIM_NUM_MAX_TEMP"
!!  write(*,*) ao_val_prim_num_max_temp
!!  print*, ''
!!
!!  provide ao_val_prim_num
!!  write(*,*) "AO_PP_PRIM_NUM"
!!  write(*,*) ao_val_prim_num
!!  print*, ''
!!
!!  provide ao_val_prim_num_max
!!  write(*,*) "AO_PP_PRIM_NUM_MAX"
!!  write(*,*) ao_val_prim_num_max
!!  print*, ''
!!
!!  provide ao_expo_transp
!!  write(*,*) "AO_EXPO_TRANSP, shape = ", shape(ao_expo_transp)
!!  do row = 1, ao_prim_num_max
!!    write(*,'(100F12.7)'), (ao_expo_transp(row,col), col = 1, ao_num)
!!  end do
!!  print*, ''
!!
!!  provide ao_val_expo_transp_temp
!!  write(*,*) "AO_PP_EXPO_TRANSP_TEMP, shape = ", shape(ao_val_expo_transp_temp)
!!  do row = 1, ao_val_prim_num_max_temp
!!    write(*,'(100F12.7)'), (ao_val_expo_transp_temp(row,col), col = 1, ao_num)
!!  end do
!!  print*, ''
!!
!!  provide ao_val_expo_transp
!!  write(*,*) "AO_PP_EXPO_TRANSP, shape = ", shape(ao_val_expo_transp)
!!  do row = 1, ao_val_prim_num_max
!!    write(*,'(100F12.7)'), (ao_val_expo_transp(row,col), col = 1, ao_num)
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
!!  do row = 1, ao_val_prim_num_max_temp
!!    write(*,'(100F12.7)'), (ao_pp_coef_transp_temp(row,col), col = 1, ao_num)
!!  end do
!!  print*, ''
!!
!!  provide ao_pp_coef_transp
!!  write(*,*) "AO_PP_COEF_TRANSP, shape = ", shape(ao_pp_coef_transp)
!!  do row = 1, ao_val_prim_num_max
!!    write(*,'(100F12.7)'), (ao_pp_coef_transp(row,col), col = 1, ao_num)
!!  end do
!!  print*, ''
!!  ! END TESTING VALENCE PRIMITIVES
 

end subroutine main
