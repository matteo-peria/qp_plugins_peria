program core_exchange_potential_test

  BEGIN_DOC
  ! TODO : Put the documentation of the program here
  END_DOC

  implicit none

!  my_grid_becke  = .True.
!  PROVIDE tc_grid1_a tc_grid1_r
!  my_n_pt_r_grid = tc_grid1_r
!  my_n_pt_a_grid = tc_grid1_a
!  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

!  call write_int(6, my_n_pt_r_grid, 'radial  external grid over')
!  call write_int(6, my_n_pt_a_grid, 'angular external grid over')

!  my_extra_grid_becke  = .True.
!  PROVIDE tc_grid2_a tc_grid2_r
!  my_n_pt_r_extra_grid = tc_grid2_r
!  my_n_pt_a_extra_grid = tc_grid2_a
!  touch my_extra_grid_becke my_n_pt_r_extra_grid my_n_pt_a_extra_grid
!
!  call write_int(6, my_n_pt_r_extra_grid, 'radial  internal grid over')
!  call write_int(6, my_n_pt_a_extra_grid, 'angular internal grid over')
!
  call main()

end program

! ---

subroutine main()
  implicit none
  double precision :: difference
  integer :: i,j,ii,jj
  
  print*,'NUMERIC'
  do i = 1, n_act_orb
   ii = list_act(i)
   write(*,'(100(F16.10,X))') core_exchange_pot_numeric(ii,n_core_orb+1:mo_num)
  enddo

  print*,'EXACT'
  do i = 1, n_act_orb
   ii = list_act(i)
   write(*,'(100(F16.10,X))') core_exchange_pot_exact(ii,n_core_orb+1:mo_num)
!   do j = 1, n_act_orb
!    jj = list_act(j)
!    numeric += dabs(x_pot_mo_prov(ii,jj) - v_tmp(ii,jj))
!   enddo
  enddo

  print*,'difference = ', sum(abs(core_exchange_pot_numeric(2:mo_num,2:mo_num)-&
                                  core_exchange_pot_exact(2:mo_num,2:mo_num)))

end subroutine
