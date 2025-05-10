program core_xpot_test_new

  BEGIN_DOC
  ! Testing best resolution of extra and adaptive grid 
  END_DOC

  implicit none

  print*, 'IM IN THE MAINNNNNNNNNNNNNNNNNNNN'

  !touch my_grid_becke !my_extra_grid_becke my_grid_adapt

  ! Need to set low number so that when these providers are initialized
  ! in the main subroutine they do not take long time to be created
  my_n_pt_r_grid = 6
  my_n_pt_a_grid = 6
  my_n_pt_r_extra_grid = 6
  my_n_pt_a_extra_grid = 6
  my_n_pt_r_extra_grid = 6
  my_n_pt_a_extra_grid = 6

  my_grid_becke = .True.
  my_extra_grid_becke = .True.
  my_grid_adapt = .True.


  print*, 'exit from THE MAINNNNNNNNNNNNNNNNNNNN'

  call main()

end program

! ---

subroutine main()
  implicit none
  double precision :: difference
  integer :: i,j,ii,jj

  integer :: n_radial_points(6) = [25, 50, 75, 100, 150, 200]
  integer :: n_angular_points(8) = [110, 146, 170, 194, 230, 266, 302, 590]

  ! Typical couples of values
  ! radial     angular 
  !     23        170
  !     50        194
  !     75        302
  !     99        590

  print*, 'LOOP'

  write(*,'(A)') repeat('-', 70)
  integer :: i_r, i_a
  integer :: j_r, j_a
  integer :: k_r, k_a

  do i_r = 1, 1 !3 !size(n_radial_points)
    do i_a = 1, 1 !3 !size(n_angular_points)
      n_points_radial_grid = n_radial_points(i_r)
      n_points_integration_angular = n_angular_points(i_a)

      touch n_points_radial_grid n_points_integration_angular

      print*, 'FIRST GRID'
      provide n_points_radial_grid
      !provide n_points_integration_angular
      !provide n_points_grid_per_atom 
      !provide n_points_final_grid

      !print*, 'N radial points                 = ', n_points_radial_grid
      !print*, 'N angular points                = ', n_points_integration_angular
      !print*, 'N total points (before pruning) = ', n_points_grid_per_atom 
      !print*, 'N total points (after pruning)  = ', n_points_final_grid
      write(*,*) ''

      write(*,'(A)') repeat('=', 70)

!!      do j_r = 1, 3 !size(n_radial_points)
!!        do j_a = 1, 3 !size(n_angular_points)
!!          ! To be changed later with another loop
!!          k_r = j_r
!!          k_a = j_a
!!
!!          !! STANDARD EXTRA GRID
!!          n_points_extra_radial_grid = n_radial_points(j_r)
!!          n_points_extra_integration_angular = n_angular_points(j_a)
!!          touch n_points_extra_radial_grid n_points_extra_integration_angular
!!          print*, 'EXTRA GRID'
!!          !print*, 'N radial points                 = ', n_points_extra_radial_grid
!!          !print*, 'N angular points                = ', n_points_extra_integration_angular
!!          !print*, 'N total points (before pruning) = ', n_points_extra_grid_per_atom 
!!          !print*, 'N total points (after pruning)  = ', n_points_extra_final_grid
!!          write(*,*) ''
!!          write(*,'(A)') repeat('-', 70)
!!
!!          !! ADAPTIVE GRID
!!          n_points_radial_grid_adapt = n_radial_points(k_r)
!!          n_points_integration_angular_adapt = n_angular_points(k_a)
!!          touch n_points_radial_grid_adapt n_points_integration_angular_adapt
!!          !print*, 'ADAPTIVE GRID'
!!          !print*, 'N radial points                 = ', n_points_radial_grid_adapt          
!!          !print*, 'N angular points                = ', n_points_integration_angular_adapt
!!          !write(*,*) ''
!!
!!
!!        end do
!!        write(*,'(A)') repeat('%', 70)
!!      end do
    end do
  end do




!!  print*, 'LOOOP'
!!
!!  write(*,'(A)') repeat('-', 70)
!!  integer :: i_r, i_a
!!  integer :: j_r, j_a
!!  integer :: k_r, k_a
!!  do i_r = 1, 1 !3 !size(n_radial_points)
!!    do i_a = 1, 1 !3 !size(n_angular_points)
!!      my_n_pt_r_grid = n_radial_points(i_r)
!!      my_n_pt_a_grid = n_angular_points(i_a)
!!      touch my_grid_becke !my_n_pt_r_grid my_n_pt_a_grid
!!      
!!      print*, 'FIRST GRID'
!!      print*, 'N radial points                 = ', n_points_radial_grid
!!      print*, 'N angular points                = ', n_points_integration_angular
!!      print*, 'N total points (before pruning) = ', n_points_grid_per_atom 
!!      print*, 'N total points (after pruning)  = ', n_points_final_grid
!!
!!      write(*,'(A)') repeat('-', 70)
!!
!!      do j_r = 1, 3 !size(n_radial_points)
!!        do j_a = 1, 3 !size(n_angular_points)
!!          ! To be changed later with another loop
!!          k_r = j_r
!!          k_a = j_a
!!
!!          ! STANDARD EXTRA GRID
!!          my_n_pt_r_extra_grid = n_radial_points(j_r)
!!          my_n_pt_a_extra_grid = n_angular_points(j_a)
!!          touch my_extra_grid_becke
!!          
!!          print*, 'EXTRA_GRID WITH PRUNING'
!!          print*, 'N radial points extra grid      = ', n_points_extra_radial_grid
!!          print*, 'N angular points extra grid     = ', n_points_extra_integration_angular 
!!          print*, 'N total points (before pruning) = ', n_points_extra_grid_per_atom 
!!          print*, 'N total points (after pruning)  = ', n_points_extra_final_grid
!!
!!          difference = sum(abs(core_xpot_numeric(2:mo_num,2:mo_num)-&
!!                             & core_xpot_exact(2:mo_num,2:mo_num)))
!!          print*, 'Difference with exact result    = ', difference
!!
!!          write(*,'(A)') repeat('=', 70)
!!          !print*, 'EXTRA_GRID WITHOUT PRUNING'
!!
!!
!!          !! STANDARD ADAPTIVE GRID
!!          !my_n_pt_r_grid_adapt = n_radial_points(k_r)
!!          !my_n_pt_a_grid_adapt = n_angular_points(k_a)
!!          !touch my_grid_adapt
!!          !
!!          !print*, 'ADAPTIVE GRID WITHOUT PRUNING'
!!          !print*, 'N radial points adaptive grid      = ', n_points_radial_grid_adapt 
!!          !print*, 'N angular points adaptive grid     = ', n_points_integration_angular_adapt
!!          !print*, 'N total points (no pruning) = ', n_total_adapt_grid 
!!
!!          !difference = sum(abs(core_xpot_numeric_adapt_old(2:mo_num,2:mo_num)-&
!!          !                   & core_xpot_exact(2:mo_num,2:mo_num)))
!!          !print*, 'Difference with exact result    = ', difference
!!        end do
!!        write(*,'(A)') repeat('=', 70)
!!      end do
!!    end do
!!  end do


end subroutine
