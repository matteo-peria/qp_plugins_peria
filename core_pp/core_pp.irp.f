program core_pp
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC

  my_grid_becke  = .True.
  my_extra_grid_becke  = .True.
  my_grid_adapt        = .True.

  my_n_pt_r_grid       = my_n_pt_r_grid
  my_n_pt_a_grid       = my_n_pt_a_grid

  my_n_pt_r_extra_grid = my_n_pt_r_extra_grid
  my_n_pt_a_extra_grid = my_n_pt_a_extra_grid
  
  my_n_pt_r_grid_adapt = my_n_pt_r_extra_grid
  my_n_pt_a_grid_adapt = my_n_pt_a_extra_grid

  touch my_extra_grid_becke 
  touch my_grid_adapt my_n_pt_r_grid_adapt my_n_pt_a_grid_adapt 

  !call main()

end program
