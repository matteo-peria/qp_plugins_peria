 ! SUBSTITUTE OF GET_VORONOI_PARTITION 
 subroutine get_adaptive_grid_becke_weights3( rp, float_grid, fixed_grid     &
                                           , float_grid_becke_boundary_shift &
                                           , becke_weights_at_float_grid     &
                                           , becke_weights_at_fixed_grid)

  implicit none
  BEGIN_DOC
  !
  ! Computes the Becke-cell weights for a set of points in a grid obtained as 
  ! the union of 'float_grid' and 'fixed_grid', the first consisting in 1 cell
  ! and the second made of nucl_num cells.
  ! Therefore, the Voronoi partition is made of (nucl_num+1) cells, 
  ! where the nucl_num corresponds to to the total number of real atom presents
  ! in the system (corresponding to the 'fixed' part of the grid), 
  ! while +1 stands for an additional spherical cell centered around 'rp', 
  ! treated as an extra ghost-atom (corresponding to the 'floating' grid).
  !
  !extra atom at rp and characterized by float_grid_becke_boundary_shift
  END_DOC
  double precision, intent(in)  :: rp(3)
  double precision, intent(in)  :: float_grid(3, n_points_ang_float_grid, n_points_rad_float_grid,1)
  double precision, intent(in)  :: fixed_grid(3, grid3_ang_size, grid3_rad_size, nucl_num)
  double precision, intent(in)  :: float_grid_becke_boundary_shift(nucl_num)
  !
  double precision, intent(out) :: becke_weights_at_float_grid(n_points_ang_float_grid, n_points_rad_float_grid, 1)
  double precision, intent(out) :: becke_weights_at_fixed_grid(grid3_ang_size, grid3_rad_size, nucl_num)
  integer :: i_nucl,k,l
  !
  double precision :: r(3)
  double precision :: weights_per_atom(nucl_num+1)


  ! Initialize
  becke_weights_at_fixed_grid(:,:,:) = 0.d0
  becke_weights_at_float_grid(:,:,:) = 0.d0
  
  ! Loop over all points belonging to real atoms grid (fixed grid)
  !print*,'fixed'
  do i_nucl = 1, nucl_num
    do k = 1, grid3_rad_size - 1
      do l = 1, grid3_ang_size 
        r(1:3) = fixed_grid(1:3,l,k,i_nucl)
        !print*, "CALL GET_BECKE_FUNCTIONS_GENERAL2(R,RP,FLOAT_GRID_BECKE_BOUNDARY_SHIFT, NUCL_NUM+1, WEIGHTS_PER_ATOM)"
        print*, "call get_becke_functions_general3 at i_nucl, k, l = ", i_nucl, k, l
        call get_becke_functions_general3(r,rp,float_grid_becke_boundary_shift, nucl_num+1, weights_per_atom)
        !becke_weights(l,k,i_nucl+1) = weights_per_atom(i_nucl+1)
        !becke_weights_at_fixed_grid(l,k,i_nucl) = weights_per_atom(i_nucl+1)
        becke_weights_at_fixed_grid(l,k,i_nucl) = weights_per_atom(i_nucl+1)
         
      enddo
    enddo
  enddo

  if (silence_floating_grid.eq..true.) then
    ! No becke weights are computed for the floating grid, which is suppressed
  else
    ! Loop over all points belonging to ghost-atom grid (floating grid)
    !print*,'floating'
    i_nucl = 1
    do k = 1, n_points_rad_float_grid - 1
      do l = 1, n_points_ang_float_grid 
        r(1:3) = float_grid(1:3,l,k,i_nucl)
        !print*, "CALL GET_BECKE_FUNCTIONS_GENERAL2(R,RP,FLOAT_GRID_BECKE_BOUNDARY_SHIFT, NUCL_NUM+1, WEIGHTS_PER_ATOM)"
        print*, "call get_becke_functions_general3 at i_nucl, k, l = ", i_nucl, k, l
        call get_becke_functions_general3(r,rp,float_grid_becke_boundary_shift, nucl_num+1, weights_per_atom)
        becke_weights_at_float_grid(l,k,i_nucl) = weights_per_atom(i_nucl)
         !write(*,'(3I4,E13.6)') i_nucl, k, l, becke_weights_at_float_grid(l,k,i_nucl)
        
        if(isnan(becke_weights_at_float_grid(l,k,i_nucl))) then
          print*,'isnan(becke_weights_at_float_grid(l,k,i_nucl))'
          print*,l,k,i_nucl
          print*, becke_weights_at_float_grid(l,k,i_nucl)
          stop                                                                                        
        endif                                                                   
      enddo
    enddo
  end if
end subroutine


subroutine get_becke_functions_general3(r, rp, slater_inter_per_input, n_weights, weights)
 implicit none
 BEGIN_DOC
 !
 ! Computes the Becke-weight associated to each atom, 'weights',
 ! at the point in real space 'r'.
 ! The 3D space has been already Voronoi-partitioned in atom-centered Becke-cells
 ! plus a new one associated to a ghost atom centered in 'rp'.
 ! The new Becke-cell is characterized by its center (rp) and its
 ! a_ij parameter defining offset Voronoi boundaries, 'slater_inter_per_input'.
 !
 !  weights(i_nucl) = P(r,i_nucl)/sum_k P(r,k) 
 !  BUT NOW \sum_k' P(r,k) = P(r,new_atom) + \sum_k P(r,k)
 !
 ! WARNING the first entry corresponds to the fake atom "rp"
 ! This other version generalize the previous one, now named 
 !
 !   get_becke_functions_general2_nosilence
 !
 ! so that we can silence the floating grid and compute only the becke weights 
 ! on the usual fixed grid
 END_DOC
 !
 double precision, intent(in) :: r(3)
 double precision, intent(in) :: rp(3)
 double precision, intent(in) :: slater_inter_per_input(nucl_num)
 integer, intent(in) :: n_weights
 !
 !double precision, intent(out):: weights(nucl_num+1) !generalized with the following:
 double precision, intent(out):: weights(n_weights)
 !
 double precision :: cell_function_becke_extra_atom3
 double precision :: cell_function_becke_general_adapt3
 integer :: i

 integer :: fixed_grid_start

 !print*, "N_WEIGHTS ===================================", n_weights

 if (silence_floating_grid.eq..true.) then
   weights(1) = 0.d0
!!   print*, "Floating grid has been silenced. Now the adaptive_grid3 should be equal to grid3"
!!   fixed_grid_start = 1
 else
!!   fixed_grid_start = 2
   ! Weight of the extra ghost-atom (located in rp) on the point r
   ! This correspond to P_1(r) where 1 is the ghost-Becke-cell (Eq. 13 of Becke paper)
   !>print*, "weights(1) = cell_function_becke_extra_atom3"
   weights(1) = cell_function_becke_extra_atom3(r, rp, slater_inter_per_input)
  !print*,'distance',  norm2(r-rp), 'weights(1)', weights(1)
   
 end if
 
 ! Weights of all the other real atoms (located in i) on the point 'r'
 do i = 1, nucl_num
 !do i = fixed_grid_start, nucl_num
   !! Compute P_i(r) for 'i' belonging to real nucleus (Eq. 13 of Becke paper)
   !>print*, "weights(i+1) = cell_function_becke_general_adapt3"
   weights(i+1) = cell_function_becke_general_adapt3(r, i, rp, slater_inter_per_input) 
   !weights(i) = cell_function_becke_general_adapt3(r, i-1, rp, slater_inter_per_input) 
 enddo

 if (any(isnan(weights))) then
   print*, 'WARNING: array weights has some NaN'
   print*, 'weights: '
   print*, weights
   print*, 'sum(weights): ', sum(weights)
   stop
 endif

 ! Normalisation computed as w_i(r) = P_i(r) / \sum_n P_n(r) (Eq. 22 Becke paper)
 weights = weights/sum(weights)
end subroutine




double precision function cell_function_becke_extra_atom3(r, rp, slater_inter_perp)
  BEGIN_DOC
  ! 
  ! Return the Becke-weight of a given point 'r', belonging to the Becke-cell 
  ! assigned to a ghost atom located in 'rp',with respect to a 
  ! Becke-partitioning based on a list of atoms supplied with the ghost atom,
  ! whose Slater(-Bragg) radius with respect all the other atoms is stored in 
  ! 'slater_inter_per_input'.
  ! 
  ! This subroutine is used also by the nosilence version of the adaptive grid.
  ! 
  ! See Becke (1988, JCP,88(4))
  ! r(1:3)                       :: x,y,z coordinantes of the current point
  END_DOC

  implicit none
  double precision, intent(in) :: r(3)
  double precision, intent(in) :: rp(3)
  double precision, intent(in) :: slater_inter_perp(nucl_num)
  !
  integer           :: j
  double precision  :: mu_ij, nu_ij
  double precision  :: r_i, r_j, R_ij
  double precision :: step_function_becke
  
  cell_function_becke_extra_atom3 = 1.d0

  ! Distance of the point 'r' from 'rp', the center of the ghost atom
  r_i = dsqrt(sum((r(1:3) - rp(1:3) ) * (r(1:3) - rp(1:3)))) 

  ! Loop over all the other real atoms
  do j = 1, nucl_num
    ! Distance of the point 'r' from the j-th atom center
    r_j = dsqrt(sum((r(1:3) - nucl_coord_transp(1:3,j) ) * (r(1:3) - nucl_coord_transp(1:3,j)))) 
    ! Internuclear distance
    R_ij = dsqrt(sum((rp(1:3) - nucl_coord_transp(1:3,j) ) * (rp(1:3) - nucl_coord_transp(1:3,j))))
    ! Elliptical coordinate (Becke Eq. 11)
    mu_ij = (r_i-r_j)/R_ij
    ! Heteronuclear correction for off-set Voronoi walls (Becke Eq. A2)
    nu_ij = mu_ij + slater_inter_perp(j) * (1.d0 - mu_ij*mu_ij)
    ! Update Becke-weight with j-th atom contribution
    cell_function_becke_extra_atom3 *= step_function_becke(nu_ij)
    !
    !write(*,'(I4,100E13.4)'), j, r_j, R_ij, mu_ij, nu_ij, cell_function_becke_extra_atom3
  enddo
  return
end function


double precision function cell_function_becke_general_adapt3( r                      &
                                                            , atom_number            &
                                                            , rp                     &
                                                            , slater_inter_per_input &
                                                            )
  BEGIN_DOC
  !
  ! Return the Becke-weight of a given point 'r', belonging to the Becke-cell 
  ! assigned to an atom whose index in the atom list is 'atom_number',
  ! with respect to a Becke-partitioning based on a list of atoms supplied with
  ! an extra atom located in 'rp' and whose Slater(-Bragg) radius with respect
  ! all the other atoms is stored in 'slater_inter_per_input'
  ! 
  ! See Becke (1988, JCP,88(4))
  ! r(1:3)                       :: x,y,z coordinantes of the current point
  END_DOC

  implicit none
  double precision, intent(in) :: r(3)
  double precision, intent(in) :: rp(3)
  double precision, intent(in) :: slater_inter_per_input(nucl_num)
  integer         , intent(in) :: atom_number
  !
  !OUTPUT
  !
  integer          :: j
  double precision :: r_i, r_j
  double precision :: mu_ij, nu_ij
  double precision :: step_function_becke
  double precision :: r_ij
  
  cell_function_becke_general_adapt3 = 1.d0

  ! Distance of the point 'r' from the nucleus belonging to the same Voronoi cell
  r_i = dsqrt(sum((r(1:3) - nucl_coord_transp(1:3,atom_number)) &
              & * (r(1:3) - nucl_coord_transp(1:3,atom_number)))) 

  ! Loop over all the real atoms for their contribution to the Becke-weight
  do j = 1, nucl_num
    ! Except the one whose Voronoi-cell to which 'r' belongs
    if (j==atom_number) cycle
    ! Distance of the point 'r' from the j-th atom center
    r_j = dsqrt(sum((r(1:3) - nucl_coord_transp(1:3,j)) * (r(1:3) - nucl_coord_transp(1:3,j)))) 
    ! Elliptical coordinate (Becke Eq. 11)
    mu_ij = (r_i - r_j) * nucl_dist_inv(atom_number,j)
    ! Heteronuclear correction for off-set Voronoi walls (Becke Eq. A2)
    nu_ij = mu_ij + slater_bragg_type_inter_distance_ua(atom_number,j) * (1.d0 - mu_ij*mu_ij)
    ! Update Becke-weight with j-th atom contribution
    cell_function_becke_general_adapt3 *= step_function_becke(nu_ij)
    !write(*,'(I4,100E13.4)'), j, r_j, mu_ij, nu_ij, cell_function_becke_general_adapt3
  enddo

  if (silence_floating_grid.eq..false.) then
    ! Contribute coming from the extra atom positioned at 'rp'
    r_j = dsqrt(sum((r(1:3) - rp(1:3)) * (r(1:3) - rp(1:3)))) 
    ! Internuclear distance
    R_ij = dsqrt(sum((rp(1:3) - nucl_coord_transp(1:3,atom_number) ) * (rp(1:3) - nucl_coord_transp(1:3,atom_number))))
    ! Elliptical coordinate (Becke Eq. 11)
    mu_ij = (r_i - r_j) / R_ij
    ! Heteronuclear correction for off-set Voronoi walls (Becke Eq. A2)
    nu_ij = mu_ij + slater_inter_per_input(atom_number) * (1.d0 - mu_ij*mu_ij)
    ! Update Becke-weight with j-th atom contribution
    cell_function_becke_general_adapt3 *= step_function_becke(nu_ij)
    !write(*,'(100(F16.6,X))'),r_j,r_i,R_ij,mu_ij,nu_ij,slater_inter_per_input(atom_number)
    !!write(*,'(100E13.4)'), r_j, mu_ij, nu_ij, cell_function_becke_general_adapt3
  end if
 
  return
end function
