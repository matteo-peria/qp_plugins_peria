 ! SUBSTITUTE OF GET_VORONOI_PARTITION 
 subroutine get_adaptive_grid_becke_weights_nosilence(rp, float_grid, fixed_grid, &
                                          & float_grid_becke_boundary_shift, &
                                          & becke_weights_at_float_grid, &
                                          & becke_weights_at_fixed_grid)
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
  double precision, intent(in) :: rp(3)
  double precision, intent(in) :: float_grid(3, n_points_ang_float_grid, n_points_rad_float_grid,1)
  double precision, intent(in) :: fixed_grid(3, n_points_ang_extra_grid, n_points_rad_extra_grid, nucl_num)
  double precision, intent(in) :: float_grid_becke_boundary_shift(nucl_num)
  !
  double precision, intent(out):: becke_weights_at_float_grid(n_points_ang_float_grid, n_points_rad_float_grid, 1)
  double precision, intent(out):: becke_weights_at_fixed_grid(n_points_ang_extra_grid, n_points_rad_extra_grid, nucl_num)
  integer :: i_nucl,k,l
  !
  double precision :: r(3)
  double precision :: weights_per_atom(nucl_num+1)
  
  ! Loop over all points belonging to real atoms grid (fixed grid)
  do i_nucl = 1, nucl_num
    do k = 1, n_points_rad_extra_grid - 1
      do l = 1, n_points_ang_extra_grid
        r(1:3) = fixed_grid(1:3,l,k,i_nucl)
        call get_becke_functions_general2(r,rp,float_grid_becke_boundary_shift,weights_per_atom)
        !becke_weights(l,k,i_nucl+1) = weights_per_atom(i_nucl+1)
        !becke_weights_at_fixed_grid(l,k,i_nucl) = weights_per_atom(i_nucl+1)
        becke_weights_at_fixed_grid(l,k,i_nucl) = weights_per_atom(i_nucl+1)
        !write(*,'(3I4,E13.6)') i_nucl, k, l, becke_weights_at_fixed_grid(l,k,i_nucl)
      enddo
    enddo
  enddo

 ! Loop over all points belonging to ghost-atom grid (floating grid)
 i_nucl = 1
 !for each radial grid attached to the "jth" atom
 do k = 1, n_points_rad_float_grid - 1
   ! for each angular point attached to the "jth" atom
   do l = 1, n_points_ang_float_grid 
     r(1:3) = float_grid(1:3,l,k,i_nucl)

     call get_becke_functions_general2(r,rp,float_grid_becke_boundary_shift,weights_per_atom)
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

end subroutine


subroutine get_becke_functions_general2_nosilence( r                      &
                                                 , r_input                &
                                                 , slater_inter_per_input &
                                                 , weights_per_atom       &
                                                 )
 implicit none
 BEGIN_DOC
 !
 ! Computes the Becke-weight associated to each atom, 'weights_per_atom',
 ! at the point in real space 'r'.
 ! The 3D space has been already Voronoi-partitioned in atom-centered Becke-cells
 ! plus a new one associated to a ghost atom centered in 'r_input'.
 ! The new Becke-cell is characterized by its center (r_input) and its
 ! a_ij parameter defining offset Voronoi boundaries, 'slater_inter_per_input'
 !
 ! weights_per_atom(i_nucl) = P(r,i_nucl)/sum_k P(r,k) 
 ! BUT NOW sum_k P(r,k) = P(r,new_atom) + \sum_k = 1, nucl_num  P(r,k)
 !
 ! function defined for each atom "i" by equation (13) and (21) with k == 3
 !
 ! in terms of nucl_num atoms + 1 additional atom characterized by "r_input" and "slater_inter_per_input"
 !
 ! returns the following partition :: weights_per_atom(i) == w(r,i) =  P(r,i) / \sum_{k in partition } P(r,k)
 !
 ! WARNING the first entry corresponds to the fake atom "r_input"
 END_DOC
 !
 double precision, intent(in) :: r(3)
 double precision, intent(in) :: r_input(3)
 double precision, intent(in) :: slater_inter_per_input(nucl_num)
 !
 double precision, intent(out):: weights_per_atom(nucl_num+1)
 !
 double precision :: cell_function_becke_extra_atom2
 double precision :: cell_function_becke_general_adapt2_nosilence
 integer :: i
 
 ! Weight of the extra ghost-atom (located in r_input) on the point r
 ! This correspond to P_1(r) where 1 is the ghost-Becke-cell (Eq. 13 of Becke paper)
 weights_per_atom(1) = cell_function_becke_extra_atom2(r, r_input, slater_inter_per_input)
 
 ! Weights of all the other real atoms (located in i) on the point 'r'
 do i = 1, nucl_num
   ! Compute P_i(r) for 'i' belonging to real nucleus (Eq. 13 of Becke paper)
   weights_per_atom(i+1) = cell_function_becke_general_adapt2_nosilence(r, i, r_input, slater_inter_per_input) 
 enddo

 if (any(isnan(weights_per_atom))) then
   print*, 'AAAAAAAH THERE IS SOME NAN IN WEIGHTS PER ATOM'
   print*, 'weights_per_atom: '
   print*, weights_per_atom
   print*, 'sum(weights_per_atom): ', sum(weights_per_atom)
   stop
 endif

 ! Normalisation computed as w_i(r) = P_i(r) / \sum_n P_n(r) (Eq. 22 Becke paper)
 !print*, 'im in get_becke_function_general2'
 !print*, 'sum(weights_per_atom) = ', sum(weights_per_atom)
 weights_per_atom = weights_per_atom/sum(weights_per_atom)

end subroutine


double precision function cell_function_becke_general_adapt2_nosilence( r                       &
                                                                       , atom_number            &
                                                                       , rp                     &
                                                                       , slater_inter_per_input &
                                                                       )
  BEGIN_DOC
  ! No silence version of
  !   cell_function_becke_general_adapt2
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
  
  cell_function_becke_general_adapt2_nosilence = 1.d0

  ! Distance of the point 'r' from the nucleus belonging to the same Voronoi cell
  r_i = dsqrt(sum((r(1:3) - nucl_coord_transp(1:3,atom_number)) &
               &* (r(1:3) - nucl_coord_transp(1:3,atom_number)))) 

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
    cell_function_becke_general_adapt2_nosilence *= step_function_becke(nu_ij)
    !write(*,'(I4,100E13.4)'), j, r_j, mu_ij, nu_ij, cell_function_becke_general_adapt2_nosilence
  enddo

  ! Contribute coming from the extra atom positioned at 'rp'
  r_j = dsqrt(sum((r(1:3) - rp(1:3)) * (r(1:3) - rp(1:3)))) 
  ! Internuclear distance
  R_ij = dsqrt(sum((rp(1:3) - nucl_coord_transp(1:3,atom_number) ) * (rp(1:3) - nucl_coord_transp(1:3,atom_number))))
  ! Elliptical coordinate (Becke Eq. 11)
  mu_ij = (r_i - r_j) / R_ij
  ! Heteronuclear correction for off-set Voronoi walls (Becke Eq. A2)
  nu_ij = mu_ij + slater_inter_per_input(atom_number) * (1.d0 - mu_ij*mu_ij)
  ! Update Becke-weight with j-th atom contribution
  cell_function_becke_general_adapt2_nosilence *= step_function_becke(nu_ij)
  !write(*,'(100E13.4)'), r_j, mu_ij, nu_ij, cell_function_becke_general_adapt2_nosilence
 
  return
end function
