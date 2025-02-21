 BEGIN_PROVIDER [ double precision, mo_core_pseudo_coef   , (ao_num, n_core_pseudo_orb) ]
&BEGIN_PROVIDER [ double precision, mo_valence_pseudo_coef, (ao_num, n_valence_pseudo_orb) ]
  implicit none
  BEGIN_DOC
  ! Core MOs' coefficients on |AO| basis set
  ! Valence MOs' coefficients on |AO| basis set
  !
  ! mo_core_pseudo_coef(i,j) = coefficient of the i-th |AO| on the jth core |MO|
  ! mo_valence_pseudo_coef(i,j) = coefficient of the i-th |AO| on the jth valence |MO|
  END_DOC
  logical                        :: exists
  mo_core_pseudo_coef = mo_coef(:,list_core_pseudo)
  mo_valence_pseudo_coef = mo_coef(:,list_valence_pseudo)  
END_PROVIDER

 BEGIN_PROVIDER [ double precision, mo_core_pseudo_coef_sphe   , (ao_cart_to_sphe_num, n_core_pseudo_orb) ]
&BEGIN_PROVIDER [ double precision, mo_valence_pseudo_coef_sphe, (ao_cart_to_sphe_num, n_valence_pseudo_orb) ]
  implicit none
  BEGIN_DOC
  ! Core MOs' coefficients on |AO| basis set
  ! Valence MOs' coefficients on |AO| basis set
  !
  ! mo_core_pseudo_coef(i,j) = coefficient of the i-th |AO| on the jth core |MO|
  ! mo_valence_pseudo_coef(i,j) = coefficient of the i-th |AO| on the jth valence |MO|
  END_DOC
  mo_core_pseudo_coef_sphe = mo_coef_ao_sphe(:,list_core_pseudo)
  mo_valence_pseudo_coef_sphe = mo_coef_ao_sphe(:,list_valence_pseudo)  
END_PROVIDER

 BEGIN_PROVIDER [ double precision, P_mo_pp_valence_sphe, (ao_cart_to_sphe_num, ao_cart_to_sphe_num) ]
&BEGIN_PROVIDER [ double precision, P_mo_pp_valence_sphe_transp, (ao_cart_to_sphe_num, ao_cart_to_sphe_num) ]
  BEGIN_DOC
  ! Porjector P restricted to the pseudo-potential valence_sphe MOs ON THE CARTESIAN AO BASIS
  !
  ! P_{ij} = \sum_a^{N_val} C_{ia}C*_{ja} 
  !        = \sum_a^{N_val} C_{ia}C^{\dag}_{aj} 
  !        = [CC^{\dag}]_{ij}
  END_DOC 
  integer :: i,j
  !
  P_mo_pp_valence_sphe(:,:) = 0.d0
  call dgemm('N',                    & 
             'T',                    & 
             ao_cart_to_sphe_num,                 & 
             ao_cart_to_sphe_num,                 & 
             n_valence_pseudo_orb,   & 
             1.d0,                   & 
             mo_valence_pseudo_coef_sphe, & 
             ao_cart_to_sphe_num,                 & 
             mo_valence_pseudo_coef_sphe, & 
             ao_cart_to_sphe_num,                 & 
             1.d0,                   & 
             P_mo_pp_valence_sphe,        & 
             ao_cart_to_sphe_num                  & 
  )
  !
  ! Compute transpose
  do i = 1,ao_cart_to_sphe_num
   do j = 1,ao_cart_to_sphe_num
     P_mo_pp_valence_sphe_transp(j,i) = P_mo_pp_valence_sphe(i,j)
   enddo
  enddo
END_PROVIDER 




 BEGIN_PROVIDER [ double precision, P_mo_pp_valence, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, P_mo_pp_valence_transp, (ao_num, ao_num) ]
  BEGIN_DOC
  ! Porjector P restricted to the pseudo-potential valence MOs ON THE CARTESIAN AO BASIS
  !
  ! P_{ij} = \sum_a^{N_val} C_{ia}C*_{ja} 
  !        = \sum_a^{N_val} C_{ia}C^{\dag}_{aj} 
  !        = [CC^{\dag}]_{ij}
  END_DOC 
  integer :: i,j
  !
  P_mo_pp_valence(:,:) = 0.d0
  call dgemm('N',                    & 
             'T',                    & 
             ao_num,                 & 
             ao_num,                 & 
             n_valence_pseudo_orb,   & 
             1.d0,                   & 
             mo_valence_pseudo_coef, & 
             ao_num,                 & 
             mo_valence_pseudo_coef, & 
             ao_num,                 & 
             1.d0,                   & 
             P_mo_pp_valence,        & 
             ao_num                  & 
  )
  !
  ! Compute transpose
  do i = 1,ao_num
   do j = 1,ao_num
     P_mo_pp_valence_transp(j,i) = P_mo_pp_valence(i,j)
   enddo
  enddo
END_PROVIDER 


 BEGIN_PROVIDER [double precision, P_mo_pp_valence_sphe_old, (ao_cart_to_sphe_num, ao_cart_to_sphe_num) ]
&BEGIN_PROVIDER [double precision, P_mo_pp_valence_sphe_old_transp, (ao_cart_to_sphe_num, ao_cart_to_sphe_num) ]
  BEGIN_DOC
  ! Projector P restricted to the pseudo-potential valence MOs ON THE SHPERICAL AO BASIS
  !
  ! P_{sph} = C^TP_{cart}C
  !
  ! where C is ao_cart_to_sphe_coef
  END_DOC
  implicit none
  integer :: i,j
  double precision, allocatable :: S(:,:)
  allocate (S(ao_cart_to_sphe_num,ao_num))
  !
  P_mo_pp_valence_sphe_old = 0.d0
  P_mo_pp_valence_sphe_old_transp = 0.d0
  !
  call ao_cart_to_ao_sphe(P_mo_pp_valence,size(P_mo_pp_valence,1), & 
                          P_mo_pp_valence_sphe_old,size(P_mo_pp_valence_sphe_old,1))
! write(*,*) "ao_cart_to_sphe_num", ao_cart_to_sphe_num
! !
! ! Compute left side C^TP_{cart}=S
! call dgemm('T', 'N',                                  &
!   ao_cart_to_sphe_num, ao_num, ao_num, 1.d0,          &
!   ao_cart_to_sphe_coef, size(ao_cart_to_sphe_coef,1), &
!   P_mo_pp_valence, size(P_mo_pp_valence,2),           &
!   0.d0, S, size(S,1)                                  &
! )
! !
! ! Compute remaining right side SC
! call dgemm('N', 'N',                                          &
!     ao_cart_to_sphe_num, ao_cart_to_sphe_num, ao_num, 1.d0,   &
!     S, size(S,1),                                             &
!     ao_cart_to_sphe_coef, size(ao_cart_to_sphe_coef,2),       &
!     0.d0, P_mo_pp_valence_sphe_old, size(P_mo_pp_valence_sphe_old,1)  &
! )
! deallocate(S)
  !
  ! Compute transpose
  do i = 1,ao_cart_to_sphe_num
    do j = 1,ao_cart_to_sphe_num
      P_mo_pp_valence_sphe_old_transp(j,i) = P_mo_pp_valence_sphe_old(i,j)
    enddo
  enddo
END_PROVIDER 


 BEGIN_PROVIDER [ double precision, ao_pp_valence_sphe_coef, (ao_cart_to_sphe_num, ao_cart_to_sphe_num) ]
&BEGIN_PROVIDER [ double precision, ao_pp_valence_sphe_coef_transp, (ao_cart_to_sphe_num, ao_cart_to_sphe_num) ]
  BEGIN_DOC
  ! AOs coefficients to express the pseudopot valence molecular orbitals
  !
  ! C^v_{ij} = [CC^{\dag}S]_{ij} =  [PS]_{ij}
  !
  ! where C is mo_core_pseudo_coef,
  !       P is P_mo_pp_valence_sphe in spherical coordinates
  !       S is ao_cart_to_sphe_overlap
  END_DOC 
  implicit none
  integer :: i_ao, j_ao, i, j
  !
  ao_pp_valence_sphe_coef(:,:) = 0.d0
  call dgemm('N','N',                                                      &
      ao_cart_to_sphe_num, ao_cart_to_sphe_num, ao_cart_to_sphe_num, 1.d0, &
      P_mo_pp_valence_sphe, ao_cart_to_sphe_num,                           &
      ao_cart_to_sphe_overlap, ao_cart_to_sphe_num,                        &
      0.d0, ao_pp_valence_sphe_coef, ao_cart_to_sphe_num                        &
  )
  !
  ! Compute transpose
  do i = 1,ao_cart_to_sphe_num
    do j = 1,ao_cart_to_sphe_num
      ao_pp_valence_sphe_coef_transp(j,i) = ao_pp_valence_sphe_coef(i,j)
    enddo
  enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_pp_overlap_sphe, (ao_cart_to_sphe_num, ao_cart_to_sphe_num) ]
  implicit none
  BEGIN_DOC
  ! AOs overlap matrix in the valence subspace
  END_DOC 
  double precision, dimension(ao_cart_to_sphe_num,ao_cart_to_sphe_num) :: SCCtS
  integer :: i
  !
  call get_AB_prod(                                                     &
      ao_cart_to_sphe_overlap,ao_cart_to_sphe_num,ao_cart_to_sphe_num,  &
      ao_pp_valence_sphe_coef,ao_cart_to_sphe_num,                           &
      SCCtS                                                             &
  )
  call get_AB_prod(                                                     & 
      ao_pp_valence_sphe_coef_transp,ao_cart_to_sphe_num,ao_cart_to_sphe_num,&
      SCCtS,ao_cart_to_sphe_num,                                        &
      ao_pp_overlap_sphe                                                     &
  )
END_PROVIDER


 BEGIN_PROVIDER [ double precision, ao_pp_overlap_sphe_eigval, (ao_cart_to_sphe_num) ]
&BEGIN_PROVIDER [ double precision, ao_pp_overlap_sphe_eigvec, (ao_cart_to_sphe_num, ao_cart_to_sphe_num) ]
  implicit none
  BEGIN_DOC
  ! AOs overlap matrix in the valence subspace
  END_DOC 
  call lapack_diagd(ao_pp_overlap_sphe_eigval,ao_pp_overlap_sphe_eigvec,ao_pp_overlap_sphe,ao_cart_to_sphe_num,ao_cart_to_sphe_num)
END_PROVIDER 


!!!OLD PROVIDERS BASED ON CARTESIAN AOs
 BEGIN_PROVIDER [ double precision, ao_pp_valence_coef, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_pp_valence_coef_transp, (ao_num, ao_num) ]
  BEGIN_DOC
  ! AOs coefficients to express the pseudopot valence molecular orbitals
  !
  ! C^v_{ij} = [CC^{\dag}S]_{ij} =  [PS]_{ij}
  !
  ! where C is mo_core_pseudo_coef,
  !       P is P_mo_pp_valence
  !       S is ao_overlap
  END_DOC 
  implicit none
  integer :: i_ao, j_ao, i, j
  !
  ao_pp_valence_coef(:,:) = 0.d0
  call dgemm('N','N',ao_num,ao_num,ao_num,1.d0,       &
      P_mo_pp_valence,ao_num,        &
      ao_overlap,ao_num,                  &
      0.d0,ao_pp_valence_coef,ao_num &
  )
  !
  ! Compute transpose
  do i = 1,ao_num
    do j = 1,ao_num
      ao_pp_valence_coef_transp(j,i) = ao_pp_valence_coef(i,j)
    enddo
  enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_pp_overlap, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! AOs overlap matrix in the valence subspace
  END_DOC 
  double precision, dimension(ao_num,ao_num) :: SCCtS
  integer :: i
  !
  call get_AB_prod(ao_overlap,ao_num,ao_num,ao_pp_valence_coef,ao_num,SCCtS)
  call get_AB_prod(ao_pp_valence_coef_transp,ao_num,ao_num,SCCtS,ao_num,ao_pp_overlap)
END_PROVIDER


 BEGIN_PROVIDER [ double precision, ao_pp_overlap_eigval, (ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_pp_overlap_eigvec, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! AOs overlap matrix in the valence subspace
  END_DOC 
  call lapack_diagd(ao_pp_overlap_eigval,ao_pp_overlap_eigvec,ao_pp_overlap,ao_num,ao_num)
END_PROVIDER 
