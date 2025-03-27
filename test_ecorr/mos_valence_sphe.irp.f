!BEGIN_PROVIDER [ double precision, SsINV_CsT_temp, (ao_sphe_num, ao_num)]
! implicit none
! call get_AB_prod(ao_cart_to_sphe_overlap_inv,              &
!                  ao_sphe_num, ao_sphe_num, &
!                  ao_cart_to_sphe_coef_transp, ao_num,      &
!                  SsINV_CsT_temp)                          
!END_PROVIDER
!
!BEGIN_PROVIDER [ double precision, Sc_Cc_temp, (ao_num, mo_num)]
!  implicit none
!  call get_AB_prod(ao_overlap,ao_num,ao_num,mo_coef,mo_num,Sc_Cc_temp)
!END_PROVIDER 


BEGIN_PROVIDER [ double precision, mo_coef_in_ao_sphe_basis, (ao_sphe_num,mo_num) ]
  implicit none
  BEGIN_DOC
  ! |MO|s coefficients in the spherical |AO| basis
  !
  ! mo_coef_in_ao_sphe_basis(i,j) = coefficient of the i-th spherical |AO| on the j-th |MO|
  !
  ! These coefficients are obtained as
  !
  ! :math:`(S^s)^{-1} (C^s)^T S^c C^c`
  !
  ! where
  !
  ! :math:`(S^s)^{-1}` is :c:data:`ao_cart_to_sphe_overlap_inv`,
  ! the overlap matrix in spherical |AO| basis
  ! :math:`   (C^s)^T` is :c:data:`ao_cart_to_sphe_coef_transp`,
  ! the transpose of the coefficients' matrix to switch from cartesian to spherical |AO|
  ! :math:`       S^c` is :c:data:`ao_overlap`,
  ! the overlap matrix in the cartesian |AO| basis
  ! :math:`       C^c` is :c:data:`mo_coef`,
  ! the matrix of coefficients to expres the |MO|s in cartesian |AO|s
  END_DOC
  double precision :: SsINV_CsT_temp(ao_sphe_num, ao_num)
  double precision :: Sc_Cc_temp(ao_num, mo_num)
  ! Matrix product between first two matrices
  call get_AB_prod(ao_cart_to_sphe_overlap_inv,         &
                   ao_sphe_num, ao_sphe_num,            &
                   ao_cart_to_sphe_coef_transp, ao_num, &
                   SsINV_CsT_temp)                          
  ! Matrix product between last two matrices
  call get_AB_prod(ao_overlap,ao_num,ao_num,mo_coef,mo_num,Sc_Cc_temp)
  ! Total matrix product
  call get_AB_prod(SsINV_CsT_temp, ao_sphe_num, ao_num, &
                   Sc_Cc_temp, mo_num,                  &
                   mo_coef_in_ao_sphe_basis)
END_PROVIDER


 BEGIN_PROVIDER [ double precision, mo_core_coef_sphe   , (ao_sphe_num, n_core_pseudo_orb) ]
&BEGIN_PROVIDER [ double precision, mo_val_coef_sphe, (ao_sphe_num, n_valence_pseudo_orb) ]
  implicit none
  BEGIN_DOC
  ! Core MOs' coefficients on |AO| basis set
  ! Valence MOs' coefficients on |AO| basis set
  !
  ! mo_core_coef(i,j) = coefficient of the i-th |AO| on the jth core |MO|
  ! mo_val_coef(i,j) = coefficient of the i-th |AO| on the jth valence |MO|
  END_DOC
  mo_core_coef_sphe = mo_coef_in_ao_sphe_basis(:,list_core_pseudo)
  mo_val_coef_sphe = mo_coef_in_ao_sphe_basis(:,list_valence_pseudo)  

  !! Normalization
  !integer :: col
  !double precision :: norm
  !do col=1,n_core_pseudo_orb
  !  norm = sqrt(dot_product(mo_core_coef_sphe(:,col),mo_core_coef_sphe(:,col)))
  !  mo_core_coef_sphe(:,col) = mo_core_coef_sphe(:,col)/norm
  !end do
  !do col=1,n_valence_pseudo_orb
  !  norm = sqrt(dot_product(mo_val_coef_sphe(:,col),mo_val_coef_sphe(:,col)))
  !  mo_val_coef_sphe(:,col) = mo_val_coef_sphe(:,col)/norm
  !end do
END_PROVIDER


 BEGIN_PROVIDER [ double precision, proj_mo_val_sphe, (ao_sphe_num, ao_sphe_num) ]
&BEGIN_PROVIDER [ double precision, proj_mo_val_sphe_transp, (ao_sphe_num, ao_sphe_num) ]
  BEGIN_DOC
  ! Projector on the valence MOs subspace in the spherical AO basis 
  !
  ! :math:`P^v_{ij} = \sum_a^{N_val} \tilde{C^s}_{ia}\tilde{C^s}*_{ja} = \sum_a^{N_val} \tilde{C^s}_{ia}\tilde{C^s}^{\dag}_{aj} = [\tilde{C^s}\tilde{C^s}^{\dag}]_{ij}`
  !
  ! where :math:`\tilde{C^s}` is the matrix of the MOs coefficients expressed in 
  ! the basis of spherical AOs, restricted to the valence orbitals
  END_DOC 
  integer :: i,j
  !
  proj_mo_val_sphe(:,:) = 0.d0
  call dgemm('N',                     & 
             'T',                     & 
             ao_sphe_num,             & 
             ao_sphe_num,             & 
             n_valence_pseudo_orb,    & 
             1.d0,                    & 
             mo_val_coef_sphe, & 
             ao_sphe_num,             & 
             mo_val_coef_sphe, & 
             ao_sphe_num,             & 
             1.d0,                    & 
             proj_mo_val_sphe, & 
             ao_sphe_num              & 
  )
  !
  ! Compute transpose
  do i = 1,ao_sphe_num
   do j = 1,ao_sphe_num
     proj_mo_val_sphe_transp(j,i) = proj_mo_val_sphe(i,j)
   enddo
  enddo
END_PROVIDER 


!! BEGIN_PROVIDER [double precision, proj_mo_val_sphe_old, (ao_sphe_num, ao_sphe_num) ]
!!&BEGIN_PROVIDER [double precision, proj_mo_val_sphe_old_transp, (ao_sphe_num, ao_sphe_num) ]
!!  BEGIN_DOC
!!  ! Projector P restricted to the pseudo-potential valence MOs ON THE SHPERICAL AO BASIS
!!  !
!!  ! P_{sph} = C^TP_{cart}C
!!  !
!!  ! where C is ao_cart_to_sphe_coef
!!  END_DOC
!!  implicit none
!!  integer :: i,j
!!  double precision, allocatable :: S(:,:)
!!  allocate (S(ao_sphe_num,ao_num))
!!  !
!!  proj_mo_val_sphe_old = 0.d0
!!  proj_mo_val_sphe_old_transp = 0.d0
!!  !
!!  call ao_cart_to_ao_sphe(proj_mo_val,size(proj_mo_val,1), & 
!!                          proj_mo_val_sphe_old,size(proj_mo_val_sphe_old,1))
!!  ! Compute transpose
!!  do i = 1,ao_sphe_num
!!    do j = 1,ao_sphe_num
!!      proj_mo_val_sphe_old_transp(j,i) = proj_mo_val_sphe_old(i,j)
!!    enddo
!!  enddo
!!END_PROVIDER 
