 BEGIN_PROVIDER [ double precision, mo_core_coef   , (ao_num, n_core_pseudo_orb) ]
&BEGIN_PROVIDER [ double precision, mo_val_coef, (ao_num, n_valence_pseudo_orb) ]
  implicit none
  BEGIN_DOC
  ! Core MOs' coefficients on |AO| basis set
  ! Valence MOs' coefficients on |AO| basis set
  !
  ! mo_core_coef(i,j) = coefficient of the i-th |AO| on the jth core |MO|
  ! mo_val_coef(i,j) = coefficient of the i-th |AO| on the jth valence |MO|
  END_DOC
  integer :: i

  !mo_core_coef = mo_coef(:,list_core_pseudo)
  !mo_val_coef = mo_coef(:,list_valence_pseudo)  

  do i=1, n_core_pseudo_orb
    mo_core_coef(:,i) = mo_coef(:,list_core_pseudo(i))  
  end do


  do i=1, n_valence_pseudo_orb
    mo_val_coef(:,i) = mo_coef(:,list_valence_pseudo(i))  
  end do

  ! Normalization
  integer :: col
  double precision :: norm
  do col=1,n_core_pseudo_orb
    norm = sqrt(dot_product(mo_core_coef(:,col),mo_core_coef(:,col)))
    mo_core_coef(:,col) = mo_core_coef(:,col)/norm
  end do
  do col=1,n_valence_pseudo_orb
    norm = sqrt(dot_product(mo_val_coef(:,col),mo_val_coef(:,col)))
    mo_val_coef(:,col) = mo_val_coef(:,col)/norm
  end do
 END_PROVIDER


 BEGIN_PROVIDER [ double precision, proj_mo_val, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, proj_mo_val_transp, (ao_num, ao_num) ]
  BEGIN_DOC
  ! Projector on the valence MOs subspace in the cartesian AO basis 
  !
  ! :math:`P^v_{ij} = \sum_a^{N_val} \tilde{C^c}_{ia}\tilde{C^c}*_{ja} = \sum_a^{N_val} \tilde{C^c}_{ia}\tilde{C^c}^{\dag}_{aj} = [\tilde{C^c}\tilde{C^c}^{\dag}]_{ij}`
  !
  ! where :math:`\tilde{C^s}` is the matrix of the MOs coefficients expressed in 
  ! the basis of cartesian AOs, restricted to the valence orbitals
  END_DOC 
  integer :: i,j
  !
  proj_mo_val(:,:) = 0.d0
  call dgemm('N',                  & 
             'T',                  & 
             ao_num,               & 
             ao_num,               & 
             n_valence_pseudo_orb, & 
             1.d0,                 & 
             mo_val_coef,   & 
             ao_num,               & 
             mo_val_coef,   & 
             ao_num,               & 
             1.d0,                 & 
             proj_mo_val,   & 
             ao_num                & 
  )
  !
  ! Compute transpose
  do i = 1,ao_num
   do j = 1,ao_num
     proj_mo_val_transp(j,i) = proj_mo_val(i,j)
   enddo
  enddo
END_PROVIDER 


 BEGIN_PROVIDER [ double precision, proj_mo_core, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, proj_mo_core_transp, (ao_num, ao_num) ]
  BEGIN_DOC
  ! Projector on the core MOs subspace in the cartesian AO basis 
  !
  ! :math:`P^v_{ij} = \sum_a^{N_val} \tilde{C^c}_{ia}\tilde{C^c}*_{ja} = \sum_a^{N_val} \tilde{C^c}_{ia}\tilde{C^c}^{\dag}_{aj} = [\tilde{C^c}\tilde{C^c}^{\dag}]_{ij}`
  !
  ! where :math:`\tilde{C^s}` is the matrix of the MOs coefficients expressed in 
  ! the basis of cartesian AOs, restricted to the core orbitals
  END_DOC 
  integer :: i,j
  !
  proj_mo_core(:,:) = 0.d0
  call dgemm('N',                  & 
             'T',                  & 
             ao_num,               & 
             ao_num,               & 
             n_core_pseudo_orb, & 
             1.d0,                 & 
             mo_core_coef,   & 
             ao_num,               & 
             mo_core_coef,   & 
             ao_num,               & 
             1.d0,                 & 
             proj_mo_core,   & 
             ao_num                & 
  )
  !
  ! Compute transpose
  do i = 1,ao_num
   do j = 1,ao_num
     proj_mo_core_transp(j,i) = proj_mo_core(i,j)
   enddo
  enddo
END_PROVIDER 
