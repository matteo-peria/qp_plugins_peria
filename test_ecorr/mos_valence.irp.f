 BEGIN_PROVIDER [ double precision, mo_core_coef_notnorm,        (ao_num, n_core_pseudo_orb) ]
&BEGIN_PROVIDER [ double precision, mo_val_coef_notnorm,         (ao_num, n_valence_pseudo_orb) ]
&BEGIN_PROVIDER [ double precision, mo_core_coef_notnorm_transp, (n_core_pseudo_orb, ao_num) ]
&BEGIN_PROVIDER [ double precision, mo_val_coef_notnorm_transp,  (n_valence_pseudo_orb, ao_num) ]

  implicit none
  BEGIN_DOC
  ! Core MOs' coefficients on |AO| basis set. NOT NORMALIZED
  ! Valence MOs' coefficients on |AO| basis set. NOT NORMALIZED
  !
  ! mo_core_nn_coef_notnorm(i,j) = coefficient of the i-th |AO| on the j-th core |MO|
  ! mo_val_nn_coef_notnorm(i,j) = coefficient of the i-th |AO| on the j-th valence |MO|
  END_DOC

  ! Array slicing solution
  mo_core_coef_notnorm = mo_coef(:,list_core_pseudo)
  mo_val_coef_notnorm  = mo_coef(:,list_valence_pseudo)  

  ! Transpose
  mo_core_coef_notnorm_transp = transpose(mo_core_coef_notnorm)
  mo_val_coef_notnorm_transp = transpose(mo_val_coef_notnorm)
 END_PROVIDER


 BEGIN_PROVIDER [ double precision, mo_core_coef,        (ao_num, n_core_pseudo_orb) ]
&BEGIN_PROVIDER [ double precision, mo_val_coef,         (ao_num, n_valence_pseudo_orb) ]
&BEGIN_PROVIDER [ double precision, mo_core_coef_transp, (n_core_pseudo_orb, ao_num) ]
&BEGIN_PROVIDER [ double precision, mo_val_coef_transp,  (n_valence_pseudo_orb, ao_num) ]

  implicit none
  BEGIN_DOC
  ! Core MOs' coefficients on |AO| basis set
  ! Valence MOs' coefficients on |AO| basis set
  !
  ! mo_core_coef(i,j) = coefficient of the i-th |AO| on the jth core |MO|
  ! mo_val_coef(i,j) = coefficient of the i-th |AO| on the jth valence |MO|
  END_DOC

  !! ARRAY SLICING SOLUTION
  !mo_core_coef = mo_coef(:,list_core_pseudo)
  !mo_val_coef = mo_coef(:,list_valence_pseudo)  
  !! Normalization
  !integer :: row
  !double precision :: norm_core_mos(n_core_pseudo_orb)
  !double precision :: norm_valence_mos(n_valence_pseudo_orb)
  !norm_core_mos(:)    = norm2(mo_core_coef_notnorm(:,:), dim=1)
  !norm_valence_mos(:) = norm2(mo_val_coef_notnorm(:,:), dim=1)
  !do row = 1, ao_num 
  !  mo_core_coef(row,:) = mo_core_coef_notnorm(row,:) / norm_core_mos(:)    
  !  mo_val_coef(row,:)  = mo_val_coef_notnorm(row,:)  / norm_valence_mos(:) 
  !end do

  ! LOOP SOLUTION
  integer :: i
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

  ! Transpose
  mo_core_coef_transp = transpose(mo_core_coef)
  mo_val_coef_transp = transpose(mo_val_coef)
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
