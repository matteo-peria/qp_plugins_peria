
use bitmasks

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



 BEGIN_PROVIDER [ integer, n_core_pseudo_orb]
&BEGIN_PROVIDER [ integer, n_valence_pseudo_orb]
  implicit none
  BEGIN_DOC
  ! Number of frozen core orbitals considered in the core pseudo potential
  END_DOC
  integer :: i
  n_core_pseudo_orb = n_core_orb
  n_valence_pseudo_orb = mo_num - n_core_orb !n_core_pseudo_orb
!  if(.False.)then
!    n_core_pseudo_orb = 0
!    do i = 1, mo_num
!      if(mo_class(i) == 'core_pseudo')then
!        n_core_pseudo_orb += 1
!      endif
!    enddo
!  endif
  call write_int(6,n_core_pseudo_orb, 'Number of core_pseudo     MOs')
  call write_int(6,n_valence_pseudo_orb, 'Number of valence_pseudo  MOs')
END_PROVIDER 


 BEGIN_PROVIDER [integer, dim_list_core_pseudo_orb]
&BEGIN_PROVIDER [integer, dim_list_valence_pseudo_orb]
  implicit none
  BEGIN_DOC
  ! dimensions for the allocation of list_core_pseudo.
  ! it is at least 1
  END_DOC
  dim_list_core_pseudo_orb = max(n_core_pseudo_orb,1)
  dim_list_valence_pseudo_orb = max(n_valence_pseudo_orb,1)
END_PROVIDER


 BEGIN_PROVIDER [ integer(bit_kind), core_pseudo_bitmask , (N_int,2) ]
&BEGIN_PROVIDER [ integer(bit_kind), valence_pseudo_bitmask , (N_int,2) ]
  implicit none
  BEGIN_DOC
  ! Bitmask identifying the core_pseudo MOs and the valence_pseudo MOs 
  END_DOC
  core_pseudo_bitmask(:,:)  = 0_bit_kind
  valence_pseudo_bitmask(:,:)  = 0_bit_kind
  if(n_core_pseudo_orb > 0)then
    call list_to_bitstring( core_pseudo_bitmask(1,1), list_core_pseudo, n_core_pseudo_orb, N_int)
    call list_to_bitstring( core_pseudo_bitmask(1,2), list_core_pseudo, n_core_pseudo_orb, N_int)
    valence_pseudo_bitmask(1,1) = not(core_pseudo_bitmask(1,1))
    valence_pseudo_bitmask(1,2) = not(core_pseudo_bitmask(1,2))
  endif
END_PROVIDER


 BEGIN_PROVIDER [ integer, list_core_pseudo           , (dim_list_core_pseudo_orb) ]
&BEGIN_PROVIDER [ integer, list_valence_pseudo        , (dim_list_valence_pseudo_orb) ]
&BEGIN_PROVIDER [ integer, list_core_pseudo_reverse   , (mo_num) ]
&BEGIN_PROVIDER [ integer, list_valence_pseudo_reverse, (mo_num) ]
   implicit none
   BEGIN_DOC
   ! List of MO indices which are in the core_pseudo and valence_pseudo.
   END_DOC
   integer :: i_mo, i_core_pseudo, i_valence_pseudo
   
   list_core_pseudo = 0
   list_core_pseudo_reverse = 0
   list_valence_pseudo = 0
   list_valence_pseudo_reverse = 0

   ! The following two lines are a provisory short cut 
   ! to define the core of the pseudo-potential
   !list_core_pseudo = list_core
   !list_core_pseudo_reverse = list_core_reverse 

   ! list_core and its reverse are built using the following lines
   ! We recycle them in order to define also the list_valence_pseudo at once 
   !if(.False.)then
   if(.true.)then
     i_core_pseudo=0
     i_valence_pseudo=0
     do i_mo = 1, mo_num
       !if(mo_class(i) == 'core_pseudo')then ! In future version there will be an appropriate label for the core of PP
       if(mo_class(i_mo) == 'Core')then
         i_core_pseudo += 1
         list_core_pseudo(i_core_pseudo) = i_mo
         list_core_pseudo_reverse(i_mo) = i_core_pseudo
       else
         i_valence_pseudo += 1
         list_valence_pseudo(i_valence_pseudo) = i_mo
         list_valence_pseudo_reverse(i_mo) = i_valence_pseudo
       endif
     enddo
     print *,  'core_pseudo MOs:'
     print *,  list_core_pseudo(1:n_core_pseudo_orb)
     print *,  'valence_pseudo MOs:'
     print *,  list_valence_pseudo(1:n_valence_pseudo_orb)
   endif
END_PROVIDER
