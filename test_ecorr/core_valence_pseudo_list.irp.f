
use bitmasks

 BEGIN_PROVIDER [ integer, n_core_pseudo_orb]
&BEGIN_PROVIDER [ integer, n_valence_pseudo_orb]
  implicit none
  integer :: i
  BEGIN_DOC
  ! Number of frozen core orbitals considered in the core pseudo potential
  END_DOC

  print*, "Computing separation between core and valence when using pseudop"
  print*, "n_core_orb", n_core_orb
  print*, "mo_num", mo_num
  n_core_pseudo_orb = n_core_orb
  n_valence_pseudo_orb = mo_num - n_core_orb !n_core_pseudo_orb
  if(.False.)then
    n_core_pseudo_orb = 0
    do i = 1, mo_num
      if(mo_class(i) == 'core_pseudo')then
        n_core_pseudo_orb += 1
      endif
    enddo
  endif
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
  print*, "IM IN dim_list_valence_pseudo_orb"
  dim_list_core_pseudo_orb = max(n_core_pseudo_orb,1)
  dim_list_valence_pseudo_orb = max(n_valence_pseudo_orb,1)
END_PROVIDER


!BEGIN_PROVIDER [integer, dim_list_valence_pseudo_orb]
!  implicit none
!  BEGIN_DOC
!  ! dimensions for the allocation of list_valence_pseudo.
!  ! it is at least 1
!  END_DOC
!   dim_list_valence_pseudo_orb = max(n_valence_pseudo_orb,1)
!END_PROVIDER


 BEGIN_PROVIDER [ integer(bit_kind), core_pseudo_bitmask , (N_int,2) ]
&BEGIN_PROVIDER [ integer(bit_kind), valence_pseudo_bitmask , (N_int,2) ]
  implicit none
  BEGIN_DOC
  ! Bitmask identifying the core_pseudo MOs and the valence_pseudo MOs 
  END_DOC
  core_pseudo_bitmask  = 0_bit_kind
  valence_pseudo_bitmask  = 0_bit_kind
  if(n_core_pseudo_orb > 0)then
    call list_to_bitstring( core_pseudo_bitmask(1,1), list_core_pseudo, n_core_pseudo_orb, N_int)
    call list_to_bitstring( core_pseudo_bitmask(1,2), list_core_pseudo, n_core_pseudo_orb, N_int)
    valence_pseudo_bitmask(1,1) = not(core_pseudo_bitmask(1,1))
    valence_pseudo_bitmask(1,2) = not(core_pseudo_bitmask(1,2))
  endif
END_PROVIDER


! BEGIN_PROVIDER [ integer, list_core_pseudo        , (dim_list_core_pseudo_orb) ]
!&BEGIN_PROVIDER [ integer, list_core_pseudo_reverse, (mo_num) ]
!   implicit none
!   BEGIN_DOC
!   ! List of MO indices which are in the core_pseudo.
!   END_DOC
!   integer                        :: i, n
!   list_core_pseudo = 0
!   list_core_pseudo_reverse = 0
!
!   list_core_pseudo = list_core
!   list_core_pseudo_reverse = list_core_reverse 
!   if(.False.)then
!    n=0
!    do i = 1, mo_num
!      if(mo_class(i) == 'core_pseudo')then
!        n += 1
!        list_core_pseudo(n) = i
!        list_core_pseudo_reverse(i) = n
!      endif
!    enddo
!    print *,  'core_pseudo MOs:'
!    print *,  list_core_pseudo(1:n_core_pseudo_orb)
!   endif   
!END_PROVIDER


 BEGIN_PROVIDER [ integer, list_core_pseudo           , (dim_list_core_pseudo_orb) ]
&BEGIN_PROVIDER [ integer, list_valence_pseudo        , (dim_list_valence_pseudo_orb) ]
&BEGIN_PROVIDER [ integer, list_core_pseudo_reverse   , (mo_num) ]
&BEGIN_PROVIDER [ integer, list_valence_pseudo_reverse, (mo_num) ]
   implicit none
   BEGIN_DOC
   ! List of MO indices which are in the core_pseudo and valence_pseudo.
   END_DOC
   integer :: i_mo, i_core_pseudo, i_valence_pseudo
   !
   print*, 'IM COMPUTING VARIOUS LISTS'
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
     print*, "MO_NUM", mo_num
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
