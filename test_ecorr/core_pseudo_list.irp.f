
use bitmasks

BEGIN_PROVIDER [ integer, n_core_pseudo_orb]
 implicit none
 integer :: i
 BEGIN_DOC
 ! number of frozen core orbitals considered in the core pseudo potential
 END_DOC
 n_core_pseudo_orb = n_core_orb
 if(.False.)then
  n_core_pseudo_orb = 0
  do i = 1, mo_num
    if(mo_class(i) == 'core_pseudo')then
      n_core_pseudo_orb += 1
    endif
  enddo
  
 endif
  call write_int(6,n_core_pseudo_orb, 'Number of core_pseudo     MOs')
END_PROVIDER 

BEGIN_PROVIDER [integer, dim_list_core_pseudo_orb]
  implicit none
  BEGIN_DOC
  ! dimensions for the allocation of list_core_pseudo.
  ! it is at least 1
  END_DOC
   dim_list_core_pseudo_orb = max(n_core_pseudo_orb,1)
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), core_pseudo_bitmask , (N_int,2) ]
   implicit none
   BEGIN_DOC
   ! Bitmask identifying the core_pseudo MOs 
   END_DOC
   core_pseudo_bitmask  = 0_bit_kind
   if(n_core_pseudo_orb > 0)then
     call list_to_bitstring( core_pseudo_bitmask(1,1), list_core_pseudo, n_core_pseudo_orb, N_int)
     call list_to_bitstring( core_pseudo_bitmask(1,2), list_core_pseudo, n_core_pseudo_orb, N_int)
   endif
 END_PROVIDER

 BEGIN_PROVIDER [ integer, list_core_pseudo        , (dim_list_core_pseudo_orb) ]
&BEGIN_PROVIDER [ integer, list_core_pseudo_reverse, (mo_num) ]
   implicit none
   BEGIN_DOC
   ! List of MO indices which are in the core_pseudo.
   END_DOC
   integer                        :: i, n
   list_core_pseudo = 0
   list_core_pseudo_reverse = 0

   list_core_pseudo = list_core
   list_core_pseudo_reverse = list_core_reverse 
   if(.False.)then
    n=0
    do i = 1, mo_num
      if(mo_class(i) == 'core_pseudo')then
        n += 1
        list_core_pseudo(n) = i
        list_core_pseudo_reverse(i) = n
      endif
    enddo
    print *,  'core_pseudo MOs:'
    print *,  list_core_pseudo(1:n_core_pseudo_orb)
   endif
   
END_PROVIDER
 
