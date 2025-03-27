!!! Test providers for aos_valence.irp.f
!!
!!BEGIN_PROVIDER [ double precision, ao_val_overlap_as_matprod, (ao_num,ao_num)]
!!  implicit none
!!  BEGIN_DOC
!!  ! Matrix of the nuclear-electron potential for valence orbitals
!!  ! obtained as matrix-base change from the cartesian AOs and the valence AOs
!!  END_DOC
!!  !
!!  call ao_cart_to_ao_val(ao_overlap,                &
!!                         ao_num,                    &
!!                         ao_val_overlap_as_matprod, &
!!                         ao_num                     &
!!  )
!!END_PROVIDER
!!
!!
!!BEGIN_PROVIDER [ double precision, ao_val_overlap_normed_as_matprod, (ao_num,ao_num)]
!!  implicit none
!!  BEGIN_DOC
!!  ! Matrix of the nuclear-electron potential for valence orbitals
!!  ! obtained as matrix-base change from the cartesian AOs and the valence AOs
!!  END_DOC
!!  !
!!  call ao_cart_to_ao_val_normed(ao_overlap,                &
!!                                ao_num,                    &
!!                                ao_val_overlap_normed_as_matprod, &
!!                                ao_num                     &
!!  )
!!END_PROVIDER
