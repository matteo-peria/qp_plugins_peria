! Test providers for aos_valence_sphe.irp.f

BEGIN_PROVIDER [ double precision, ao_val_overlap_sphe_as_matprod, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Matrix of the nuclear-electron potential for valence orbitals
  ! obtained as matrix-base change from the spherical AOs and the valence AOs
  END_DOC
  !
  call ao_sphe_to_ao_val_eigvec_sphe(ao_cart_to_sphe_overlap,                &
                                     ao_num,                         &
                                     ao_val_overlap_sphe_as_matprod, &
                                     ao_num                          &
  )
END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_val_overlap_sphe_normed_as_matprod, (ao_sphe_num,ao_sphe_num)]
  implicit none
  BEGIN_DOC
  ! Matrix of the nuclear-electron potential for valence orbitals
  ! obtained as matrix-base change from the spherical AOs and the valence AOs
  END_DOC
  !
  call ao_sphe_to_ao_val_sphe_normed(ao_cart_to_sphe_overlap,                  &
                                ao_sphe_num,                           &
                                ao_val_overlap_sphe_normed_as_matprod, &
                                ao_sphe_num                            &
  )
END_PROVIDER


