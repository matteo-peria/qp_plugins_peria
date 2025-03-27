! Test providers for aos_coreence.irp.f

BEGIN_PROVIDER [ double precision, ao_core_overlap_as_matprod, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Matrix of the nuclear-electron potential for core orbitals
  ! obtained as matrix-base change from the cartesian AOs and the core AOs
  END_DOC
  !
  !call ao_cart_to_ao_core(ao_overlap,                 &
  call ao_cart_to_ao_core_eigvec(ao_overlap,                 &
                          ao_num,                     &
                          ao_core_overlap_as_matprod, &
                          ao_num                      &
  )
END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_core_overlap_normed_as_matprod, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Matrix of the nuclear-electron potential for core orbitals
  ! obtained as matrix-base change from the cartesian AOs and the core AOs
  END_DOC
  !
  call ao_cart_to_ao_core_normed(ao_overlap,                &
                                 ao_num,                    &
                                 ao_core_overlap_normed_as_matprod, &
                                 ao_num                     &
  )
END_PROVIDER


