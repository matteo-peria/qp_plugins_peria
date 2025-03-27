! Spherical AOs used to define the valence-subspace

 BEGIN_PROVIDER [ double precision, ao_val_coef_sphe, (ao_sphe_num, ao_sphe_num) ]
&BEGIN_PROVIDER [ double precision, ao_val_coef_sphe_transp, (ao_sphe_num, ao_sphe_num) ]
  BEGIN_DOC
  ! AOs coefficients to express the pseudopot valence molecular orbitals
  !
  ! C^v_{ij} = [CC^{\dag}S]_{ij} =  [PS]_{ij}
  !
  ! where C is mo_core_coef,
  !       P is proj_mo_val_sphe in spherical coordinates
  !       S is ao_cart_to_sphe_overlap
  END_DOC 
  implicit none
  integer :: i_ao, j_ao, i, j
  !
  ao_val_coef_sphe(:,:) = 0.d0
  call dgemm('N','N',                              &
      ao_sphe_num, ao_sphe_num, ao_sphe_num, 1.d0, &
      proj_mo_val_sphe, ao_sphe_num,        &
      ao_cart_to_sphe_overlap, ao_sphe_num,        &
      0.d0, ao_val_coef_sphe, ao_sphe_num   &
  )
  !
  ! Compute transpose
  do i = 1,ao_sphe_num
    do j = 1,ao_sphe_num
      ao_val_coef_sphe_transp(j,i) = ao_val_coef_sphe(i,j)
    enddo
  enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_val_overlap_sphe, (ao_sphe_num, ao_sphe_num) ]
  implicit none
  BEGIN_DOC
  ! Spherical AOs overlap matrix in the valence subspace
  END_DOC 
  double precision, dimension(ao_sphe_num,ao_sphe_num) :: SSCCtS
  integer :: i
  !
  call get_AB_prod(                                     &
      ao_cart_to_sphe_overlap,ao_sphe_num,ao_sphe_num,  &
      ao_val_coef_sphe,ao_sphe_num,                     &
      SSCCtS                                            &
  )
  call get_AB_prod(                                    & 
      ao_val_coef_sphe_transp,ao_sphe_num,ao_sphe_num, &
      SSCCtS,ao_sphe_num,                              &
      ao_val_overlap_sphe                              &
  )
END_PROVIDER


 BEGIN_PROVIDER [ double precision, ao_val_overlap_sphe_eigval, (ao_sphe_num) ]
&BEGIN_PROVIDER [ double precision, ao_val_overlap_sphe_eigvec, (ao_sphe_num, ao_sphe_num) ]
  implicit none
  BEGIN_DOC
  ! AOs overlap matrix in the valence subspace expressed in spherical AOs
  END_DOC 
  call lapack_diagd(ao_val_overlap_sphe_eigval, ao_val_overlap_sphe_eigvec, &
                    ao_val_overlap_sphe, ao_sphe_num, ao_sphe_num           &
  )
 END_PROVIDER 


 BEGIN_PROVIDER [ double precision, ao_val_overlap_sphe_eigvec_coef, (ao_sphe_num, ao_sphe_num) ]
  implicit none
  BEGIN_DOC
  ! Coefficients of the eigenvectors of the overlap matrix (on the basis of 
  ! valence-projected AOs) developed on the spherical basis 
  END_DOC
  call get_AB_prod(ao_val_coef_sphe, ao_sphe_num, ao_sphe_num, &
                   ao_val_overlap_sphe_eigvec, ao_sphe_num,    &
                   ao_val_overlap_sphe_eigvec_coef             &
  )
 END_PROVIDER 


 BEGIN_PROVIDER [ double precision, ao_val_overlap_sphe_evec_overlap, (ao_sphe_num, ao_sphe_num) ]
  implicit none
  BEGIN_DOC
  ! Overlap between the eigenvectors of the overlap matrix on the spherical basis
  END_DOC
  call ao_sphe_to_ao_val_eigvec_sphe(ao_cart_to_sphe_overlap,          &
                                     ao_sphe_num,                      &
                                     ao_val_overlap_sphe_evec_overlap, &
                                     ao_sphe_num                       &
  )
 END_PROVIDER 


subroutine ao_sphe_to_ao_val_eigvec_sphe(A_ao_sphe,LDA_ao_sphe,A_ao_val_eigvec_sphe,LDA_ao_val_eigvec_sphe)
  implicit none
  BEGIN_DOC
  ! Transform matrix A from the |AO| spherical basis to the valence |AO| basis 
  ! obtained as the eigenvectors of the overlap matrix of the 
  ! (spherical) valence basis
  !
  ! :math:`A^{eig,sphe}=B^T.A^{sphe}.B`
  !
  ! where :math:`B` is :c:data:`ao_val_coef`, 
  ! the matrix of coefficients from the cartesian AO basis to valence AO one,
  ! and :math:`B` is :c:data:`ao_val_coef_transp`, its transpose.
  END_DOC
  integer, intent(in)           :: LDA_ao_sphe,LDA_ao_val_eigvec_sphe
  double precision, intent(in)  :: A_ao_sphe(LDA_ao_sphe,ao_sphe_num)
  double precision, intent(out) :: A_ao_val_eigvec_sphe(LDA_ao_val_eigvec_sphe,ao_sphe_num)
  double precision, allocatable :: T(:,:)
  !
  allocate (T(ao_sphe_num,ao_sphe_num) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T

  call dgemm('N','N', ao_sphe_num, ao_sphe_num, ao_sphe_num,       &
      1.d0, A_ao_sphe, LDA_ao_sphe,                                &
      ao_val_overlap_sphe_eigvec_coef, size(ao_val_coef,1), &
      0.d0, T, size(T,1)                                           &
  )
  ! Notice that for the following dgemm we could have used
  ! ao_val_coef_transp, but instead we transposed with the 'T' argument
  call dgemm('T','N', ao_sphe_num, ao_sphe_num, ao_sphe_num,             &
      1.d0, ao_val_overlap_sphe_eigvec_coef, size(ao_val_coef,1), &
      T, ao_sphe_num,                                                    &
      0.d0, A_ao_val_eigvec_sphe, size(A_ao_val_eigvec_sphe,1)           &
  )
  deallocate(T)
end subroutine


 BEGIN_PROVIDER [ double precision, ao_val_sphe_coef_normed, (ao_sphe_num, ao_sphe_num) ]
&BEGIN_PROVIDER [ double precision, ao_val_sphe_coef_normed_transp, (ao_sphe_num, ao_sphe_num) ]
  implicit none
  integer :: i,j
  !
  ao_val_sphe_coef_normed(:,:) = 0.d0
  do i = 1,ao_sphe_num
    if(ao_val_overlap_sphe_eigval(i).gt.1.d-6)then
      do j = 1,ao_sphe_num
        ao_val_sphe_coef_normed(j,i) =  ao_val_overlap_sphe_eigvec_coef(j,i) / dsqrt(ao_val_overlap_sphe_eigval(i))
      enddo
    end if
  enddo
  !
  ! Compute transpose
  ao_val_sphe_coef_normed_transp(:,:) = 0.d0
  do i = 1,ao_sphe_num
    do j = 1,ao_sphe_num
      ao_val_sphe_coef_normed_transp(j,i) = ao_val_sphe_coef_normed(i,j)
    enddo
  enddo
END_PROVIDER


subroutine ao_sphe_to_ao_val_sphe_normed(A_ao_sphe, LDA_ao_sphe, A_ao_val_normed, LDA_ao_val_normed)
  implicit none
  integer, intent(in)           :: LDA_ao_sphe, LDA_ao_val_normed
  double precision, intent(in)  :: A_ao_sphe(LDA_ao_sphe, ao_sphe_num)
  double precision, intent(out) :: A_ao_val_normed(LDA_ao_val_normed, ao_sphe_num)
  double precision, allocatable :: T(:,:)
  !
  allocate (T(ao_sphe_num, ao_sphe_num))
  call get_AB_prod(A_ao_sphe,          &
                   ao_sphe_num,        &
                   ao_sphe_num,        &
                   ao_val_coef_normed, &
                   ao_sphe_num,        &
                   T                   &
  )
  call get_AB_prod(ao_val_coef_normed_transp, &
                   ao_sphe_num,               &
                   ao_sphe_num,               &
                   T, ao_sphe_num,            &
                   A_ao_val_normed            &
  )
end subroutine 


 BEGIN_PROVIDER [ double precision, ao_val_sphe_normed_overlap, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! AOs overlap matrix in the valence subspace using normed eigenvectors
  END_DOC 
  !
  call ao_sphe_to_ao_val_sphe_normed(ao_cart_to_sphe_overlap, &
                                     ao_sphe_num,             &
                                     ao_val_sphe_normed_overlap,   &
                                     ao_sphe_num              &
  )
 END_PROVIDER
