
 BEGIN_PROVIDER [ double precision, ao_val_coef_normed, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_val_coef_normed_transp, (ao_num, ao_num) ]
  implicit none
  integer :: i,j
  BEGIN_DOC
 ! Coefficients, on the AO basis, of the valence projected AO basis 
 !
 ! these functions are orthonormalized, and span a vector space of dimension mo_num - n_core_orb
  END_DOC
  ! 
  ao_val_coef_normed(:,:) = 0.d0
  do i = 1,ao_num
    if(ao_val_overlap_eigval(i).gt.1.d-6)then
      do j = 1,ao_num
        ao_val_coef_normed(j,i) =  ao_val_overlap_eigvec_coef(j,i) / dsqrt(ao_val_overlap_eigval(i))
      enddo
    end if
  enddo
  !
  ! Compute transpose
  ao_val_coef_normed_transp(:,:) = 0.d0
  do i = 1,ao_num
    do j = 1,ao_num
      ao_val_coef_normed_transp(j,i) = ao_val_coef_normed(i,j)
    enddo
  enddo
END_PROVIDER

subroutine ao_cart_to_ao_val_normed(A_ao_cart, LDA_ao_cart, A_ao_val_normed, LDA_ao_val_normed)
  implicit none
 BEGIN_DOC
 ! Transfrom a matrix from the usual AO basis to the valence projected, orthonormal, AO basis 
 !
 ! use the coefficients ao_val_coef_normed
 END_DOC
  integer, intent(in)           :: LDA_ao_cart, LDA_ao_val_normed
  double precision, intent(in)  :: A_ao_cart(LDA_ao_cart, ao_num)
  double precision, intent(out) :: A_ao_val_normed(LDA_ao_val_normed,ao_num)
  double precision, allocatable :: T(:,:)
  !
  allocate (T(ao_num,ao_num) )
  call get_AB_prod(A_ao_cart, ao_num, ao_num, ao_val_coef_normed, ao_num, T)
  call get_AB_prod(ao_val_coef_normed_transp, ao_num, ao_num, T, ao_num, A_ao_val_normed)
end subroutine 


 BEGIN_PROVIDER [ double precision, ao_val_normed_overlap, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! AOs overlap matrix in the valence subspace using normed eigenvectors
  ! 
  ! should be the identity in the valence space and zero in the core space 
  END_DOC 
  !
  call ao_cart_to_ao_val_normed(ao_overlap,            &
                                ao_num,                &
                                ao_val_normed_overlap, &
                                ao_num                 &
  )
 END_PROVIDER

!use tester

! Cartesian AOs used to define the valence-subspace

 BEGIN_PROVIDER [ double precision, ao_val_coef, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_val_coef_transp, (ao_num, ao_num) ]
  BEGIN_DOC
  ! AOs coefficients to express the valence projected AOS
  !
  ! P^val|ket{AO_i}> = \sum_{j} ao_val_coef(j,i) \ket{AO_j}
  !
  ! ao_val_coef(j,i) =  [PS]_{ji}
  !
  ! where P is proj_mo_val
  !       S is ao_overlap
  END_DOC 
  implicit none
  integer :: i_ao, j_ao, i, j
  !
  ao_val_coef(:,:) = 0.d0
  call dgemm('N','N',ao_num,ao_num,ao_num,1.d0, &
      proj_mo_val,ao_num,                       &
      ao_overlap,ao_num,                        &
      0.d0,ao_val_coef,ao_num                   &
  )
  !
  ! Compute transpose
  do i = 1,ao_num
    do j = 1,ao_num
      ao_val_coef_transp(j,i) = ao_val_coef(i,j)
    enddo
  enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_val_overlap, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! Overlap between the valence projected AO basis: if P^val|ket{AO_i}> = \ket{tilde{AO}_i} 
  ! 
  ! ao_val_overlap(j,i) = \bra{tilde{AO}_j} \ket{tilde{AO}_i} = ao_val_coef^T S ao_val_coef
  END_DOC 
  double precision, dimension(ao_num,ao_num) :: SCCtS
  integer :: i,j
  call get_AB_prod(ao_overlap,ao_num,ao_num,ao_val_coef,ao_num,SCCtS)
  call get_AB_prod(ao_val_coef_transp,ao_num,ao_num,SCCtS,ao_num,ao_val_overlap)
  !
  ! Avoid orbital mixing by setting epsilons to hard zeros
  do i = 1, ao_num
    do j = 1, ao_num
      ! remember why we did it ?
      if(dabs(ao_val_overlap(j,i)).lt.1.d-10)then
        ao_val_overlap(j,i) = 0.d0
      endif
    enddo
  enddo
END_PROVIDER


 BEGIN_PROVIDER [ double precision, ao_val_overlap_eigval, (ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_val_overlap_eigvec, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! Eigenvectors and eigenvalues of the Overlap matrix between the valence projected AOs
  ! 
  ! In practice, we use an SVD (more stable) which turns out to be the same as Overlap matrix are positive definite
  END_DOC 
  double precision :: Vt(ao_num,ao_num)
  call svd(ao_val_overlap,ao_num,ao_val_overlap_eigvec,ao_num,ao_val_overlap_eigval,Vt,& 
           ao_num,ao_num,ao_num)
 END_PROVIDER 


 BEGIN_PROVIDER [ double precision, ao_val_overlap_eigvec_coef, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! Coefficients, on the AO basis, of the eigenvectors of the overlap matrix of valence projected AOs 
  ! 
  END_DOC
  call get_AB_prod(ao_val_coef, ao_num, ao_num,   &
                   ao_val_overlap_eigvec, ao_num, &
                   ao_val_overlap_eigvec_coef     &
  )
 END_PROVIDER 


 BEGIN_PROVIDER [ double precision, ao_val_overlap_diag, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! Overlap between the eigenvectors of the overlap matrix on the cartesian basis
  END_DOC
  call ao_cart_to_ao_val_eigvec(ao_overlap,                 &
                               ao_num,                      &
                               ao_val_overlap_diag, &
                               ao_num                       &
  )
 END_PROVIDER 

!---





subroutine ao_cart_to_ao_val_eigvec(A_ao_cart,LDA_ao_cart,A_ao_val,LDA_ao_val)
  implicit none
  BEGIN_DOC
  ! Transform matrix A from the |AO| cartesian basis to the projected valence |AO| basis 
  !
  ! :math:`B^T.A^c.B`
  !
  ! where :math:`B` is :c:data:`ao_val_coef`, 
  ! the matrix of coefficients from the cartesian AO basis to valence AO one,
  ! and :math:`B` is :c:data:`ao_val_coef_transp`, its transpose.
  END_DOC
  integer, intent(in)            :: LDA_ao_cart,LDA_ao_val
  double precision, intent(in)   :: A_ao_cart(LDA_ao_cart,ao_num)
  double precision, intent(out)  :: A_ao_val(LDA_ao_val,ao_num)
  !
  double precision :: T(ao_num,ao_num)
  !
  !allocate (T(ao_num,ao_num) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
  !
  call dgemm('N','N', ao_num, ao_num, ao_num,  &
      1.d0, A_ao_cart, LDA_ao_cart,            &
      ao_val_overlap_eigvec_coef, size(ao_val_coef,1), &
      0.d0, T, size(T,1)                       &
  )
  ! Notice that for the following dgemm we could have used
  ! ao_val_coef_transp, but instead we transposed with the 'T' argument
  !call dgemm('N','N', ao_num, ao_num, ao_num,        &
  call dgemm('T','N', ao_num, ao_num, ao_num,        &
      1.d0, ao_val_overlap_eigvec_coef, size(ao_val_coef,1), &
      T, ao_num,                                     &
      0.d0, A_ao_val, size(A_ao_val,1)               &
  )
  !deallocate(T)
end subroutine
