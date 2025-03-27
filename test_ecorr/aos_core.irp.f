! Cartesian AOs used to define the coreence-subspace

 BEGIN_PROVIDER [ double precision, ao_core_coef, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_core_coef_transp, (ao_num, ao_num) ]
  BEGIN_DOC
  ! AOs coefficients to express the pseudopot core molecular orbitals
  !
  ! C^v_{ij} = [CC^{\dag}S]_{ij} =  [PS]_{ij}
  !
  ! where C is mo_core_coef,
  !       P is proj_mo_core
  !       S is ao_overlap
  END_DOC 
  implicit none
  integer :: i_ao, j_ao, i, j
  !
  ao_core_coef(:,:) = 0.d0
  call dgemm('N','N',ao_num,ao_num,ao_num,1.d0, &
      proj_mo_core,ao_num,                   &
      ao_overlap,ao_num,                        &
      0.d0,ao_core_coef,ao_num                  &
  )
  !
  ! Compute transpose
  do i = 1,ao_num
    do j = 1,ao_num
      ao_core_coef_transp(j,i) = ao_core_coef(i,j)
    enddo
  enddo
END_PROVIDER


!! BEGIN_PROVIDER [ double precision, ao_core_coef_normalized, (ao_num, ao_num) ]
!!&BEGIN_PROVIDER [ double precision, ao_core_coef_normalized_transp, (ao_num, ao_num) ]
!!  ! FAILED ATTEMPT?
!!  implicit none
!!  integer :: i 
!!  double precision :: norm
!!  ao_core_coef_normalized = ao_core_coef
!!  ao_core_coef_normalized_transp = ao_core_coef_transp
!!  do i=1,ao_num
!!    !norm = sqrt(dot_product(ao_core_coef(:,i),ao_core_coef(:,i)))
!!    norm = sqrt(dot_product(ao_core_coef(i,:),ao_core_coef(i,:)))
!!    !ao_core_coef_normalized(i,:) = ao_core_coef_normalized(i,:)/norm
!!    !ao_core_coef_normalized_transp(:,i) =  ao_core_coef_normalized_transp(:,i)/norm
!!  end do
!!END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_core_overlap, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! AOs overlap matrix in the coreence subspace
  END_DOC 
  double precision, dimension(ao_num, ao_num) :: SCCtS
  integer :: i
  !
  ! NOTICE THAT WE CAN USE ao_cart_to_ao_core transformation subroutine
  call get_AB_prod(ao_overlap,ao_num,ao_num,ao_core_coef,ao_num,SCCtS)
  call get_AB_prod(ao_core_coef_transp,ao_num,ao_num,SCCtS,ao_num,ao_core_overlap)
END_PROVIDER


 BEGIN_PROVIDER [ double precision, ao_core_overlap_eigval, (ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_core_overlap_eigvec, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! AOs overlap matrix in the coreence subspace expressed in cartesian AOs
  END_DOC 
  call lapack_diagd(ao_core_overlap_eigval, ao_core_overlap_eigvec, &
                    ao_core_overlap, ao_num, ao_num                 &
  )
END_PROVIDER 


 BEGIN_PROVIDER [ double precision, ao_core_overlap_eigvec_coef, (ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! Coefficients of the eigenvectors of the overlap matrix (on the basis of 
 ! core-projected AOs) developed on the cartesian basis 
 END_DOC
 call get_AB_prod(ao_core_coef, ao_num, ao_num,   &
                  ao_core_overlap_eigvec, ao_num, &
                  ao_core_overlap_eigvec_coef     &
 )
 END_PROVIDER 


 BEGIN_PROVIDER [ double precision, ao_core_overlap_evec_overlap, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! Overlap between the eigenvectors of the overlap matrix on the cartesian basis
  END_DOC
  call ao_cart_to_ao_core_eigvec(ao_overlap,                 &
                               ao_num,                      &
                               ao_core_overlap_evec_overlap, &
                               ao_num                       &
 )
 END_PROVIDER 


subroutine ao_cart_to_ao_core_eigvec(A_ao_cart,LDA_ao_cart,A_ao_core_eigvec,LDA_ao_core_eigvec)
  implicit none
  BEGIN_DOC
  ! Transform matrix A from the |AO| cartesian basis to the core |AO| basis 
  ! obtained as the eigenvectors of the overlap matrix of the (cartesian) core basis
  !
  ! :math:`B^T.A^c.B`
  !
  ! where :math:`B` is :c:data:`ao_core_overlap_eigvec_coef`, 
  ! the matrix of coefficients from the cartesian AO basis to core AO one,
  ! and :math:`B` is :c:data:`ao_core_coef_transp`, its transpose.
  END_DOC
  integer, intent(in)           :: LDA_ao_cart, LDA_ao_core_eigvec
  double precision, intent(in)  :: A_ao_cart(LDA_ao_cart, ao_num)
  double precision, intent(out) :: A_ao_core_eigvec(LDA_ao_core_eigvec,ao_num)
  double precision :: T(ao_num,ao_num)
  !
  !allocate (T(ao_num,ao_num) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
  !
  call dgemm('N','N', ao_num, ao_num, ao_num,       &
      1.d0, A_ao_cart, LDA_ao_cart,                                &
      ao_core_overlap_eigvec_coef, size(ao_core_coef,1), &
      0.d0, T, size(T,1)                                           &
  )
  ! Notice that for the following dgemm we could have used
  ! ao_core_coef_transp, but instead we transposed with the 'T' argument
  call dgemm('T','N', ao_num, ao_num, ao_num,             &
      1.d0, ao_core_overlap_eigvec_coef, size(ao_core_coef,1), &
      T, ao_num,                                                    &
      0.d0, A_ao_core_eigvec, size(A_ao_core_eigvec,1)           &
  )
  !deallocate(T)
end subroutine


 BEGIN_PROVIDER [ double precision, ao_core_coef_normed, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_core_coef_normed_transp, (ao_num, ao_num) ]
  implicit none
  integer :: i,j
  !
  ao_core_coef_normed(:,:) = 0.d0
  do i = 1,ao_num
    if(ao_core_overlap_eigval(i).gt.1.d-6)then
      do j = 1,ao_num
        ao_core_coef_normed(j,i) =  ao_core_overlap_eigvec_coef(j,i) / dsqrt(ao_core_overlap_eigval(i))
      enddo
    end if
  enddo
  !
  ! Compute transpose
  ao_core_coef_normed_transp(:,:) = 0.d0
  do i = 1,ao_num
    do j = 1,ao_num
      ao_core_coef_normed_transp(j,i) = ao_core_coef_normed(i,j)
    enddo
  enddo
END_PROVIDER


subroutine ao_cart_to_ao_core_normed(A_ao_cart, LDA_ao_cart, A_ao_core_normed, LDA_ao_core_normed)
  implicit none
  integer, intent(in)           :: LDA_ao_cart, LDA_ao_core_normed
  double precision, intent(in)  :: A_ao_cart(LDA_ao_cart, ao_num)
  double precision, intent(out) :: A_ao_core_normed(LDA_ao_core_normed,ao_num)
  double precision :: T(ao_num,ao_num)
  !
  !allocate (T(ao_num,ao_num) )
  call get_AB_prod(A_ao_cart, ao_num, ao_num, ao_core_coef_normed, ao_num, T)
  call get_AB_prod(ao_core_coef_normed_transp, ao_num, ao_num, T, ao_num, A_ao_core_normed)
  !deallocate(T)
end subroutine 


 BEGIN_PROVIDER [ double precision, ao_core_normed_overlap, (ao_num, ao_num) ]
  implicit none
  BEGIN_DOC
  ! AOs overlap matrix in the coreence subspace using normed eigenvectors
  END_DOC 
  !
  call ao_cart_to_ao_core_normed(ao_overlap, ao_num, ao_core_normed_overlap, ao_num)
 END_PROVIDER


subroutine ao_cart_to_ao_core(A_ao_cart,LDA_ao_cart,A_ao_core,LDA_ao_core)
  implicit none
  BEGIN_DOC
  ! Transform matrix A from the |AO| cartesian basis to the core |AO| basis
  !
  ! :math:`(B^{sc})^T.A^c.B^{sc}`
  !
  ! where :math:`B^{sc}` is :c:data:`ao_core_coef`, 
  ! the matrix of coefficients from the cartesian AO basis to core AO one,
  ! and :math:`B^{sc}` is :c:data:`ao_core_coef_transp`, its transpose.
  END_DOC
  integer, intent(in)            :: LDA_ao_cart,LDA_ao_core
  double precision, intent(in)   :: A_ao_cart(LDA_ao_cart,ao_num)
  double precision, intent(out)  :: A_ao_core(LDA_ao_core,ao_num)
  double precision  :: T(ao_num,ao_num)
  !
  !allocate (T(ao_num,ao_num) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
  !
  call dgemm('N','T', ao_num, ao_num, ao_num,  &
      1.d0, A_ao_cart, LDA_ao_cart,            &
      ao_core_coef_normed, size(ao_core_coef,1), &
      0.d0, T, size(T,1)                       &
  )
  ! Notice that for the following dgemm we could have used
  ! ao_core_coef_transp, but instead we transposed with the 'T' argument
  call dgemm('N','N', ao_num, ao_num, ao_num,        &
      1.d0, ao_core_coef_normed, size(ao_core_coef,1), &
      T, ao_num,                                     &
      0.d0, A_ao_core, size(A_ao_core,1)               &
  )
  !deallocate(T)
end subroutine
