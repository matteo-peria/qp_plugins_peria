use bitmasks ! you need to include the bitmasks_module.f90 features

BEGIN_PROVIDER [ integer, n_det_cisd ]
  BEGIN_DOC
  ! Number of determinants which have a degree of excitation <= to 2
  END_DOC
  implicit none
  integer :: i, degree
  double precision :: norm
  n_det_cisd = 0
  norm = 0.d0
  do i = 1, N_det
    call get_excitation_degree(ref_bitmask,psi_det(1,1,i),degree,N_int)
    norm += psi_coef(i,1) * psi_coef(i,1)
    if (degree.le.2) then
       n_det_cisd += 1
    endif
  enddo
END_PROVIDER 


BEGIN_PROVIDER [ integer(bit_kind), psi_det_cisd_proj, (N_int,2,n_det_cisd) ]
&BEGIN_PROVIDER [ integer, psi_core_val_decomposition, (n_det_cisd) ]
&BEGIN_PROVIDER [ double precision, psi_coef_cisd_proj, (n_det_cisd) ]
  BEGIN_DOC
  ! psi_det_cisd_proj: CISD determinants in bit_kind format
  ! psi_coef_cisd_proj: CISD determinants respective coefficients
  ! psi_core_val_decomposition: indices of CISD determinants involving either core or valence excitations
  END_DOC
  implicit none
  integer :: i, degree, i_cisd, core_fill
  ! First determinant is always the HF one
  psi_det_cisd_proj(1:N_int,1:2,1) = psi_det(1:N_int,1:2,1)
  ! First coefficient is the HF one
  psi_coef_cisd_proj(1) = psi_coef(1,1)
  ! Initialise the decomposition to zero
  psi_core_val_decomposition = 0
  ! Iterator on the CISD determinants
  i_cisd = 1
  ! Skip first determinant because it is the HF reference
  do i = 2, N_det
    call get_excitation_degree(ref_bitmask,psi_det(1,1,i),degree,N_int)
    if (degree.le.2) then
       i_cisd += 1
       psi_det_cisd_proj(1:N_int,1:2,i_cisd) = psi_det(1:N_int,1:2,i)
       ! ibits(n,i,l) get the value(s) starting in position i and ending in i+l
       core_fill = ibits(psi_det(1,1,i), 0, 1) + ibits(psi_det(1,2,i), 0, 1)
       ! core_fill = 0 --> empty core
       ! core_fill = 1 --> half empty core
       ! core_fill = 2 --> full core
       psi_core_val_decomposition(i_cisd) = core_fill+3*(degree-1)
       ! core_fill | degree | result
       !         0 |      1 | 0 = impossible
       !         1 |      1 | 1 = single excitation (core)
       !         2 |      1 | 2 = single excitation (valence)
       !         0 |      2 | 3 = double excitation (core)
       !         1 |      2 | 4 = double excitation (core+valence)
       !         2 |      2 | 5 = double excitation (valence)
       psi_coef_cisd_proj(i_cisd) = psi_coef(i,1)
    endif
  enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, e_corr_cisd_proj]
&BEGIN_PROVIDER [ double precision, e_corr_cisd_proj_cc]
&BEGIN_PROVIDER [ double precision, e_corr_cisd_proj_vv]
&BEGIN_PROVIDER [ double precision, e_corr_cisd_proj_cv]
  BEGIN_DOC
  ! e_corr_cisd_proj: CISD correlation energy computed as projection of FCI only on S and D determinants
  END_DOC
  implicit none
  integer :: i, degree
  integer :: mask_cc(n_det_cisd), mask_vv(n_det_cisd), mask_cv(n_det_cisd)
  double precision :: h0i(n_det_cisd)
  double precision :: amplitude(n_det_cisd)
  e_corr_cisd_proj = 0.d0
  !
  ! Skip first determinant because it is the HF one
  do i = 2, n_det_cisd
  ! htilde_mu_mat_opt_bi_ortho_tot(psi_det(1,1,j), psi_det(1,1,i), N_int, htot)
    call i_H_j(ref_bitmask,psi_det_cisd_proj(1,1,i),N_int,h0i(i))
  enddo
  amplitude = psi_coef_cisd_proj/psi_coef_cisd_proj(1) 
  ! Mask to separate core, valence and core+valence contributions
  ! Fortran conversion from logical to integer is: .false.-->0, .true.-->-1, hence the prefactor -1
  mask_cc = -1*((psi_core_val_decomposition .eq. 1).or.(psi_core_val_decomposition .eq. 3))
  mask_vv = -1*((psi_core_val_decomposition .eq. 2).or.(psi_core_val_decomposition .eq. 5))
  mask_cv = -1*(psi_core_val_decomposition .eq. 4)
  ! CISD correlation energy decomposition and total
  e_corr_cisd_proj_cc = sum(amplitude*h0i*mask_cc)
  e_corr_cisd_proj_vv = sum(amplitude*h0i*mask_vv)
  e_corr_cisd_proj_cv = sum(amplitude*h0i*mask_cv)
  e_corr_cisd_proj = e_corr_cisd_proj_cc + e_corr_cisd_proj_vv + e_corr_cisd_proj_cv
END_PROVIDER 
