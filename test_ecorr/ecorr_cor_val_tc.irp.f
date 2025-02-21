use bitmasks ! you need to include the bitmasks_module.f90 features

BEGIN_PROVIDER [ integer, n_det_cisd_tc ]
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
    norm += psi_l_coef_bi_ortho(i,1) * psi_r_coef_bi_ortho(i,1)
    if (degree.le.2) then 
       n_det_cisd += 1
    endif
  enddo
END_PROVIDER 


BEGIN_PROVIDER [ integer(bit_kind), psi_det_cisd_proj_tc, (N_int,2,n_det_cisd) ]
&BEGIN_PROVIDER [ integer, psi_core_val_decomposition_tc, (n_det_cisd) ]
&BEGIN_PROVIDER [ double precision, psi_coef_cisd_proj_tc_r, (n_det_cisd) ]
&BEGIN_PROVIDER [ double precision, psi_coef_cisd_proj_tc_l, (n_det_cisd) ]
  BEGIN_DOC
  ! psi_det_cisd_proj: CISD determinants in bit_kind format
  ! psi_coef_cisd_proj_tc_r: RIGHT coefficients of CISD determinants found in TC computations
  ! psi_coef_cisd_proj_tc_l: LEFT coefficients of CISD determinants found in TC computations
  ! psi_core_val_decomposition: indices of CISD determinants involving either core or valence excitations
  END_DOC
  implicit none
  integer :: i, degree, i_cisd, core_fill, core_test
  double precision :: cisd_norm
  ! First determinant is always the HF one 
  psi_det_cisd_proj_tc(1:N_int,1:2,1) = psi_det(1:N_int,1:2,1)
  ! Left/Right HF has different coefficients
  psi_coef_cisd_proj_tc_r(1) = psi_r_coef_bi_ortho(1,1)
  psi_coef_cisd_proj_tc_l(1) = psi_l_coef_bi_ortho(1,1)
  ! Norm is initialised with the HF contribution 
  cisd_norm =  psi_r_coef_bi_ortho(1,1)*psi_l_coef_bi_ortho(1,1)
  ! Initialise the decomposition to zero
  psi_core_val_decomposition_tc = 0
  ! Iterator on the CISD determinants
  i_cisd = 1
  ! Skip first determinant because it is the HF reference
  do i = 2, N_det
    call get_excitation_degree(ref_bitmask,psi_det(1,1,i),degree,N_int)
    cisd_norm += psi_l_coef_bi_ortho(i,1) * psi_r_coef_bi_ortho(i,1)
    if (degree.le.2) then
       i_cisd += 1
       psi_det_cisd_proj_tc(1:N_int,1:2,i_cisd) = psi_det(1:N_int,1:2,i)
       ! ibits(n,i,l) get the value(s) starting in position i and ending in i+l
       core_fill = ibits(psi_det(1,1,i), 0, 1) + ibits(psi_det(1,2,i), 0, 1)
     ! core_test = popcnt(iand(core_bitmask(1,1),psi_det(1,1,i)))
     ! core_test += popcnt(iand(core_bitmask(1,2),psi_det(1,2,i)))
     ! if(core_test.ne.core_fill)then
     !  print*,'core_test    core_fill = ',core_test  , core_fill 
     !  call debug_det(psi_det(1,1,i),N_int)
     ! endif
       ! core_fill = 0 --> empty core
       ! core_fill = 1 --> half empty core
       ! core_fill = 2 --> full core
       psi_core_val_decomposition_tc(i_cisd) = core_fill+3*(degree-1)
       ! core_fill | degree | result of core_fill+3*(degree-1)
       !         0 |      1 | 0 = impossible
       !         1 |      1 | 1 = single excitation (core)
       !         2 |      1 | 2 = single excitation (valence)
       !         0 |      2 | 3 = double excitation (core)
       !         1 |      2 | 4 = double excitation (core+valence)
       !         2 |      2 | 5 = double excitation (valence)
       psi_coef_cisd_proj_tc_r(i_cisd) = psi_r_coef_bi_ortho(i,1)
       psi_coef_cisd_proj_tc_l(i_cisd) = psi_l_coef_bi_ortho(i,1)
    endif
  enddo
  print*,'norm total = ', cisd_norm
END_PROVIDER


BEGIN_PROVIDER [ double precision, e_corr_cisd_proj_tc_r]
&BEGIN_PROVIDER [ double precision, e_corr_cisd_proj_tc_l]
&BEGIN_PROVIDER [ double precision, e_corr_cisd_proj_tc_r_cc]
&BEGIN_PROVIDER [ double precision, e_corr_cisd_proj_tc_r_vv]
&BEGIN_PROVIDER [ double precision, e_corr_cisd_proj_tc_r_cv]
&BEGIN_PROVIDER [ double precision, e_corr_cisd_proj_tc_l_cc]
&BEGIN_PROVIDER [ double precision, e_corr_cisd_proj_tc_l_vv]
&BEGIN_PROVIDER [ double precision, e_corr_cisd_proj_tc_l_cv]
  BEGIN_DOC
  ! e_corr_cisd_proj_tc_r: CISD correlation energy computed as projection of FCI only on S and D RIGHT determinants
  ! e_corr_cisd_proj_tc_l: CISD correlation energy computed as projection of FCI only on S and D LEFT determinants
  ! Decomposition:
  !   cc : core contributions only 
  !   vv : valence contributions only 
  !   cv : core and valence contributions
  END_DOC
  implicit none
  integer :: i, degree
  integer :: mask_cc(n_det_cisd), mask_vv(n_det_cisd), mask_cv(n_det_cisd)
  double precision, allocatable :: delta_e(:)
  allocate(delta_e(n_det_cisd))
  double precision :: h0i(n_det_cisd),hii
  double precision :: hi0(n_det_cisd)
  double precision :: amplitude_r(n_det_cisd)
  double precision :: amplitude_l(n_det_cisd)
  e_corr_cisd_proj_tc_r = 0.d0
  e_corr_cisd_proj_tc_l = 0.d0
  e_corr_cisd_proj_tc_r_cc = 0.d0
  e_corr_cisd_proj_tc_l_cc = 0.d0
  e_corr_cisd_proj_tc_r_cv = 0.d0
  e_corr_cisd_proj_tc_l_cv = 0.d0
  e_corr_cisd_proj_tc_r_vv = 0.d0
  e_corr_cisd_proj_tc_l_vv = 0.d0
  !
  ! Skip first determinant because it is the HF one
  do i = 2, n_det_cisd
    call htilde_mu_mat_opt_bi_ortho_tot(psi_det_cisd_proj_tc(1,1,i),psi_det_cisd_proj_tc(1,1,i),N_int,hii)
    delta_e(i) = ref_tc_energy_tot - hii
    call htilde_mu_mat_opt_bi_ortho_tot(ref_bitmask,psi_det_cisd_proj_tc(1,1,i),N_int,h0i(i))
    call htilde_mu_mat_opt_bi_ortho_tot(psi_det_cisd_proj_tc(1,1,i),ref_bitmask,N_int,hi0(i))
    !OLD call i_H_j(ref_bitmask,psi_det_cisd_proj_tc(1,1,i),N_int,h0i(i))
  enddo
  amplitude_r = psi_coef_cisd_proj_tc_r/psi_coef_cisd_proj_tc_r(1) 
  amplitude_l = psi_coef_cisd_proj_tc_l/psi_coef_cisd_proj_tc_l(1) 
  ! Mask to separate core, valence and core+valence contributions
  ! Fortran conversion from logical to integer is: .false.-->0, .true.-->-1, hence the prefactor -1
  mask_cc = -1*((psi_core_val_decomposition_tc .eq. 1).or.(psi_core_val_decomposition_tc .eq. 3))
  mask_vv = -1*((psi_core_val_decomposition_tc .eq. 2).or.(psi_core_val_decomposition_tc .eq. 5))
  mask_cv = -1*(psi_core_val_decomposition_tc .eq. 4)
  ! CISD correlation energy decomposition and total
  ! Right
  e_corr_cisd_proj_tc_r_cc = sum(amplitude_r*h0i*mask_cc)
  e_corr_cisd_proj_tc_r_vv = sum(amplitude_r*h0i*mask_vv)
  e_corr_cisd_proj_tc_r_cv = sum(amplitude_r*h0i*mask_cv)
  e_corr_cisd_proj_tc_r = e_corr_cisd_proj_tc_r_cc + e_corr_cisd_proj_tc_r_vv + e_corr_cisd_proj_tc_r_cv
  ! Left
  e_corr_cisd_proj_tc_l_cc = sum(amplitude_l*hi0*mask_cc)
  e_corr_cisd_proj_tc_l_vv = sum(amplitude_l*hi0*mask_vv)
  e_corr_cisd_proj_tc_l_cv = sum(amplitude_l*hi0*mask_cv)
  e_corr_cisd_proj_tc_l = e_corr_cisd_proj_tc_l_cc + e_corr_cisd_proj_tc_l_vv + e_corr_cisd_proj_tc_l_cv
  write(*,*) "sum over masck_cc", sum(mask_cc)
  write(*,*) "sum over masck_cv", sum(mask_cv)
  write(*,*) "sum over masck_vv", sum(mask_vv)
  double precision :: cv_corr_per_orb(5),cv_corr_per_orb_pert(5)
  integer                        :: h1,p1,s1,h2,p2,s2
  integer                        :: exc(0:2,2,2)
  double precision               :: phase
  cv_corr_per_orb = 0.d0
  cv_corr_per_orb_pert = 0.d0
  do i = 2, n_det_cisd
   if(psi_core_val_decomposition_tc(i) .eq. 4)then
    call get_excitation(ref_bitmask,psi_det_cisd_proj_tc(1,1,i),exc,degree,phase,N_int)
    call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
    if(h1.ne.1)then
     cv_corr_per_orb(h1) += amplitude_r(i)*h0i(i)
     cv_corr_per_orb_pert(h1) += hi0(i)/delta_e(i)*h0i(i)
    else
     cv_corr_per_orb(h2) += amplitude_r(i)*h0i(i)
     cv_corr_per_orb_pert(h2) += hi0(i)/delta_e(i)*h0i(i)
    endif
   endif
  enddo
  double precision :: accu
  accu = 0.d0
  do i = 2,5
   accu += cv_corr_per_orb(i)
   print*,'i',i
   print*,cv_corr_per_orb(i),cv_corr_per_orb_pert(i)
  enddo
  print*,'accu = ',accu
END_PROVIDER 
