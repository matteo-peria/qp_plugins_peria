program test_ecorr
  implicit none
  BEGIN_DOC
  ! Test of the plugin on correlation energy 
  END_DOC
  read_wf = .true.      ! Setting by default |true|: read the wave function from the |EZFIO| file 
  touch read_wf         ! Update all the subroutines using read_wf on its new value
  call routine_pouet
end

subroutine routine_pouet
  implicit none
  integer :: i
  double precision :: norm
  norm = 0.d0
  do i = 1,n_det_cisd
!   print*,i
!   call debug_det(psi_det_cisd_proj(1,1,i),N_int)
!   print*,psi_coef_cisd_proj(i)
    norm += psi_coef_cisd_proj(i)*psi_coef_cisd_proj(i)
  enddo
  print*,'norm = ', norm

  print*, ""
  print*, "Correlation energy from S and D contributions decomposition:"
  write(*,'(A25,A15)') 'From', 'Value'
  write(*,'(A25,F15.7)') 'Core only (S and D)',        e_corr_cisd_proj_cc
  write(*,'(A25,F15.7)') 'Valence only (S and D)',     e_corr_cisd_proj_vv
  write(*,'(A25,F15.7)') 'Core and valence (S)',       e_corr_cisd_proj_cv
  write(*,'(A25,F15.7)') 'Total correlation',          e_corr_cisd_proj   
end subroutine

