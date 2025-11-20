program core_tcxc_test

  BEGIN_DOC
  ! Testing the core exchange potential integral in TC
  END_DOC

  implicit none

  write(*,'(A)') repeat('=', 70)
  write(*,*) "TEST 0"
  call test_ao_overlap_exact_vs_numeric

end program core_tcxc_test


subroutine test_ao_overlap_exact_vs_numeric
  use iso_fortran_env, only: out_unit => output_unit
  implicit none
  BEGIN_DOC
  ! Check that the AO overlap matrix computed numerically is reasonable
  END_DOC
  double precision :: difference


  write(*,*) "ao_overlap_grid1 VS ao_overlap"
  
  write(*,*) "Depending on ao_num..."
  call write_int(out_unit, ao_num, 'ao_num')                   
  write(*,*) "... Expected sizes:"
  call write_int(out_unit, ao_num*ao_num   , 'ao_num^2')                   
  call write_int(out_unit, size(ao_overlap), 'size(ao_overlap)')                   

  write(*,*) 

  difference = sum(abs(ao_overlap_grid1(:,:) - ao_overlap(:,:)))

  write(*,*) "Difference =           ", difference
  write(*,*) "Difference/n_entries = ", difference/size(ao_overlap)
  
  write(*,*) 

end subroutine test_ao_overlap_exact_vs_numeric
