 BEGIN_PROVIDER [double precision, ao_two_e_tc_tot_numeric, (ao_num, ao_num, ao_num, ao_num)]
  BEGIN_DOC
  ! Temporary copy of PROVIDER ao_two_e_tc_tot_numeric computed in its numeric version
  END_DOC
  implicit none
  integer                       :: i, j, k, l, m, ipoint
  double precision              :: weight1, ao_k_r, ao_i_r
  double precision              :: der_envsq_x, der_envsq_y, der_envsq_z, lap_envsq
  double precision              :: time0, time1
  double precision, allocatable :: c_mat(:,:,:)
  logical,          external    :: ao_two_e_integral_zero
  double precision, external    :: get_ao_two_e_integral
  double precision, external    :: ao_two_e_integral

  ! if tc_integ_type .eq. numeric it will compute
  PROVIDE tc_integ_type
  PROVIDE env_type
  PROVIDE j2e_type
  PROVIDE j1e_type

  call wall_time(time0)

  print *, ' providing ao_two_e_tc_tot_numeric ...'
  print*, ' j2e_type: ', j2e_type
  print*, ' j1e_type: ', j1e_type
  print*, ' env_type: ', env_type


  PROVIDE tc_integ_type
  print*, ' approach for integrals: ', tc_integ_type

  ! ---

  PROVIDE int2_grad1_u12_sqrd_ao_numeric

  if(tc_save_mem_loops) then

    print*, ' LOOPS are used to evaluate Hermitian part of ao_two_e_tc_tot_numeric ...'

    !$OMP PARALLEL                                              &
    !$OMP DEFAULT (NONE)                                        &
    !$OMP PRIVATE (i, j, k, l, ipoint, ao_i_r, ao_k_r, weight1) &
    !$OMP SHARED (ao_num, n_points_final_grid, ao_two_e_tc_tot_numeric, &
    !$OMP         aos_in_r_array_transp, final_weight_at_r_vector, int2_grad1_u12_sqrd_ao_numeric)
    !$OMP DO COLLAPSE(3)
    do i = 1, ao_num
      do k = 1, ao_num
        do l = 1, ao_num
          do j = 1, ao_num
            ao_two_e_tc_tot_numeric(j,l,k,i) = 0.d0
            do ipoint = 1, n_points_final_grid
              weight1 = final_weight_at_r_vector(ipoint)
              ao_i_r = aos_in_r_array_transp(ipoint,i)
              ao_k_r = aos_in_r_array_transp(ipoint,k)
              ao_two_e_tc_tot_numeric(j,l,k,i) = ao_two_e_tc_tot_numeric(j,l,k,i) + int2_grad1_u12_sqrd_ao_numeric(j,l,ipoint) * weight1 * ao_i_r * ao_k_r
            enddo
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  else

    print*, ' DGEMM are used to evaluate Hermitian part of ao_two_e_tc_tot_numeric ...'

    allocate(c_mat(n_points_final_grid,ao_num,ao_num))
    !$OMP PARALLEL               &
    !$OMP DEFAULT (NONE)         &
    !$OMP PRIVATE (i, k, ipoint) &
    !$OMP SHARED (aos_in_r_array_transp, c_mat, ao_num, n_points_final_grid, final_weight_at_r_vector)
    !$OMP DO SCHEDULE (static)
    do i = 1, ao_num
      do k = 1, ao_num
        do ipoint = 1, n_points_final_grid
          c_mat(ipoint,k,i) = final_weight_at_r_vector(ipoint) * aos_in_r_array_transp(ipoint,i) * aos_in_r_array_transp(ipoint,k)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, 1.d0                 &
              , int2_grad1_u12_sqrd_ao_numeric(1,1,1), ao_num*ao_num, c_mat(1,1,1), n_points_final_grid &
              , 0.d0, ao_two_e_tc_tot_numeric(1,1,1,1), ao_num*ao_num)
    deallocate(c_mat)
  endif

  FREE int2_grad1_u12_sqrd_ao_numeric

  call wall_time(time1)
  print*, ' done with Hermitian part after (min) ', (time1 - time0) / 60.d0
  call print_memory_usage()

  ! ---

  if(.not. var_tc) then

    PROVIDE int2_grad1_u12_vect_ao_numeric

    if(tc_save_mem_loops) then

      print*, ' LOOPS are used to evaluate non-Hermitian part of ao_two_e_tc_tot_numeric ...'
      print*, ' TODO : this can be improved by doing 3 separated step'

      !$OMP PARALLEL                                                 &
      !$OMP DEFAULT (NONE)                                           &
      !$OMP PRIVATE (i, j, k, l, ipoint, ao_i_r, ao_k_r, weight1)    &
      !$OMP SHARED (ao_num, n_points_final_grid, ao_two_e_tc_tot_numeric,    &
      !$OMP         aos_in_r_array_transp, final_weight_at_r_vector, &
      !$OMP         int2_grad1_u12_vect_ao_numeric, aos_grad_in_r_array_transp_bis)
      !$OMP DO COLLAPSE(3)
      do i = 1, ao_num
        do k = 1, ao_num
          do l = 1, ao_num
            do j = 1, ao_num
              do ipoint = 1, n_points_final_grid
                weight1 = 0.5d0 * final_weight_at_r_vector(ipoint)
                ao_i_r  = aos_in_r_array_transp(ipoint,i)
                ao_k_r  = aos_in_r_array_transp(ipoint,k)
                ao_two_e_tc_tot_numeric(j,l,k,i) = ao_two_e_tc_tot_numeric(j,l,k,i) &
                                         - weight1 * int2_grad1_u12_vect_ao_numeric(j,l,ipoint,1) * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,1) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,1)) &
                                         - weight1 * int2_grad1_u12_vect_ao_numeric(j,l,ipoint,2) * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,2) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,2)) &
                                         - weight1 * int2_grad1_u12_vect_ao_numeric(j,l,ipoint,3) * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,3) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,3))
              enddo
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

    else

      print*, ' DGEMM are used to evaluate non-Hermitian part of ao_two_e_tc_tot_numeric ...'

      allocate(c_mat(n_points_final_grid,ao_num,ao_num))
      do m = 1, 3
        !$OMP PARALLEL                                                              &
        !$OMP DEFAULT (NONE)                                                        &
        !$OMP PRIVATE (i, k, ipoint, weight1, ao_i_r, ao_k_r)                       & 
        !$OMP SHARED (aos_in_r_array_transp, aos_grad_in_r_array_transp_bis, c_mat, & 
        !$OMP         ao_num, n_points_final_grid, final_weight_at_r_vector, m)
        !$OMP DO SCHEDULE (static)
        do i = 1, ao_num
          do k = 1, ao_num
            do ipoint = 1, n_points_final_grid

              weight1 = 0.5d0 * final_weight_at_r_vector(ipoint)
              ao_i_r  = aos_in_r_array_transp(ipoint,i)
              ao_k_r  = aos_in_r_array_transp(ipoint,k)

              c_mat(ipoint,k,i) = weight1 * (ao_k_r * aos_grad_in_r_array_transp_bis(ipoint,i,m) - ao_i_r * aos_grad_in_r_array_transp_bis(ipoint,k,m))
            enddo
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid, -1.d0           &
                  , int2_grad1_u12_vect_ao_numeric(1,1,1,m), ao_num*ao_num, c_mat(1,1,1), n_points_final_grid &
                  , 1.d0, ao_two_e_tc_tot_numeric(1,1,1,1), ao_num*ao_num)
      enddo
      deallocate(c_mat)

    end if

  endif ! var_tc

  call wall_time(time1)
  print*, ' done with non-Hermitian part after (min) ', (time1 - time0) / 60.d0
  call print_memory_usage()

  ! ---

  call sum_A_At(ao_two_e_tc_tot_numeric(1,1,1,1), ao_num*ao_num)

  ! ---

  double precision :: integ_val

  print*, ' adding ERI to ao_two_e_tc_tot_numeric ...'

  if(tc_save_mem) then
    print*, ' ao_integrals_map will not be used'
    !$OMP PARALLEL DEFAULT(NONE)                     &
    !$OMP PRIVATE(i, j, k, l, integ_val) & 
    !$OMP SHARED(ao_num, ao_two_e_tc_tot_numeric)
    !$OMP DO COLLAPSE(3)
    do j = 1, ao_num
      do l = 1, ao_num
        do i = 1, ao_num
          do k = 1, ao_num
            if(.not. ao_two_e_integral_zero(i,j,k,l)) then
                          ! i,k : r1    j,l : r2
              integ_val = ao_two_e_integral(i,k,j,l)
              ao_two_e_tc_tot_numeric(k,i,l,j) = ao_two_e_tc_tot_numeric(k,i,l,j) + integ_val
            endif
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  else
    !print*, ' ao_integrals_map will be used'
    !PROVIDE ao_integrals_map
    print*,'Cholesky vectors will be used '
    double precision :: get_ao_integ_chol,eri
    eri = get_ao_integ_chol(1,1,1,1) ! FOR OPENMP 
    !$OMP PARALLEL DEFAULT(NONE)                            &
    !!$OMP SHARED(ao_num, ao_two_e_tc_tot_numeric, ao_integrals_map) &
    !$OMP SHARED(ao_num, ao_two_e_tc_tot_numeric) &
    !$OMP PRIVATE(i, j, k, l,eri)
    !$OMP DO COLLAPSE(3)
    do j = 1, ao_num
      do l = 1, ao_num
        do i = 1, ao_num
          do k = 1, ao_num
             !< 1:i, 2:j | 1:k, 2:l > 
             !eri =  get_ao_two_e_integral(i, j, k, l, ao_integrals_map)
             eri = get_ao_integ_chol(i,k,j,l)
            ao_two_e_tc_tot_numeric(k,i,l,j) = ao_two_e_tc_tot_numeric(k,i,l,j) + eri
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    !FREE ao_integrals_map
  endif

  ! Why free this stuff?
  if((tc_integ_type .eq. "numeric") .and. (.not. tc_save_mem)) then
    FREE int2_grad1_u12_vect_ao_numeric int2_grad1_u12_sqrd_ao_numeric
  endif

  call wall_time(time1)
  print*, ' Wall time for ao_two_e_tc_tot_numeric (min) = ', (time1 - time0) / 60.d0
  call print_memory_usage()

END_PROVIDER 


 BEGIN_PROVIDER [double precision, int2_grad1_u12_vect_ao_numeric, (ao_num,ao_num,n_points_final_grid,3)]
&BEGIN_PROVIDER [double precision, int2_grad1_u12_sqrd_ao_numeric, (ao_num,ao_num,n_points_final_grid)  ]
  BEGIN_DOC
  ! Temporary copy of
  ! PROVIDER int2_grad1_u12_ao_num
  ! PROVIDER int2_grad1_u12_square_ao_num
  ! in the specific case of numerical evaluation
  END_DOC
  implicit none
  integer                       :: ipoint, i, j, m, jpoint
  integer                       :: n_blocks, n_rest, n_pass
  integer                       :: i_blocks, i_rest, i_pass, ii
  double precision              :: time0, time1
  double precision              :: mem, n_double
  double precision, allocatable :: tmp(:,:,:)
  double precision, allocatable :: tmp_grad1_u12(:,:,:), tmp_grad1_u12_squared(:,:)

  print*, ' providing int2_grad1_u12_vect_ao_numeric & int2_grad1_u12_sqrd_ao_numeric ...'
  call wall_time(time0)

  PROVIDE final_weight_at_r_vector_extra aos_in_r_array_extra

  allocate(tmp(n_points_extra_final_grid,ao_num,ao_num))
  !$OMP PARALLEL               &
  !$OMP DEFAULT (NONE)         &
  !$OMP PRIVATE (j, i, jpoint) &
  !$OMP SHARED (tmp, ao_num, n_points_extra_final_grid, final_weight_at_r_vector_extra, aos_in_r_array_extra_transp)
  !$OMP DO SCHEDULE (static)
  do j = 1, ao_num
    do i = 1, ao_num
      do jpoint = 1, n_points_extra_final_grid
        tmp(jpoint,i,j) = final_weight_at_r_vector_extra(jpoint) * aos_in_r_array_extra_transp(jpoint,i) * aos_in_r_array_extra_transp(jpoint,j)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call total_memory(mem)
  mem      = max(1.d0, qp_max_mem - mem)
  n_double = mem * 1.d8
  n_blocks = int(min(n_double / (n_points_extra_final_grid * 4.d0), 1.d0*n_points_final_grid))
  n_rest   = int(mod(n_points_final_grid, n_blocks))
  n_pass   = int((n_points_final_grid - n_rest) / n_blocks)

  call write_int(6, n_pass, 'Number of passes')
  call write_int(6, n_blocks, 'Size of the blocks')
  call write_int(6, n_rest, 'Size of the last block')

  
  allocate(tmp_grad1_u12_squared(n_points_extra_final_grid,n_blocks))
  allocate(tmp_grad1_u12(n_points_extra_final_grid,n_blocks,3))
  
  do i_pass = 1, n_pass
    ii = (i_pass-1)*n_blocks + 1
  
    !$OMP PARALLEL                   &
    !$OMP DEFAULT (NONE)             &
    !$OMP PRIVATE (i_blocks, ipoint) &
    !$OMP SHARED (n_blocks, n_points_extra_final_grid, ii, final_grid_points, tmp_grad1_u12, tmp_grad1_u12_squared)
    !$OMP DO 
    do i_blocks = 1, n_blocks
      ipoint = ii - 1 + i_blocks ! r1
      call get_grad1_u12_withsq_r1_seq(ipoint, n_points_extra_final_grid, tmp_grad1_u12(1,i_blocks,1) &
                                                                        , tmp_grad1_u12(1,i_blocks,2) &
                                                                        , tmp_grad1_u12(1,i_blocks,3) &
                                                                        , tmp_grad1_u12_squared(1,i_blocks))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    do m = 1, 3
      call dgemm( "T", "N", ao_num*ao_num, n_blocks, n_points_extra_final_grid, 1.d0                     &
                , tmp(1,1,1), n_points_extra_final_grid, tmp_grad1_u12(1,1,m), n_points_extra_final_grid &
                , 0.d0, int2_grad1_u12_vect_ao_numeric(1,1,ii,m), ao_num*ao_num) 
    enddo
    call dgemm( "T", "N", ao_num*ao_num, n_blocks, n_points_extra_final_grid, -0.5d0                         &
              , tmp(1,1,1), n_points_extra_final_grid, tmp_grad1_u12_squared(1,1), n_points_extra_final_grid &
              , 0.d0, int2_grad1_u12_sqrd_ao_numeric(1,1,ii), ao_num*ao_num) 
  enddo
  
  deallocate(tmp_grad1_u12, tmp_grad1_u12_squared)
  
  if(n_rest .gt. 0) then
  
    allocate(tmp_grad1_u12_squared(n_points_extra_final_grid,n_rest))
    allocate(tmp_grad1_u12(n_points_extra_final_grid,n_rest,3))
  
    ii = n_pass*n_blocks + 1

    !$OMP PARALLEL                 &
    !$OMP DEFAULT (NONE)           &
    !$OMP PRIVATE (i_rest, ipoint) &
    !$OMP SHARED (n_rest, n_points_extra_final_grid, ii, final_grid_points, tmp_grad1_u12, tmp_grad1_u12_squared)
    !$OMP DO 
    do i_rest = 1, n_rest
      ipoint = ii - 1 + i_rest ! r1
      call get_grad1_u12_withsq_r1_seq(ipoint, n_points_extra_final_grid, tmp_grad1_u12(1,i_rest,1) &
                                                                        , tmp_grad1_u12(1,i_rest,2) &
                                                                        , tmp_grad1_u12(1,i_rest,3) &
                                                                        , tmp_grad1_u12_squared(1,i_rest))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
    do m = 1, 3
      call dgemm( "T", "N", ao_num*ao_num, n_rest, n_points_extra_final_grid, 1.d0                       &
                , tmp(1,1,1), n_points_extra_final_grid, tmp_grad1_u12(1,1,m), n_points_extra_final_grid &
                , 0.d0, int2_grad1_u12_vect_ao_numeric(1,1,ii,m), ao_num*ao_num) 
    enddo
    call dgemm( "T", "N", ao_num*ao_num, n_rest, n_points_extra_final_grid, -0.5d0                           &
              , tmp(1,1,1), n_points_extra_final_grid, tmp_grad1_u12_squared(1,1), n_points_extra_final_grid &
              , 0.d0, int2_grad1_u12_sqrd_ao_numeric(1,1,ii), ao_num*ao_num) 

    deallocate(tmp_grad1_u12, tmp_grad1_u12_squared)
  endif

  deallocate(tmp)

  call wall_time(time1)
  print*, ' wall time for int2_grad1_u12_vect_ao_numeric & int2_grad1_u12_sqrd_ao_numeric = (min)', (time1-time0) / 60.d0
  call print_memory_usage()

END_PROVIDER
