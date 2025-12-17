BEGIN_PROVIDER [double precision, ao_two_e_tc_tot_adapt, (ao_num, ao_num, ao_num, ao_num)]
  BEGIN_DOC
  ! Compute TC 2-electrons integrals comprising the K-integrals (both the 
  ! hermitian and the non-hermitian integrals) and the usual Coulomb integrals.
  !
  ! Chemist notation: (ki|lj)
  ! Physicist notation: <lk|ij>
  !
  ! ao_two_e_tc_tot_adapt(k,i,l,j) = (ki| V^TC(r_12) |lj) = <lk| V^TC(r_12) |ji> 
  !
  ! where V^TC(r_12) is the total TC operator. 
  ! The result is computed as
  !
  ! ao_two_e_tc_tot_adapt(k,i,l,j) = tc_grad_and_lapl_ao(k,i,l,j) 
  !                          + tc_grad_square_ao(k,i,l,j) 
  !                          + ao_two_e_coul(k,i,l,j)
  ! where:
  !
  ! tc_grad_and_lapl_ao(k,i,l,j) = 
  !   = < kl | -1/2 \Delta_1 u(r1,r2) - \grad_1 u(r1,r2) . \grad_1 | ij >
  !  = -1/2 \int dr1 (phi_k(r1) \grad_r1 phi_i(r1) - phi_i(r1) \grad_r1 phi_k(r1)) . \int dr2\grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2) 
  !  =  1/2 \int dr1 (phi_k(r1) \grad_r1 phi_i(r1) - phi_i(r1) \grad_r1 phi_k(r1)) . \int dr2 (-1) \grad_r1 u(r1,r2) \phi_l(r2) \phi_j(r2) 
  !
  ! tc_grad_square_ao(k,i,l,j) = -1/2 <kl| |\grad_1 u(r1,r2)|^2 + |\grad_2 u(r1,r2)|^2 |ij>
  !
  ! ao_two_e_coul(k,i,l,j) = <kl| 1/r12 |ij> = (ki| 1/r12 |lj)
  !
  END_DOC
  implicit none
  integer                       :: i, j, k, l, m, i1
  double precision              :: w1, ao_k_r1, ao_i_r1
  double precision              :: time0, time1
  double precision, allocatable :: c_mat_grad_sqrd(:,:,:)
  double precision, allocatable :: c_mat_grad_vect(:,:,:,:)
  logical,          external    :: ao_two_e_integral_zero
  double precision, external    :: get_ao_two_e_integral
  double precision, external    :: ao_two_e_integral

  !double precision :: int2_grad1_u12_vect_ao_at_r1(ao_num, ao_num, 3)
  !double precision :: int2_grad1_u12_sqrd_ao_at_r1(ao_num, ao_num)
  !!
  !double precision :: int2_grad1_u12_vect_ao(ao_num, ao_num, n_points_final_grid, 3)
  !double precision :: int2_grad1_u12_sqrd_ao(ao_num, ao_num, n_points_final_grid)

  double precision, allocatable :: int2_grad1_u12_vect_ao_at_r1(:,:,:)
  double precision, allocatable :: int2_grad1_u12_sqrd_ao_at_r1(:,:)
  !
  double precision, allocatable :: int2_grad1_u12_vect_ao(:,:,:,:)
  double precision, allocatable :: int2_grad1_u12_sqrd_ao(:,:,:)

  ! Allocate once for all the floating grid and the weights on fixed and float
  double precision, allocatable :: grid_float_points(:,:,:,:)
  double precision, allocatable :: grid_fixed_weights(:,:,:)
  double precision, allocatable :: grid_float_weights(:,:,:)
  ! Dummy variables necessary for pruning on floating grid (not coded yet)
  integer :: n_fixed_pts_effective(nucl_num)
  integer :: n_float_pts_effective
  integer :: n_pts_effective_max


  double precision :: r1(3)
  double precision :: ao_i_grad1(3)
  double precision :: ao_k_grad1(3)

  print *, ' providing AO_TWO_E_TC_TOT_ADAPT ...'
  ao_two_e_tc_tot_adapt = 0.d0

  allocate(int2_grad1_u12_vect_ao_at_r1(ao_num, ao_num, 3))
  allocate(int2_grad1_u12_sqrd_ao_at_r1(ao_num, ao_num))

  allocate(int2_grad1_u12_vect_ao(ao_num, ao_num, n_points_final_grid, 3))
  allocate(int2_grad1_u12_sqrd_ao(ao_num, ao_num, n_points_final_grid))

  ! New
  allocate(grid_float_points(3, n_points_ang_float_grid, n_points_rad_float_grid, 1))
  allocate(grid_fixed_weights(n_points_ang_extra_grid, n_points_rad_extra_grid, nucl_num))
  allocate(grid_float_weights(n_points_ang_float_grid, n_points_rad_float_grid, 1))


  call wall_time(time0)

  print *, ' providing AO_TWO_E_TC_TOT_ADAPT ...'
  !print*, ' j2e_type: ', j2e_type
  !print*, ' env_type: ', env_type
  

  ! start: added this temporarily
  read_tc_integ = .false.
  ! end: added this temporarily
  if (read_tc_integ) then
    print*, ' Reading ao_two_e_tc_tot_adapt from ', trim(ezfio_filename) // '/work/ao_two_e_tc_tot_adapt'
    open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/ao_two_e_tc_tot_adapt', action="read")
    do i = 1, ao_num
      read(11) ao_two_e_tc_tot_adapt(:,:,:,i)
    enddo
    close(11)

  else

    ! In a first version of the subroutine, when it was based on a fixed grid 
    ! in r2, the non-Hermitian and the Hermitian contributions to the integral 
    ! were computed in separate loops, for code readability and code testing.
    ! Now, using an adaptive grid in r2, which depends on the value of r1, 
    ! implies computing each time the new floating grid and the weights of all
    ! fixed+float grid in r2, doubling the computational load if done twice.
    ! Therefore, non-hermitian and hermitian part are no longer separated and
    ! are treated instead in the same loop.

    !PROVIDE tc_integ_type
    !print*, ' approach for integrals: ', tc_integ_type
    !PROVIDE int2_grad1_u12_square_ao

    int2_grad1_u12_vect_ao_at_r1(:,:,:) = 0.d0 
    int2_grad1_u12_sqrd_ao_at_r1(:,:) = 0.d0

!    if(tc_save_mem_loops) then
      print*, ' LOOPS are used to evaluate Hermitian part of ao_two_e_tc_tot_adapt ...'
      !!$OMP PARALLEL                                                    &
      !!$OMP DEFAULT (NONE)                                              &
      !!$OMP PRIVATE (i, j, k, l, ao_i_r1, ao_k_r1, ao_i_grad1, ao_k_grad1)  &
      !!$OMP SHARED  ( n_points_final_grid, i1, r1, w1, ao_num           &
      !!$OMP         , final_grid_points, final_weight_at_r_vector       &
      !!$OMP         , int2_grad1_u12_vect_ao_at_r1                      &
      !!$OMP         , int2_grad1_u12_sqrd_ao_at_r1                      &
      !!$OMP         , aos_in_r_array_transp                             &
      !!$OMP         , aos_grad_in_r_array_transp_bis                    &
      !!$OMP         , ao_two_e_tc_tot_adapt                             )
      do i1 = 1, n_points_final_grid
        r1(1:3) = final_grid_points(1:3,i1)
        w1 = final_weight_at_r_vector(i1)
        ! old function 
        !call get_int2_grad1_u12_ao_adapt(i1, int2_grad1_u12_vect_ao_at_r1, int2_grad1_u12_sqrd_ao_at_r1)
        ! new function need the floating grid to be created before
        call get_adaptive_grid(r1                                            &
                             , grid_points_extra_per_atom, grid_float_points &                      
                             , grid_fixed_weights, grid_float_weights        &         
                             , n_fixed_pts_effective, n_float_pts_effective  &
                             , n_pts_effective_max)
        !print*, "GET_ADAPTIVE_GRID: floating grid:"
        !do i = 1, size(grid_float_points, 3)-1
        !  do j = 1, size(grid_float_points, 2)
        !    print*, i, j, grid_float_points(1:3,j,i,1)
        !  end do
        !end do
        !call get_int2_grad1_u12_ao( i1, ao_num                   &
        call get_int2_grad1_u12_ao_full_adaptive( i1, ao_num                   &
                                  , n_points_ang_extra_grid      &
                                  , n_points_rad_extra_grid      & 
                                  , 1                            &
                                  , n_points_ang_float_grid      &
                                  , n_points_rad_float_grid      &
                                  , nucl_num                     &
                                  , grid_points_extra_per_atom   & 
                                  , grid_float_points            &
                                  , grid_fixed_weights           &
                                  , grid_float_weights           &
                                  , int2_grad1_u12_vect_ao_at_r1 &
                                  , int2_grad1_u12_sqrd_ao_at_r1)
        !!$OMP DO COLLAPSE(3)
        do i = 1, ao_num
          do k = 1, ao_num
            do l = 1, ao_num
              do j = 1, ao_num

                !print*, l, j, int2_grad1_u12_vect_ao_at_r1(j,l,:), int2_grad1_u12_sqrd_ao_at_r1(j,l)

                ao_two_e_tc_tot_adapt(j,l,k,i) = 0.d0
                ! AOs values in space
                ao_i_r1 = aos_in_r_array_transp(i1,i)
                ao_k_r1 = aos_in_r_array_transp(i1,k)
                ! AOs gradients in space
                ao_i_grad1(1:3) = aos_grad_in_r_array_transp_bis(i1,i,1:3)
                ao_k_grad1(1:3) = aos_grad_in_r_array_transp_bis(i1,k,1:3)
                ! HERMITIAN CONTRIBUTION
                !print*,  "Hermitian components"
                !print*, "w1 * int2_grad1_u12_sqrd_ao_at_r1(j,l) * ao_i_r1 * ao_k_r1"
                !print*, w1, int2_grad1_u12_sqrd_ao_at_r1(j,l), ao_i_r1, ao_k_r1
                ao_two_e_tc_tot_adapt(j,l,k,i) = ao_two_e_tc_tot_adapt(j,l,k,i) &
                    & + w1 * int2_grad1_u12_sqrd_ao_at_r1(j,l) * ao_i_r1 * ao_k_r1
                ! NON-HERMITIAN CONTRIBUTION
                !print*,  "NON Hermitian components"
                !print*, "w1, int2_grad1_u12_vect_ao_at_r1(j,l,1:3), ao_k_r1, ao_i_grad1(1:3), ao_i_r1, ao_k_grad1(1:3)"
                !print*, w1, int2_grad1_u12_vect_ao_at_r1(j,l,1:3), ao_k_r1, ao_i_grad1(1:3), ao_i_r1, ao_k_grad1(1:3)
                w1 = 0.5d0 * w1
                ao_two_e_tc_tot_adapt(j,l,k,i) = ao_two_e_tc_tot_adapt(j,l,k,i) &
                    & - w1 * sum(int2_grad1_u12_vect_ao_at_r1(j,l,1:3) * (ao_k_r1 * ao_i_grad1(1:3) - ao_i_r1 * ao_k_grad1(1:3) ) )
              enddo
            enddo
          enddo
        enddo
        !!$OMP END DO
      enddo
      !!$OMP END PARALLEL

!    else
!      print*, ' DGEMM are used to evaluate Hermitian part of ao_two_e_tc_tot_adapt ...'
!
!      int2_grad1_u12_vect_ao(:,:,:,:) = 0.d0 
!      int2_grad1_u12_sqrd_ao(:,:,:) = 0.d0
!
!      allocate(c_mat_grad_sqrd(n_points_final_grid,ao_num,ao_num))
!      allocate(c_mat_grad_vect(n_points_final_grid,ao_num,ao_num,3))
!
!      do i1 = 1, n_points_final_grid
!        r1(1:3) = final_grid_points(1:3,i1)
!        w1 = final_weight_at_r_vector(i1)
!        ! Get integrals in r2 for the current value of r1
!        call get_int2_grad1_u12_ao(i1, int2_grad1_u12_vect_ao_at_r1, int2_grad1_u12_sqrd_ao_at_r1)
!        ! Fill the array keeping track of the integrals for different values of r1
!        int2_grad1_u12_vect_ao(:,:,i1,:) = int2_grad1_u12_vect_ao_at_r1(:,:,:)
!        int2_grad1_u12_sqrd_ao(:,:,i1)   = int2_grad1_u12_sqrd_ao_at_r1(:,:)  
!        !$OMP PARALLEL               &
!        !$OMP DEFAULT (NONE)         &
!        !$OMP PRIVATE (i, k, ao_i_r1, ao_k_r1, ao_i_grad1, ao_k_grad1)  &
!        !$OMP SHARED  ( i1, r1, w1, n_points_final_grid  &
!        !$OMP         , final_grid_points                &
!        !$OMP         , final_weight_at_r_vector         &
!        !$OMP         , int2_grad1_u12_vect_ao_at_r1     &
!        !$OMP         , int2_grad1_u12_sqrd_ao_at_r1     &
!        !$OMP         , aos_in_r_array_transp            &
!        !$OMP         , aos_grad_in_r_array_transp_bis   &
!        !$OMP         , ao_num                           &       
!        !$OMP         , c_mat_grad_vect, c_mat_grad_sqrd )
!        !$OMP DO
!        do i = 1, ao_num
!          do k = 1, ao_num
!            ! AOs values in space
!            ao_i_r1 = aos_in_r_array_transp(i1,i)
!            ao_k_r1 = aos_in_r_array_transp(i1,k)
!            ! AOs gradients in space
!            ao_i_grad1(1:3) = aos_grad_in_r_array_transp_bis(i1,i,1:3)
!            ao_k_grad1(1:3) = aos_grad_in_r_array_transp_bis(i1,k,1:3)
!            ! \int dr_1 (\chi_k(1)\grad\chi_i(1) - \chi_i(1)\grad\chi_k(1) )
!            c_mat_grad_vect(i1,k,i,1:3) = 0.5 * w1 * (ao_k_r1*ao_i_grad1(1:3) - ao_i_r1*ao_k_grad1(1:3))
!            ! \int dr_1 \chi_k(1)\chi_i(1)
!            c_mat_grad_sqrd(i1,k,i) = w1 * ao_i_r1 * ao_k_r1
!          enddo
!        enddo
!        !$OMP END DO
!        !$OMP END PARALLEL
!      enddo
!
!      ! HERMITIAN CONTRIBUTION
!      call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid &
!                , 1.d0, int2_grad1_u12_square_ao(1,1,1), ao_num*ao_num        &
!                , c_mat_grad_sqrd(1,1,1), n_points_final_grid                 &
!                , 0.d0, ao_two_e_tc_tot_adapt(1,1,1,1), ao_num*ao_num)
!
!      ! NON-HERMITIAN CONTRIBUTION
!      call dgemm( "N", "N", ao_num*ao_num, ao_num*ao_num, n_points_final_grid &
!                , -1.d0, int2_grad1_u12_ao(1,1,1,1), ao_num*ao_num            &
!                , c_mat_grad_vect(1,1,1,1), n_points_final_grid               &
!                , 1.d0, ao_two_e_tc_tot_adapt(1,1,1,1), ao_num*ao_num)
!
!      deallocate(c_mat_grad_sqrd)
!      deallocate(c_mat_grad_vect)
!
!    endif !tc_save_mem_loops
 
    call wall_time(time1)
    print*, 'TC two electrons K integrals computed in (min) ', (time1 - time0) / 60.d0
    call print_memory_usage()

    ! Symmetrize result
    call sum_A_At(ao_two_e_tc_tot_adapt(1,1,1,1), ao_num*ao_num)

    ! ELECTRON REPULSION INTEGRAL
    print*, ' adding ERI to ao_two_e_tc_tot_adapt ...'

    if(tc_save_mem) then
      print*, ' ao_integrals_map will not be used'
      !$OMP PARALLEL DEFAULT(NONE)                     &
      !$OMP PRIVATE(i, j, k, l) & 
      !$OMP SHARED(ao_num, ao_two_e_tc_tot_adapt)
      !$OMP DO COLLAPSE(3)
      do j = 1, ao_num
        do l = 1, ao_num
          do i = 1, ao_num
            do k = 1, ao_num
              if(.not. ao_two_e_integral_zero(i,j,k,l)) then
                ! Orbitals i,k : r1
                ! Orbitals j,l : r2     
                ao_two_e_tc_tot_adapt(k,i,l,j) = ao_two_e_tc_tot_adapt(k,i,l,j) + ao_two_e_integral(i,k,j,l)
              endif
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    else
      print*,'Cholesky vectors will be used '
      double precision :: get_ao_integ_chol, eri
      eri = get_ao_integ_chol(1,1,1,1) ! FOR OPENMP 
      !$OMP PARALLEL DEFAULT(NONE)                            &
      !$OMP SHARED(ao_num, ao_two_e_tc_tot_adapt) &
      !$OMP PRIVATE(i, j, k, l,eri)
      !$OMP DO COLLAPSE(3)
      do j = 1, ao_num
        do l = 1, ao_num
          do i = 1, ao_num
            do k = 1, ao_num
              eri = get_ao_integ_chol(i,k,j,l)
              ao_two_e_tc_tot_adapt(k,i,l,j) = ao_two_e_tc_tot_adapt(k,i,l,j) + eri
            enddo
          enddo
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    endif

  endif ! read_tc_integ


  ! Not sure if it necessary to deallocate everything

  deallocate(int2_grad1_u12_vect_ao_at_r1)
  deallocate(int2_grad1_u12_sqrd_ao_at_r1)

  deallocate(int2_grad1_u12_vect_ao)
  deallocate(int2_grad1_u12_sqrd_ao)

  deallocate(grid_float_points)
  deallocate(grid_fixed_weights)
  deallocate(grid_float_weights)


  if(write_tc_integ .and. mpi_master) then
    print*, ' Saving ao_two_e_tc_tot_adapt in ', trim(ezfio_filename) // '/work/ao_two_e_tc_tot_adapt'
    open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/ao_two_e_tc_tot_adapt', action="write")
    call ezfio_set_work_empty(.False.)
    do i = 1, ao_num
      write(11) ao_two_e_tc_tot_adapt(:,:,:,i)
    enddo
    close(11)
    call ezfio_set_tc_keywords_io_tc_integ('Read')
  endif

  call wall_time(time1)
  print*, ' Wall time for ao_two_e_tc_tot_adapt (min) = ', (time1 - time0) / 60.d0
  call print_memory_usage()

END_PROVIDER 
