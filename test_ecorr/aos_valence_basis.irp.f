! Construction of the primitives and of the basis of the cartesian AOs
! defining the valence-subspace
! ao_val_coef_normed_normed,

 BEGIN_PROVIDER [ double precision, ao_coef_transp, (ao_prim_num_max, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_expo_transp, (ao_prim_num_max, ao_num) ]
  ! Necessary for other providers
  implicit none
  integer :: row, col
  do row=1,ao_prim_num_max
    do col=1,ao_num
      ao_coef_transp(row,col) = ao_coef(col,row)
      ao_expo_transp(row,col) = ao_expo(col,row)
    end do
  end do
END_PROVIDER


 BEGIN_PROVIDER [ integer, ao_pp_prim_num_temp, (ao_num)]
  BEGIN_DOC
  ! temporary number of primitives defining each pseudo-potential AO 
  END_DOC
  implicit none
  integer, dimension(:), allocatable :: nonzero_indx
  integer :: col, n_nonzero

  ao_pp_prim_num_temp = 0
  do col = 1, ao_num
    if (abs(ao_val_overlap_eigval(col)) > 1.d-8) then
      print*, "ao_val_overlap_eigval(col) is DIFFERENT FROM ZERO: ", ao_val_overlap_eigval(col)
      ! nonzero_indx are the indices where the pp-|AO| coefficients are non-zero
      !nonzero_indx = get_nonzero_indices_1d(ao_val_coef_normed(:,col))
      ! EXPLICIT INTERFACE NEEDED
      n_nonzero = count(ao_val_coef_normed(:,col) /= 0.d0)
      allocate(nonzero_indx(n_nonzero))
      call get_nonzero_indices_1d(ao_val_coef_normed(:,col), &
                                  ao_num,                    &
                                  nonzero_indx,                   &
                                  n_nonzero)
      ! Sum of the number of primitives involved in defining the j-th pp-AO
      ao_pp_prim_num_temp(col) = sum(ao_prim_num(nonzero_indx))
      deallocate(nonzero_indx)
    else
      print*, "ao_val_overlap_eigval(col) is ZERO: ", ao_val_overlap_eigval(col)
      cycle
    end if
  end do
END_PROVIDER


 BEGIN_PROVIDER [ integer, ao_pp_prim_num_max_temp ]
  BEGIN_DOC
  ! Temporary maximum number of primitives employed in the definition of 
  ! the pseudo-potential |AO| basis
  END_DOC
  implicit none
  ao_pp_prim_num_max_temp = maxval(ao_pp_prim_num_temp)
  print*,'ao_pp_prim_num_max_temp = ',ao_pp_prim_num_max_temp
END_PROVIDER


 BEGIN_PROVIDER [ integer, ao_pp_prim_num_max ]
&BEGIN_PROVIDER [ integer, ao_pp_prim_num, (ao_num)]
&BEGIN_PROVIDER [ double precision, ao_pp_coef_transp_temp, (ao_pp_prim_num_max_temp, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_pp_expo_transp_temp, (ao_pp_prim_num_max_temp, ao_num) ]
  BEGIN_DOC
  ! Sorted primitives-to-pp-|AO|
  END_DOC
  implicit none
  integer, dimension(:), allocatable :: nonzero_indx
  double precision, dimension(:), allocatable :: nonzero_array
  integer, dimension(ao_pp_prim_num_max_temp) :: iorder
  double precision, dimension(ao_pp_prim_num_max_temp,2) :: sorter
  integer, dimension(ao_pp_prim_num_max_temp) :: partitioning
  integer :: i, row, col, n_nonzero

  ao_pp_prim_num(:) = 0
  ao_pp_prim_num_max = 0 
  ao_pp_expo_transp_temp(:,:) = 0.d0
  ao_pp_coef_transp_temp(:,:) = 0.d0
 
  ! Compute the new |AO| basis coefficients
  do col = 1, ao_num
    if (abs(ao_val_overlap_eigval(col)) > 1.d-8) then
      !nonzero_indx = get_nonzero_indices_1d(ao_val_coef_normed(:,col))
      ! EXPLICIT INTERFACE NEEDED
      n_nonzero = count(ao_val_coef_normed(:,col) /= 0.d0)
      allocate(nonzero_indx(n_nonzero))
      call get_nonzero_indices_1d(ao_val_coef_normed(:,col), &
                                  ao_num,                    &
                                  nonzero_indx,                   &
                                  n_nonzero)
      integer :: lambda, i_ao_pp_prim, i_ao_prim
      ! Loop over all the non-zero rows in the column vector ao_val_coef_normed(:,col)
      i_ao_pp_prim = 1
      do row = 1, size(nonzero_indx)
        ! row index of ao_val_coef_normed matrix
        lambda = nonzero_indx(row)
        ! Loop over the total number of primitives corresponding to row-lambda 
        do i_ao_prim = 1, ao_prim_num(lambda)
          ao_pp_expo_transp_temp(i_ao_pp_prim,col) = ao_expo_transp(i_ao_prim,lambda)
          ao_pp_coef_transp_temp(i_ao_pp_prim,col) = &
                        & ao_coef_transp(i_ao_prim,lambda)*ao_val_coef_normed(lambda,col)
          i_ao_pp_prim += 1
        end do
      end do
      ! Sort the exponents of all the primitives defining the col-th pp-AO
      ! and their coefficients accordingly
      sorter(:,1) = ao_pp_expo_transp_temp(:,col) !(1:ao_pp_prim_num_temp(col),col)
      sorter(:,2) = ao_pp_coef_transp_temp(:,col) !(1:ao_pp_prim_num_temp(col),col)
      iorder = [(i, i=1,ao_pp_prim_num_temp(col))]
      call dsort(sorter(1,1),iorder,ao_pp_prim_num_temp(col))
      call dset_order(sorter(1,2),iorder,ao_pp_prim_num_temp(col))
      ao_pp_expo_transp_temp(:,col) = sorter(:,1)
      ao_pp_coef_transp_temp(:,col) = sorter(:,2)
      ! Remove exponent repetition among all ao_pp_expo_transp_temp
      ! First locate nonzero exponents
      n_nonzero = count(ao_pp_expo_transp_temp(:,col) /= 0.d0)
      print*,'n_nonzero = ',n_nonzero
      integer :: iii
!     do iii = 1,size(ao_pp_expo_transp_temp,1)
!      print*,ao_pp_expo_transp_temp(iii,col)
!     enddo
!     pause
      deallocate(nonzero_indx)
      allocate(nonzero_indx(n_nonzero))
      call get_nonzero_indices_1d(ao_pp_expo_transp_temp(:,col), &
                                  ao_pp_prim_num_temp(col),              &
                                  nonzero_indx,                          &
                                  n_nonzero)
      allocate(nonzero_array(n_nonzero))
      nonzero_array(:) = 0.d0
      nonzero_array =  ao_pp_expo_transp_temp(nonzero_indx,col)
      partitioning = 0
      ! Sum the coefficients with same exponent in ao_pp_coef_transp_temp
      !partitioning = get_equivalence_partition_1d_dp(ao_pp_expo_transp_temp(:,col),tolerance=1d-7)
      ! EXPLICIT INTERFACE NEEDED
      call get_equivalence_partition_1d_dp(nonzero_array, & 
                                          n_nonzero,      & 
                                          1d-7,           & 
                                          partitioning)
      ! The total number of true different primitives is given by the total number
      ! of different non-zero partitions, which are numbered with increasing order.
      ! So the total number of true different primitives is the last partition number
      ao_pp_prim_num(col) = partitioning(n_nonzero)
      double precision :: avg_expo = 0.d0
      double precision :: sum_coef = 0.d0
      do row = 1, ao_pp_prim_num(col)
        ! Compute average exponent between those that are equal according to partitioning
        ! (these exponent are expected to actually be the same, 
        write(*,*) "avg_expo = ", sum(merge(ao_pp_expo_transp_temp(:,col), 0.d0, partitioning==row)),"/", count(partitioning==row)

        avg_expo = sum(merge(ao_pp_expo_transp_temp(:,col), 0.d0, partitioning==row))/count(partitioning==row)
        ! Set to zero all the equal exponents belonging to the same partition
        ao_pp_expo_transp_temp(:,col) = merge(0.d0, ao_pp_expo_transp_temp(:,col), partitioning==row)
        ! Set to the average exponent as the final exponent
        ao_pp_expo_transp_temp(row,col) = avg_expo
        ! Compute the new coefficient as the sum of coefficients of primitives having the same exponent
        sum_coef = sum(merge(ao_pp_coef_transp_temp(:,col), 0.d0, partitioning==row))
        ! Set to zero all the entries corresponding to same-exponent primitives
        ao_pp_coef_transp_temp(:,col) = merge(0.d0, ao_pp_coef_transp_temp(:,col), partitioning==row)
        ! And set the new coefficient  
        ao_pp_coef_transp_temp(row,col) = sum_coef
      end do
      deallocate(nonzero_indx)
      deallocate(nonzero_array)
    else
      cycle
   end if
  end do
  ao_pp_prim_num_max = maxval(ao_pp_prim_num)
END_PROVIDER


 BEGIN_PROVIDER [ double precision, ao_pp_coef_transp, (ao_pp_prim_num_max, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_pp_expo_transp, (ao_pp_prim_num_max, ao_num) ]
  implicit none
  ao_pp_coef_transp = ao_pp_coef_transp_temp(1:ao_pp_prim_num_max,1:ao_num)
  ao_pp_expo_transp = ao_pp_expo_transp_temp(1:ao_pp_prim_num_max,1:ao_num)
END_PROVIDER


 BEGIN_PROVIDER [ double precision, ao_pp_coef, (ao_num, ao_pp_prim_num_max) ]
&BEGIN_PROVIDER [ double precision, ao_pp_expo, (ao_num, ao_pp_prim_num_max) ]
  implicit none
  integer :: row, col
  do row = 1, ao_pp_prim_num_max
    do col = 1, ao_num
      ao_pp_coef(col,row) = ao_pp_coef_transp(row,col)
      ao_pp_expo(col,row) = ao_pp_expo_transp(row,col)
    end do
  end do
END_PROVIDER


 BEGIN_PROVIDER [ double precision, ao_pp_coef_normalized, (ao_num, ao_pp_prim_num_max) ]
&BEGIN_PROVIDER [ double precision, ao_pp_coef_normalization_factor, (ao_num) ]    
  implicit none                                                                 
  BEGIN_DOC                                                                     
  ! Coefficients including the |AO| normalization                               
  END_DOC                                                                       
  double precision               :: norm,overlap_x,overlap_y,overlap_z,C_A(3), c
  integer                        :: l, powA(3), nz                              
  integer                        :: i,j,k                                       
  !
  nz=100                                                                        
  C_A(1:3) = 0.d0                                                                 
  ao_pp_coef_normalized = 0.d0                                                     
                                                                                
  ao_pp_coef_normalized = ao_pp_coef

  do i=1,ao_num                                                                 

    if (abs(ao_val_overlap_eigval(i)) > 1.d-8) then
      powA(1:3) = ao_power(i,1:3)                                                     
      ! Normalization of the primitives                                           
      if (primitives_normalized) then
        do j=1,ao_pp_prim_num(i)                                                     
          call overlap_gaussian_xyz(C_A,C_A,ao_pp_expo(i,j),ao_pp_expo(i,j), &          
             powA,powA,overlap_x,overlap_y,overlap_z,norm,nz)                     
          write(*,*) "ao_pp_coef_normalized(i,j) = ", ao_pp_coef(i,j), "/", dsqrt(norm)                      
          ao_pp_coef_normalized(i,j) = ao_pp_coef(i,j)/dsqrt(norm)                      
        enddo                                                                     
      else                                                                        
        do j=1,ao_pp_prim_num(i)                                                     
          ao_pp_coef_normalized(i,j) = ao_pp_coef(i,j)                                  
        enddo                                                                     
      endif                                                                       
                                                                                  
      ! Normalization of the contracted basis functions                           
      norm = 0.d0                                                                 
      do j=1,ao_pp_prim_num(i)                                                       
        do k=1,ao_pp_prim_num(i)                                                     
          call overlap_gaussian_xyz(C_A,C_A,ao_pp_expo(i,j),ao_pp_expo(i,k),powA,powA,overlap_x,overlap_y,overlap_z,c,nz)
          norm = norm+c*ao_pp_coef_normalized(i,j)*ao_pp_coef_normalized(i,k)           
        enddo                                                                     
      enddo                                                                       
      write(*,*) "ao_pp_coef_normalization_factor(i) = ", 1.d0, "/", dsqrt(norm)                          
      ao_pp_coef_normalization_factor(i) = 1.d0/dsqrt(norm)                          
                                                                                  
      if (ao_normalized) then                                                     
        do j=1,ao_pp_prim_num(i)                                                     
          ao_pp_coef_normalized(i,j) = ao_pp_coef_normalized(i,j) * ao_pp_coef_normalization_factor(i)
        enddo                                                                     
      else                                                                        
        ao_pp_coef_normalization_factor(i) = 1.d0                                    
      endif                                                                       
    else 
      cycle
    end if
  enddo                                            
END_PROVIDER


 BEGIN_PROVIDER [ double precision, ao_pp_coef_normalized_transp, (ao_pp_prim_num_max, ao_num) ]
  implicit none
  integer :: row, col
  do row = 1,ao_num
    do col = 1,ao_pp_prim_num_max      
      ao_pp_coef_normalized_transp(col,row) = ao_pp_coef_normalized(row,col)
    end do
  end do
END_PROVIDER


BEGIN_PROVIDER [double precision, ao_pp_expo_transp_per_nucl, (ao_pp_prim_num_max,N_AOs_max,nucl_num) ]
 implicit none
 integer :: i,j,k,l
 do i = 1, nucl_num
  do j = 1,Nucl_N_Aos(i)
   k = Nucl_Aos_transposed(j,i)
   do l = 1, ao_pp_prim_num(k)
    ao_pp_expo_transp_per_nucl(l,j,i) = ao_pp_expo_transp(l,k)
   enddo
  enddo
 enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_pp_coef_normalized_transp_per_nucl, (ao_pp_prim_num_max,N_AOs_max,nucl_num) ]
 implicit none
 integer :: i,j,k,l
 do i = 1, nucl_num
  do j = 1,Nucl_N_Aos(i)
   k = Nucl_Aos_transposed(j,i)
   do l = 1, ao_pp_prim_num(k)
    ao_pp_coef_normalized_transp_per_nucl(l,j,i) = ao_pp_coef_normalized_transp(l,k)
   enddo
  enddo
 enddo
END_PROVIDER

