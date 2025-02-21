BEGIN_PROVIDER [ double precision, overlap_cart, (ao_cart_to_sphe_num,ao_cart_to_sphe_num  )]
 implicit none
 integer :: m,n,k,l
 overlap_cart = 0.d0
!do m = 1,ao_cart_to_sphe_num
! do n = 1,ao_cart_to_sphe_num

!  do k = 1, ao_num
!   do l = 1, ao_num
!     overlap_cart(m,n) += ao_cart_to_sphe_coef_2(l,m)  * ao_overlap(l,k) * ao_cart_to_sphe_coef_2(k,n)
!   enddo
!  enddo

! enddo
!enddo

 call ao_cart_to_ao_sphe(ao_overlap,ao_num,overlap_cart,ao_cart_to_sphe_num)
 double precision :: accu
 accu = 0.d0
 do m = 1,ao_cart_to_sphe_num
  do n = 1,ao_cart_to_sphe_num
    accu += dabs(overlap_cart(m,n) - ao_cart_to_sphe_overlap(m,n))
    print*,overlap_cart(m,n) , ao_cart_to_sphe_overlap(m,n),dabs(overlap_cart(m,n) - ao_cart_to_sphe_overlap(m,n))
  enddo
 enddo
 print*,'ao_cart_to_sphe_num = ',ao_cart_to_sphe_num
 print*,'accu = ',accu
END_PROVIDER 

BEGIN_PROVIDER [ double precision, mo_overlap_sphe, (mo_num, mo_num)]
 implicit none
 integer :: i,j,mu,nu
 mo_overlap_sphe = 0.d0
 do i = 1, mo_num
  do j = 1, mo_num
   do mu = 1, ao_cart_to_sphe_num
    do nu = 1, ao_cart_to_sphe_num
     mo_overlap_sphe(j,i) += mo_coef_ao_sphe(nu,i) * mo_coef_ao_sphe(mu,j) * ao_cart_to_sphe_overlap(mu,nu)
    enddo
   enddo
  enddo
 enddo
 double precision :: accu
 accu = 0.d0
 print*,'MO OVERLAP'
 do i = 1, mo_num
  write(*,'(100(F16.10,X))') mo_overlap_sphe(:,i)
  do j = 1, mo_num
   accu += dabs( mo_overlap_sphe(j,i) - mo_overlap(j,i))
  enddo
 enddo
 print*,'accu MO = ',accu
END_PROVIDER 


BEGIN_PROVIDER [ double precision, SsINV_CsT_prov, (ao_cart_to_sphe_num, ao_num)]
 implicit none
 call get_AB_prod(ao_cart_to_sphe_overlap_inv,ao_cart_to_sphe_num,ao_cart_to_sphe_num, &
        ao_cart_to_sphe_coef_transp, ao_num,SsINV_CsT_prov)
END_PROVIDER


BEGIN_PROVIDER [ double precision, Sc_Cc_prov, (ao_num, mo_num)]
  implicit none
  call get_AB_prod(ao_overlap,ao_num,ao_num,mo_coef,mo_num,Sc_Cc_prov)
END_PROVIDER 
