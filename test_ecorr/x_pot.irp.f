double precision function exchange_pot_kernel(r,rp)
 implicit none
 double precision,intent(in) :: r(3),rp(3)
 BEGIN_DOC
 ! returns the kernel of the Core exchange operator
 !
 ! = \sum_{i\in core_pseudo} phi_i(r) phi_i(rp)/|r-rp|
 END_DOC
 double precision :: distance
 distance = (r(1) - rp(1))*(r(1) - rp(1))
 distance+= (r(2) - rp(2))*(r(2) - rp(2))
 distance+= (r(3) - rp(3))*(r(3) - rp(3))
 distance = dsqrt(distance)
 integer :: i,ii
 double precision :: mos_array_r(mo_num),mos_array_rp(mo_num)
 call give_all_mos_at_r(r,mos_array_r)
 call give_all_mos_at_r(rp,mos_array_rp)
 exchange_pot_kernel = 0.D0
 do ii = 1, n_core_pseudo_orb
  i = list_core_pseudo(ii)
  exchange_pot_kernel += mos_array_rp(i) * mos_array_r(i)
 enddo
 if(distance.gt.1.d-10)then
  exchange_pot_kernel = - exchange_pot_kernel/distance 
 else
  exchange_pot_kernel = 0.d0
 endif

end

BEGIN_PROVIDER [ double precision, x_pot_mo_prov, (mo_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! <phi_i | V_x^{core_pseudo} | phi_j >
 END_DOC
 integer :: i,j,ii,k
 double precision :: get_two_e_integral
 x_pot_mo_prov = 0.d0
 do k = 1, mo_num
  do j = 1, mo_num
   do ii = 1, n_core_pseudo_orb
    i = list_core_pseudo(ii)
    x_pot_mo_prov(j,k) -= get_two_e_integral(i,j,k,i,mo_integrals_map) ! <ij|ki>
   enddo
  enddo
 enddo

END_PROVIDER 
