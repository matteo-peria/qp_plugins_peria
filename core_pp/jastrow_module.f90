module jastrow_module
  implicit none
  contains

pure function j12(r1, r2, mu)
  !BEGIN_DOC
  ! j(mu,r12) = 1/2 r12 (1 - erf(mu r12)) - 1/2 (sqrt(pi) * mu) e^{-(mu*r12)^2}
  !END_DOC
  implicit none
  include 'constants.include.F'
  double precision, intent(in) :: r1(3), r2(3), mu
  !
  double precision :: j12
  !
  double precision :: r12

  r12 = norm2(r1(1:3) - r2(1:3))
  j12 = 0.5d0 * r12 * (1.d0 - derf(mu*r12)) - inv_sq_pi_2 * dexp(-mu*r12*mu*r12) / mu
end function j12


!double precision function deriv_j12_r12(r12,mu) result(deriv)
pure function deriv_j12_r12(r1,r2,mu) result(deriv)
  !BEGIN_DOC
  ! d/dr12 j(mu,r12)
  !END_DOC
  ! INPUT
  double precision, intent(in) :: r1(3), r2(3), mu
  ! RETURN
  double precision :: deriv
  !
  double precision :: r12
  
  r12 = norm2(r1(1:3) - r2(1:3))
  deriv = 0.5d0 * (1.d0 - derf(mu * r12))
end function


!pure function grad1_j12(r1,r2,mu) result(grad)
function grad1_j12(r1,r2,mu) result(grad)
  !BEGIN_DOC
  ! grad j(mu,r12)
  !END_DOC
  implicit none
  ! INPUT
  double precision, intent(in) :: r1(3), r2(3), mu
  ! RETURN
  double precision :: grad(3)
  !
  double precision :: r12_vect(3)
  double precision :: r12
  double precision :: deriv
  !double precision, external :: deriv_r12_j12
 
  !print*, "IM IN GRAD_J12"
  grad(:) = 0.d0

  r12_vect = r1(1:3) - r2(1:3)
  r12 = norm2(r12_vect)
  if(r12 .lt. 1d-10) return
  deriv = deriv_j12_r12(r1,r2,mu)
  grad(1:3) = deriv * r12_vect(1:3)
  grad = grad/r12
end function


end module jastrow_module

