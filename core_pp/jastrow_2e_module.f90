module jastrow_module
  implicit none

  include 'constants.include.F'

  contains


pure function j12_mu(r1, r2, mu)
  !BEGIN_DOC
  ! $$  j(r1,r2;mu) 
  !   = u(r12;mu)
  !   = 1/2 r12 (1 - erf(mu r12)) - 1/2 (sqrt(pi) * mu) e^{-(mu*r12)^2}$$
  !END_DOC
  implicit none
  double precision, intent(in) :: r1(3)
  double precision, intent(in) :: r2(3)
  double precision, intent(in) :: mu
  !
  double precision :: j12_mu
  !
  double precision :: r12

  r12 = norm2(r1(1:3) - r2(1:3))
  j12_mu = u12_mu(r12,mu)
  !j12_mu = 0.5d0 * r12 * (1.d0 - derf(mu*r12)) - inv_sq_pi_2 * dexp(-mu*r12*mu*r12) / mu
end function j12_mu


pure function u12_mu(x, mu)
  !BEGIN_DOC
  ! Originally from: pot_j_gauss.irp.f:13:
  ! double precision function j_mu_simple(x,mu)
  ! $$ u_{\mu}(x) = 0.5  (1 - \erf(mu x)) - 1/[2 sqrt(pi)mu] exp(-(x*mu)^2) $$
  !END_DOC
  implicit none
  double precision, intent(in):: x
  double precision, intent(in):: mu
  ! OUTPUT
  double precision :: u12_mu
  u12_mu = 0.5d0 * x * (1.d0 - derf(mu*x)) - 0.5d0 * inv_sq_pi/mu *  dexp(-x*mu*x*mu)
end function


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

