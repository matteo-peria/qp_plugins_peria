! Couldn't put this stuff in a standalone module because it depends on too many
! providers that would need to be passed as argument to the functions
!
! This file contains some copies with alternative names of important functions
! defined in the plugin non_h_ints_mu, which we can't put in NEED to avoid
! circular dependencies

!double precision function envelope_nucl(r, env_type)
! now env_type is a provider so no need to be passed as argument
double precision function envelope_nucl(r)
  ! TEMPORARY TESTING FUNCTION
  ! duplicate of:
  !   double precision function env_nucl(r) 
  ! from:
  !   jast_deriv_utils.irp.f:197
  implicit none
  double precision, intent(in) :: r(3)
  !character(:), allocatable :: env_type
  !
  integer :: i
  double precision :: a, d, e

  envelope_nucl = 1.d0

  select case(trim(env_type))

  case("Sum_Slat")
    do i = 1, nucl_num
      a = env_expo(i)
      d = sum( (r(:)-nucl_coord(i,:))*(r(:)-nucl_coord(i,:)) )
      envelope_nucl = envelope_nucl - env_coef(i) * exp(-a*dsqrt(d))
    enddo

  case("Prod_Gauss")
    do i = 1, nucl_num
      a = env_expo(i)
      d = sum( (r(:)-nucl_coord(i,:))*(r(:)-nucl_coord(i,:)) )
      e = 1.d0 - exp(-a*d)
      envelope_nucl = envelope_nucl * e
    enddo

  case("Sum_Gauss")
    do i = 1, nucl_num
      a = env_expo(i)
      d = sum( (r(:)-nucl_coord(i,:))*(r(:)-nucl_coord(i,:)) )
      envelope_nucl = envelope_nucl - env_coef(i) * exp(-a*d)
    enddo

  case("Sum_Quartic")
    do i = 1, nucl_num
      a = env_expo(i)
      d = sum( (r(:)-nucl_coord(i,:))*(r(:)-nucl_coord(i,:)) )
      envelope_nucl = envelope_nucl - env_coef(i) * exp(-a*d*d)
    enddo

  case default 
    print *, ' Error in envelope_nucl: Unknown env_type = ', env_type
    stop

  end select

end function

! can't have a pure function in irp.f files
!pure function j_mu_env(r1,r2,mu,env_type)
double precision function j_mu_env(r1,r2,mu)
  use jastrow_2e_module, only : j12_mu
  implicit none
  ! INPUT
  double precision, intent(in) :: r1(3)
  double precision, intent(in) :: r2(3)
  double precision, intent(in) :: mu
  !character(:), allocatable :: env_type
  double precision, external :: envelope_nucl
  !! OUTPUT
  !double precision :: j_mu_env

  j_mu_env = envelope_nucl(r1) * envelope_nucl(r2) * j12_mu(r1,r2,mu)
end function
