!#####################################################################
! Catani-Seymour P and K kernels for a two-coloured-parton (q q~) Born, where
! T_a.T_b/T_a^2 = -1. Transcribed from DYN's PK_functions_mod.
!
! Every kernel is written k(z) = k0(z) + L*k1(z), L = log(muF^2/shat), and
! combines the P, Kbar and Ktilde pieces. The plus prescription is applied by the
! caller to the whole of qq_plus.
!#####################################################################
module dy_kp_kernels_mod

  use iso_fortran_env, only : real64
  implicit none
  private
  public :: qq_plus, qq_regular, qq_delta, gq_regular

  real(real64), parameter :: PI_VALUE = 3.141592653589793238462643d0
  real(real64), parameter :: CF = 4d0/3d0, TR = 0.5d0

contains

  subroutine qq_plus(z, k0, k1)
  real(real64), intent(in)  :: z
  real(real64), intent(out) :: k0, k1
  k0 = CF*2d0/(1d0 - z)*(log((1d0 - z)/z) + log(1d0 - z))
  k1 = -CF*(1d0 + z*z)/(1d0 - z)
  end subroutine qq_plus

  subroutine qq_regular(z, k0, k1)
  real(real64), intent(in)  :: z
  real(real64), intent(out) :: k0, k1
  k0 = CF*(-(1d0 + z)*log((1d0 - z)/z) + (1d0 - z)) - CF*(1d0 + z)*log(1d0 - z)
  k1 = 0d0
  end subroutine qq_regular

! The delta term has no muF dependence.
  real(real64) function qq_delta() result(k0)
  k0 = -CF*(5d0 - PI_VALUE*PI_VALUE) - CF*PI_VALUE*PI_VALUE/3d0
  end function qq_delta

  subroutine gq_regular(z, k0, k1)
  real(real64), intent(in)  :: z
  real(real64), intent(out) :: k0, k1
  real(real64) :: splitting
  splitting = TR*(z*z + (1d0 - z)**2)
  k0 = splitting*log((1d0 - z)/z) + 2d0*TR*z*(1d0 - z) + splitting*log(1d0 - z)
  k1 = -splitting
  end subroutine gq_regular

end module dy_kp_kernels_mod
