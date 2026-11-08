!#####################################################################
! Initial-initial Catani-Seymour dipoles for p p -> e+ e-.
!
! A two-coloured-parton Born has T1.T2 = -CF, so no colour-correlation call is
! needed and the colour factor cancels against 1/T_ai^2. g -> q q~ splittings
! into a quark Born carry no spin correlation.
!
! A q~ q Born is the q q~ amplitude with beam legs 1 and 2 exchanged.
!#####################################################################
module dipole_mod

  use iso_fortran_env, only : real64
  use olp_mod,         only : olp_born, NFLAVOUR
  implicit none
  private
  public :: dipoles_evaluate

  real(real64), parameter :: PI_VALUE = 3.141592653589793238462643d0
  real(real64), parameter :: CF = 4d0/3d0, TR = 0.5d0

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! momenta      the real kinematics
! ptilde       reduced kinematics per emitter leg, from ii_dipole_map
! dipole_qq(f, leg)   q -> q g, on the q q~ Born, for each emitter leg
! dipole_qqbar_born(f, leg)   g -> q q~, gluon on beam `leg`, reduced Born q q~
! dipole_qbarq_born(f, leg)   g -> q q~, gluon on beam `leg`, reduced Born q~ q
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine dipoles_evaluate(alphas, momenta, ptilde, dipole_qq,             &
                              dipole_qqbar_born, dipole_qbarq_born)
  real(real64), intent(in)  :: alphas, momenta(0:3,5), ptilde(0:3,4,2)
  real(real64), intent(out) :: dipole_qq(NFLAVOUR,2)
  real(real64), intent(out) :: dipole_qqbar_born(NFLAVOUR,2)
  real(real64), intent(out) :: dipole_qbarq_born(NFLAVOUR,2)

  real(real64) :: dot12, dot_emission(2), x_ab, kernel_q, kernel_g, prefactor
  real(real64) :: born_qqbar, born_qbarq, swapped(0:3,4), born_distinct(2,2)
  integer :: flavour, leg

  dot12 = dot(momenta(:,1), momenta(:,2))
  dot_emission(1) = dot(momenta(:,1), momenta(:,5))
  dot_emission(2) = dot(momenta(:,2), momenta(:,5))
  x_ab = (dot12 - dot_emission(1) - dot_emission(2))/dot12

  kernel_q = CF*(2d0/(1d0 - x_ab) - (1d0 + x_ab))
  kernel_g = TR*(1d0 - 2d0*x_ab*(1d0 - x_ab))

  do leg = 1, 2
     prefactor = 8d0*PI_VALUE*alphas/(2d0*x_ab*dot_emission(leg))
     swapped(:,1)   = ptilde(:,2,leg)
     swapped(:,2)   = ptilde(:,1,leg)
     swapped(:,3:4) = ptilde(:,3:4,leg)
! With massless quarks s, b equal d and c equals u, so the Borns of the two
! distinct flavours are reused for the rest.
     do flavour = 1, NFLAVOUR
        if (flavour <= 2) then
           call olp_born(flavour, ptilde(:,:,leg), born_distinct(flavour,1))
           call olp_born(flavour, swapped,         born_distinct(flavour,2))
        endif
        born_qqbar = born_distinct(mod(flavour-1,2)+1, 1)
        born_qbarq = born_distinct(mod(flavour-1,2)+1, 2)
        dipole_qq(flavour,leg)         = prefactor*kernel_q*born_qqbar
        dipole_qqbar_born(flavour,leg) = prefactor*kernel_g*born_qqbar
        dipole_qbarq_born(flavour,leg) = prefactor*kernel_g*born_qbarq
     end do
  end do
  end subroutine dipoles_evaluate

  real(real64) function dot(p, q) result(d)
  real(real64), intent(in) :: p(0:3), q(0:3)
  d = p(0)*q(0) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
  end function dot

end module dipole_mod
