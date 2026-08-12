!#################################################################################!
! The leading-order integrand for  p p -> e+ e-.
!
! Four VEGAS dimensions:
!   yy(1), yy(2)   the beam momentum fractions x1, x2, thrown flat
!   yy(3), yy(4)   handed to PSGen as the final-state phase-space variables
!
! All four are ordinary VEGAS dimensions and all four get adapted. That is the
! point of driving PSGen rather than letting a generator own the loop: an extra
! non-phase-space variable -- a Mellin contour, say -- would simply be yy(5).
!
! NOTHING HERE ASSUMES TWO FINAL-STATE PARTICLES. x1 and x2 are handed over as
! they come, and the mass window is applied to m(e+e-) rebuilt from the
! generated momenta. Both statements stay true when a gluon is added: only the
! leg count and psgen_ndim change.
!#################################################################################!
  double precision function flo2_LO(yy, vwgt)

! xmur is deliberately absent: the tree amplitude carries no power of gs, so muR
! cannot enter it. Nothing here is allowed to depend on the renormalization scale.
  use globals_mod,         only : s, mee_min2, mee_max2, xmuf
  use constants_mod,       only : hbarc2
  use misc_utils_mod,      only : dot
  use matrix_elements_mod, only : evaluate_born
  use psgen_interface,     only : psgen_event, psgen_ndim

  implicit none

  real(8), intent(in) :: yy(10)
  real(8), intent(in) :: vwgt

  real(8) :: xa, xb
  real(8) :: xr(psgen_ndim), plab(0:3,4), flux, dlips
  real(8) :: q(0:3), mee2
  real(8) :: f1(-6:6), f2(-6:6), amp(5), sig
  integer :: i

    flo2_LO = 0d0

!--- Beam momentum fractions -----------------------------------------------!
! Flat, handed straight to PSGen, which takes x1 and x2 as arguments and never
! generates them. No change of variables, so no jacobian.
!
! Deliberately kept free of any importance sampling. A tau = x1*x2 mapping
! measures 25x better on the MC error at equal statistics, but the sampling is
! then something the integrand has to get right for each multiplicity, and
! keeping x1, x2 raw is what lets this integrand carry over to the real
! emission unchanged. VEGAS still adapts these two dimensions like any other.
    xa = yy(1)
    xb = yy(2)
    if (xa <= 0d0 .or. xb <= 0d0) return

!--- Phase space -----------------------------------------------------------!
    xr(1) = yy(3)
    xr(2) = yy(4)

    call psgen_event(sqrt(s), xa, xb, xr, plab, flux, dlips)
    if (dlips <= 0d0) return

!--- Mass window -----------------------------------------------------------!
! Rebuilt from the momenta, NOT from x1*x2*s. The two agree only at LO. Doing
! it this way is what lets the same line survive the addition of a gluon, where
! shat and m(ee)^2 part company.
    q    = plab(:,3) + plab(:,4)
    mee2 = dot(q, q)
    if (mee2 < mee_min2 .or. mee2 > mee_max2) return

!--- Parton luminosity and matrix element ----------------------------------!
    call pdf(xa, xmuf, f1)
    call pdf(xb, xmuf, f2)

    call evaluate_born(plab(:,1), plab(:,2), plab(:,3), plab(:,4), amp)

! Flavour sum over d, u, s, c, b. Both beam orderings share one amp(i): the
! amplitude is NOT symmetric under swapping the incoming quark and antiquark --
! that asymmetry is the forward-backward asymmetry -- but its antisymmetric part
! integrates to zero over the full angular range, so for a total with no angular
! cut this is exact. It stops being exact the moment an eta or cos(theta) cut
! is introduced; at that point the two orderings need separate amplitudes.
    sig = 0d0
    do i = 1, 5
       sig = sig + (f1(-i)*f2(i) + f1(i)*f2(-i))*amp(i)
    end do

!--- Weight ----------------------------------------------------------------!
! dsigma = hbarc2/flux * |M|^2 * dLIPS, with flux = 2*shat from PSGen and
! dLIPS = cps/dlinv = 1/(8 pi) for a massless two-body final state. That is
! hbarc2/(16 pi shat) -- identical to the xnorm this code used with the
! hand-written kinvar2 kinematics, which is what makes the switch checkable.
    flo2_LO = hbarc2/flux*sig*dlips

  end function flo2_LO
!#################################################################################!
