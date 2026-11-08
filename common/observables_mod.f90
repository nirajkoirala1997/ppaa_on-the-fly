!#####################################################################
! Kinematic observables of a single four-vector p(0:3) in the
! (E,px,py,pz) convention used throughout this codebase.
!
! Everything here is a plain function of one four-vector and nothing
! else -- no dispatch, no observable "kinds", no leg lists. A composite
! object is built by the caller, which adds the four-vectors it wants
! and passes the sum:
!
!   m(e+e-)     obs_mass(p(:,3) + p(:,4))
!   pT(ll)      obs_pt  (p(:,3) + p(:,4))
!   MET         obs_pt  (p(:,5) + p(:,6))     the neutrino-pair pT
!   m4l         obs_mass(sum(p(:,3:6), dim=2))
!
! so an invariant mass needs no function of its own: it is obs_mass of
! the summed legs. The card name -> observable mapping lives in
! analysis_mod, next to the rest of the run-card plumbing.
!
! Legs are counted over ALL particles, the two incoming partons being
! 1 and 2, so leg 3 is the first final-state particle.
!#####################################################################

module observables_mod

  use iso_fortran_env, only: real64
  implicit none
  private
  public :: obs_pt, obs_rapidity, obs_pseudorapidity, obs_mass
  public :: obs_delta_phi, obs_delta_r, obs_mt

  real(real64), parameter :: TWOPI = 6.283185307179586476925287d0
  real(real64), parameter :: PIVAL = 3.141592653589793238462643d0

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Lorentz scalar product. Private: obs_mass is its only user in the
! tree, so it lives with the formula that needs it rather than in a
! shared bag of helpers. cuts_mod deliberately keeps its own copy of
! this algebra -- see the note above its private routines.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function dot(p, q) result(d)
  real(real64), intent(in) :: p(0:3), q(0:3)
  d = p(0)*q(0) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
  end function dot

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Transverse momentum.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_pt(p) result(val)
  real(real64), intent(in) :: p(0:3)
  val = sqrt(p(1)**2 + p(2)**2)
  end function obs_pt

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Mass, sqrt(p.p). Also the invariant mass of a system, when p is the
! sum of that system's four-vectors. Clipped at zero: a massless leg
! can come back with p.p slightly negative from rounding.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_mass(p) result(val)
  real(real64), intent(in) :: p(0:3)
  val = sqrt(max(0d0, dot(p,p)))
  end function obs_mass

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Rapidity y = 0.5*log((E+pz)/(E-pz)).
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_rapidity(p) result(val)
  real(real64), intent(in) :: p(0:3)
  real(real64) :: num, den
  num = p(0) + p(3)
  den = p(0) - p(3)
  if (num <= 0d0 .or. den <= 0d0) then
     val = 0d0
     return
  endif
  val = 0.5d0*log(num/den)
  end function obs_rapidity

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Pseudorapidity eta = 0.5*log((|p|+pz)/(|p|-pz)).
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_pseudorapidity(p) result(val)
  real(real64), intent(in) :: p(0:3)
  real(real64) :: pabs, num, den
  pabs = sqrt(p(1)**2 + p(2)**2 + p(3)**2)
  num = pabs + p(3)
  den = pabs - p(3)
  if (num <= 0d0 .or. den <= 0d0) then
     val = 0d0
     return
  endif
  val = 0.5d0*log(num/den)
  end function obs_pseudorapidity

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Azimuthal separation, wrapped to [0,pi].
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_delta_phi(p, q) result(val)
  real(real64), intent(in) :: p(0:3), q(0:3)
  val = atan2(p(2), p(1)) - atan2(q(2), q(1))
  do while (val >   PIVAL) ; val = val - TWOPI ; end do
  do while (val <= -PIVAL) ; val = val + TWOPI ; end do
  val = abs(val)
  end function obs_delta_phi

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Delta R = sqrt(delta_eta^2 + delta_phi^2).
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_delta_r(p, q) result(val)
  real(real64), intent(in) :: p(0:3), q(0:3)
  real(real64) :: delta_eta
  delta_eta = obs_pseudorapidity(p) - obs_pseudorapidity(q)
  val = sqrt(delta_eta**2 + obs_delta_phi(p, q)**2)
  end function obs_delta_r

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Massless transverse mass of two objects,
!   mT^2 = 2*pT(a)*pT(b)*(1 - cos(dphi))
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_mt(p, q) result(val)
  real(real64), intent(in) :: p(0:3), q(0:3)
  val = sqrt(max(0d0, 2d0*obs_pt(p)*obs_pt(q)*(1d0 - cos(obs_delta_phi(p,q)))))
  end function obs_mt

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module observables_mod

!#####################################################################
