!#####################################################################
! Kinematic observables computed from raw four-vectors p(0:3) in the
! (E,px,py,pz) convention used throughout this codebase.
!
! These are the quantities a histogram is binned in and the quantities
! the cut is built from, so a new distribution is a one-line addition
! wherever events are processed.
!
! OBSERVABLES TAKE A LEG *SET*, NOT A LEG.
!   A four-lepton final state is mostly described by composite objects:
!   m(e+e-), the transverse momentum of the dilepton system, the missing
!   transverse momentum (which is nothing but the pT of the neutrino
!   pair). Binding an observable to a single leg cannot express any of
!   them, so obs_evaluate takes a list of legs, sums their four-vectors
!   and evaluates the observable on that sum. A one-element list
!   reproduces the per-leg case exactly, so `leg = 4` in a run card is
!   just `legs = 4` written shorter.
!
!   Legs are counted over ALL particles, the two incoming partons being
!   1 and 2, so leg 3 is the first final-state particle.
!#####################################################################

module observables_mod

  use iso_fortran_env, only: real64
  use misc_utils_mod, only: dot
  implicit none
  private
  public :: obs_pt, obs_rapidity, obs_pseudorapidity, obs_phi
  public :: obs_energy, obs_mass, obs_inv_mass, obs_delta_r, obs_delta_phi
  public :: OBSK_INV_MASS, OBSK_PT, OBSK_RAPIDITY, OBSK_PSEUDORAP, OBSK_ENERGY
  public :: OBSK_MASS, OBSK_DPHI, OBSK_DR, OBSK_MT
  public :: obs_kind_from_name, obs_nlegs_required, obs_evaluate

  ! Observable kinds, so a run card can name a distribution and the
  ! generator can dispatch on it without a compile-time list.
  integer, parameter :: OBSK_INV_MASS  = 1   ! whole final state, no legs
  integer, parameter :: OBSK_PT        = 2
  integer, parameter :: OBSK_RAPIDITY  = 3
  integer, parameter :: OBSK_PSEUDORAP = 4
  integer, parameter :: OBSK_ENERGY    = 5
  integer, parameter :: OBSK_MASS      = 6   ! mass of the summed leg set
  integer, parameter :: OBSK_DPHI      = 7   ! two objects (see the half-split rule)
  integer, parameter :: OBSK_DR        = 8   ! two objects
  integer, parameter :: OBSK_MT        = 9   ! two objects

  real(real64), parameter :: TWOPI = 6.283185307179586476925287d0
  real(real64), parameter :: PIVAL = 3.141592653589793238462643d0

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Map a card name to an observable kind. Returns 0 if unrecognised.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer function obs_kind_from_name(name) result(kind)
  character(len=*), intent(in) :: name
  select case (trim(adjustl(name)))
    case ('inv_mass','invmass')         ; kind = OBSK_INV_MASS
    case ('mass','m')                   ; kind = OBSK_MASS
    case ('pt','pT')                    ; kind = OBSK_PT
    case ('rapidity','y')               ; kind = OBSK_RAPIDITY
    case ('pseudorapidity','eta')       ; kind = OBSK_PSEUDORAP
    case ('energy','e')                 ; kind = OBSK_ENERGY
    case ('dphi','delta_phi')           ; kind = OBSK_DPHI
    case ('dr','delta_r')               ; kind = OBSK_DR
    case ('mt','transverse_mass')       ; kind = OBSK_MT
    case default                        ; kind = 0
  end select
  end function obs_kind_from_name

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! What this observable requires of the leg list:
!   0  no legs at all -- it applies to the final state as a whole
!  -1  one or more
!  -2  an even number, at least two (the half-split kinds below)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer function obs_nlegs_required(kind) result(nreq)
  integer, intent(in) :: kind
  select case (kind)
    case (OBSK_INV_MASS)              ; nreq =  0
    case (OBSK_DPHI, OBSK_DR, OBSK_MT); nreq = -2
    case default                      ; nreq = -1
  end select
  end function obs_nlegs_required

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Evaluate an observable on one event.
!
! ONE-OBJECT KINDS sum every listed leg and evaluate on the sum, so
!   legs = 4        is the pT of leg 4
!   legs = 3,4      is the pT of the leg-3 + leg-4 system
!   legs = 5,6      is, for a neutrino pair, the missing transverse
!                   momentum -- no separate MET observable is needed
!
! TWO-OBJECT KINDS (dphi, dr, mt) split the list in half, sum each half
! and evaluate between the two resulting objects. A separation is a
! property of a pair, not of a sum, so it cannot use the rule above; the
! half-split is what lets the same syntax express both
!   legs = 3,4      Delta phi between the two charged leptons, and
!   legs = 3,4,5,6  Delta phi between the dilepton and the neutrino pair.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_evaluate(kind, legs, nlegs, qmom, n) result(val)
  integer,      intent(in) :: kind, nlegs, n
  integer,      intent(in) :: legs(:)
  real(real64), intent(in) :: qmom(0:3,n)

  real(real64) :: psum(0:3), pa(0:3), pb(0:3)
  integer :: i, nh

  select case (kind)

    case (OBSK_INV_MASS)
       val = obs_inv_mass(qmom(:,3:n))

    case (OBSK_DPHI, OBSK_DR, OBSK_MT)
       nh = nlegs/2
       pa = 0d0
       pb = 0d0
       do i = 1, nh
          pa = pa + qmom(:,legs(i))
       end do
       do i = nh+1, nlegs
          pb = pb + qmom(:,legs(i))
       end do
       select case (kind)
         case (OBSK_DPHI) ; val = obs_delta_phi(pa, pb)
         case (OBSK_DR)   ; val = obs_delta_r(pa, pb)
         case default     ; val = obs_mt(pa, pb)
       end select

    case default
       psum = 0d0
       do i = 1, nlegs
          psum = psum + qmom(:,legs(i))
       end do
       select case (kind)
         case (OBSK_PT)        ; val = obs_pt(psum)
         case (OBSK_RAPIDITY)  ; val = obs_rapidity(psum)
         case (OBSK_PSEUDORAP) ; val = obs_pseudorapidity(psum)
         case (OBSK_ENERGY)    ; val = obs_energy(psum)
         case (OBSK_MASS)      ; val = obs_mass(psum)
         case default          ; val = 0d0
       end select

  end select

  end function obs_evaluate

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Transverse momentum.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_pt(p) result(val)
  real(real64), intent(in) :: p(0:3)
  val = sqrt(p(1)**2 + p(2)**2)
  end function obs_pt

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Energy.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_energy(p) result(val)
  real(real64), intent(in) :: p(0:3)
  val = p(0)
  end function obs_energy

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
! Azimuthal angle in (-pi,pi].
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_phi(p) result(val)
  real(real64), intent(in) :: p(0:3)
  val = atan2(p(2), p(1))
  end function obs_phi

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Mass of a single four-vector.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_mass(p) result(val)
  real(real64), intent(in) :: p(0:3)
  val = sqrt(max(0d0, dot(p,p)))
  end function obs_mass

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Invariant mass of a set of four-vectors, p(0:3,1:n).
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_inv_mass(p) result(val)
  real(real64), intent(in) :: p(:,:)          ! (0:3, n)
  real(real64) :: psum(0:3)
  integer :: i
  psum = 0d0
  do i = 1, size(p,2)
     psum = psum + p(:,i)
  end do
  val = sqrt(max(0d0, dot(psum,psum)))
  end function obs_inv_mass

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Azimuthal separation, wrapped to [0,pi].
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_delta_phi(p, q) result(val)
  real(real64), intent(in) :: p(0:3), q(0:3)
  val = obs_phi(p) - obs_phi(q)
  do while (val >   PIVAL) ; val = val - TWOPI ; end do
  do while (val <= -PIVAL) ; val = val + TWOPI ; end do
  val = abs(val)
  end function obs_delta_phi

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Delta R = sqrt(delta_eta^2 + delta_phi^2).
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_delta_r(p, q) result(val)
  real(real64), intent(in) :: p(0:3), q(0:3)
  real(real64) :: deta, dphi
  deta = obs_pseudorapidity(p) - obs_pseudorapidity(q)
  dphi = obs_delta_phi(p, q)
  val = sqrt(deta**2 + dphi**2)
  end function obs_delta_r

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Transverse mass of two objects,
!
!   mT^2 = 2*pT(a)*pT(b)*(1 - cos(dphi))
!
! the massless form, which is the one used for a system whose second
! object is reconstructed from missing transverse momentum and therefore
! has no longitudinal information.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_mt(p, q) result(val)
  real(real64), intent(in) :: p(0:3), q(0:3)
  val = sqrt(max(0d0, 2d0*obs_pt(p)*obs_pt(q)*(1d0 - cos(obs_delta_phi(p,q)))))
  end function obs_mt

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module observables_mod

!#####################################################################
