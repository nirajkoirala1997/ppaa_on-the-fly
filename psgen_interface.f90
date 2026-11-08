!#################################################################################!
! Thin wrapper over PSGen (Kolodziej), the generated 'u u~ -> e+ e-' kinematics.
!
! PSGen is a phase-space MAPPING, not an integrator: psgen() never draws a random
! number, it consumes x(1:ndim) handed to it and returns the momenta, the initial
! flux factor and the phase-space normalization. That is the whole reason it is
! used here -- VEGAS (vsup.f) keeps the integration loop and keeps seeing every
! integration variable, phase-space or not.
!
! PSGen's own multichannel integrator, carlos.f, is deliberately not linked. At
! this multiplicity nkin = 1 and it would have nothing to adapt anyway.
!
! Everything this module exists to absorb is a sharp edge in psgen.f itself; each
! one is marked TRAP below.
!#################################################################################!
module psgen_interface

  use inprms_ps, only : dp, psnac, im_test, ip_test, dm_test, dp_test
  use kinparams, only : nexp, nfsp, nkin, m, aw, xmc, iwkin, aw_w
  use kincuts,   only : icuts
  use drivekk,   only : ikin, icall
  use psnrmtest, only : ipsnrm, ik_sch, ik_tch, ik_gch
  use ntdisc,    only : ntd0, ntd1, ntd2
  use fourmom,   only : p

  implicit none
  private

  public :: psgen_init, psgen_event
  public :: psgen_ndim, psgen_nexp, psgen_nkin

! Number of random numbers the final-state kinematics consumes, and the number
! of external legs. Both compile-time in PSGen: ndim = 2 and nexp = 4 here.
  integer, parameter :: psgen_ndim = 3*nfsp - 4
  integer, parameter :: psgen_nexp = nexp
  integer, parameter :: psgen_nkin = nkin

! TRAP: psgen assigns `flux` only when the partonic energy differs from the one
! it saw last (psgen.f:48, ecm0 is saved), so on a repeat it leaves the caller's
! variable untouched. Upstream cs_sect.f copes with a `save flux`; same here.
  real(dp) :: flux_saved = 0.0_dp

contains

!=================================================================!
! Call once, before the first psgen_event. `unit` receives PSGen's banner,
! the cut summary and the mass table.
  subroutine psgen_init(unit)
  implicit none
  integer, intent(in) :: unit
  integer :: i
  external :: param_trans

! No cuts inside PSGen. Two reasons. The mass window is closed by the tau
! sampling instead, so there is nothing left to cut on here; and PSGen's
! define_cuts.f ships with ictlb = 1 (5 degree lepton-beam) and ienl = 1
! (E_lep > 5 GeV) ENABLED, which would quietly bias the comparison against the
! KALEU reference. Set before param_trans -- define_cuts branches on icuts.
    icuts = 0

! Phase-space normalization cross-check and the s/t/g channel filters: off.
    ipsnrm = 0
    ik_sch = 0
    ik_tch = 0
    ik_gch = 0
    psnac  = 1.0d-3

! Per-event mass and momentum-conservation tests. Cheap to switch on when the
! kinematics are under suspicion, noisy otherwise.
    im_test = 0
    ip_test = 0
    dm_test = 1.0d-5
    dp_test = 1.0d-11

! Discarded-call counters and the channel index.
    icall = 0
    ntd0  = 0
    ntd1  = 0
    ntd2  = 0
    ikin  = 1

! Masses, particle names, cuts, particle indices, banner.
    call param_trans(unit)

! Massless external legs. The generated kinparams.f sets m = (mu0, mu0, me, me),
! i.e. 5 MeV quarks and a 0.511 MeV electron, but OpenLoops is registered with
! all of them massless -- momenta must sit on the mass shell it expects.
! Overridden here rather than in PSGen/kinparams.f, which is generated code and
! would lose the edit on regeneration.
    m = 0.0_dp

! param_trans printed the generated masses above, so say plainly that they no
! longer hold -- otherwise that table is the most misleading thing in the log.
    write(unit,'(a)') '#-------------------------------------------------'//   &
                      '----------------------------#'
    write(unit,'(a)') '# All external masses overridden to ZERO after the '//  &
                      'table above:'
    write(unit,'(a)') '# Recola is registered with massless quarks and '//  &
                      'leptons.'
    write(unit,'(a)') '#-------------------------------------------------'//   &
                      '----------------------------#'

! Multichannel weights. With nkin = 1 this is just aw(1) = 1, but it is NOT
! optional: psncls.f builds the phase-space weight as dlinv = sum_j aw(j)/dlips_j,
! so an uninitialised aw makes every weight garbage.
    aw     = 1.0_dp/real(nkin, dp)
    iwkin  = 0
    aw_w   = 0.0_dp
    xmc(0) = 0.0_dp
    do i = 1, nkin
       xmc(i) = xmc(i-1) + aw(i)
    end do

  end subroutine psgen_init

!=================================================================!
! One phase-space point.
!
!   ecm_had  hadronic CMS energy
!   xa, xb   beam momentum fractions (this module does not generate them --
!            PSGen takes them as input, which is why the tau/y sampling lives
!            in the integrand)
!   xr       psgen_ndim random numbers in (0,1)
!
!   plab     the nexp momenta, boosted to the LAB frame
!   flux     initial-state flux factor, 2*shat for massless beams
!   dlips    phase-space normalization; 0 means the point is unusable
  subroutine psgen_event(ecm_had, xa, xb, channel, xr, plab, flux, dlips)
  implicit none
  real(dp), intent(in)  :: ecm_had, xa, xb
  integer,  intent(in)  :: channel
  real(dp), intent(in)  :: xr(psgen_ndim)
  real(dp), intent(out) :: plab(0:3, psgen_nexp)
  real(dp), intent(out) :: flux, dlips

  real(dp) :: pcm(0:3, psgen_nexp)
  real(dp) :: vz, gz
  integer  :: i
  external :: psgen

! TRAP: when kincls fails, psgen skips its entire body (psgen.f:60) and leaves
! `dlips` with whatever the caller had in it. Zero it first, always.
    dlips = 0.0_dp
    flux  = 0.0_dp
    plab  = 0.0_dp

    ikin = channel
    call psgen(channel, ecm_had, xa, xb, xr, psgen_ndim, flux, dlips)

! The saved-flux workaround for the trap noted at the top of the module.
    if (flux > 0.0_dp) flux_saved = flux
    flux = flux_saved

    if (dlips <= 0.0_dp) return

! TRAP: fourmom::p is a global array that the next call overwrites. Copy out now.
!
! It also holds momenta in the PARTONIC CM frame. Boost to the lab with the same
! transform psgen uses internally when it imposes cuts (psgen.f:64-73), because
! everything downstream -- the dipole reduction, the PK terms, the ISR collinear
! structure -- is written for lab-frame momenta.
    pcm = p

    vz = (xa - xb)/(xa + xb)
    gz = 0.5_dp*(xa + xb)/sqrt(xa*xb)

    do i = 1, psgen_nexp
       plab(0,i) = gz*(pcm(0,i) + vz*pcm(3,i))
       plab(1,i) = pcm(1,i)
       plab(2,i) = pcm(2,i)
       plab(3,i) = gz*(pcm(3,i) + vz*pcm(0,i))
    end do

  end subroutine psgen_event

end module psgen_interface
!#################################################################################!
