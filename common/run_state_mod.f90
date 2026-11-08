!#####################################################################
! Run-wide LIVE STATE for runLO -- what main has decided, in the form the
! integrand consumes it. Kept deliberately small: every number the card
! supplies -- masses, widths, GF, pi, hbarc2 -- lives in run_card_mod and
! is read from there directly by whoever needs it, and the selection is in
! cuts_mod. What is left here is what main hands the integrand and
! nothing else.
!
! Note the difference from run_card_mod, which is what the two names are
! for: `ecm` there is a number the user typed, `ecm_had` and `s` here are
! what the run is actually using.
!#####################################################################

module run_state_mod

  use iso_fortran_env, only: real64
  implicit none

! HADRONIC CMS ENERGY, AND ITS SQUARE.
! The integrand needs sqrt(s) every point -- PSGen takes the hadronic
! energy, not its square.
  real(real64) :: ecm_had
  real(real64) :: s

! RENORMALIZATION AND FACTORIZATION SCALES.
!
! xmur IS LIVE HERE, and that is the sharpest difference from the q q~
! skeleton this build was converted from. There the tree amplitude had
! four electroweak vertices and no power of gs, muR could not enter, and
! the old comment in this very spot forbade the integrand from reading
! xmur at all -- muR-independence was a free correctness check.
!
! That check is GONE. The gluon-initiated amplitude is a quark box
! carrying two powers of gs, so |M_loop|^2 ~ alpha_s(muR)^2 and the
! answer genuinely moves with muR. A run that does NOT move when
! mur_fixed moves now signals a bug (a frozen alpha_s), which is the
! exact opposite of the old rule. See olp_set_alphas in olp_mod.f90.
  real(real64) :: xmur, xmuf
  logical      :: scale_dynamic

! Re-centres the 7-point (muR,muF) grid (scales_mod) without moving its
! shape -- see run_card_mod::central_factor, which this mirrors exactly
! as xmur/xmuf mirror mur_fixed/muf_fixed.
  real(real64) :: xcentral_factor

! RECOLA PROCESS ID. One process, not six.
!
! The q q~ skeleton needed six: u/c shared one amplitude, d/s another, b
! stood alone (its charged-current partner is the TOP, exchanged in the
! t channel of the W+W- topologies), and each of those three needed both
! beam orderings because the amplitude is not symmetric under swapping
! the incoming quark and antiquark.
!
! None of that applies to two gluons. There is no flavour to sum over,
! and the two beam orderings are the SAME process -- the gluon is its own
! antiparticle, so `g g` is symmetric under swapping which beam supplies
! which gluon, and the single luminosity term f1(0)*f2(0) already covers
! both. (Note this is an initial state, so there is no identical-particle
! symmetry factor either: those belong to final states.)
  integer, parameter :: NPROC_LO = 1
  integer :: id_LO(NPROC_LO)

end module run_state_mod

!#####################################################################
