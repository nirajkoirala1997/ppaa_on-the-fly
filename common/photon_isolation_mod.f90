!#################################################################################!
! Frixione (smooth-cone) photon isolation, process-independent, process
! MULTIPLICITY-independent -- ported from MATRIX's src/event/event.photon.
! isolation.cpp (standalone Fortran port and its own 6-case verification
! driver were prepared separately, in the MATRIX exploration notes; this file
! is that port, renamed to match this repository's file=module convention and
! wired to read its parameters from the card the way every other module here
! does, in place of that port's own default-valued derived type).
!
! Reference: S. Frixione, "Isolated photons in perturbative QCD",
! Phys.Lett. B429 (1998) 369, hep-ph/9801442.
!
! PHYSICS: a photon is isolated if, for every sub-cone radius Delta swept out
! by the partons found within the isolation cone (ordered by increasing
! Delta), the cumulative transverse energy of those partons stays below a
! bound H(Delta) that vanishes as Delta -> 0. This kills the photon
! fragmentation/collinear singularity while still allowing arbitrary soft
! radiation collinear to the photon (unlike a fixed hard cone cut) -- see
! niraj_is_working_towards_nlo_ggza.md section 2.4 for why this is REQUIRED
! here, not optional: the q g -> Z gamma q channel has a genuine final-state
! QED collinear singularity (the photon on an open quark line, via the axial
! Z g g triangle with one gluon off shell) that no pt cut on the photon alone
! regulates, since the photon and the quark can both be hard and collinear to
! each other.
!
! DELIBERATELY NOT PART OF cuts_mod.f90 (onshell/cuts_mod.f90, the Z+photon
! acceptance cut every piece of this study reads). That module only ever
! looks at legs 3 and 4 (Z, photon) -- exactly the two-body Born kinematics
! every piece except rs_sa's real term has. Isolation instead needs the
! radiated parton (leg 5), which only exists in a real-emission event, so it
! is kept as an independent module with its own call, not folded into the
! shared acceptance cut.
!
! ONE CALL PER CANDIDATE PHOTON, against however many QCD partons the event
! actually has (never the photon itself, never a lepton) -- npart may be 0,
! which returns .true. unconditionally (nothing to isolate against, exactly
! the case every Born-kinematics piece of this study is in).
!#################################################################################!
module photon_isolation_mod

  use iso_fortran_env,  only : real64
  use run_card_mod,     only : iso_mode, iso_n, iso_epsilon, iso_et_max,    &
                               iso_delta_0, iso_m_hard

  implicit none
  private
  public :: photon_isolation_init, photon_isolated

  integer,  parameter :: dp = real64

! [isolation] card values -- see photon_isolation_init. Defaults match
! MATRIX's own ppeexa03 defaults (mode 1 = static cone, (1-cos)/(1-cos0)
! ratio, eps*ET_gamma normalisation -- the "ATLAS-style" convention per the
! MATRIX manual).
  integer  :: mode      = 1
  real(dp) :: frix_n     = 1.0_dp
  real(dp) :: epsilon    = 0.5_dp
  real(dp) :: ET_max     = 5.0_dp
  real(dp) :: delta_0    = 0.4_dp
  real(dp) :: M_hard     = 91.188_dp

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! Call once, after run_card_mod's read_input_file. [isolation] mode/n/epsilon/
! et_max/delta_0/m_hard are read and validated THERE, not here: all four
! programs share one card and only rs_sa links this module, so run_card_mod
! (linked into all four) has to be the one that claims the keys -- see its
! [isolation] block for the defaults. See the module header for what each
! one means (frixione_isolated's own comments give the closed forms).
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  subroutine photon_isolation_init()
  mode     = iso_mode
  frix_n   = iso_n
  epsilon  = iso_epsilon
  ET_max   = iso_et_max
  delta_0  = iso_delta_0
  M_hard   = iso_m_hard
  end subroutine photon_isolation_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! p_gamma        photon 4-momentum, (E,px,py,pz)
! p_parton(:,i)  the npart QCD parton 4-momenta of THIS event (exclude the
!                photon and the Z/leptons) -- npart may be 0
! isolated       .true. if the photon passes the cumulative sub-cone test
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  logical function photon_isolated(p_gamma, p_parton, npart) result(isolated)
  real(dp), intent(in) :: p_gamma(0:3)
  integer,  intent(in) :: npart
  real(dp), intent(in) :: p_parton(0:3, max(npart,1))

  integer  :: i, j, k, nsel
  real(dp) :: dR(max(npart,1)), ET(max(npart,1))
  integer  :: order(max(npart,1))
  real(dp) :: delta0_eff, norm, ET_gamma, E_sum, H_i, tmp
  integer  :: itmp

  isolated = .true.
  if (npart <= 0) return

  ET_gamma = transverse_pT(p_gamma)

  select case (mode)
  case (3)
     delta0_eff = M_hard / (ET_gamma*sqrt(epsilon))
  case default
     delta0_eff = delta_0
  end select

  select case (mode)
  case (2)
     norm = ET_max
  case default
     norm = epsilon*ET_gamma
  end select

  nsel = 0
  do i = 1, npart
     tmp = delta_R(p_gamma, p_parton(:,i))
     if (tmp < delta0_eff) then
        nsel = nsel + 1
        dR(nsel) = tmp
        ET(nsel) = transverse_pT(p_parton(:,i))
        order(nsel) = nsel
     endif
  end do
  if (nsel == 0) return

! Insertion sort by increasing Delta -- nsel is small (O(1) partons at this
! study's multiplicity), so this need not be anything fancier.
  do i = 2, nsel
     do j = i, 2, -1
        if (dR(order(j)) < dR(order(j-1))) then
           itmp = order(j) ; order(j) = order(j-1) ; order(j-1) = itmp
        else
           exit
        endif
     end do
  end do

  E_sum = 0.0_dp
  do k = 1, nsel
     i = order(k)
     H_i = frixione_bound(dR(i), delta0_eff, frix_n, norm)
     E_sum = E_sum + ET(i)
     if (E_sum > H_i) then
        isolated = .false.
        return
     endif
  end do

  end function photon_isolated

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! H(Delta) = norm * [(1-cos Delta)/(1-cos Delta0)]^n
! Evaluated as [sin(Delta/2)/sin(Delta0/2)]^2 (1-cos x = 2 sin^2(x/2)): same
! value, but no cancellation as Delta -> 0, where 1-cos(Delta) rounds to 0.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  real(dp) function frixione_bound(delta, delta0, n, norm) result(H)
  real(dp), intent(in) :: delta, delta0, n, norm
  real(dp) :: ratio
  ratio = (sin(0.5_dp*delta)/sin(0.5_dp*delta0))**2
  if (abs(n - 1.0_dp) < 1.0d-14) then
     H = norm*ratio
  else
     H = norm*ratio**n
  endif
  end function frixione_bound

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! DeltaR^2 = Deta^2 + Dphi^2, the ordinary eta-phi measure.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  real(dp) function delta_R(pgam, ppar) result(dR)
  real(dp), intent(in) :: pgam(0:3), ppar(0:3)
  real(dp) :: eta1, eta2, phi1, phi2, deta, dphi
  real(dp), parameter :: pi = 3.141592653589793238462643_dp

  eta1 = pseudorapidity(pgam)
  eta2 = pseudorapidity(ppar)
  phi1 = atan2(pgam(2), pgam(1))
  phi2 = atan2(ppar(2), ppar(1))

  deta = eta1 - eta2
  dphi = abs(phi1 - phi2)
  if (dphi > pi) dphi = 2.0_dp*pi - dphi

  dR = sqrt(deta*deta + dphi*dphi)
  end function delta_R

  real(dp) function pseudorapidity(p) result(eta)
  real(dp), intent(in) :: p(0:3)
  real(dp) :: pmod, pz, pt
  pz = p(3)
  pt = sqrt(p(1)**2 + p(2)**2)
  pmod = sqrt(p(1)**2 + p(2)**2 + p(3)**2)
  if (pt < 1.0d-14) then
     eta = sign(1.0d30, pz)
  else
     eta = 0.5_dp*log((pmod+pz)/(pmod-pz))
  endif
  end function pseudorapidity

  real(dp) function transverse_pT(p) result(v)
  real(dp), intent(in) :: p(0:3)
  v = sqrt(p(1)**2 + p(2)**2)
  end function transverse_pT

end module photon_isolation_mod
!#################################################################################!
