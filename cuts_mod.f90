!#####################################################################
! [cuts] for p p -> gamma gamma: photon pT and pseudorapidity, and a window on
! the diphoton invariant mass.
!
! Shared by lo_sa, kp_sa, virt_sa and rs_sa through an explicit ../ rule in
! each makefile, so the acceptance cannot differ between the pieces. Legs 3 and
! 4 are the photons; the radiated parton (leg 5) is never cut here. Frixione
! isolation needs that parton, so rs_sa applies it separately.
!
! pta_min and maa_min ARE REQUIRED AND MUST BE > 0: they regulate the t- and
! u-channel quark poles and the region m(gamma gamma) -> 0.
!#####################################################################

module cuts_mod

  use iso_fortran_env, only: real64
  use card_parser_mod, only: card_get_real
  implicit none
  private
  public :: cuts_init, event_cut, cuts_label, cuts_maa_window

  real(real64), save :: pta_min = 0d0, etaa_max = 0d0
  real(real64), save :: maa_min = 0d0, maa_max = 0d0

contains

  subroutine cuts_init(ecm)
  real(real64), intent(in) :: ecm

  pta_min  = card_get_real('cuts','pta_min',  -1d0)
  etaa_max = card_get_real('cuts','etaa_max', -1d0)
  maa_min  = card_get_real('cuts','maa_min',  -1d0)
  maa_max  = card_get_real('cuts','maa_max',   0d0)

  if (pta_min <= 0d0)  call fail('pta_min is REQUIRED and must be > 0')
  if (etaa_max <= 0d0) call fail('etaa_max is REQUIRED and must be > 0')
  if (maa_min <= 0d0)  call fail('maa_min is REQUIRED and must be > 0')
  if (maa_max <= 0d0) maa_max = ecm
  if (maa_max <= maa_min) call fail('maa_max must exceed maa_min')
  end subroutine cuts_init

! momenta(0:3, n), n >= 4. Only legs 3 and 4 are read.
  logical function event_cut(momenta) result(keep)
  real(real64), intent(in) :: momenta(0:,:)
  real(real64) :: pair(0:3), maa_squared

  keep = photon_accepted(momenta(:,3)) .and. photon_accepted(momenta(:,4))
  if (.not. keep) return

  pair = momenta(:,3) + momenta(:,4)
  maa_squared = pair(0)*pair(0) - pair(1)*pair(1) - pair(2)*pair(2) - pair(3)*pair(3)
  keep = maa_squared >= maa_min*maa_min .and. maa_squared <= maa_max*maa_max
  end function event_cut

! The pT cut comes first, so the pseudorapidity is only taken at pt > 0.
  logical function photon_accepted(photon) result(accepted)
  real(real64), intent(in) :: photon(0:3)
  real(real64) :: pt

  pt = sqrt(photon(1)*photon(1) + photon(2)*photon(2))
  accepted = pt >= pta_min
  if (accepted) accepted = abs(asinh(photon(3)/pt)) <= etaa_max
  end function photon_accepted

! The m(gamma gamma) window, for a phase-space generator to sample inside.
  subroutine cuts_maa_window(window_low, window_high)
  real(real64), intent(out) :: window_low, window_high
  window_low  = maa_min
  window_high = maa_max
  end subroutine cuts_maa_window

  function cuts_label() result(lab)
  character(len=160) :: lab
  write(lab,'(A,F0.3,A,F0.3,A,F0.3,A,F0.3)') 'pta>', pta_min, ',etaa<', etaa_max, &
       ',maa=', maa_min, ':', maa_max
  end function cuts_label

  subroutine fail(msg)
  character(len=*), intent(in) :: msg
  print '(A)','ERROR input card [cuts]: '//trim(msg)
  stop 1
  end subroutine fail

end module cuts_mod
