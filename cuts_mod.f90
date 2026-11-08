!#####################################################################
! [cuts] for p p -> e+ e-: a window on the dilepton invariant mass.
!
! Shared by lo_sa, kp_sa, virt_sa and rs_sa through an explicit ../ rule in
! each makefile, so the acceptance cannot differ between the pieces. Legs 3 and
! 4 are the leptons; the radiated parton (leg 5) is never cut.
!
! mee_min IS REQUIRED AND MUST BE > 0: the photon propagator diverges as
! m(ee) -> 0.
!
! There are deliberately no lepton rapidity cuts. One q q~ amplitude serves
! both beam orderings, which is only valid for observables symmetric under
! the swap of the two leptons' roles.
!#####################################################################

module cuts_mod

  use iso_fortran_env, only: real64
  use card_parser_mod, only: card_get_real
  implicit none
  private
  public :: cuts_init, event_cut, cuts_label, cuts_mee_window

  real(real64), save :: mee_min = 0d0, mee_max = 0d0

contains

  subroutine cuts_init(ecm)
  real(real64), intent(in) :: ecm

  mee_min = card_get_real('cuts','mee_min', -1d0)
  mee_max = card_get_real('cuts','mee_max',  0d0)

  if (mee_min <= 0d0) &
     call fail('mee_min is REQUIRED and must be > 0 -- the photon pole '// &
               'at m(ee) = 0 is not integrable')
  if (mee_max <= 0d0) mee_max = ecm
  if (mee_max <= mee_min) call fail('mee_max must exceed mee_min')
  end subroutine cuts_init

! momenta(0:3, n), n >= 4. Only legs 3 and 4 are read.
  logical function event_cut(momenta) result(keep)
  real(real64), intent(in) :: momenta(0:,:)
  real(real64) :: pair(0:3), mee_squared

  pair = momenta(:,3) + momenta(:,4)
  mee_squared = pair(0)*pair(0) - pair(1)*pair(1) - pair(2)*pair(2) - pair(3)*pair(3)
  keep = mee_squared >= mee_min*mee_min .and. mee_squared <= mee_max*mee_max
  end function event_cut

! The m(ee) window, for a phase-space generator to sample inside.
  subroutine cuts_mee_window(window_low, window_high)
  real(real64), intent(out) :: window_low, window_high
  window_low  = mee_min
  window_high = mee_max
  end subroutine cuts_mee_window

  function cuts_label() result(lab)
  character(len=160) :: lab
  write(lab,'(A,F0.3,A,F0.3)') 'mee=', mee_min, ':', mee_max
  end function cuts_label

  subroutine fail(msg)
  character(len=*), intent(in) :: msg
  print '(A)','ERROR input card [cuts]: '//trim(msg)
  stop 1
  end subroutine fail

end module cuts_mod
