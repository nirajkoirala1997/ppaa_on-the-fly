!#####################################################################
! Booking, filling and reporting of the run's distributions.
!
! This is the layer between the run card and histogram_mod: it turns the
! card's [histograms] rows into booked histograms, evaluates the named
! observable on each accepted event, and writes the report. Nothing here
! knows about the process -- it works off the lab-frame momenta and a
! leg list, so adding a gluon changes only the leg count it is handed.
!
! WHY THE FILL IS DRIVEN FROM vsup. VEGAS runs two passes: an adaptation
! pass that refines the grid, and a production pass with the grid frozen
! that produces the quoted number. Only the second may be histogrammed --
! events from the adaptation pass are drawn from a different (and still
! moving) density, so mixing them in would bias every bin. vsup therefore
! calls
!
!     histo_begin_production      once, when the frozen-grid pass starts
!     histo_end_iteration         after each production iteration
!
! and until the first of those the fill below is a no-op. That is the
! whole reason the integrand can call analysis_fill unconditionally.
!
! The two routines at the BOTTOM OF THIS FILE, outside the module, are
! plain external subroutines so that vsup.f -- fixed-form F77 -- can call
! them without depending on the compiler's module name mangling.
!#####################################################################

module analysis_mod

  use iso_fortran_env, only: real64
  use histogram_mod,   only: histo1d_type, histo1d_init, histo1d_reset,       &
                             histo1d_fill,                                    &
                             histo1d_end_iteration, histo1d_rescale_ceiling,  &
                             histo1d_write,                                   &
                             histo1d_write_state, histo1d_read_state,         &
                             histo1d_mean, histo1d_error,                     &
                             write_report_header, write_total_xsec
  use scales_mod,      only: icentral
  use pdfvar_mod,      only: nvar
  use observables_mod, only: obs_pt, obs_rapidity, obs_pseudorapidity, obs_mass, &
                             obs_delta_phi, obs_delta_r, obs_mt
  use card_parser_mod, only: card_nhist, card_hist_name, card_hist_int,       &
                             card_hist_real, card_hist_str, card_hist_intlist,&
                             card_dump

  implicit none
  private
  public :: analysis_book_from_card, analysis_nhists
  public :: analysis_begin_production, analysis_fill, analysis_end_iteration
  public :: analysis_total, analysis_total_error, analysis_niter
  public :: analysis_variation, analysis_variation_error
  public :: analysis_rescale_for_ceiling
  public :: analysis_write_report
  public :: analysis_write_state, analysis_read_state

  integer, parameter :: MAXHIST = 32
  integer, parameter :: MAXLEGS = 8

  ! The observables a run card may name. Kept here rather than in
  ! observables_mod because this is card plumbing: observables_mod holds
  ! the formulas, this holds the spelling the card uses to reach them.
  ! Each one is a function of a single four-vector, so the leg list is
  ! summed first and the observable applied to the sum -- which is why
  ! `legs = 5,6` with `obs = pt` IS the missing transverse momentum, and
  ! why an invariant mass is just `mass` over every final-state leg.
  !
  ! Named IOBS_ rather than OBS_: Fortran is case-insensitive, so an
  ! OBS_MASS parameter would collide with the obs_mass function itself.
  integer, parameter :: IOBS_MASS = 1, IOBS_PT = 2, IOBS_Y = 3, IOBS_ETA = 4
  ! Two-object kinds: the leg list is split in half, each half summed, and
  ! the observable evaluated between the two resulting four-vectors.
  integer, parameter :: IOBS_DPHI = 5, IOBS_DR = 6, IOBS_MT = 7

  integer, save :: nhists = 0
  type(histo1d_type), save :: hists(MAXHIST)
  integer, save :: hobs(MAXHIST)          = 0    ! observable kind
  integer, save :: hnleg(MAXHIST)         = 0    ! how many legs it names
  integer, save :: hlegs(MAXLEGS,MAXHIST) = 0    ! which ones

  ! The total cross section, accumulated by exactly the machinery the
  ! bins use: a one-bin histogram of unit width that every accepted
  ! event falls into. Reusing histo1d_type here is not a trick to save
  ! code -- it is what guarantees the quoted total and the bin contents
  ! cannot drift apart through two copies of the same algebra.
  type(histo1d_type), save :: htot

  logical, save :: filling = .false.

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Book every distribution the run card asks for. Errors are fatal: a
! misspelled observable or a bad range means the run would produce a
! report that silently lacks a distribution somebody asked for.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine analysis_book_from_card()

  integer :: ih, ik, nb, nleg, i
  real(real64) :: xlo, xhi
  character(len=32)  :: nm
  character(len=128) :: obsname
  integer :: legs(MAXLEGS)

  nhists = 0

  do ih = 1, card_nhist()

     nm      = card_hist_name(ih)
     obsname = card_hist_str(ih, 'obs', '')

     ik = obs_kind(obsname)
     if (ik == 0) then
        print '(A)','ERROR histogram "'//trim(nm)//'": unknown observable "'// &
                    trim(obsname)//'"'
        stop 1
     endif

     ! `leg = 3` is just `legs = 3` written shorter; accept both, and
     ! refuse a row that gives neither when the observable needs legs.
     call card_hist_intlist(ih, 'legs', legs, nleg)
     if (nleg == 0) then
        legs(1) = card_hist_int(ih, 'leg', 0)
        if (legs(1) > 0) nleg = 1
     endif

     ! Every observable is a function of one four-vector applied to the
     ! summed legs, so they all want the same thing: at least one leg.
     if (nleg < 1) then
        print '(A)','ERROR histogram "'//trim(nm)//'": observable "'// &
                    trim(obsname)//'" needs at least one leg'
        stop 1
     endif

     if (ik >= IOBS_DPHI .and. (nleg < 2 .or. mod(nleg,2) /= 0)) then
        print '(A)','ERROR histogram "'//trim(nm)//'": observable "'// &
                    trim(obsname)//'" needs an even number of legs (at least two)'
        stop 1
     endif

     do i = 1, nleg
        if (legs(i) < 1) then
           print '(A)','ERROR histogram "'//trim(nm)//'": leg indices start at 1'
           stop 1
        endif
     end do

     nb  = card_hist_int (ih, 'nbins', 0)
     xlo = card_hist_real(ih, 'xmin', 0d0)
     xhi = card_hist_real(ih, 'xmax', 0d0)

     if (nb <= 0) then
        print '(A)','ERROR histogram "'//trim(nm)//'": nbins must be positive'
        stop 1
     endif
     if (xhi <= xlo) then
        print '(A)','ERROR histogram "'//trim(nm)//'": xmax must exceed xmin'
        stop 1
     endif

     if (nhists >= MAXHIST) then
        print '(A)','ERROR analysis_book_from_card: increase MAXHIST'
        stop 1
     endif

     nhists = nhists + 1
     hobs(nhists)          = ik
     hnleg(nhists)         = nleg
     hlegs(:,nhists)       = 0
     hlegs(1:nleg,nhists)  = legs(1:nleg)

     ! The label carries the leg list, so the report says 'mass(3,4)' rather
     ! than 'mass' -- otherwise mee and mvv are indistinguishable in a file
     ! somebody reads six months from now. Same shape as the KALEU build's,
     ! which makes the two reports comparable by eye as well as by script.
     call histo1d_init(hists(nhists), trim(nm), trim(leg_label(obsname, legs, nleg)), &
                       nb, xlo, xhi)

  end do

  end subroutine analysis_book_from_card

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer function analysis_nhists() result(n)
  n = nhists
  end function analysis_nhists

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 'mass' + legs (3,4)  ->  'mass(3,4)'.  An observable that takes no
! legs keeps its bare name.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function leg_label(obsname, legs, nleg) result(lab)
  character(len=*), intent(in) :: obsname
  integer,          intent(in) :: legs(:), nleg
  character(len=64) :: lab
  character(len=8)  :: num
  integer :: i

  lab = trim(adjustl(obsname))
  if (nleg <= 0) return

  lab = trim(lab)//'('
  do i = 1, nleg
     write(num,'(I0)') legs(i)
     if (i > 1) lab = trim(lab)//','
     lab = trim(lab)//trim(num)
  end do
  lab = trim(lab)//')'

  end function leg_label

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Card name -> observable. Returns 0 for anything unrecognised, which
! the caller turns into a fatal error.
!
! 'inv_mass' is accepted as a synonym for 'mass': an invariant mass is
! the mass of the summed legs, so an old card that says
! `obs = inv_mass  legs = 3,4,5,6` still means what it always meant.
! (It must now name its legs, though -- there is no bare form.)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer function obs_kind(name) result(kind)
  character(len=*), intent(in) :: name
  select case (trim(adjustl(name)))
    case ('mass','m','inv_mass','invmass') ; kind = IOBS_MASS
    case ('pt','pT')                       ; kind = IOBS_PT
    case ('rapidity','y')                  ; kind = IOBS_Y
    case ('pseudorapidity','eta')          ; kind = IOBS_ETA
    case ('dphi','delta_phi')              ; kind = IOBS_DPHI
    case ('dr','delta_r')                  ; kind = IOBS_DR
    case ('mt','transverse_mass')          ; kind = IOBS_MT
    case default                           ; kind = 0
  end select
  end function obs_kind

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Evaluate an observable on one four-vector -- the sum of the legs the
! run card named. Dispatch is on the integer kind rather than the name:
! this runs once per histogram per accepted event.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function obs_value(kind, p) result(val)
  integer,      intent(in) :: kind
  real(real64), intent(in) :: p(0:3)
  select case (kind)
    case (IOBS_MASS) ; val = obs_mass(p)
    case (IOBS_PT)   ; val = obs_pt(p)
    case (IOBS_Y)    ; val = obs_rapidity(p)
    case (IOBS_ETA)  ; val = obs_pseudorapidity(p)
    case default     ; val = 0d0
  end select
  end function obs_value

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Called by vsup when the frozen-grid production pass begins. Every
! accumulator is zeroed here rather than at booking time, so that
! anything the adaptation pass happened to reach is discarded.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine analysis_begin_production()
  integer :: ih
  do ih = 1, nhists
     call histo1d_reset(hists(ih))
  end do
  call histo1d_init(htot, 'total', 'sigma', 1, 0d0, 1d0)
  filling = .true.
  end subroutine analysis_begin_production

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Checkpoint round-trip for every accumulator this module owns.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine analysis_write_state(iunit)
  integer, intent(in) :: iunit
  integer :: ih

  write(iunit,'(I12)') nhists
  do ih = 1, nhists
     call histo1d_write_state(iunit, hists(ih))
  end do
  call histo1d_write_state(iunit, htot)

  end subroutine analysis_write_state

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Requires every histogram already booked -- it restores into an
! existing binning and never creates one, so the caller must have run
! analysis_book_from_card first.
!
! The final `filling = .true.` is the whole reason this is not just a
! loop over histo1d_read_state: it is the resume-mode stand-in for
! analysis_begin_production, which vsup never calls on a resume because
! it would zero exactly what was just restored.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine analysis_read_state(iunit)
  integer, intent(in) :: iunit
  integer :: ih, rnhists

  read(iunit,'(I12)') rnhists
  if (rnhists /= nhists) then
     print '(A,I0,A,I0)','ERROR analysis_read_state: '//                &
           'the card books ', nhists, ' histograms, the checkpoint holds ', &
           rnhists
     stop 1
  endif

  do ih = 1, nhists
     call histo1d_read_state(iunit, hists(ih))
  end do

  ! htot is booked by analysis_begin_production, which a resume never
  ! calls -- so it has to be booked here, with the same binning, before
  ! its accumulators can be read back into it.
  call histo1d_init(htot, 'total', 'sigma', 1, 0d0, 1d0)
  call histo1d_read_state(iunit, htot)

  filling = .true.

  end subroutine analysis_read_state

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! One accepted event. `weight` is the integrand value times the VEGAS
! vwgt; see histogram_mod's header for why that normalisation is the
! one that makes the accumulators come out as means.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine analysis_fill(qmom, n, weight)
  integer,      intent(in) :: n
  real(real64), intent(in) :: qmom(0:3,n)
  real(real64), intent(in) :: weight(nvar)

  integer :: ih, i
  real(real64) :: xval, psum(0:3), first_object(0:3), second_object(0:3)
  integer :: half_nleg

  if (.not. filling) return

  ! Fills the unit bin, so the total picks up every accepted event
  ! whether or not it lands inside any booked range.
  call histo1d_fill(htot, 0.5d0, weight)

  do ih = 1, nhists
! A histogram naming a leg beyond n cannot be filled from this qmom -- e.g.
! a real-emission caller filling a dipole's REDUCED kinematics (n=4, no
! gluon leg) against a histogram booked on leg 5. Skipping it here, rather
! than reading out of bounds, is also the physically right answer: that
! term has no such leg, so it should not contribute to that histogram, not
! be truncated into one that doesn't exist. Never triggered before this,
! since every prior caller filled with n covering every booked leg.
     if (maxval(hlegs(1:hnleg(ih),ih)) > n) cycle
     if (hobs(ih) >= IOBS_DPHI) then
        half_nleg = hnleg(ih)/2
        first_object  = sum(qmom(:,hlegs(1:half_nleg,ih)), dim=2)
        second_object = sum(qmom(:,hlegs(half_nleg+1:hnleg(ih),ih)), dim=2)
        select case (hobs(ih))
          case (IOBS_DPHI) ; xval = obs_delta_phi(first_object, second_object)
          case (IOBS_DR)   ; xval = obs_delta_r(first_object, second_object)
          case default     ; xval = obs_mt(first_object, second_object)
        end select
     else
        psum = 0d0
        do i = 1, hnleg(ih)
           psum = psum + qmom(:,hlegs(i,ih))
        end do
        xval = obs_value(hobs(ih), psum)
     endif
     call histo1d_fill(hists(ih), xval, weight)
  end do

  end subroutine analysis_fill

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine analysis_end_iteration()
  integer :: ih
  if (.not. filling) return
  do ih = 1, nhists
     call histo1d_end_iteration(hists(ih))
  end do
  call histo1d_end_iteration(htot)
  end subroutine analysis_end_iteration

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! The histogram-side total and its error, in fb. Same estimator as the
! bins, so comparing it against what VEGAS quotes is a genuine check.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function analysis_total() result(val)
  val = 0d0
  if (htot%nbins == 0) return
  val = histo1d_mean(htot, 1, icentral)
  end function analysis_total

  real(real64) function analysis_total_error() result(err)
  err = 0d0
  if (htot%nbins == 0) return
  err = histo1d_error(htot, 1, icentral)
  end function analysis_total_error

  integer function analysis_niter() result(nit)
  nit = htot%nit
  end function analysis_niter

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Call once, right after vsup returns, before ANY read of htot or hists
! (analysis_variation/_error, analysis_total/_error, analysis_write_report
! -- every one of them is built on histo1d_mean/histo1d_error). Corrects
! for an epsrel early exit; see histo1d_rescale_ceiling's own header for
! the derivation. it_prod is the card's ceiling, the same value already
! passed to checkpoint_configure and vsup by every caller.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine analysis_rescale_for_ceiling(it_prod)
  integer, intent(in) :: it_prod
  integer :: ih
  call histo1d_rescale_ceiling(htot, it_prod)
  do ih = 1, nhists
     call histo1d_rescale_ceiling(hists(ih), it_prod)
  end do
  end subroutine analysis_rescale_for_ceiling

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! The histogram-side total and its error at any flat variation `iv`, not
! only the central one. This is the ONLY estimator available away from
! the central column: vsup's own grid adapts on flo_LO's return value
! alone, so there is no independent VEGAS number for the other 6 scale
! points or any PDF member. A caller comparing a non-central point
! against a reference uses this, not analysis_total/analysis_total_error.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function analysis_variation(iv) result(val)
  integer, intent(in) :: iv
  val = 0d0
  if (htot%nbins == 0) return
  val = histo1d_mean(htot, 1, iv)
  end function analysis_variation

  real(real64) function analysis_variation_error(iv) result(err)
  integer, intent(in) :: iv
  err = 0d0
  if (htot%nbins == 0) return
  err = histo1d_error(htot, 1, iv)
  end function analysis_variation_error

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Combined report: run header, the echoed card, the total, then each
! distribution.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine analysis_write_report(filename, label, seed, n_prod,              &
                                   ecm_val, pdf_name, cuts, scales, amplitude, &
                                   ans, sd, chi2, quiet)
  character(len=*), intent(in) :: filename, label, pdf_name, cuts, scales
  character(len=*), intent(in) :: amplitude
  integer,          intent(in) :: seed, n_prod
  real(real64),     intent(in) :: ecm_val, ans, sd, chi2
! quiet: skip the "report written to" line -- for a caller that already shows
! the output path itself (e.g. lo_sa's PDF-member summary), so the path is
! not echoed twice. Defaults to printing it, unchanged for every other caller.
  logical,          intent(in), optional :: quiet

  integer, parameter :: iunit = 78
  integer :: ih, ios, iv
  real(real64) :: hmean(nvar), herr(nvar)

  open(unit=iunit, file=filename, status='replace', action='write', iostat=ios)
  if (ios /= 0) then
     print '(A)','ERROR analysis_write_report: cannot open '//trim(filename)
     return
  endif

  ! htot%nit is the number of iterations actually accumulated, which with an
  ! epsrel early exit is not the card's it_prod. The same value is already
  ! used a few lines below for write_total_xsec.
  call write_report_header(iunit, label, seed, n_prod, htot%nit, ecm_val, &
                           pdf_name, cuts, scales, amplitude)
  call card_dump(iunit)

  ! The histogram estimator for EVERY variation column -- ans/sd/chi2
  ! above are the true VEGAS estimate, valid only for the central column;
  ! write_total_xsec substitutes these in for every other one.
  do iv = 1, nvar
     hmean(iv) = histo1d_mean(htot, 1, iv)
     herr(iv)  = histo1d_error(htot, 1, iv)
  end do

  call write_total_xsec(iunit, label, ans, sd, chi2, hmean, herr, &
                        htot%hw(1,:), htot%hww(1,:), htot%nit)
  do ih = 1, nhists
     call histo1d_write(iunit, hists(ih), label, n_prod)
  end do

  close(iunit)

  if (.not. (present(quiet) .and. quiet)) &
     print '(A)','  report written to '//trim(filename)

  end subroutine analysis_write_report

end module analysis_mod

!#####################################################################
! F77-callable hooks for vsup.f. Plain external subroutines: the fixed-
! form integrator must not have to know how this compiler mangles module
! procedure names.
!#####################################################################

subroutine histo_begin_production()
  use analysis_mod, only : analysis_begin_production
  implicit none
  call analysis_begin_production()
end subroutine histo_begin_production

subroutine histo_end_iteration()
  use analysis_mod, only : analysis_end_iteration
  implicit none
  call analysis_end_iteration()
end subroutine histo_end_iteration

!#####################################################################
