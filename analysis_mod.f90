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
! The three routines at the BOTTOM OF THIS FILE, outside the module, are
! plain external subroutines so that vsup.f -- fixed-form F77 -- can call
! them without depending on the compiler's module name mangling.
!#####################################################################

module analysis_mod

  use iso_fortran_env, only: real64
  use histogram_mod,   only: histo1d_type, histo1d_init, histo1d_reset,       &
                             histo1d_fill,                                    &
                             histo1d_end_iteration, histo1d_write,            &
                             histo1d_mean, histo1d_error,                     &
                             write_report_header, write_total_xsec
  use observables_mod, only: obs_kind_from_name, obs_nlegs_required, obs_evaluate
  use input_card_mod,  only: card_nhist, card_hist_name, card_hist_int,       &
                             card_hist_real, card_hist_str, card_hist_intlist,&
                             card_dump

  implicit none
  private
  public :: analysis_book_from_card, analysis_nhists
  public :: analysis_begin_production, analysis_fill, analysis_end_iteration
  public :: analysis_total, analysis_total_error, analysis_niter
  public :: analysis_write_report

  integer, parameter :: MAXHIST = 32
  integer, parameter :: MAXLEGS = 8

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

  integer :: ih, ik, nb, nleg, nreq, i
  real(real64) :: xlo, xhi
  character(len=32)  :: nm
  character(len=128) :: obsname
  integer :: legs(MAXLEGS)

  nhists = 0

  do ih = 1, card_nhist()

     nm      = card_hist_name(ih)
     obsname = card_hist_str(ih, 'obs', '')

     ik = obs_kind_from_name(obsname)
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

     nreq = obs_nlegs_required(ik)
     select case (nreq)
       case (0)
          if (nleg /= 0) then
             print '(A)','ERROR histogram "'//trim(nm)//'": observable "'// &
                         trim(obsname)//'" takes no leg list'
             stop 1
          endif
       case (-1)
          if (nleg < 1) then
             print '(A)','ERROR histogram "'//trim(nm)//'": observable "'// &
                         trim(obsname)//'" needs at least one leg'
             stop 1
          endif
       case (-2)
          if (nleg < 2 .or. mod(nleg,2) /= 0) then
             print '(A)','ERROR histogram "'//trim(nm)//'": observable "'// &
                         trim(obsname)//'" needs an even number of legs (at least two)'
             stop 1
          endif
     end select

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

     call histo1d_init(hists(nhists), trim(nm), trim(obsname), nb, xlo, xhi)

  end do

  end subroutine analysis_book_from_card

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer function analysis_nhists() result(n)
  n = nhists
  end function analysis_nhists

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
! One accepted event. `weight` is the integrand value times the VEGAS
! vwgt; see histogram_mod's header for why that normalisation is the
! one that makes the accumulators come out as means.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine analysis_fill(qmom, n, weight)
  integer,      intent(in) :: n
  real(real64), intent(in) :: qmom(0:3,n)
  real(real64), intent(in) :: weight

  integer :: ih
  real(real64) :: xval

  if (.not. filling) return

  ! Fills the unit bin, so the total picks up every accepted event
  ! whether or not it lands inside any booked range.
  call histo1d_fill(htot, 0.5d0, weight)

  do ih = 1, nhists
     xval = obs_evaluate(hobs(ih), hlegs(1:hnleg(ih),ih), hnleg(ih), qmom, n)
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
  val = histo1d_mean(htot, 1)
  end function analysis_total

  real(real64) function analysis_total_error() result(err)
  err = 0d0
  if (htot%nbins == 0) return
  err = histo1d_error(htot, 1)
  end function analysis_total_error

  integer function analysis_niter() result(nit)
  nit = htot%nit
  end function analysis_niter

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Combined report: run header, the echoed card, the total, then each
! distribution.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine analysis_write_report(filename, label, seed, n_prod, it_prod,     &
                                   ecm_val, pdf_name, cuts, scales,            &
                                   xmur, xmuf, cross, cerr, chi2)
  character(len=*), intent(in) :: filename, label, pdf_name, cuts, scales
  integer,          intent(in) :: seed, n_prod, it_prod
  real(real64),     intent(in) :: ecm_val, xmur, xmuf, cross, cerr, chi2

  integer, parameter :: iunit = 78
  integer :: ih, ios

  open(unit=iunit, file=filename, status='replace', action='write', iostat=ios)
  if (ios /= 0) then
     print '(A)','ERROR analysis_write_report: cannot open '//trim(filename)
     return
  endif

  call write_report_header(iunit, label, seed, n_prod, it_prod, ecm_val, &
                           pdf_name, cuts, scales)
  call card_dump(iunit)
  call write_total_xsec(iunit, label, xmur, xmuf, cross, cerr, chi2,      &
                        analysis_total(), analysis_total_error(),         &
                        htot%hw(1), htot%hww(1), htot%nit)
  do ih = 1, nhists
     call histo1d_write(iunit, hists(ih), label, n_prod)
  end do

  close(iunit)

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
