!#####################################################################
! On-the-fly histogramming for a VEGAS run.
!
! A generic 1D histogram filled inline in the production loop that
! produces the quoted cross section -- one run gives the full
! differential distribution instead of needing one run per bin.
!
! Ported from the KALEU Drell-Yan code's histogram_mod, with one
! deliberate difference kept below.
!
! nvar VARIATIONS, ONE SAMPLE ESTIMATOR. hw/hww/hcur carry a second
! (variation) index, exactly as the KALEU code's do, and the column
! layout below -- bin columns, then nvar (value,error) pairs, then nvar
! (sum_w,sum_w2) pairs -- matches its report format so the same merge
! tooling reads both. What does NOT change is the estimator itself (see
! next paragraph): each variation column is folded per ITERATION, not
! per event, because that is what this build's VEGAS integrator needs.
!
! THE SAMPLE UNIT IS A VEGAS ITERATION, NOT AN EVENT. KALEU draws
! independent points, so it can take the per-point spread as the error.
! VEGAS cannot: with MDS = 1 the points inside one iteration are
! STRATIFIED, deliberately anti-correlated, and the per-point spread
! overestimates the true variance of their sum. What IS independent is
! one iteration from the next -- which is exactly the estimator vsup
! uses for the total. So a bin accumulates into hcur during an
! iteration, and histo1d_end_iteration folds that one number into
! hw/hww. n below is therefore the number of PRODUCTION ITERATIONS.
!
! With
!     u_i = the bin's contribution from iteration i, in units where
!           the integrand is called with vwgt = WGT/ITMX
!     hw  = sum_i u_i          (already the mean, the 1/ITMX is in vwgt)
!     hww = sum_i u_i^2
! the per-iteration estimates are v_i = n*u_i, and
!     mean(v) = hw
!     err     = sqrt( (n*hww - hw^2) / (n-1) )
!
! RAW ACCUMULATORS. Each row carries hw, hww and n alongside the
! physical d(sigma)/d(obs), so several independent (seed) runs can be
! merged exactly: sum the accumulators, recompute the mean once. That is
! only valid between runs using the SAME points-per-iteration, since the
! sample unit carries that in its variance -- the header records it.
!#####################################################################

module histogram_mod

  use iso_fortran_env, only: real64
  use pdfvar_mod,      only: nvar
  implicit none
  private
  public :: histo1d_type
  public :: histo1d_init, histo1d_reset, histo1d_fill, histo1d_end_iteration
  public :: histo1d_rescale_ceiling
  public :: histo1d_write
  public :: histo1d_write_state, histo1d_read_state
  public :: histo1d_mean, histo1d_error
  public :: write_report_header, write_total_xsec

  type :: histo1d_type
     character(len=64) :: name   = ''
     character(len=32) :: xlabel = ''
     integer :: nbins = 0
     real(real64) :: xmin = 0d0, xmax = 0d0, binwidth = 0d0
     ! hcur is the current iteration; hw/hww are folded once per
     ! iteration by histo1d_end_iteration. nit counts the folds.
     ! Second index is the flat variation, sized nvar at booking time --
     ! see pdfvar_mod for what that index means.
     real(real64), allocatable :: hcur(:,:), hw(:,:), hww(:,:)
     integer :: nit = 0
  end type histo1d_type

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine histo1d_init(h, name, xlabel, nbins, xmin, xmax)
  type(histo1d_type), intent(inout) :: h
  character(len=*),   intent(in)    :: name, xlabel
  integer,            intent(in)    :: nbins
  real(real64),       intent(in)    :: xmin, xmax

  h%name     = name
  h%xlabel   = xlabel
  h%nbins    = nbins
  h%xmin     = xmin
  h%xmax     = xmax
  h%binwidth = (xmax - xmin)/real(nbins, real64)

  if (allocated(h%hcur)) deallocate(h%hcur)
  if (allocated(h%hw))   deallocate(h%hw)
  if (allocated(h%hww))  deallocate(h%hww)
  allocate(h%hcur(nbins,nvar), h%hw(nbins,nvar), h%hww(nbins,nvar))
  h%hcur = 0d0
  h%hw   = 0d0
  h%hww  = 0d0
  h%nit  = 0

  end subroutine histo1d_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Zero the accumulators, keeping the binning. Separate from
! histo1d_init because the caller that needs it -- the start of the
! production pass -- is re-zeroing histograms it booked earlier, and
! rebinning them from their own name/range fields would mean passing a
! variable's components back into itself as arguments, which Fortran
! does not allow to alias.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine histo1d_reset(h)
  type(histo1d_type), intent(inout) :: h
  if (h%nbins == 0) return
  h%hcur = 0d0
  h%hw   = 0d0
  h%hww  = 0d0
  h%nit  = 0
  end subroutine histo1d_reset

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! One event. `weight(nvar)` is the integrand's variation vector times
! the VEGAS vwgt, i.e. f_is(x)*WGT/ITMX per column -- see the module
! header for why that particular normalisation is the one that makes hw
! come out as the mean.
!
! Events outside [xmin,xmax] are simply not binned. They still count in
! the total cross section, which is what makes a histogram range a
! plotting choice rather than a cut.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine histo1d_fill(h, xval, weight)
  type(histo1d_type), intent(inout) :: h
  real(real64),       intent(in)    :: xval, weight(nvar)
  integer :: ibin

  if (h%nbins == 0) return
  if (xval < h%xmin .or. xval >= h%xmax) return

  ibin = int((xval - h%xmin)/h%binwidth) + 1
  if (ibin < 1 .or. ibin > h%nbins) return

  h%hcur(ibin,:) = h%hcur(ibin,:) + weight

  end subroutine histo1d_fill

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Close one VEGAS iteration: fold the running bin sums into the
! accumulators and start the next iteration from zero.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine histo1d_end_iteration(h)
  type(histo1d_type), intent(inout) :: h

  if (h%nbins == 0) return
  h%hw   = h%hw  + h%hcur
  h%hww  = h%hww + h%hcur**2
  h%hcur = 0d0
  h%nit  = h%nit + 1

  end subroutine histo1d_end_iteration

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Correct hw/hww for an epsrel early exit.
!
! Every fill carries weight = integrand * vwgt, with vwgt = WGT/ITMX fixed
! by vsup BEFORE the production loop starts -- it cannot know in advance
! how many iterations will actually run, only the card's ceiling ITMX
! (it_prod). Folding nit < ITMX such contributions therefore leaves
! hw = (nit/it_prod) * true_mean, not the mean itself: the constant
! implicitly assumed nit == it_prod, which a converged early exit breaks.
!
! Called once, right after vsup returns, before anything reads hw/hww
! (histo1d_mean/histo1d_error, and everything built on them: the console
! scale-variation table, the written report). Multiplying hw by
! c = it_prod/nit and hww by c**2 restores exactly the same mean/error
! that folding with vwgt = WGT/nit from the start would have produced --
! see niraj_is_working_towards_nlo_ggza.md for the derivation. A no-op
! when nit >= it_prod (the ceiling-reached and pre-epsrel cases), so
! every already-validated run is untouched.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine histo1d_rescale_ceiling(h, it_prod)
  type(histo1d_type), intent(inout) :: h
  integer,             intent(in)    :: it_prod
  real(real64) :: c
  if (h%nbins == 0) return
  if (h%nit <= 0 .or. h%nit >= it_prod) return
  c = real(it_prod, real64) / real(h%nit, real64)
  h%hw  = h%hw  * c
  h%hww = h%hww * c * c
  end subroutine histo1d_rescale_ceiling

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Bin content and its error, both as d(sigma)/d(obs), for variation
! column `iv`. The two are kept as functions so the report writer and
! anything else that wants the numbers (the console summary, a
! cross-check against the total) agree on the estimator by construction
! rather than by copied algebra.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function histo1d_mean(h, ib, iv) result(val)
  type(histo1d_type), intent(in) :: h
  integer,            intent(in) :: ib, iv
  val = h%hw(ib,iv)/h%binwidth
  end function histo1d_mean

  real(real64) function histo1d_error(h, ib, iv) result(err)
  type(histo1d_type), intent(in) :: h
  integer,            intent(in) :: ib, iv
  real(real64) :: dn
  err = 0d0
  if (h%nit < 2) return
  dn  = real(h%nit, real64)
  err = sqrt(max(0d0, dn*h%hww(ib,iv) - h%hw(ib,iv)**2)/(dn - 1d0))/h%binwidth
  end function histo1d_error

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine histo1d_write(iunit, h, contrib, npoints_per_iter)
  use pdfvar_mod, only: pdfvar_label
  use scales_mod, only: scale_label
  integer,            intent(in) :: iunit
  type(histo1d_type), intent(in) :: h
  character(len=*),   intent(in) :: contrib
  integer,            intent(in) :: npoints_per_iter

  integer :: ib, iv
  real(real64) :: lo, hi, mid

  if (h%nbins == 0) return

  write(iunit,'(A)')     '# ==== HISTOGRAM '//trim(contrib)//' '//trim(h%name)//' ===='
  write(iunit,'(A,I0)')  '# nbins = ', h%nbins
  write(iunit,'(A,I0)')  '# niter = ', h%nit
  write(iunit,'(A,I0)')  '# npoints_per_iter = ', npoints_per_iter
  write(iunit,'(A)')     '# xlabel = '//trim(h%xlabel)
  write(iunit,'(A)')     '# scale_order = '//scale_label
  write(iunit,'(A)')     '# variation_order = '//trim(pdfvar_label())
  write(iunit,'(A,I0,A)') '# columns 1-3: bin_low bin_high bin_center ; ' // &
                          'columns 4-', 3+2*nvar, ': (dsigma_dx error) per variation'
  write(iunit,'(A,I0,A,I0,A)') '# columns ', 4+2*nvar, '-', 3+4*nvar, &
                          ': (sum_w sum_w2) per variation -- raw per-iteration ' // &
                          'accumulators (NOT divided by the bin width)'
  write(iunit,'(A)',advance='no') '# columns: bin_low bin_high bin_center'
  do iv = 1, nvar
     write(iunit,'(A)',advance='no') ' dsigma_dx[fb] error[fb]'
  end do
  do iv = 1, nvar
     write(iunit,'(A)',advance='no') ' sum_w sum_w2'
  end do
  write(iunit,*)

  ! Physical values first, raw accumulators in a trailing block, so a
  ! plotting script keeps the column positions it expects for the
  ! central variation and the rest can be appended without disturbing
  ! them.
  do ib = 1, h%nbins
     lo  = h%xmin + real(ib-1,real64)*h%binwidth
     hi  = lo + h%binwidth
     mid = 0.5d0*(lo+hi)
     write(iunit,'(3E24.16)',advance='no') lo, hi, mid
     do iv = 1, nvar
        write(iunit,'(2E24.16)',advance='no') histo1d_mean(h,ib,iv), histo1d_error(h,ib,iv)
     end do
     do iv = 1, nvar
        write(iunit,'(2E24.16)',advance='no') h%hw(ib,iv), h%hww(ib,iv)
     end do
     write(iunit,*)
  end do

  write(iunit,'(A)') '# ==== END HISTOGRAM '//trim(contrib)//' '//trim(h%name)//' ===='
  write(iunit,*)

  end subroutine histo1d_write

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Checkpoint round-trip. Deliberately NOT histo1d_write above: that one
! writes physical bin content (divided by the bin width, errors folded)
! for a human or a plotting script. This one writes the raw state and
! nothing else, because it has to come back bit-identical.
!
! hcur is absent on purpose. A checkpoint is only ever taken at an
! iteration boundary, and histo1d_end_iteration zeroes hcur as its
! second-to-last act, so hcur is provably zero at every point this can
! be called from.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine histo1d_write_state(iunit, h)
  integer,            intent(in) :: iunit
  type(histo1d_type), intent(in) :: h

  write(iunit,'(A)')    trim(h%name)
  write(iunit,'(3I12)') h%nbins, nvar, h%nit
  if (h%nbins > 0) then
     write(iunit,'(E26.17)') h%hw
     write(iunit,'(E26.17)') h%hww
  endif

  end subroutine histo1d_write_state

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! The mirror. `h` must already be booked -- this restores accumulators
! into an existing binning, it never creates one. name and nbins are
! checked rather than assigned for exactly that reason: a [histograms]
! row edited between the kill and the resume would otherwise pour one
! histogram's sums into another's bins and produce a plausible,
! silently wrong plot.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine histo1d_read_state(iunit, h)
  integer,            intent(in)    :: iunit
  type(histo1d_type), intent(inout) :: h

  character(len=64) :: rname
  integer :: rbins, rnvar, rnit

  read(iunit,'(A)') rname
  read(iunit,'(3I12)') rbins, rnvar, rnit

  if (trim(rname) /= trim(h%name)) then
     print '(A)','ERROR histo1d_read_state: checkpoint holds histogram "'// &
                 trim(rname)//'" where "'//trim(h%name)//'" was booked'
     stop 1
  endif
  if (rbins /= h%nbins) then
     print '(A,I0,A,I0)','ERROR histo1d_read_state: histogram "'//   &
                 trim(h%name)//'" is booked with nbins = ', h%nbins, &
                 ' but the checkpoint holds ', rbins
     stop 1
  endif
  ! A checkpoint written under a different variation count (an older
  ! build, or [pdf] variations toggled between the kill and the resume)
  ! cannot be poured into this run's arrays -- fail loudly rather than
  ! silently misplace half the columns.
  if (rnvar /= nvar) then
     print '(A,I0,A,I0)','ERROR histo1d_read_state: histogram "'//    &
                 trim(h%name)//'" is booked with nvar = ', nvar,      &
                 ' but the checkpoint holds ', rnvar
     stop 1
  endif

  if (h%nbins > 0) then
     read(iunit,'(E26.17)') h%hw
     read(iunit,'(E26.17)') h%hww
     h%hcur = 0d0
  endif
  h%nit = rnit

  end subroutine histo1d_read_state

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Report header: everything a merge step needs to check that two files
! describe the same calculation.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine write_report_header(iunit, mode, seed, n_prod, niter, ecm_val, &
                                 pdf_name, cuts, scales, amplitude)
  use scales_mod, only: nscale, scale_label
  use pdfvar_mod, only: pdfvar_active, pdfvar_nmem, pdfvar_central, pdfvar_label
  integer,          intent(in) :: iunit, seed, n_prod, niter
  character(len=*), intent(in) :: mode, pdf_name, cuts, scales, amplitude
  real(real64),     intent(in) :: ecm_val

  write(iunit,'(A)')        '# ==== RUN ===='
  write(iunit,'(A)')        '# mode = '//trim(mode)
  write(iunit,'(A)')        '# generator = PSGen + VEGAS (vsup)'
  write(iunit,'(A,I0)')     '# seed = ', seed
  ! Points per production iteration and the iteration count, separately:
  ! the iteration IS the sample unit for every error in this file, so a
  ! merge step has to see both rather than only their product.
  write(iunit,'(A,I0)')     '# npoints_per_iter = ', n_prod
  write(iunit,'(A,I0)')     '# niter = ', niter
  write(iunit,'(A,I0)')     '# npoints_total = ', n_prod*niter
  write(iunit,'(A,E24.16)') '# ecm = ', ecm_val
  write(iunit,'(A)')        '# pdf = '//trim(pdf_name)
  ! WHICH PROVIDER AND WHICH PROCESS STRING ACTUALLY PRODUCED THIS FILE.
  ! It was hardcoded to 'Recola' until a MadGraph campaign wrote ten
  ! reports claiming to be Recola's. A report that misdescribes its own
  ! amplitude is worse than one that omits it: the numbers survive on
  ! disk long after the run, and there is nothing else in the file to
  ! contradict the label. It is passed in for the same reason cuts and
  ! scales are -- these are the fields a merge step checks, and a merge
  ! of two providers must be refusable.
  write(iunit,'(A)')        '# amplitude = '//trim(amplitude)
  ! Runs made under different selections or different scales describe
  ! different quantities and must never be merged, so both travel with
  ! the things a merge step checks, not only in the echoed card below.
  write(iunit,'(A)')        '# cuts = '//trim(cuts)
  write(iunit,'(A)')        '# scales = '//trim(scales)
  write(iunit,'(A,I0)')     '# nscale = ', nscale
  write(iunit,'(A)')        '# scale_order = '//scale_label
  ! The PDF-member axis. npdfvar is 0 when [pdf] variations is off, which
  ! keeps nvar = nscale and the file byte-compatible with a pre-PDF-axis
  ! run apart from these two extra lines. Runs spanning a different
  ! member count or a different central member are different quantities
  ! (tools/merge_runs.py's CONSISTENCY check), same as a different cut.
  write(iunit,'(A,I0)')     '# npdfvar = ', merge(pdfvar_nmem()+1, 0, pdfvar_active())
  write(iunit,'(A,I0)')     '# pdf_central_member = ', pdfvar_central()
  write(iunit,'(A,I0)')     '# nvar = ', nvar
  write(iunit,'(A)')        '# variation_order = '//trim(pdfvar_label())
  write(iunit,'(A)')        '# ==== END RUN ===='
  write(iunit,*)

  end subroutine write_report_header

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Total cross section block, one row per flat variation.
!
! Two estimators for the CENTRAL variation, on purpose: VEGAS quotes a
! variance-WEIGHTED average over iterations (chi2_per_it says whether
! that weighting was legitimate); the histogram quotes a plain
! unweighted mean, because a weighted one cannot be merged bin by bin
! across runs. They must agree well within the error, and a disagreement
! is the cheapest available signal that the run has not converged --
! see main.f90's pull check.
!
! EVERY OTHER VARIATION HAS ONLY THE HISTOGRAM ESTIMATOR. vsup's own
! grid adapts on the central column alone (see the integrand), so there
! is no independent VEGAS number for the other 6 scale points or any PDF
! member -- their `vegas`/`vegas_err` columns are the histogram values
! too, and chi2_per_it is written as 0 there (not applicable outside the
! one column vsup actually tracked). This still gives merge_runs.py a
! genuine per-run error to weight seeds by in every column, since the
! histogram accumulator is filled once per iteration for every
! variation regardless of which one VEGAS adapts on.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine write_total_xsec(iunit, contrib, ans, sd, chi2, &
                              htot, herr, sum_w, sum_w2, nit)
  use scales_mod, only: nscale, icentral, mur_factors, muf_factors
  use pdfvar_mod, only: pdfvar_central
  integer,          intent(in) :: iunit, nit
  character(len=*), intent(in) :: contrib
  real(real64),     intent(in) :: ans, sd, chi2
  real(real64),     intent(in) :: htot(nvar), herr(nvar), sum_w(nvar), sum_w2(nvar)

  integer :: is, imem
  real(real64) :: fmur, fmuf, v, verr, c2

  write(iunit,'(A)')      '# ==== TOTAL '//trim(contrib)//' ===='
  write(iunit,'(A,I0)')   '# niter = ', nit
  write(iunit,'(A)')      '# columns: mur[factor] muf[factor] pdf_member vegas[fb] ' // &
                          'vegas_err[fb] chi2_per_it histo[fb] histo_err[fb] sum_w sum_w2'

  do is = 1, nvar
     if (is <= nscale) then
        fmur = mur_factors(is)
        fmuf = muf_factors(is)
        imem = pdfvar_central()
     else
        fmur = 1d0
        fmuf = 1d0
        imem = is - nscale - 1
     endif

     if (is == icentral) then
        v = ans ; verr = sd ; c2 = chi2
     else
        v = htot(is) ; verr = herr(is) ; c2 = 0d0
     endif

     write(iunit,'(2F10.4,I6,7E24.16)') fmur, fmuf, imem, v, verr, c2, &
                                       htot(is), herr(is), sum_w(is), sum_w2(is)
  end do

  write(iunit,'(A)') '# ==== END TOTAL '//trim(contrib)//' ===='
  write(iunit,*)

  end subroutine write_total_xsec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module histogram_mod

!#####################################################################
