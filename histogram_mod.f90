!#####################################################################
! On-the-fly histogramming for a VEGAS run.
!
! A generic 1D histogram filled inline in the production loop that
! produces the quoted cross section -- one run gives the full
! differential distribution instead of needing one run per bin.
!
! Ported from the KALEU Drell-Yan code's histogram_mod, with two
! deliberate differences.
!
! ONE VARIATION, NOT nvar. That code carries a bin grid per scale/PDF
! variation. Here every accumulator is a plain 1D array: this build
! produces the CENTRAL SCALE ONLY. Widening it later means giving hw,
! hww and hcur a second (variation) index and looping the fill -- the
! file format below already leaves the column block where it would go.
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
  implicit none
  private
  public :: histo1d_type
  public :: histo1d_init, histo1d_reset, histo1d_fill, histo1d_end_iteration
  public :: histo1d_write
  public :: histo1d_mean, histo1d_error
  public :: write_report_header, write_total_xsec

  type :: histo1d_type
     character(len=64) :: name   = ''
     character(len=32) :: xlabel = ''
     integer :: nbins = 0
     real(real64) :: xmin = 0d0, xmax = 0d0, binwidth = 0d0
     ! hcur is the current iteration; hw/hww are folded once per
     ! iteration by histo1d_end_iteration. nit counts the folds.
     real(real64), allocatable :: hcur(:), hw(:), hww(:)
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
  allocate(h%hcur(nbins), h%hw(nbins), h%hww(nbins))
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
! One event. `weight` is the integrand value times the VEGAS vwgt, i.e.
! f(x)*WGT/ITMX -- see the module header for why that particular
! normalisation is the one that makes hw come out as the mean.
!
! Events outside [xmin,xmax] are simply not binned. They still count in
! the total cross section, which is what makes a histogram range a
! plotting choice rather than a cut.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine histo1d_fill(h, xval, weight)
  type(histo1d_type), intent(inout) :: h
  real(real64),       intent(in)    :: xval, weight
  integer :: ibin

  if (h%nbins == 0) return
  if (xval < h%xmin .or. xval >= h%xmax) return

  ibin = int((xval - h%xmin)/h%binwidth) + 1
  if (ibin < 1 .or. ibin > h%nbins) return

  h%hcur(ibin) = h%hcur(ibin) + weight

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
! Bin content and its error, both as d(sigma)/d(obs). The two are kept
! as functions so the report writer and anything else that wants the
! numbers (the console summary, a cross-check against the total) agree
! on the estimator by construction rather than by copied algebra.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(real64) function histo1d_mean(h, ib) result(val)
  type(histo1d_type), intent(in) :: h
  integer,            intent(in) :: ib
  val = h%hw(ib)/h%binwidth
  end function histo1d_mean

  real(real64) function histo1d_error(h, ib) result(err)
  type(histo1d_type), intent(in) :: h
  integer,            intent(in) :: ib
  real(real64) :: dn
  err = 0d0
  if (h%nit < 2) return
  dn  = real(h%nit, real64)
  err = sqrt(max(0d0, dn*h%hww(ib) - h%hw(ib)**2)/(dn - 1d0))/h%binwidth
  end function histo1d_error

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine histo1d_write(iunit, h, contrib, npoints_per_iter)
  integer,            intent(in) :: iunit
  type(histo1d_type), intent(in) :: h
  character(len=*),   intent(in) :: contrib
  integer,            intent(in) :: npoints_per_iter

  integer :: ib
  real(real64) :: lo, hi, mid

  if (h%nbins == 0) return

  write(iunit,'(A)')     '# ==== HISTOGRAM '//trim(contrib)//' '//trim(h%name)//' ===='
  write(iunit,'(A,I0)')  '# nbins = ', h%nbins
  write(iunit,'(A,I0)')  '# niter = ', h%nit
  write(iunit,'(A,I0)')  '# npoints_per_iter = ', npoints_per_iter
  write(iunit,'(A)')     '# xlabel = '//trim(h%xlabel)
  write(iunit,'(A)')     '# columns: bin_low bin_high bin_center dsigma_dx error sum_w sum_w2'
  write(iunit,'(A)')     '# sum_w, sum_w2 are the raw per-iteration accumulators ' // &
                         '(NOT divided by the bin width) -- see the header block'

  ! Physical values first, raw accumulators trailing, so a plotting
  ! script keeps the column positions it expects and the accumulators
  ! can be appended without disturbing them.
  do ib = 1, h%nbins
     lo  = h%xmin + real(ib-1,real64)*h%binwidth
     hi  = lo + h%binwidth
     mid = 0.5d0*(lo+hi)
     write(iunit,'(7E24.16)') lo, hi, mid, &
                              histo1d_mean(h,ib), histo1d_error(h,ib), &
                              h%hw(ib), h%hww(ib)
  end do

  write(iunit,'(A)') '# ==== END HISTOGRAM '//trim(contrib)//' '//trim(h%name)//' ===='
  write(iunit,*)

  end subroutine histo1d_write

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Report header: everything a merge step needs to check that two files
! describe the same calculation.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine write_report_header(iunit, mode, seed, n_prod, it_prod, ecm_val, &
                                 pdf_name, cuts, scales)
  integer,          intent(in) :: iunit, seed, n_prod, it_prod
  character(len=*), intent(in) :: mode, pdf_name, cuts, scales
  real(real64),     intent(in) :: ecm_val

  write(iunit,'(A)')        '# ==== RUN ===='
  write(iunit,'(A)')        '# mode = '//trim(mode)
  write(iunit,'(A)')        '# generator = PSGen + VEGAS (vsup)'
  write(iunit,'(A,I0)')     '# seed = ', seed
  ! Points per production iteration and the iteration count, separately:
  ! the iteration IS the sample unit for every error in this file, so a
  ! merge step has to see both rather than only their product.
  write(iunit,'(A,I0)')     '# npoints_per_iter = ', n_prod
  write(iunit,'(A,I0)')     '# niter = ', it_prod
  write(iunit,'(A,I0)')     '# npoints_total = ', n_prod*it_prod
  write(iunit,'(A,E24.16)') '# ecm = ', ecm_val
  write(iunit,'(A)')        '# pdf = '//trim(pdf_name)
  write(iunit,'(A)')        '# amplitude = OpenLoops'
  ! Runs made under different selections or different scales describe
  ! different quantities and must never be merged, so both travel with
  ! the things a merge step checks, not only in the echoed card below.
  write(iunit,'(A)')        '# cuts = '//trim(cuts)
  write(iunit,'(A)')        '# scales = '//trim(scales)
  write(iunit,'(A)')        '# nvar = 1'
  write(iunit,'(A)')        '# variation_order = central'
  write(iunit,'(A)')        '# ==== END RUN ===='
  write(iunit,*)

  end subroutine write_report_header

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Total cross section block.
!
! Two numbers for the same quantity, on purpose. VEGAS quotes a
! variance-WEIGHTED average over iterations (and a chi2/it to say
! whether that weighting was legitimate); the histograms quote a plain
! unweighted mean, because a weighted one cannot be merged bin by bin
! across runs. They must agree well within the error, and `histo` here
! is the same estimator the bins use, applied to every accepted event.
! A disagreement between the two is the cheapest available signal that
! the run has not converged.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine write_total_xsec(iunit, contrib, xmur, xmuf, cross, cerr, chi2, &
                              htot, herr, sum_w, sum_w2, nit)
  integer,          intent(in) :: iunit, nit
  character(len=*), intent(in) :: contrib
  real(real64),     intent(in) :: xmur, xmuf, cross, cerr, chi2
  real(real64),     intent(in) :: htot, herr, sum_w, sum_w2

  write(iunit,'(A)')      '# ==== TOTAL '//trim(contrib)//' ===='
  write(iunit,'(A,I0)')   '# niter = ', nit
  ! mur/muf are ABSOLUTE scales in GeV, not the factors the KALEU report
  ! carries: this build has one scale point, so there is no grid for a
  ! factor to index into. Under a dynamic scale they hold the fixed
  ! fallback and the '# scales =' line above is what actually applies.
  write(iunit,'(A)')      '# columns: mur[GeV] muf[GeV] vegas[fb] vegas_err[fb] ' // &
                          'chi2_per_it histo[fb] histo_err[fb] sum_w sum_w2'
  write(iunit,'(2F10.4,7E24.16)') xmur, xmuf, cross, cerr, chi2, &
                                  htot, herr, sum_w, sum_w2
  write(iunit,'(A)') '# ==== END TOTAL '//trim(contrib)//' ===='
  write(iunit,*)

  end subroutine write_total_xsec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module histogram_mod

!#####################################################################
