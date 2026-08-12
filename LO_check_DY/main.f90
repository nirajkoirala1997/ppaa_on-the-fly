program DY_at_LO
!#################################################################################!
! Leading-order  p p -> e+ e-  (Drell-Yan).
!
!     ./runLO [seed] [outdir]
!
! seed    positive integer RNG seed      (default: [run] seed from the card)
! outdir  where the report is written    (default: [output] dir from the card)
!
! Everything else comes from ../input.input, which fully determines the run: the
! same card, code and PDF set reproduce the numbers exactly. The card is echoed
! into the report.
!
! Phase space from PSGen, integration by VEGAS (vsup). PSGen supplies the mapping
! only; VEGAS keeps the loop and adapts all four dimensions.
!
! Distributions are filled ON THE FLY, inside the same production pass that
! produces the quoted cross section, so one run gives the full differential
! result. Which distributions is a run-card question, not a compile-time one --
! see the [histograms] section of input.input.
!
! ONE RUN IS ONE SEED. Independent statistics come from running several seeds in
! parallel and merging their reports afterwards -- tools/run_seeds.sh, then
! tools/merge_runs.py, which sums the raw accumulators every report carries. So
! every report here is seed-tagged, LO_<seed>.dat, and the untagged LO.dat is
! reserved for the merged one; see the naming note further down.
!
! WHAT THIS BUILD PRODUCES. At the card as shipped -- seed 40, 300k x 10
! adaptation then 3M x 20 production, mee > 50, fixed muF = mZ, PDF4LHC21_40
! member 0 --
!
!     sigma_LO = 1.86226247e+06 +- 1.522e+03 fb,  chi2/it = 0.806
!
! reproduced to every printed digit of every iteration before and after the
! histogramming was added, so that number is the regression target. It is a
! SINGLE-SEED number and stays one: campaign statistics belong in the merged
! report, never in this comment.
!
! IT IS NOT ref_data/LO_total.txt. That file, SETUP A of the KALEU Drell-Yan
! code, records 1.8520052e+06 +- 1.593e+02 fb -- 0.55% low, about 7 sigma, which
! is far too large to be statistical. The stored file is known to have drifted
! from the live KALEU code, so it is not a reference this build should be tuned
! against; compare against a KALEU RUN, not against that file.
!#################################################################################!
    use iso_fortran_env, only : real64, output_unit
    use globals_mod
    use constants_mod,   only : init_mass, set_masses_from_card, GF, alpha_gfermi
    use psgen_interface, only : psgen_init, psgen_ndim
    use input_card_mod,  only : card_report_unused
    use analysis_mod,    only : analysis_book_from_card, analysis_nhists,       &
                                analysis_total, analysis_total_error,           &
                                analysis_niter, analysis_write_report
    use misc_utils_mod,  only : read_input_file, cuts_label, scales_label,      &
                                seed, n_adapt, it_adapt, n_prod, it_prod,       &
                                out_dir, out_file, ecm, pdf_set, pdf_member,    &
                                scale_dynamic_in => scale_dynamic,              &
                                mur_fixed, muf_fixed,                           &
                                mZ_in, wZ_in, mW_in, wW_in, mH_in, mtop_in,     &
                                gfermi, mee_min, mee_max

    implicit none

    integer            :: ndim, ios
    real(real64)       :: ans, sd, chi2, htot, herr, pull
    character(len=256) :: outfile, arg

    external flo2_LO

!--- Run card ---------------------------------------------------------------!
    call read_input_file('../input.input')

!--- Command line -----------------------------------------------------------!
! Two overrides, and only two: the seed and where the report goes. Everything
! that changes the PHYSICS stays in the card, so a report plus its echoed card
! is still a complete description of the run.
!
! Both exist for the same reason -- a campaign runs one process per seed
! (tools/run_seeds.sh) and needs to tell each job who it is and where to write.
    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg)
       read(arg,*,iostat=ios) seed
       if (ios /= 0) then
          print '(A)','ERROR: seed argument is not an integer: '//trim(arg)
          stop 1
       endif
       ! Same rule as the card's. BRM48I takes IJKL >= 0 and would silently
       ! accept 0, giving every job that passed it the identical stream --
       ! which is exactly what a seed campaign exists to avoid.
       if (seed <= 0) then
          print '(A)','ERROR: seed must be a positive integer'
          stop 1
       endif
    endif

    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg)
       ! Silent truncation here would send a batch job's report somewhere
       ! other than where the launcher will look for it, and the launcher
       ! judges success by whether that file appeared.
       if (len_trim(arg) > len(out_dir)) then
          print '(A,I0,A)','ERROR: outdir argument is longer than ', len(out_dir), &
                           ' characters: '//trim(arg)
          stop 1
       endif
       out_dir = arg(1:len(out_dir))   ! length checked just above
    endif

    s             = ecm*ecm
    xmur          = mur_fixed
    xmuf          = muf_fixed
    scale_dynamic = scale_dynamic_in

! The mass window, squared. Applied to m(e+e-) from the momenta, not to x1*x2.
    mee_min2 = mee_min**2
    mee_max2 = mee_max**2

!--- Physics setup ----------------------------------------------------------!
    call init_mass
    call set_masses_from_card(mZ_in, wZ_in, mW_in, wW_in, mH_in, mtop_in)
    GF = gfermi
    ge = alpha_gfermi()

    call initpdfsetbyname(pdf_set)
    call initPDF(pdf_member)

    call ol_LO_init

!--- Distributions ----------------------------------------------------------!
! Booked before anything runs, so a bad [histograms] row fails in the first
! second of a run rather than after the integration.
    call analysis_book_from_card

! Everything has now had its turn at the card, so anything still unclaimed is a
! typo. Fatal on purpose: silently ignoring a misspelled key is exactly the
! failure the keyed format exists to prevent.
    call card_report_unused

!--- Banner -----------------------------------------------------------------!
    print *, " "
    print *, "===================================================="
    print *, "  p p -> e+ e-   at LO      [ PSGen + VEGAS ]"
    print *, "===================================================="
    write(*,'(a,i0)')         "  seed           = ", seed
    write(*,'(a,f12.2,a)')    "  ecm            = ", ecm, " GeV"
    write(*,'(a,a,a,i0,a)')   "  PDF            = ", trim(pdf_set),            &
                              "  (member ", pdf_member, ")"
    write(*,'(a,a)')          "  scales         = ", trim(scales_label())
    write(*,'(a,f9.2,a,f9.2,a)')                                              &
                              "  m(e+e-) window = [", mee_min, ",", mee_max, " ] GeV"
    write(*,'(a,es22.15)')    "  alpha (Gmu)    = ", ge
    write(*,'(a,es22.15)')    "  1/alpha        = ", 1.0_real64/ge
    write(*,'(a,i3)')         "  VEGAS dims     = ", 2 + psgen_ndim
    write(*,'(a,i0,a,i0)')    "  adaptation     = ", n_adapt, " x ", it_adapt
    write(*,'(a,i0,a,i0)')    "  production     = ", n_prod,  " x ", it_prod
    write(*,'(a,i0)')         "  histograms     = ", analysis_nhists()
    print *, "===================================================="
    print *, " "
! Redirected to a log by tools/run_seeds.sh, stdout is block-buffered, so the
! banner would otherwise sit in the buffer for the whole run. Flushed here and
! after every VEGAS iteration (vsup.f) the log tracks the run live.
    flush(output_unit)

!--- Phase space ------------------------------------------------------------!
    call psgen_init(6)

!--- Integrate --------------------------------------------------------------!
    ndim = 2 + psgen_ndim

    call brm48i(seed, 0, 0)
    call vsup(ndim, n_adapt, it_adapt, n_prod, it_prod, flo2_LO, ans, sd, chi2)

!--- Result -----------------------------------------------------------------!
! Two estimators for the same number. VEGAS quotes a variance-WEIGHTED average
! over the production iterations; the histograms accumulate a plain unweighted
! mean, because a weighted one cannot be merged bin by bin across runs. They
! must agree well inside the error -- a pull of a few sigma here is the
! cheapest available signal that the run has not converged, and it costs
! nothing, since the histogram total is filled by the same events.
    htot = analysis_total()
    herr = analysis_total_error()

    pull = 0.0_real64
    if (sd > 0.0_real64 .or. herr > 0.0_real64) &
       pull = (htot - ans)/sqrt(sd**2 + herr**2)

    print *, " "
    print *, "===================================================="
    write(*,'(a,es16.8,a,es10.3,a)') "  sigma_LO   = ", ans,  "  +-  ", sd,  "  fb  [VEGAS]"
    write(*,'(a,es16.8,a,es10.3,a)') "  sigma_hist = ", htot, "  +-  ", herr, "  fb  [histograms]"
    write(*,'(a,f10.4)')             "  chi^2/it   = ", chi2
    write(*,'(a,f10.4,a,i0,a)')      "  pull       = ", pull,                  &
                                     "        (", analysis_niter(), " production iterations)"
    print *, "===================================================="
    print *, " "

!--- Report -----------------------------------------------------------------!
! EVERY report is seed-tagged: [output] file = LO.dat is written as LO_<seed>.dat.
!
! That is not decoration. tools/merge_runs.py combines a campaign by globbing
! LO_*.dat, and the merged result lands on the untagged LO.dat -- which the glob
! cannot match, so a merged report can never be fed back into its own merge and
! have its points counted twice. The naming is the guard.
    call execute_command_line('mkdir -p '//trim(out_dir))
    outfile = trim(out_dir)//'/'//trim(seed_tagged(out_file, seed))
    call analysis_write_report(trim(outfile), 'LO', seed, n_prod, it_prod,     &
                               ecm, trim(pdf_set), trim(cuts_label()),         &
                               trim(scales_label()), xmur, xmuf,               &
                               ans, sd, chi2)

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! 'LO.dat', 1002  ->  'LO_1002.dat'
!
! The tag goes before the extension rather than after the whole name, so the
! file keeps whatever suffix the card chose and stays openable by anything that
! dispatches on it. Split at the LAST dot: a card naming its report
! 'LO.central.dat' means the extension is 'dat', not 'central.dat'. A name with
! no dot at all just gets the tag appended.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

  function seed_tagged(fname, iseed) result(tagged)
  character(len=*), intent(in) :: fname
  integer,          intent(in) :: iseed
  character(len=len(fname)+16) :: tagged
  character(len=16) :: tag
  integer :: idot

  write(tag,'(A,I0)') '_', iseed

  idot = index(trim(fname), '.', back=.true.)
  if (idot > 0) then
     tagged = fname(1:idot-1)//trim(tag)//trim(fname(idot:))
  else
     tagged = trim(fname)//trim(tag)
  endif

  end function seed_tagged

end program DY_at_LO
