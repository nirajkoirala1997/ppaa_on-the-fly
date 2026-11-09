program drell_yan_rs
!#################################################################################!
! p p -> gamma gamma at NLO: real emission minus Catani-Seymour dipoles.
!
!     ./runRS [seed] [outdir]
!
! Everything else comes from ../input.input. One run is one seed; reports are
! seed-tagged and the untagged name is reserved for a merged one.
!#################################################################################!
    use iso_fortran_env, only : real64, output_unit
    use run_state_mod,   only : ecm_had, s, xmur, xmuf, scale_dynamic, xcentral_factor
    use cuts_mod,        only : cuts_init, cuts_label, cuts_maa_window
    use photon_isolation_mod, only : photon_isolation_init
    use psgen_interface, only : psgen_init, psgen_ndim
    use card_parser_mod, only : card_report_unused, card_stamp
    use olp_mod,         only : olp_init, olp_set_collier_dir
    use analysis_mod,    only : analysis_book_from_card, analysis_total,          &
                                analysis_total_error, analysis_write_report,      &
                                analysis_rescale_for_ceiling
    use checkpoint_mod,  only : checkpoint_configure, checkpoint_exists, checkpoint_resume
    use pdfvar_mod,      only : pdfvar_init
    use run_card_mod,    only : read_input_file, scales_label, seed, n_adapt,     &
                                it_adapt, n_prod, it_prod, epsrel, resume_run,    &
                                out_dir, ecm, pdf_set, pdf_member, pdf_variations,&
                                scale_dynamic_in => scale_dynamic, mur_fixed,     &
                                muf_fixed, central_factor, mZ_in, sqsc_in

    implicit none

    character(len=*), parameter :: contribution = 'rs'
    character(len=*), parameter :: amplitude_label = 'recola: q q~, q g, g q -> A A + parton'

    integer            :: ndim, ios, start_it
    logical            :: converged
    real(real64)       :: ans, sd, chi2, alphasPDF, window_low, window_high
    character(len=256) :: outfile, arg, collier_dir

    external frs, alphasPDF

    call read_input_file('../input.input', 'rs')

    if (command_argument_count() >= 1) then
       call get_command_argument(1, arg)
       read(arg,*,iostat=ios) seed
       if (ios /= 0 .or. seed <= 0) then
          print '(A)','ERROR: seed must be a positive integer: '//trim(arg)
          stop 1
       endif
    endif
    if (command_argument_count() >= 2) then
       call get_command_argument(2, arg)
       if (len_trim(arg) > len(out_dir)) then
          print '(A)','ERROR: outdir argument is too long: '//trim(arg)
          stop 1
       endif
       out_dir = arg(1:len(out_dir))
    endif

    call card_stamp(trim(out_dir))
    out_dir = trim(out_dir)//'/rs'

    ecm_had         = ecm
    s               = ecm*ecm
    xmur            = mur_fixed
    xmuf            = muf_fixed
    scale_dynamic   = scale_dynamic_in
    xcentral_factor = central_factor

    call cuts_init(ecm)
    call photon_isolation_init()
    call pdfvar_init(pdf_set, pdf_member, pdf_variations)

    write(collier_dir,'(A,I0)') '.collier/rs_', seed
    call olp_set_collier_dir(collier_dir)
    call olp_init(alphasPDF(mZ_in), mZ_in, .false., .true.)

    call analysis_book_from_card
    call card_report_unused

    ndim = 2 + psgen_ndim + 1
    write(*,'(a,i0)')      '  seed       = ', seed
    write(*,'(a,a,a,i0)')  '  PDF        = ', trim(pdf_set), ' member ', pdf_member
    write(*,'(a,a)')       '  scales     = ', trim(scales_label())
    write(*,'(a,a)')       '  cuts       = ', trim(cuts_label())
    write(*,'(a,i0,a,i0)') '  adaptation = ', n_adapt, ' x ', it_adapt
    write(*,'(a,i0,a,i0)') '  production = ', n_prod,  ' x ', it_prod
    flush(output_unit)

    call cuts_maa_window(window_low, window_high)
    call psgen_init(6, sqsc_in, window_low, window_high)
    flush(output_unit)

    call execute_command_line('mkdir -p '//trim(out_dir))
    call checkpoint_configure(resume_run, seed, out_dir, ndim, n_prod, it_prod)

    if (checkpoint_exists()) then
       call checkpoint_resume(start_it)
    else
       call brm48i(seed, 0, 0)
       start_it = 0
    endif

    call vsup(ndim, n_adapt, it_adapt, n_prod, it_prod, epsrel, frs,           &
              ans, sd, chi2, converged, start_it > 0, start_it)

! An epsrel early exit otherwise biases the histogram means low.
    call analysis_rescale_for_ceiling(it_prod)

    write(*,'(a,es16.8,a,es10.3,a)') '  sigma_rs   = ', ans, ' +- ', sd, ' fb [VEGAS]'
    write(*,'(a,es16.8,a,es10.3,a)') '  sigma_hist = ', analysis_total(), ' +- ', &
                                     analysis_total_error(), ' fb [histograms]'

    write(outfile,'(A,A,I0,A)') trim(out_dir), '/rs_', seed, '.dat'
    call analysis_write_report(trim(outfile), contribution, seed, n_prod, ecm,    &
                               trim(pdf_set), trim(cuts_label()),                 &
                               trim(scales_label()), amplitude_label, ans, sd,    &
                               chi2, quiet=.true.)
    write(*,'(a)') '  report written to '//trim(outfile)

end program drell_yan_rs
!#################################################################################!
