program DY_at_LO
!#################################################################################!
! Leading-order  p p -> e+ e-  (Drell-Yan).
!
! Phase space from PSGen, integration by VEGAS (vsup). PSGen supplies the mapping
! only; VEGAS keeps the loop and adapts all four dimensions.
!
! Target for this build: ref_data/LO_total.txt SETUP A of the KALEU Drell-Yan
! code, sigma = 1.8520052e+06 +- 1.593e+02 fb.
!#################################################################################!
    use iso_fortran_env, only : real64
    use globals_mod
    use constants_mod,   only : init_mass, hbarc2, GF, alpha_gfermi
    use psgen_interface, only : psgen_init, psgen_ndim
    use misc_utils_mod,  only : read_input_file, ecm, pdf_LO, pt_LO, its_LO,   &
                                run_tag, iprint, LO_file,                      &
                                mee_min, mee_max, mur_fixed, muf_fixed, gfermi

    implicit none

    integer            :: npt, its, ndim
    real(real64)       :: ans, sd, chi2
    character(len=200) :: outfile

    external flo2_LO

!--- Run card ---------------------------------------------------------------!
    call read_input_file('../input.input')

    s    = ecm*ecm
    xmur = mur_fixed
    xmuf = muf_fixed

! The mass window, squared. Applied to m(e+e-) from the momenta, not to x1*x2.
    mee_min2 = mee_min**2
    mee_max2 = mee_max**2

    npt = int(pt_LO)
    its = its_LO

!--- Physics setup ----------------------------------------------------------!
    call init_mass
    GF = gfermi
    ge = alpha_gfermi()

    call initpdfsetbyname(pdf_LO)
    call initPDF(0)

    call ol_LO_init

!--- Banner -----------------------------------------------------------------!
    print *, " "
    print *, "===================================================="
    print *, "  p p -> e+ e-   at LO      [ PSGen + VEGAS ]"
    print *, "===================================================="
    write(*,'(a,f12.2,a)')    "  ecm            = ", ecm, " GeV"
    write(*,'(a,a)')          "  PDF            = ", trim(pdf_LO)
    write(*,'(a,f12.4,a)')    "  muR = muF      = ", xmur, " GeV"
    write(*,'(a,f9.2,a,f9.2,a)')                                              &
                              "  m(e+e-) window = [", mee_min, ",", mee_max, " ] GeV"
    write(*,'(a,es22.15)')    "  alpha (Gmu)    = ", ge
    write(*,'(a,es22.15)')    "  1/alpha        = ", 1.0_real64/ge
    write(*,'(a,i3)')         "  VEGAS dims     = ", 2 + psgen_ndim
    print *, "===================================================="
    print *, " "

!--- Phase space ------------------------------------------------------------!
    call psgen_init(6)

!--- Integrate --------------------------------------------------------------!
    ndim = 2 + psgen_ndim

    call brm48i(40, 0, 0)
    call vsup(ndim, npt, its, flo2_LO, ans, sd, chi2)

!--- Result -----------------------------------------------------------------!
    print *, " "
    print *, "===================================================="
    write(*,'(a,es16.8,a,es10.3,a)') "  sigma_LO = ", ans, "  +-  ", sd, "  fb"
    write(*,'(a,f10.4)')             "  chi^2/it = ", chi2
    print *, "===================================================="
    print *, " "

    if (iprint == 1) then
       call execute_command_line('mkdir -p ../summary/'//trim(run_tag))
       outfile = '../summary/'//trim(run_tag)//'/'//trim(LO_file)
       open(unit=20, file=trim(outfile), status='unknown')
       write(20,'(a)') '# sigma_LO [fb]        error [fb]         chi2'
       write(20,'(3es24.15)') ans, sd, chi2
       close(20)
       print *, " Written: ", trim(outfile)
    end if

end program DY_at_LO
