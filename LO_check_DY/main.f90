program DY_at_LO
    use iso_fortran_env
    use globals_mod
    use printframes
    use constants_mod
    use kinparams
! LOCAL PARAMETERS
    integer            :: iprint 
    real(real64)       :: ans, sd, chi2
    character(len=100) :: run_tag, filename, mode
    real(real64), allocatable :: ai_lo2(:), err(:)

! EXTERNAL FUNCTIONS
    external flo2_LO

    ! Read vegas input
    open(unit=10, file='../run.vegas.dat', status='unknown')
    do i = 1, 9
        read(10, *)
    end do
    read(10, *) pt_LO
    read(10, *) its_LO
    close(10)
    npt_LO = pt_LO

    ! Read other parameter
    open(unit=10, file='../param_card.dat', status='unknown')
    read(10, *) ge
    close(10)

    ! Machine and run parameters
    open(unit=15, file='../run.machine.dat', status='unknown')
    read(15, *) mid
    read(15, *) ecm
    read(15, *) pdf_name_LO
    read(15, *) it_max
    read(15, *) xq_initial
    read(15, *) step_size
    read(15, *) run_tag
    read(15, *) iprint
    pdf_name_LO = 'MMHT2014lo68cl'
    close(15)

    ! Output filename
    open(unit=20, file='../output_files.dat', status='unknown')
    do i = 1, 3
        read(20, *)
    end do
    read(20, *) filename
    close(20)

    if (iprint == 1) call output(run_tag, filename)

    ! Initialization
    call init_mass
    call initpdfsetbyname(pdf_name_LO)
    call initPDF(0)
    call ol_LO_init

    s   = ecm * ecm

    print *, " "
    print *, "____________________________________"
    print *, " Calculating LO_DY"
    print *, "____________________________________"
    print *, " "

    call printframe1(pt_LO, its_LO)

    allocate(ai_lo2(it_max), err(it_max))

    xq = xq_initial

    do j = 1, it_max
        call printframe2(xq)

        call brm48i(40, 0, 0)
        call vsup(5, npt_LO, its_LO, flo2_LO, ans, sd, chi2)

        ai_lo2(j) = ans
        err(j)    = sd

        mode = "Leading Order"
        call printframe3(mode, ans, sd, chi2)

        xq = xq + step_size
    end do

    call printframe4(mode)

    xq = xq_initial
    do j = 1, it_max
        write(*,'(i7,3e27.15)') int(xq), ai_lo2(j), err(j)
        xq = xq + step_size
    end do

    if (iprint == 1) then
        open(unit=20, file='../summary/'//trim(run_tag)//'/'//trim(filename), status='unknown')
        xq = xq_initial
        do i = 1, it_max
            write(20, *) xq, ai_lo2(i), err(i)
            xq = xq + step_size
        end do
        close(20)
    end if

    deallocate(ai_lo2, err)

end program DY_at_LO 
