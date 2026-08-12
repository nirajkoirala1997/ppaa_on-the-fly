program DY_at_1loop
    use iso_fortran_env
    use globals_mod
    use printframes
    use constants_mod
! LOCAL PARAMETERS
    integer            :: iprint,mid 
    real(real64)       :: ans, sd, chi2
    character(len=100) :: run_tag, filename, mode
    real(real64), allocatable :: ai_nlo2(:), err(:)

! EXTERNAL FUNCTIONS
    external fnlo2

    ! Read vegas input
    open(unit=10, file='../run.vegas.dat', status='unknown')
    do i = 1, 3
        read(10, *)
    end do
    read(10, *) pt_1loop
    read(10, *) its_1loop
    close(10)
    npt_1loop = pt_1loop

    ! Read other parameter
    open(unit=10, file='../param_card.dat', status='unknown')
    read(10, *) ge
    close(10)

    ! Machine and run parameters
    open(unit=15, file='../run.machine.dat', status='unknown')
    read(15, *) mid
    read(15, *) ecm
    read(15, *) pdf_name_NLO
    read(15, *) it_max
    read(15, *) xq_initial
    read(15, *) step_size
    read(15, *) run_tag
    read(15, *) iprint
    close(15)

    ! Output filename
    open(unit=20, file='../output_files.dat', status='unknown')
        read(20, *)
    read(20, *) filename
    close(20)

    if (iprint == 1) call output(run_tag, filename)

    ! Initialization
    call init_mass
    call initpdfsetbyname(pdf_name_NLO)
    call initPDF(0)
    call ol_NLO_init_1loop

    s   = ecm * ecm

    print *, " "
    print *, "____________________________________"
    print *, " Calculating virtual contribution   "
    print *, "____________________________________"
    print *, " "

    call printframe1(pt_1loop, its_1loop)

    allocate(ai_nlo2(it_max), err(it_max))

    xq = xq_initial

    do j = 1, it_max
        call printframe2(xq)

        call brm48i(40, 0, 0)
        call vsup(3, npt_1loop, its_1loop, fnlo2, ans, sd, chi2)

        ai_nlo2(j) = ans
        err(j)    = sd

        mode = "virtual contribution"
        call printframe3(mode, ans, sd, chi2)

        xq = xq + step_size
    end do

    call printframe0(mode)

    xq = xq_initial
    do j = 1, it_max
        write(*,'(i7,3e27.15)') int(xq), ai_nlo2(j), err(j)
        xq = xq + step_size
    end do

    if (iprint == 1) then
        open(unit=20, file='../summary/'//trim(run_tag)//'/'//trim(filename), status='unknown')
        xq = xq_initial
        do i = 1, it_max
            write(20, *) xq, ai_nlo2(i), err(i)
            xq = xq + step_size
        end do
        close(20)
    end if

    deallocate(ai_nlo2, err)

end program DY_at_1loop 
