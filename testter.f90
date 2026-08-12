program test_input_reader
    use iso_fortran_env, only: real64
    implicit none

    call read_input_file('input.input')
end program test_input_reader

!-------------------------------------------
subroutine read_input_file(filename)
    use iso_fortran_env, only: real64
    implicit none

    character(len=*), intent(in) :: filename
    character(len=256) :: line
    integer :: iunit, ios
    character(len=100) :: run_tag, pdf_LO, pdf_NLO
    real(real64) :: ecm, invmass_base, step_size
    integer :: iprint
    real(real64) :: pt_real, pt_virtual, pt_reg1, pt_reg2, pt_LO, pt_plus1, pt_plus2
    integer :: its_real, its_virtual, its_reg1, its_reg2, its_LO, its_plus1, its_plus2
    character(len=100) :: LO_file, real_file, virtual_file, PK_file
    character(len=100) :: PK_plus_file, PK_reg_file, PK_delta_file
    character(len=100) :: LO_all_file, LO_ref_file
    character(len=100) :: real_all_file, real_ref_file
    character(len=100) :: virtual_all_file, virtual_ref_file
    character(len=100) :: PK_all_file, PK_ref_file
    character(len=100) :: PK_plus_all_file, PK_reg_all_file, PK_delta_all_file
    character(len=100) :: NLO_ref_file

    iunit = 99
    open(unit=iunit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, "Error opening file:", trim(filename)
        stop
    end if

    call next_data_line(iunit, line); read(line,*) ecm
    call next_data_line(iunit, line); read(line,*) invmass_base
    call next_data_line(iunit, line); read(line,*) step_size

    call next_data_line(iunit, line); read(line,*) run_tag
    call next_data_line(iunit, line); read(line,*) iprint

    call next_data_line(iunit, line); read(line,*) pdf_LO
    call next_data_line(iunit, line); read(line,*) pdf_NLO

    call next_data_line(iunit, line); read(line,*) pt_real
    call next_data_line(iunit, line); read(line,*) its_real

    call next_data_line(iunit, line); read(line,*) pt_virtual
    call next_data_line(iunit, line); read(line,*) its_virtual

    call next_data_line(iunit, line); read(line,*) pt_reg1, pt_reg2
    call next_data_line(iunit, line); read(line,*) its_reg1, its_reg2

    call next_data_line(iunit, line); read(line,*) pt_LO
    call next_data_line(iunit, line); read(line,*) its_LO

    call next_data_line(iunit, line); read(line,*) pt_plus1, pt_plus2
    call next_data_line(iunit, line); read(line,*) its_plus1, its_plus2

    call next_data_line(iunit, line); read(line,*) LO_file
    call next_data_line(iunit, line); read(line,*) real_file
    call next_data_line(iunit, line); read(line,*) virtual_file
    call next_data_line(iunit, line); read(line,*) PK_file
    call next_data_line(iunit, line); read(line,*) PK_plus_file
    call next_data_line(iunit, line); read(line,*) PK_reg_file
    call next_data_line(iunit, line); read(line,*) PK_delta_file

    call next_data_line(iunit, line); read(line,*) LO_all_file
    call next_data_line(iunit, line); read(line,*) LO_ref_file
    call next_data_line(iunit, line); read(line,*) real_all_file
    call next_data_line(iunit, line); read(line,*) real_ref_file
    call next_data_line(iunit, line); read(line,*) virtual_all_file
    call next_data_line(iunit, line); read(line,*) virtual_ref_file
    call next_data_line(iunit, line); read(line,*) PK_all_file
    call next_data_line(iunit, line); read(line,*) PK_ref_file
    call next_data_line(iunit, line); read(line,*) PK_plus_all_file
    call next_data_line(iunit, line); read(line,*) PK_reg_all_file
    call next_data_line(iunit, line); read(line,*) PK_delta_all_file
    call next_data_line(iunit, line); read(line,*) NLO_ref_file

    close(iunit)

    print *, "✔ Input read summary:"
    print *, "  ECM =", ecm
    print *, "  invmass_base =", invmass_base
    print *, "  step_size =", step_size
    print *, "  run_tag =", trim(run_tag)
    print *, "  pdf_LO =", trim(pdf_LO)
    print *, "  pdf_NLO =", trim(pdf_NLO)
    print *, "  VEGAS (real):", pt_real, "pts,", its_real, "iterations"
    print *, "  VEGAS (virtual):", pt_virtual, "pts,", its_virtual, "iterations"
end subroutine read_input_file

subroutine next_data_line(unit, line)
    implicit none
    integer, intent(in) :: unit
    character(len=100), intent(out) :: line
    character(len=100) :: temp
    do
        read(unit,'(A)',end=100) temp
        if (trim(temp) == "") cycle                          ! skip empty lines
        if (index(adjustl(temp), "!") == 1) cycle           ! skip lines starting with "cc"
        line = temp
        return
    end do
100 continue
    line = ""   ! End of file or no valid line found
end subroutine next_data_line
