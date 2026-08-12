module printframes
  use iso_fortran_env, only: real64
  implicit none
contains

  subroutine printframe4(name)
    character(len=*), intent(in) :: name

     print*,"  "
     print*,"  "
     print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
     print*,"    Calculating "//trim(name)
     print*,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
     print*,"  "
  end subroutine printframe4

  subroutine printframe1(pt1, its1)
    real(real64), intent(in) :: pt1
    integer, intent(in) :: its1

    print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,'(A,ES25.0)') "Using Vegas points:", pt1
    write(*,'(A,I10)') "       Iteration:              ", its1
    print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  end subroutine printframe1

  subroutine printframe2(xq)
    real(real64), intent(in) :: xq
    print *, " "
    write(*,'(A,I5,A)') achar(27)//'[1;32m'//"For xq = ", int(xq), achar(27)//'[0m'
    print *, " "
  end subroutine printframe2

  subroutine printframe3(name, xintegral, error, chi2)
    character(len=*), intent(in) :: name
    real(real64), intent(in) :: xintegral, error, chi2

    print *, "  "
    print *, "  "
    write(*,'(A,ES15.8,2A,ES15.8)') achar(27)//'[1;32m'//"Integral "//trim(name)//": ",xintegral, achar(27)//'[0m','  +- ', error
    write(*,'(A,ES15.8)') "with chi^2    = ", chi2
    print *, " "
    print *, " "
  end subroutine printframe3

  subroutine printframe0(name)
    character(len=*), intent(in) :: name

    print *, " "
    write(*,'(A,A,A,A)') achar(27)//'[1;32m'//"     xq         Integral ", &
         trim(name), "                 error", achar(27)//'[0m'
  end subroutine printframe0

end module printframes

