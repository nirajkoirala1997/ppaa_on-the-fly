!#################################################
module globals_mod
    use iso_fortran_env
    implicit none

    public :: xq,ge,s
 
! OPENLOOOPS
    integer :: id_LO(5), id_NLO_1r(5), id_NLO_1loop(5)
    real(real64) :: ge,s,xq
end module globals_mod
!#################################################
