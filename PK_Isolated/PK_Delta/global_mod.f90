!=================================================
module globals_mod
    use iso_fortran_env
    implicit none
    real(real64) :: s
    real(real64) :: ge

! OPENLOOOPS
    integer :: id_LO(5), id_NLO_1r(5), id_NLO_1loop(5)

! VEGAS PARAMETERS
    real(real64) :: pt_delta
    integer :: its_delta,npt_delta

! MACHINE PARAMETERS
    real(real64) :: xq, xq_initial, step_size, ecm
    character(len=50)  :: pdf_name_NLO
    integer :: j, it_max
! CUTS USED 
    real(real64) :: delta

! External functions
  
end module globals_mod
!=================================================
