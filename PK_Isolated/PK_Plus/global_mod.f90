!=================================================
module globals_mod
    use iso_fortran_env
    implicit none
    real(real64) :: s
    real(real64) :: ge

! OPENLOOOPS
    integer :: id_LO(5), id_NLO_1r(5), id_NLO_1loop(5)

! VEGAS PARAMETERS
    real(real64) :: pt_Plus
    integer :: its_Plus,npt_Plus
    real(real64) :: pt_Plus1
    integer :: its_Plus1,npt_Plus1
    real(real64) :: pt_Plus2
    integer :: its_Plus2,npt_Plus2

! MACHINE PARAMETERS
    real(real64) :: xq, xq_initial, step_size, ecm
    character(len=50)  :: pdf_name_NLO
    integer :: i,j,k,it_max

! External functions
  
end module globals_mod
!=================================================
