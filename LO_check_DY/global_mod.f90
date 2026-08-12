!#################################################
module globals_mod
    use iso_fortran_env
    implicit none
    real(real64) :: s
    real(real64) :: AL, ge
    real(real64) :: xq

!MASS OF PARTICLES
    real(real64) :: m1, m2, m3, m4, m5

! OPENLOOOPS
    integer :: id_LO(5), id_NLO_1r(5), id_NLO_1loop(5)

! VEGAS PARAMETERS
    real(real64) :: pt_LO 
    integer :: its_LO,npt_LO

! MACHINE PARAMETERS
    real(real64) :: xq_initial, step_size, ecm
    character(len=50)  :: pdf_name_LO
    integer :: it_max
  
end module globals_mod
!#################################################
