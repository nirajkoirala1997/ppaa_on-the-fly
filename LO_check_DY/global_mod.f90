!#################################################
! Run-wide state for runLO. Kept deliberately small: everything the run card
! supplies lives in misc_utils_mod, everything physical in constants_mod. What
! is left here is what the integrand needs and what the shared modules expect
! to find (mat.amp.f90 wants ge and the id_* arrays, phasespace.f90 wants s).
module globals_mod
    use iso_fortran_env
    implicit none

! HADRONIC CMS ENERGY SQUARED
    real(real64) :: s

! ELECTROWEAK COUPLING (alpha, in the Gfermi scheme)
    real(real64) :: ge

! INVARIANT-MASS WINDOW, SQUARED.
! Applied as a cut on m(e+e-) reconstructed from the generated momenta, never
! on x1*x2. Those are the same thing only at LO, where shat = m(ee)^2 exactly;
! with a gluon in the final state shat >= m(ee)^2 and a window on x1*x2 would
! throw away real-emission configurations that belong in the cross section.
    real(real64) :: mee_min2, mee_max2

! FIXED FACTORIZATION AND RENORMALIZATION SCALES
    real(real64) :: xmur, xmuf

! OPENLOOPS PROCESS IDS
    integer :: id_LO(5), id_NLO_1r(5), id_NLO_1loop(5)

end module globals_mod
!#################################################
