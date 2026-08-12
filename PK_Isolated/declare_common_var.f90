    ! Declare locals   [ common for all plus functions ]
    real(real64) :: xa, xb, xt, x                           ! Integration variables
    real(real64) :: almin,almax,alpha
    real(real64) :: xmin,xmax

    real(real64) :: AL, ALP
    real(real64) :: s12,invmass, scale, xmur, xmuf, xmuf2  ! Scale variables
    real(real64) :: sp, rsp
    real(real64) :: xlow, xhigh                             ! Binning Variables
    real(real64) :: amp(5), SumPK(2),sig                    ! amplitudes
    real(real64) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3)      ! Phase space points
    real(real64) :: xnorm, pf, flux, ps2, azmth, pin        ! Overall Constants
    real(real64) :: xjac,xjac4,xtjac        !jacobians
    real(real64) :: f1(-6:6), f2(-6:6)
    real(real64) :: Pterm,Kterm,Tfactors

! External functions
    real(real64) :: alphasPDF


