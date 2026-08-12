!-------------------------------------------------------------!
double precision function fnlo2_PKDel(yy,vwgt)
    use globals_mod, only : s,xq,delta
    use kinematics_mod
    use misc_utils_mod
    use constants_mod
    use matrix_elements_mod
    use PK_functions_mod

    ! Input/output
    real(8), intent(in)  :: yy(10)
    real(8), intent(in)  :: vwgt
    real(8)              :: wgt_result

    real(8) ::  AL

    ! Declare locals
    real(8) :: s12, sp, pin,xmuf2, rs, xa, xb, xt, rsp, xlow, xhigh
    real(8) :: xinvmass, scale_local, xmu2, xmur, xmuf
    real(8) :: sig, xnorm, wgt,alphasPDF
    real(8) :: SumPlus,x 
    integer :: ipass

    real(8) :: f1(-6:6), f2(-6:6)
    real(8) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), SumPK(2)
    real(8) :: amp(5)
    real(8) :: Pterm,Kterm,Tfactors
    real(8) :: ALP,scale,pf,flux,ps2,azmth, Qmass 
    real(real64) :: almin,almax,alpha,xjac4,xtjac

    xlow = xq - eps
    xhigh = xq + eps

   
    xa  = yy(1)
    xb  = yy(2)
    xt  = yy(3)

    xjac = 2.0
    xnorm = hbarc2

    rs  = sqrt(s)
    sp  = xa * xb * s
    rsp = sqrt(sp)


        call kinvar2_PK(xa, xb, xt, Qmass, p1, p2, p3, p4)
        scale  = Qmass

        if ( scale >= xlow .and. scale <= xhigh ) then 
         
    
         s12  = 2.0*dot(p1,p2) 
         rsp  = sqrt(sp)
         pin  = 0.5 * rsp
         pf   = 0.5 * rsp
         flux = 4.0 * pin * rsp
    
         xmuf   = Qmass
         xmur   = xmuf
         xmuf2  = xmuf * xmuf
         AL     = alphasPDF(xmur)
         ALP    = AL / (2.0 * Pi)
    
         azmth  = 2.00 * Pi
         ps2    = 1.00 / (4.00 * Pi * Pi) * (pf / (4.0 * rsp)) * azmth
    
         call pdf(xa, xmuf, f1)
         call pdf(xb, xmuf, f2)
    
         Tfactors = -1.0
         Pterm = ALP * PK('P','qq','Delta', x) * Tfactors * log(xmuf2 / s12)
         Kterm = ALP * PK('Kb','qq','Delta', x) - ALP * Tfactors * PK('Kt','qq','Delta',x)
    
         call evaluate_born(p1, p2, p3, p4, amp)
         
         !==================================================!
         sig = 0.0
         do i = 1, Nf 
         sig = sig + (( f1(-i)*f2(+i)+f1(+i)*f2(-i) ) * amp(i)) 
         enddo
         sig = 2.0 * sig * (Pterm + Kterm)
         !==================================================!

      fnlo2_PKDel = xnorm * (sig / flux * ps2 * xjac * vwgt) / vwgt / 2.0/ eps
         else
       fnlo2_PKDel = 0d0
         endif

end function fnlo2_PKDel!-------------------------------------------------------------!
!double precision function fnlo2_PKDel(yy,vwgt)
!    use globals_mod, only : s,xq,delta
!    use kinematics_mod
!    use misc_utils_mod
!    use constants_mod
!    use matrix_elements_mod
!    use PK_functions_mod
!
!    ! Input/output
!    real(8), intent(in)  :: yy(10)
!    real(8), intent(in)  :: vwgt
!    real(8)              :: wgt_result
!
!    real(8) ::  AL
!
!    ! Declare locals
!    real(8) :: s12, sp, pin,xmuf2, rs, xa, xb, xt, rsp, eps, xlow, xhigh
!    real(8) :: xinvmass, scale_local, xmu2, xmur, xmuf
!    real(8) :: sig, xnorm, wgt,alphasPDF
!    real(8) :: SumPlus,x 
!    integer :: ipass
!
!    real(8) :: f1(-6:6), f2(-6:6)
!    real(8) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), SumPK(2)
!    real(8) :: amp(5)
!    real(8) :: Pterm,Kterm,Tfactors
!    real(8) :: ALP,scale,pf,flux,ps2,azmth, Qmass 
!    real(real64) :: almin,almax,alpha,xjac4,xtjac
!
!    ! Assign input values
!    tau = xq*xq/s
!   
!    xa  = (1.0-tau)*yy(1)+tau
!    xb  = tau/xa 
!    xt  = yy(2)
!
!    xjac = 2.0*(1d0 - tau)
!    xnorm = hbarc2*2.0*xq/s/xa
!
!    rs  = sqrt(s)
!    sp  = xa * xb * s
!    rsp = sqrt(sp)
!
!
!        call kinvar2_PK(xa, xb, xt, Qmass, p1, p2, p3, p4)
!    
!         s12  = 2.0*dot(p1,p2) 
!         rsp  = sqrt(sp)
!         pin  = 0.5 * rsp
!         pf   = 0.5 * rsp
!         flux = 4.0 * pin * rsp
!    
!         xmuf   = Qmass
!         xmur   = xmuf
!         xmuf2  = xmuf * xmuf
!         AL     = alphasPDF(xmur)
!         ALP    = AL / (2.0 * Pi)
!    
!         azmth  = 2.00 * Pi
!         ps2    = 1.00 / (4.00 * Pi * Pi) * (pf / (4.0 * rsp)) * azmth
!    
!         call pdf(xa, xmuf, f1)
!         call pdf(xb, xmuf, f2)
!    
!         Tfactors = -1.0
!         Pterm = ALP * PK('P','qq','Delta', x) * Tfactors * log(xmuf2 / s12)
!         Kterm = ALP * PK('Kb','qq','Delta', x) - ALP * Tfactors * PK('Kt','qq','Delta',x)
!    
!         call evaluate_born(p1, p2, p3, p4, amp)
!         
!         !==================================================!
!         sig = 0.0
!         do i = 1, Nf 
!         sig = sig + (( f1(-i)*f2(+i)+f1(+i)*f2(-i) ) * amp(i)) 
!         enddo
!         sig = 2.0 * sig * (Pterm + Kterm)
!         !==================================================!
!
!      fnlo2_PKDel = xnorm * (sig / flux * ps2 * xjac * vwgt) / vwgt
!
!end function fnlo2_PKDel
!
!
!
!
!
!
!
!
!
!
