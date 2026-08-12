!-------------------------------------------------------------!
double precision function fnlo2_PKReg(yy,vwgt)
    use globals_mod, only : s,xq,delta
    use kinematics_mod
    use misc_utils_mod
    use constants_mod
    use matrix_elements_mod
    use PK_functions_mod

    ! Input/output
    real(8), intent(in)  :: yy(10)
    real(8), intent(in)  :: vwgt
    real(8) ::  AL

    ! Declare locals
    real(8) :: s12, sp, pin,xmuf2, rs, xa, xb, xt, rsp, xlow, xhigh
    real(8) :: xinvmass, scale_local, xmu2, xmur, xmuf
    real(8) :: sig_qq,sig_gq, xnorm, alphasPDF
    real(8) :: SumPlus,x 
    integer :: ipass

    real(8) :: f1(-6:6), f2(-6:6)
    real(8) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), SumPK(2)
    real(8) :: amp(5)
    real(8) :: Pterm_qq,Kterm_qq,Tfactors
    real(8) :: Pterm_gq,Kterm_gq
    real(8) :: ALP,scale,pf,flux,ps2,azmth, Qmass 
    real(real64) :: almin,almax,alpha,xjac4,xtjac

    ! Assign input values
    xa  = yy(1)
    xb  = yy(2)
    xt  = yy(3)
    x   = yy(4)

    xjac = 2.0 
    xnorm = hbarc2

    rs  = sqrt(s)
    rsp = sqrt(xa * xb * s)

    xlow  = xq - eps
    xhigh = xq + eps

    do k = 1, 2
        if (k == 1) then
          call kinvar2_PK(x * xa, xb, xt, Qmass, p1, p2, p3, p4)
        elseif (k ==2) then
          call kinvar2_PK(xa, x * xb, xt, Qmass, p1, p2, p3, p4)
        end if
        scale = Qmass
    
        if (scale >= xlow .and. scale <= xhigh) then
    
         sp   = 2.0 * dot(p1, p2)
         s12  = sp
         rsp  = sqrt(sp)
         pin  = 0.5 * rsp
         pf   = 0.5 * rsp
         flux = 4.0 * pin * rsp
    
         xmuf   = scale
         xmur   = xmuf
         xmuf2  = xmuf * xmuf
         AL     = alphasPDF(xmur)
         ALP    = AL / (2.0 * Pi)
    
         azmth  = 2.00 * Pi
         ps2    = 1.00 / (4.00 * Pi * Pi) * (pf / (4.0 * rsp)) * azmth
    
         call pdf(xa, xmuf, f1)
         call pdf(xb, xmuf, f2)
         Tfactors = -1.0
         Pterm_qq = ALP * PK('P','qq','Regular', x) * Tfactors * log(xmuf2 / s12)
         Kterm_qq = ALP * PK('Kb','qq','Regular', x) - ALP * Tfactors * PK('Kt','qq','Regular',x)


         Pterm_gq = ALP * PK('P','gq','Regular', x) * Tfactors * log(xmuf2 / s12)
         Kterm_gq = ALP * PK('Kb','gq','Regular', x) - ALP * Tfactors * PK('Kt','gq','Regular',x)
    
         call evaluate_born(p1, p2, p3, p4, amp)
         

         !==================================================!
         sig_qq = 0.0
         do i = 1, Nf 
         sig_qq = sig_qq + ( f1(-i)*f2(+i)+f1(+i)*f2(-i) ) * amp(i)
         enddo
         sig_qq = sig_qq * (Pterm_qq + Kterm_qq)
         !==================================================!
         !==================================================!
         sig_gq = 0.0
         do i = 1, Nf 
         sig_gq = sig_gq + ( f1(0)*f2(+i)+f1(0)*f2(-i) ) * amp(i)
         enddo
         sig_gq = sig_gq * (Pterm_gq + Kterm_gq)
         !==================================================!


         SumPK(k)= xnorm * (sig_qq + 1d0*sig_gq)  * ps2 * xjac /flux/2.0/eps 
!         SumPK(k)= xnorm * sig_qq  * ps2 * xjac /flux/2.0/eps 
        else
        SumPK = 0.0
        end if
      
      end do
    
      fnlo2_PKReg = SumPK(1) + SumPK(2)
end function fnlo2_PKReg
