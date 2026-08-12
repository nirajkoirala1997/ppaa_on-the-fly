double precision function xint_PlusA(yy,vwgt)
    use globals_mod
    use kinematics_mod
    use misc_utils_mod
    use constants_mod
    use matrix_elements_mod
    use PK_functions_mod

    ! Input/output
    real(real64), intent(in)  :: yy(10)
    real(real64), intent(in)  :: vwgt

    include '../declare_common_var.f90'

    ! Assign input values
    xa  = yy(1)
    xb  = yy(2)
    xt  = yy(3)

!...linear maaping 
!    xmin = 0d0 
!    xmax = 1d0 - delta_cut
!    xjac4 = xmax - xmin
!    x = xjac4 * yy(4) + xmin

! exponential mapping
    almin = delta_cut
    almax = 1.0d0
    alpha = almin*(almax/almin)**yy(4)
    xjac4 = alpha*dlog(almax/almin)
    x = 1.0d0 - alpha

    xtjac = 2.0
    xjac = xjac4*xtjac 
    xnorm = hbarc2

    xlow  = xq - eps
    xhigh = xq + eps

    do k = 1, 2
        if (k == 1) then
          call kinvar2_PK(x * xa, xb, xt, invmass, p1, p2, p3, p4)
        else
          call kinvar2_PK(xa, x * xb, xt, invmass, p1, p2, p3, p4)
        end if
        scale = invmass
    
        if ( scale >= xlow .and. scale <= xhigh) then
    
         s12  = 2.0 * dot(p1, p2)
         sp   = s12 
         rsp  = sqrt(sp)

    
         xmuf   = scale
         xmur   = xmuf
         xmuf2  = xmuf * xmuf
         AL     = alphasPDF(xmur)
         ALP    = AL / (2.0 * Pi)
    
    
         call pdf(xa, xmuf, f1)
         call pdf(xb, xmuf, f2)
    
         Tfactors = -1.0
         Pterm = ALP * PK('P','qq','Plus', x) * Tfactors * dlog(xmuf2 / s12)
         Kterm = ALP * PK('Kb','qq','Plus', x) - ALP * Tfactors * PK('Kt','qq','Plus',x)
    
         call evaluate_born(p1, p2, p3, p4, amp)

         !==================================================!
         sig = 0.0
         do i = 1, Nf 
         sig = sig + (( f1(-i)*f2(+i)+f1(+i)*f2(-i) ) * amp(i)) 
         enddo
         sig = sig * (Pterm + Kterm)
         !==================================================!


         pin  = 0.5 * rsp
         pf   = 0.5 * rsp
         flux = 4.0 * pin * rsp

         azmth  = 2.00 * Pi
         ps2    = 1.00 / (4.00 * Pi * Pi) * (pf / (4.0 * rsp)) * azmth

         SumPK(k)= xjac*xnorm*sig*ps2/flux/2.0/eps

        else
         SumPK = 0.0
        end if
      end do
    
      xint_PlusA = SumPK(1) + SumPK(2)
end function xint_PlusA

!===================================================!
double precision function xint_Plus_h(yy, vwgt) 
  use globals_mod
  use kinematics_mod
  use misc_utils_mod
  use constants_mod
  use matrix_elements_mod
  use PK_functions_mod
  implicit none

  ! Arguments
  real(real64), intent(in)  :: yy(10), vwgt

  include '../declare_common_var.f90'

  ! Local variables
  real(real64) :: Pplus, xKbar_h, xKt_h, SumPlus
  real(real64) :: zero, one, one_m_del

  xlow  = xq - eps
  xhigh = xq + eps

  xa    = yy(1) 
  xb    = yy(2) 
  xt    = yy(3)
  xjac  = 2.0

  sp  = xa * xb * S
  rsp = sqrt(sp)

  call kinvar2_PK(xa, xb, xt, invmass, p1, p2, p3, p4)
  scale = invmass

 if (scale >= xlow .and. scale <= xhigh ) then 

  s12 = 2.0*dot(p1,p2)
  xmuf = scale 
  xmuf2 = xmuf**2
  xmur = xmuf 
  AL = alphasPDF(xmur)
  ALP = AL/2d0/Pi

  call pdf(xa,xmuf,f1)
  call pdf(xb,xmuf,f2)

  zero       = 0d0
  one        = 1d0
  one_m_del  = one - delta_cut

  Pplus   = -dlog(xmuf2 / s12)

  Pplus   = Pplus * (xint_h(zero,'Pqq')          - xint_h(one_m_del,'Pqq'))
  xKbar_h =          xint_h(zero,'Kbarqq')       - xint_h(one_m_del,'Kbarqq')
  xKt_h   =          xint_h(zero,'Ktqq')         - xint_h(one_m_del,'Ktqq')

  SumPlus = ALP*(Pplus + xKbar_h + xKt_h)

  call evaluate_born(p1, p2, p3, p4, amp)
         
  !==================================================!
  sig = 0.0
  do i = 1, Nf 
  sig = sig + (( f1(-i)*f2(+i)+f1(+i)*f2(-i) ) * amp(i)*SumPlus) 
  enddo
  sig = sig*2.0 ! for 2 legs 
  !==================================================!


  azmth = 2.0d0*pi
  pf    = 0.5d0*rsp
  pin   = 0.5d0*rsp
  ps2   = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth
  xnorm = hbarc2

  flux = 4.0d0*pin*rsp

  xint_Plus_h = xnorm*sig/flux*ps2*xjac/2.0/eps
 else
 xint_Plus_h = 0d0
 endif

end function xint_Plus_h

 ! this routine can be turned on for specific Q integration
!!===================================================!
!double precision function xint_Plus_h(yy, vwgt) 
!  use globals_mod
!  use kinematics_mod
!  use misc_utils_mod
!  use constants_mod
!  use matrix_elements_mod
!  use PK_functions_mod
!  implicit none
!
!  ! Arguments
!  real(real64), intent(in)  :: yy(10), vwgt
!  real(real64)              :: res
!
!  ! Variables
!  real(real64) :: xa, xb, xajac, tau, xjac
!  real(real64) :: sp, rsp, xmuf2
!  real(real64) :: AL, ALP, s12, coef, flux
!  real(real64) :: Pplus, xKbar_h, xKt_h, SumPlus
!  real(real64) :: sig, pi_1, flux_1, xnorm
!  real(real64) :: zero, one, one_m_del
!
!  ! Particle momenta and PDF-related variables
!  real(real64) :: p1(0:3), p2(0:3), p3(0:3),p4(0:3)
!  real(real64) :: f1(-6:6), f2(-6:6), amp(5) 
!  real(real64) :: ps2, pin, p(0:3,4), invmass, wgt1, xc, azmth, pf
!
!  ! Common variables from modules
!  real(real64) :: xmuf, xmur, eps, scale
!  double precision alphasPDF
!  external alphasPDF
!
!  ! --- Phase space setup ---
!  tau   = xq*xq / S
!  xajac = 1d0 - tau
!  xa    = tau + xajac * yy(1)
!  xb    = tau / xa
!  xc    = yy(2)
!  xjac  = xajac*2.0
!
!  sp  = xa * xb * S
!  rsp = sqrt(sp)
!
!  call kinvar2_PK(xa, xb, xc, invmass, p1, p2, p3, p4)
!
!  s12 = 2.0*dot(p1,p2)
!  xmuf = invmass 
!  xmuf2 = xmuf**2
!  xmur = xmuf 
!  AL = alphasPDF(xmur)
!  ALP = AL/2d0/Pi
!
!  call pdf(xa,xmuf,f1)
!  call pdf(xb,xmuf,f2)
!
!  zero       = 0d0
!  one        = 1d0
!  one_m_del  = one - delta_cut
!
!  Pplus   = -dlog(xmuf2 / s12)
!
!  Pplus   = Pplus * (xint_h(zero,'Pqq')          - xint_h(one_m_del,'Pqq'))
!  xKbar_h =          xint_h(zero,'Kbarqq')       - xint_h(one_m_del,'Kbarqq')
!  xKt_h   =          xint_h(zero,'Ktqq')         - xint_h(one_m_del,'Ktqq')
!
!  SumPlus = ALP*(Pplus + xKbar_h + xKt_h)
!
!  call evaluate_born(p1, p2, p3, p4, amp)
!         
!  !==================================================!
!  sig = 0.0
!  do i = 1, Nf 
!  sig = sig + (( f1(-i)*f2(+i)+f1(+i)*f2(-i) ) * amp(i)*SumPlus) 
!  enddo
!  sig = sig*2.0 ! for 2 legs 
!  !==================================================!
!
!
!  azmth = 2.0d0*pi
!  pf = 0.5d0*rsp
!  ps2 = 1.0d0/(4.d0*pi*pi)*(pf/4.d0/rsp)*azmth
!  xnorm = hbarc2*2.0*xq/xa/s
!
!  pin  = 0.5d0*rsp
!  flux_1 = 4.0d0*pin*rsp
!
!  xint_Plus_h = xnorm*sig/flux_1*ps2*xjac
!
!end function xint_Plus_h

  !=============[  end  ]============!

double precision function flo2_PlusA(yy, vwgt)
  use globals_mod
  use kinematics_mod
  use misc_utils_mod
  use constants_mod
  use matrix_elements_mod
  use PK_functions_mod
  implicit none

  ! Inputs
  real(8), intent(in) :: yy(10)
  real(8), intent(in) :: vwgt

  include '../declare_common_var.f90'

! Local variables
  real(8) :: PKplus_x,wgt

   xlow = xq - eps
  xhigh = xq + eps

  xa = yy(1)
  xb = yy(2)
  xt = yy(3)

  almin = delta_cut
  almax = 1.0d0
  alpha = almin * (almax / almin) ** yy(4)
  xjac4 = alpha * log(almax / almin)
  x = 1.0d0 - alpha
  xtjac = 2.0d0
  xjac = xtjac * xjac4

  xnorm = hbarc2
  PKplus_x = 0d0
  SumPK = 0.0

  do k = 1, 2
    if (k == 1) then
      call kinvar2_PK(x * xa, xb, xt, invmass, p1, p2, p3, p4)
    else
      call kinvar2_PK(xa, x * xb, xt, invmass, p1, p2, p3, p4)
    end if
    scale = invmass

    if (scale >= xlow .and. scale <= xhigh) then
      sp = 2.0d0 * dot(p1, p2)
      rsp = sqrt(sp)
      s12 = sp
      pin = 0.5d0 * rsp
      pf = 0.5d0 * rsp
      flux = 4.0d0 * pin * rsp

      xmuf = scale
      xmur = scale
      xmuf2 = xmuf * xmuf
      AL = alphasPDF(xmur)
      ALP = AL / (2d0 * pi)

      azmth = 2.0d0 * pi
      ps2 = 1.0d0 / (4d0 * pi * pi) * (pf / (4d0 * rsp)) * azmth

      call pdf(xa, xmuf, f1)
      call pdf(xb, xmuf, f2)

      Tfactors = -1.0d0
      Pterm = ALP * PK('P', 'qq', 'Plus', x) * Tfactors * dlog(xmuf2 / s12)
      Kterm = ALP * PK('Kb', 'qq', 'Plus', x) - ALP * Tfactors * PK('Kt', 'qq', 'Plus', x)

      call evaluate_born(p1, p2, p3, p4, amp)

      sig = 0d0
      do i = 1, Nf
        sig = sig + (f1(-i) * f2(+i) + f1(+i) * f2(-i)) * amp(i)
      end do
      sig = sig * (Pterm + Kterm)

      wgt = sig / flux * ps2 * xjac * vwgt
      PKplus_x = xnorm * wgt / vwgt / 2.0 / eps

      SumPK(k) = PKplus_x
    end if
  end do

  flo2_PlusA = (SumPK(1) + SumPK(2))
end function flo2_PlusA

double precision function flo2_PlusB(yy, vwgt)
  use globals_mod
  use kinematics_mod
  use misc_utils_mod
  use constants_mod
  use matrix_elements_mod
  use PK_functions_mod
  implicit none

  ! Inputs
  real(8), intent(in) :: yy(10)
  real(8), intent(in) :: vwgt

  include '../declare_common_var.f90'

! Local variables
  real(8) :: PKplus_x, wgt

  xlow = xq - eps
  xhigh = xq + eps

  xa = yy(1)
  xb = yy(2)
  xt = yy(3)

  almin = delta_cut
  almax = 1.0d0
  alpha = almin * (almax / almin) ** yy(4)
  xjac4 = alpha * log(almax / almin)
  x = 1.0d0 - alpha
  xtjac = 2.0d0
  xjac = xtjac * xjac4

  xnorm = hbarc2
  PKplus_x = 0d0
  SumPK = 0.0

  do k = 1, 2
    call kinvar2_PK(xa, xb, xt, invmass, p1, p2, p3, p4)
    scale = invmass

    if (scale >= xlow .and. scale <= xhigh) then
      sp = 2.0d0 * dot(p1, p2)
      rsp = sqrt(sp)
      s12 = sp
      pin = 0.5d0 * rsp
      pf = 0.5d0 * rsp
      flux = 4.0d0 * pin * rsp

      xmuf = scale
      xmur = scale
      xmuf2 = xmuf * xmuf
      AL = alphasPDF(xmur)
      ALP = AL / (2d0 * pi)

      azmth = 2.0d0 * pi
      ps2 = 1.0d0 / (4d0 * pi * pi) * (pf / (4d0 * rsp)) * azmth

      call pdf(xa, xmuf, f1)
      call pdf(xb, xmuf, f2)

      Tfactors = -1.0d0
      Pterm = ALP * PK('P', 'qq', 'Plus', x) * Tfactors * dlog(xmuf2 / s12)
      Kterm = ALP * PK('Kb', 'qq', 'Plus', x) - ALP * Tfactors * PK('Kt', 'qq', 'Plus', x)

      call evaluate_born(p1, p2, p3, p4, amp)

      sig = 0d0
      do i = 1, Nf
        sig = sig + (f1(-i) * f2(+i) + f1(+i) * f2(-i)) * amp(i)
      end do
      sig = sig * (Pterm + Kterm)

      wgt = sig / flux * ps2 * xjac * vwgt
      PKplus_x = xnorm * wgt / vwgt / 2.0 / eps

      SumPK(k) = PKplus_x
    end if
  end do

  flo2_PlusB = (SumPK(1) + SumPK(2))
end function flo2_PlusB
