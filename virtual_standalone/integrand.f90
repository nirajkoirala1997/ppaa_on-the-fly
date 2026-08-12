  double precision function fnlo2(yy, vwgt) 
  use globals_mod, only : s,xq
  use kinematics_mod
  use misc_utils_mod
  use constants_mod
  use matrix_elements_mod

    ! Input/output
    real(8), intent(in)  :: yy(10)
    real(8), intent(in)  :: vwgt
    real(8)              :: wgt_result

    real(8) :: scale, AL

    ! Declare locals
    real(8) :: rs, xa, xb, rsp, xlow, xhigh
    real(8) :: xinvmass, scale_local, xmu2, xmur, xmuf
    real(8) :: sig, xnorm, wgt,alphasPDF
    real(8) :: s12
    integer :: ipass

    real(8) :: f1(-6:6), f2(-6:6)
    real(8) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3)
    real(8) :: amp(5)

!    external alphasPDF

    ! Assign input values
    xa  = yy(1)
    xb  = yy(2)
    rs  = sqrt(s)
    rsp = sqrt(xa * xb * s)

    ipass = 0
    xlow  = xq - eps
    xhigh = xq + eps

    ! Kinematics
    call kinvar2(yy,xinvmass, p1, p2, p3, p4)
    
    s12 = 2d0 * dot(p1, p2)
    scale_local = xinvmass
    scale = scale_local

    if (scale >= xlow .and. scale <= xhigh) then
      xmuf = scale
      xmur = scale
      xmu2 = xmuf ** 2

      call pdf(xa, xmuf, f1)
      call pdf(xb, xmuf, f2)
      AL = alphasPDF(xmur)
      call evaluate_1loop(AL,xmur,p1,p2,p3,p4,amp)

      !==================================================!
      sig = 0.0
      do i = 1, Nf 
      sig = sig + (( f1(-i)*f2(+i)+f1(+i)*f2(-i) ) * amp(i)) 
      enddo
      !==================================================!

      xnorm = hbarc2 / (16d0 * pi * xa * xb * s)
      wgt = xnorm * sig * vwgt
      wgt_result = wgt / vwgt / 2d0 / eps
    else
      wgt_result = 0d0
    end if
      fnlo2 = wgt_result

  end function fnlo2
!
!
!
!
!c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      function fnlo2(yy,vwgt)
!      use openloops
!      implicit double precision (a-h,o-z)
!      dimension yy(10)
!      dimension f1(-6:6),f2(-6:6),xl(15)
!      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
!     .          ,p_ex(0:3,1:4),c(1:2),Born(1:2),coef(2,-2:0),SumI(-2:0)
!     .          ,am211(0:2),am211_ir(0:2),amp2cc(6)
!c... this is for the evaluate_poles
!     .          ,xm2bare(0:2),xm2ct(0:2),xm2ir(0:2),xm2sum(0:2)
!       parameter (pi=3.14159265358979d0)
!c      parameter (hbarc2=0.3894d9)  ! in pb
!      parameter (hbarc2=0.3894d12)  ! in fb
!      parameter (EulerGamma = 0.5772156649015328606065120d0)
!
!      character*50 name
!      common/pdfname/name
!      common/leg_choice/leg
!      common/energy/s
!      common/distribution/xq
!      common/renor_scale/scale
!      common/prc_id/id_LO,id_NLO_1r,id_NLO_1loop
!      common/usedalpha/AL,ge
!      external Born_uU2eE
!       
!      rs  = dsqrt(s)
!      xa     = yy(1)
!      xb     = yy(2)
!
!      rsp = dsqrt(xa*xb*s)
!        
!      ipass = 0
!        eps = 0.5d0
!       xlow = xq - eps
!      xhigh = xq + eps
!
!      xcut = xq - 10.0d0
!
!c      if (rsp .gt. xcut) then
!
!
!        call kinvar2(yy,xinvmass,p1,p2,p3,p4)
!        s12 = 2d0*dot(p1,p2)
!         scale  = xinvmass
!c      if (scale .gt. 0d0 ) then
!        if ( scale .ge. xlow .and. scale .le. xhigh) then 
!             
!              xmuf=scale
!              xmur=scale
!              xmu2=xmuf**2
!
!              call pdf(xa,xmuf,f1)
!              call pdf(xb,xmuf,f2)
!              call setlum(f1,f2,xl)
!              AL = alphasPDF(xmur)
!              sig= xl(1)*Born_uU2eE(0,AL,scale,p1,p2,p3,p4)
!
!              xnorm=hbarc2/16d0/pi/(xa*xb*s)
!              wgt=xnorm*sig*vwgt
!              fnlo2=wgt/vwgt/2d0/eps
!       else                  
!        fnlo2=0d0
!       endif
!      return
!      end
!
!c---------------------------------------------------------------------
!
!

