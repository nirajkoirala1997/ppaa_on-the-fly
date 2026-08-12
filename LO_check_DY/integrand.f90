  double precision function flo2_LO(yy, vwgt) 
  use globals_mod, only : s,xq
  use kinematics_mod
  use misc_utils_mod
  use constants_mod, pi_mod => PI
  use matrix_elements_mod

  use drivekk
  use kinparams
  use psnrmtest
  use ntdisc
  use fourmom

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
    double precision, dimension(3) :: x
    real(8) :: flux,dlips

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
    x(1) = yy(1)
    x(2) = yy(2)
    x(3) = yy(5)
!    x(3) = yy(3)
!    x(4) = yy(4)
    xa  = 1d0 !yy(1)
    xb  = 1d0 !yy(2)
    xa  = yy(3)
    xb  = yy(4)
    call psgen(1,sqrt(s),xa,xb,x,3,flux,dlips)
!if (dlips > 0.0 ) print*,flux,dlips
                  print*,flux,dlips
	do j=1,4
	print*,p(0,j),p(1,j),p(2,j),p(3,j)
	enddo
	stop
    
    s12 = 2d0 * dot(p1, p2)
    scale_local = xinvmass
    scale = scale_local

    if (scale >= xlow .and. scale <= xhigh) then
      scale = xq 
      xmuf = scale
      xmur = scale
      xmu2 = xmuf ** 2

      call pdf(xa, xmuf, f1)
      call pdf(xb, xmuf, f2)
      AL = alphasPDF(xmur)
      call evaluate_born(p1,p2,p3,p4,amp)

     sig =     (f1(-1)*f2(+1)+f1(+1)*f2(-1))*amp(1)  & 
!     sig =        + (f1(-2)*f2(+2)+f1(+2)*f2(-2))*amp(2)  &
             + (f1(-3)*f2(+3)+f1(+3)*f2(-3))*amp(3)  
!             + (f1(-4)*f2(+4)+f1(+4)*f2(-4))*amp(4)   
!             + (f1(-5)*f2(+5)+f1(+5)*f2(-5))*amp(5)    

      xnorm = hbarc2 / (16d0 * pi_mod * xa * xb * s)
      wgt = xnorm * sig * vwgt
!      wgt_result = wgt / vwgt / 2d0 / eps
      wgt_result = wgt / vwgt    ! in matrix overview result is not divided by bin size. 
    else
      wgt_result = 0d0
    end if
      flo2_LO = wgt_result

  end function flo2_LO
!
!
!
!
!c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      function flo2_LO(yy,vwgt)
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
!              flo2_LO=wgt/vwgt/2d0/eps
!       else                  
!        flo2_LO=0d0
!       endif
!      return
!      end
!
!c---------------------------------------------------------------------
!
!

