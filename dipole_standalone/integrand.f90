      double precision function fnlo3(xx,weight)
      use globals_mod
      use kinematics_mod
      use misc_utils_mod, only : dot
      use constants_mod
      use matrix_elements_mod
      use dipoles

     !Input / Output
      real(8), intent(in) :: xx(10)
      real(8), intent(in) :: weight

     ! Declare locals
      real(8) :: rs, xa, xb, rsp,sp, xlow, xhigh
      real(8) :: xinvmass, scale_local, xmur, xmuf
      real(8) :: sig, xnorm, wgt,alphasPDF
      real(8) :: xxjac,pi_1,flux
      integer :: unphy
  
      real(8) :: f1(-6:6), f2(-6:6)
      real(8) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3), ptilde(0:3,4,2),p_ex(0:3,5)
      real(8) :: amp(5),AL,dip(5,2)

      xa = xx(1)
      xb = xx(2)
      rs = dsqrt(s)
      sp = xa*xb*s
      rsp = dsqrt(sp)

      xlow = xq - eps
      xhigh = xq + eps
      fnlo3 = 0.0_real64

       call kinvar3(xx,xxjac,xinvmass,p1,p2,p3,p4,p5,unphy)
       call reduceps(p1,p2,p3,p4,p5,ptilde)

       if (unphy == 0) then ! with zero unphysical PS points proceed
     
        scale_local = xinvmass

        if ( scale_local >= xlow .and. scale_local <= xhigh) then

          xmuf = scale_local
          xmur = scale_local
           AL  = alphasPDF(xmur)


          call pdf(xa,xmuf,f1)
          call pdf(xb,xmuf,f2)
          call evaluate_1real(AL,p1,p2,p3,p4,p5,amp)
          call p1d_to_p2d_5(p1,p2,p3,p4,p5,p_ex)
          call dipole_initial_initial(AL,p_ex ,ptilde,dip)

            !==================================================!
      sig = 0.0_real64
      do i=1 ,NF
      sig = sig + ( f1(-i)*f2(+i)+f1(+i)*f2(-i) )*( amp(i) - dip(i,1) - dip(i,2) ) 
      enddo
            !==================================================!

          pi_1 = 0.5d0*rsp
          flux = 4d0*pi_1*rsp
          xnorm=hbarc2/8d0/(2d0*Pi)**4/flux
          wgt=xxjac*xnorm*sig*weight
          fnlo3=wgt/weight/2.0/eps

         endif
       endif
      end
!---------------------------------------------------------------------
