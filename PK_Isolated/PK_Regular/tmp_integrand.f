c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Regular Terms]
      function flo2_PKReg(yy,vwgt)
      implicit double precision (a-h,o-z)
      dimension yy(10)
     &  ,f1(-6:6),f2(-6:6),xl(15)
     &  ,p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),xp1(0:3),xp2(0:3)
     &  ,p(0:3,1:4),Born(1:2),AllReg(1:2)
     &  ,SumPK1(1:2),SumPK2(1:2),SumPK3(1:2)


      double precision Pterm1,Kterm1
      double precision Pterm2,Kterm2
      double precision Pterm3,Kterm3
      parameter (pi=3.14159265358979d0)
      parameter (hbarc2=389.3856741D+9)

      character*10 P_,Kb_,Kt_
      character*10 gg_,qg_,gq_ 
      character*10 Plus_,Regular_,Delta_ 

      common/energy/s
c      common/factscale/xmuf
      common/usedalpha/AL,ge
      common/distribution/xq
      common/bin_size/eps
      common/amass/am1,am2,amH,am4,am5
      common/scales/xmuf,xmur
      external PK

      P_       ='P'
      Kb_      ='Kb'
      Kt_      ='Kt'
      gg_      ='gg'
      qg_      ='qg'
      gq_      ='gq'
      Plus_    ='Plus'
      Regular_ ='Regular'
      Delta_   ='Delta'

       tau = amH**2/S
c... for regular terms I don't need to use the upper cut as it is a flat integral.
      delta = 0d0 
c      xmin = 0d0 
      xmin = tau 

      xmax = 1.0d0 - delta
      xjac4 = (xmax - xmin)
      x = xjac4*yy(1) + xmin

c      xamin = 0d0 
      xamin = tau/x
      xamax = 1d0
      xajac =  (xamax - xamin)
      xa = xamin + xajac*yy(2)
      xb = tau/x/xa
      flo2_PKReg = 0d0
!      if (xb .ge. 1d0 ) goto 151
c... till here linear mapping is added
!
      Pterm1 = 0d0
      Kterm1 = 0d0
      sp = xa*xb*S
      rsp = dsqrt(sp)
      do k=1,2

        if ( k .eq. 1) call kinvar1(xa*x,xb,p1,p2,p3)
        if ( k .eq. 2) call kinvar1(xa,xb*x,p1,p2,p3)

      sp = 2d0*dot(p1,p2)
      rsp = dsqrt(sp)
c... same reasoning as plusA. refer there.

      AL = alphasPDF(xmur)
      ALP = AL/2d0/Pi
      CA = 3d0
      CF = 4d0/3d0

      Tb_dot_Ta_by_Ta2 = -CA/CA
      Tb_dot_Ta_by_Ta2_qg = -(CA - CF)/CA
      Tb_dot_Ta_by_Ta2_gq = -(CA - CF)/CA
      xmuf2 = xmuf*xmuf
      s12   = 2d0*dot(p1,p2)
      coef = Born_gg2h_(0,AL,p1,p2,p3)

c... gg contribution
       Tfactors = Tb_dot_Ta_by_Ta2
      Tfactors2 = Tb_dot_Ta_by_Ta2_qg
      Tfactors3 = Tb_dot_Ta_by_Ta2_gq

      Pterm1 =ALP*PK(P_,gg_,Regular_,x)*Tfactors*dlog(xmuf2/s12)
      Kterm1 =ALP*PK(Kb_,gg_,Regular_,x) -
     &  ALP*Tfactors*PK(Kt_,gg_,Regular_,x)

c... qg contribution
      Pterm2 =ALP*PK(P_,qg_,Regular_,x)*Tfactors2*dlog(xmuf2/s12)
      Kterm2 =ALP*PK(Kb_,qg_,Regular_,x) -
     &  ALP*Tfactors2*PK(Kt_,qg_,Regular_,x)

c... gq contribution
      Pterm3 =ALP*PK(P_,gq_,Regular_,x)*Tfactors3*dlog(xmuf2/s12)
      Kterm3 =ALP*PK(Kb_,gq_,Regular_,x) -
     &  ALP*Tfactors3*PK(Kt_,gq_,Regular_,x)

      SumPK1(k) = coef*(Pterm1 + Kterm1)
      SumPK2(k) = coef*(Pterm2 + Kterm2)
      SumPK3(k) = coef*(Pterm3 + Kterm3)
      enddo     

      call pdf(xa,xmuf,f1)
      call pdf(xb,xmuf,f2)
      call setlum(f1,f2,xl)
  
      sig_gg = xl(2)*(SumPK1(1) + SumPK1(2))
      sig_qg = xl(7)*(SumPK2(1) + SumPK2(2))
      sig_gq = xl(8)*(SumPK3(1) + SumPK3(2))

      sig = sig_gg + sig_qg + sig_gq

      pi_1 = PI/amH 
      flux = 2d0*sp
      xnorm = hbarc2*xajac*xjac4*pi_1/flux
      pf_factor = 2d0*amH/xa/x/S

      flo2_PKReg = xnorm*sig*pf_factor
151   return
      end
