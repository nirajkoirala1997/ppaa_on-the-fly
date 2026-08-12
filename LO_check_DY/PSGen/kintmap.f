      subroutine kintmap(imom,ctc,ecm,mt,m1,m3,m4,x,p1,p3,p4,dlips)
c----------------------------------------------------------------------c
c Returns cartezian components of the final state four momenta in the 
c c.m.s. and the phase space normalization for a 2 --> 2 t-channel 
c process. Mapping of ~1/t behaviour of the matrix element is performed.
c The factor (2*pi)^-2 is not included in dlips.
c
c Input: imom=1/(else) => four momenta p3,p4 are/(not) calculated,
c        ctc - angular cut on cosine of (p1,p3) angle,
c        ecm - c.m.s. energy,
c        mt - mas of a particle in the t-channel propagator,
c        m1,m3,m4 - external particle masses,
c        x(i) - random variables, 0 < x(i) < 1, i=1,2,
c        p1 - initial state particle four momentum.
c Output: p3,p4 - final state particle four momenta.
c
c Karol Kolodziej, University of Silesia, kolodzie@us.edu.pl
c 21 May, 2009; modified on 24 May, 2009
c----------------------------------------------------------------------c

      implicit none

      integer, parameter :: dp=kind(1.d0)
      integer :: imom
      real(dp) :: ctc,ecm,mt,m1,m3,m4,x(2),p1(0:3),p3(0:3),p4(0:3),dlips
      real(dp), parameter :: pi=3.141592653589793238d0,
     &                       dphi=2.d0*pi,dct=2.d0,cps=dphi/8.d0
      real(dp) :: delme,pp,ct,pst,phi,e3,e4,ap1,ep1,ep3,mep1,mep3,
     &            delc,tc,t1,t2,dt,t,mt2

      delme=m3+m4

      if (ecm < delme) then
       dlips=0.d0
      else
       delme=(m3-m4)*delme/ecm
       if (delme == 0.d0) then
        e3=0.5d0*ecm
        e4=e3
        pp=sqrt((e3-m3)*(e3+m3))
       else
        e3=0.5d0*(ecm+delme)
        e4=e3-delme
        if (e3 < m3 .or. e4 < m4) then 
         dlips=0.d0
         return
        end if
        if (delme > 0.d0) then
         pp=sqrt((e3-m3)*(e3+m3))
        else
         pp=sqrt((e4-m4)*(e4+m4))
        end if
       end if

       ap1=abs(p1(3))
       ep1=p1(0)+ap1
       ep3=e3+pp
       mep1=m1**2/ep1
       mep3=m3**2/ep3
c Limits in t:
       t1=(mep1-mep3)*(ep1-ep3)   ! Upper limit
       delc=1.d0-ctc
       tc=(p1(0)-e3-ap1+pp)*(p1(0)-e3+ap1-pp)-2.d0*ap1*pp*delc
       if (tc > 0.d0) then
        write(*,*)'kintmap: tc > 0, tc =',tc,' put 0'
        tc=0.d0
       end if
       if (t1 > 0.d0) then
        t1=tc
       else if (t1 > tc) then
        t1=tc
       end if
       t2=-(ep3-mep1)*(ep1-mep3) ! Lower limit

       mt2=mt*mt

       if (imom == 1) then
c Mapping of ~1/t:
c        t=t2*(t1/t2)**x(1)
c        dt=t*log(t1/t2)
        t=mt2+(t2-mt2)*((t1-mt2)/(t2-mt2))**x(1)
        dt=(t-mt2)*log((t1-mt2)/(t2-mt2))
c        dt=t2-t1
c        t=dt*x(1)+t1
        ct=(t+2.d0*p1(0)*e3-m1**2-m3**2)/(2.d0*p1(3)*pp)
        pst=pp*sqrt(1.d0-ct**2)
        phi=dphi*x(2)
        p3=(/e3,pst*cos(phi),pst*sin(phi),pp*ct/)
        p4=(/e4,-p3(1),-p3(2),-p3(3)/)
       else
        t=(p1(0)-p3(0))**2-(p1(1)-p3(1))**2-(p1(2)-p3(2))**2
     &   -(p1(3)-p3(3))**2
c        dt=t*log(t1/t2)
        dt=(t-mt2)*log((t1-mt2)/(t2-mt2))
c        dt=t2-t1
       end if

       dlips=cps*abs(dt/p1(3))/ecm

      end if

      end


