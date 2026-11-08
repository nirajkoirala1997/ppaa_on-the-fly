      subroutine kintt3(imom,ecm,sqsp,it1,mt1,ctc3,m1,m3,
     &                                it2,mt2,ctc4,m2,m4,
     &                                rand,p1,p2,p3,p4,q,dlips)
c----------------------------------------------------------------------c
c Returns cartezian components of the final state four momenta in the 
c c.m.s. and the phase space normalization of a 2 --> 3 subbrocess 
c with two t-channel poles ~1/[(p1-p3)^-mt1^2] and ~1/[(p2-p4)^-mt2^2].
c The factor (2*pi)^-5 is not included in dlips.
c 
c Input: imom=1/else => four momenta p3,p4 are/are not calculated,
c        ecm - c.m.s. energy,
c        sqsp - invariant mass of the subsystem of particles 5,6,...
c        it1=1/else the pole ~1/[(p1-p3)^-mt1^2] is mapped/not mapped
c        it2=1/else the pole ~1/[(p2-p4)^-mt2^2] is mapped/not mapped
c        ctc3 - angular cut on cosine of (p1,p3) angle,
c        ctc4 - angular cut on cosine of (p1,p4) angle,
c        mt1,mt2 - masses of particles in the t-channel propagators,
c        m1,m2,m3,m4 - external particle masses,
c        rand(i) - random variables, 0 < rand(i) < 1, i=1,2,...,5
c        p1,p2 - initial state particle four momenta.
c Output: p3,p4 - final state particle four momenta,
c         q - for mometum of the subsystem of particles 5,6,...
c         dlips - phase space normalization.
c         npsp - number of solutions for y=2*p4(0)/ecm,
c
c Karol Kolodziej, University of Silesia, kolodzie@us.edu.pl
c 29 Nov.- 2 Dec., 2010; modified on 
c 13 Jan. 2022: it1 option for x commented
c----------------------------------------------------------------------c

      use ntdisc

      implicit none

      integer, parameter :: dp=kind(1.d0),qp=kind(1.q0)
      integer :: imom,it1,it2,iy,npsp
      real(dp) :: ctc3,ctc4,mt1,mt2,rand(5),dlips
      real(dp) :: ecm,sqsp,m1,m2,m3,m4,
     &            p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),
     &            pp1
      real(dp), parameter :: pi=3.141592653589793238d0,
     &                       dphi=2.d0*pi,dct=2.d0,cps=dphi**2/8.d0,
     &                       dprec=1.d-6,dprec1=1.d-16
      real(dp) :: mu30,mu40,mu0,mu3,mu4,mu,mm,xl,dxu,xul,yl,yu,x,dx,
     &            ph3,cph3,sph3,ph4,cph4,sph4,
     &            pt3,pp3,db3,beta3,rb3,ct3,dct3,st3,
     &            e4w,pt4,pp4,db4,beta4,a4,da4,rb4,ct4,dct4,st4,ct34,yp,
     &            ym,y,x2,x4,x24,a40,act3,act4!,e3max,e4max

      real :: randy

      if (ecm < m3+m4+sqsp) then
       dlips=0.d0
       return
      end if
c Generate E3=x*ecm/2:
      mu30=m3/ecm
      mu40=m4/ecm
      mu0=sqsp/ecm
      xl=2.d0*mu30
      mu3=mu30**2
c      xu=1.d0+mu3-(mu40+mu0)**2
      dxu=(mu0+mu40-mu30)*(mu0+mu40+mu30)
c      e3max=0.5d0*ecm*(1.d0-dxu)
c      e4max=0.5d0*ecm*(1.d0-(mu0+mu30-mu40)*(mu0+mu30+mu40))
      mu4=mu40**2
      mu=mu0**2
      mm=mu3+mu4+(1.d0-mu0)*(1.d0+mu0)
c      if (it1 == 1) then                       ! 13 Jan. 2022
c Generate x according to ~1/(1-x):
      xul=dxu/(1.d0-xl)
      if (imom == 1) then
       x=1.d0-(1.d0-xl)*xul**rand(1)
       p3(0)=0.5d0*ecm*x
      else
       x=2.d0*p3(0)/ecm
      end if
      dx=-(1.d0-x)*log(xul)
c      else                                     ! 13 Jan. 2022
c Generate x according to flat distribution:
cc       xu=1.d0-dxu
c       dx=1.d0-dxu-xl
c       if (imom == 1) then
c        x=dx*rand(1)+xl
c        p3(0)=0.5d0*ecm*x
c       else
c        x=2.d0*p3(0)/ecm
c       end if
c      end if                                   ! 13 Jan. 2022

      if (p3(0) > m3) then
       pp3=sqrt((p3(0)-m3)*(p3(0)+m3))
      else
       dlips=0.d0
       return
      end if

c Generate phi3 according to flat distribution:
      if (imom == 1) then
       ph3=dphi*rand(2)
       cph3=cos(ph3)
       sph3=sin(ph3)
      else
       pt3=sqrt(p3(1)**2+p3(2)**2)
c       pp3=sqrt(pt3**2+p3(3)**2)
       if (pt3 > 0.d0) then
        cph3=p3(1)/pt3
        sph3=p3(2)/pt3
       else
        cph3=1.d0
        sph3=0.d0
       end if
      end if

c Generate ct3 in [-ctc3,ctc3]:
      if (it1 == 1) then
c according to ~1/(1-beta3*ct3):
       pp1=sqrt((p1(0)-m1)*(p1(0)+m1))
       beta3=2.d0*pp1*pp3
     &                /(mt1**2-m1**2-m3**2+2.d0*p1(0)*p3(0))
       db3=1.d0-beta3
       if (abs(db3) < dprec) then
        db3=-((m1+m3)**2-mt1**2)/(2.d0*p1(0)*p3(0))
     &                          +0.5d0*(m1/p1(0)+m3/p3(0))**2
       end if
       rb3=(1.d0-ctc3+db3*ctc3)/(1.d0+ctc3-db3*ctc3)
       if (imom == 1) then
        act3=(1.d0+ctc3-db3*ctc3)*rb3**rand(3)
        ct3=(1.d0-act3)/beta3
       else
        ct3=p3(3)/pp3
        act3=1.d0-beta3*ct3
       end if
       dct3=-(1.d0-ct3+db3*ct3)*log(rb3)/beta3
      else
c according to flat distribution:
       dct3=2.d0*ctc3
       if (imom == 1) then
        ct3=dct3*rand(3)-ctc3
       else
        ct3=p3(3)/pp3
       end if
      end if

c Generate phi4 according to flat distribution:
      if (imom == 1) then
       ph4=dphi*rand(4)
       cph4=cos(ph4)
       sph4=sin(ph4)
      else
       pt4=sqrt(p4(1)**2+p4(2)**2)
       pp4=sqrt(pt4**2+p4(3)**2)
       if (pt4 > 0.d0) then
        cph4=p4(1)/pt4
        sph4=p4(2)/pt4
       else
        cph4=-1.d0
        sph4=0.d0
       end if
      end if
c Generate ct4 in [-ctc4,ctc4]:
      if (it2 == 1) then
       ct34=-1.d0
       call findy
c       call findy1
       if (iy > 0) then
        e4w=0.5d0*ecm*yp
        if (e4w > m4) then
c according to ~1/(1-beta4*ct4):
c         beta4=2.d0*sqrt((p2(0)-m2)*(e4w-m4)*(p2(0)+m2)*(e4w+m4))
c     &        /(mt2**2-m2**2-m4**2+2.d0*p2(0)*e4w)
c         db4=1.d0-beta4
c according to ~1/(a4+ct4):
c         a4=(mt2**2-m2**2-m4**2+2.d0*p2(0)*e4w)
c     &       /(2.d0*sqrt((p2(0)-m2)*(e4w-m4)*(p2(0)+m2)*(e4w+m4)))
         x2=(m2/p2(0))**2
         x4=(m4/e4w)**2
         x24=(m2**2+m4**2-mt2**2)/(p2(0)*e4w)
         if (x2 < dprec .and. x4 < dprec) then
          a40=1.d0
          da4=0.5d0*(x2+x4-x24)
     &       +0.25d0*(1.5d0*(x2**2+x4**2)+x2*x4-x24*(x2+x4))
c         else if (x2 < dprec) then
c          a40=1.d0/sqrt(1.d0-x4)
c          da4=(x2-x24+0.5d0*x2*(1.5d0*x2-x24))/(2.d0*sqrt(1-x4))
c         else if (x4 < dprec) then
c          a40=1.d0/sqrt(1.d0-x2)
c          da4=(x4-x24+0.5d0*x4*(1.5d0*x4-x24))/(2.d0*sqrt(1-x2))
         else
          a40=1.d0
          da4=(1.d0-0.5d0*x24)/sqrt((1.d0-x2)*(1.d0-x4))-1.d0
         end if
c In the small mass limit m2<<E2, m4<<E4, mt2^2<<E2*E4:
c        db4=0.5d0*((m2**2/p2(0)+m4**2/e4w)*(1.d0/p2(0)-1.d0/e4w)
c     &             +mt2**2/(p2(0)*e4w))
c        beta4=1.d0-db4
        else
         dlips=0.d0
         return
        end if
        if (da4 <= 0.d0) then
         dlips=0.d0
         return
        end if
        a4=a40+da4
c        rb4=(a4+ctc4)/(a40-ctc4+da4)
        rb4=(a4+ctc4)/da4/(1.d0+(a40-ctc4)/da4)
        if (imom == 1) then
c         ct4=((1.d0-(1.d0+ctc4-db4*ctc4)*rb4**rand(5)))/beta4
         if (ctc4 == 1.d0 .and. a40 == 1.d0) then
          act4=da4*rb4**rand(5)
          ct4=act4-a4
         else
          act4=(a4-ctc4)*rb4**rand(5)
          ct4=act4-a4
         end if
c         ct4=(a40-ctc4+da4)*rb4**rand(5)-a4
c         if (((p2(0)-p4(0))**2                 ! added 2 Feb. 2021
c     &       -(-pp1-sqrt((p4(0)-m4)*(p4(0)+m4))*ct4)**2)/m2**2 >
c     &       -(p4(0)/p2(0)+p2(0)/p4(0))+2.d0) then
c          dlips=0.d0
c          return
c         end if
        else
         ct4=p4(3)/pp4
         act4=a4+ct4
        end if
c        dct4=-(1.d0-ct4+db4*ct4)*log(rb4)/beta4
c        
        if (1.d0+ct4 < dprec1) then 
         dct4=act4*log(rb4)
        else
         dct4=(a4+ct4)*log(rb4)
        end if
       else
c according to flat distribution:
        dct4=2.d0*ctc4
        if (imom == 1) then
         ct4=dct4*rand(5)-ctc4
        else
         ct4=p4(3)/pp4
        end if
       end if
      else
       dct4=2.d0*ctc4
       if (imom == 1) then
        ct4=dct4*rand(5)-ctc4
       else
        ct4=p4(3)/pp4
       end if
      end if

      if (imom == 1) then
       if (1.d0-ct3 < dprec1) then 
        st3=sqrt(2.d0*(act3-db3))
       else
        st3=sqrt((1.d0-ct3)*(1.d0+ct3))
       end if
       if (1.d0+ct4 < dprec1) then 
        st4=sqrt(2.d0*(act4-da4))
       else
        st4=sqrt((1.d0-ct4)*(1.d0+ct4))
       end if
      else 
       st3=sqrt((1.d0-ct3)*(1.d0+ct3))
       st4=sqrt((1.d0-ct4)*(1.d0+ct4))
      end if

      ct34=ct3*ct4+st3*st4*(cph3*cph4+sph3*sph4)

      if (imom == 1) then
       call findy
       npsp=iy
       if (iy > 0) then
        if (iy == 1) then
         p4(0)=0.5d0*ecm*yp
        else if (iy == 2) then
         call random_number(randy)
c         call ranlux(randy,1)
         if (randy < 0.5) then
          p4(0)=0.5d0*ecm*yp
         else
          p4(0)=0.5d0*ecm*ym
         end if
        end if
        if (p4(0) > m4) then
         pp4=sqrt((p4(0)-m4)*(p4(0)+m4))
        else
         dlips=0.d0
         return
        end if
c        dlips=cps*pp3*pp4/(2.d0-x*(1.d0-ct34))*dx*dct3*dct4
         dlips=cps*pp3*pp4/abs(2.d0-x+(2.d0*p4(0))/ecm*pp3/pp4*ct34)
     &                                                     *dx*dct3*dct4
c         y=2.d0*p4(0)/ecm
c         dlips=(ecm/2.d0)**2*sqrt((y**2-4.d0*mu4)**3*(x**2-4.d0*mu3))
c     &   /abs(2.d0*(2.d0*mu4*(2.d0-x)-y*(mm-x)))*dx*dct3*dct4
       else
        dlips=0.d0
        return
       end if
       if (dlips > 0.d0) then
        p3(1:3)=pp3*(/st3*cph3,st3*sph3,ct3/)
        p4(1:3)=pp4*(/st4*cph4,st4*sph4,ct4/)
c The conditions below eliminate important part of the phase space
        if (((p1(0)-p3(0))**2-(p1(1)-p3(1))**2-(p1(2)-p3(2))**2
     &     -(p1(3)-p3(3))**2)/m1**2 >= 
     &                   1.d0*(-(p3(0)/p1(0)+p1(0)/p3(0))+2.d0)) then 
c     &                      -(e3max/p1(0)+p1(0)/e3max)+2.d0) then 
         ntd1=ntd1+1
         dlips=0.d0
         return
        end if
        if (((p2(0)-p4(0))**2-(p2(1)-p4(1))**2-(p2(2)-p4(2))**2
     &     -(p2(3)-p4(3))**2)/m2**2 >= 
     &                   1.d0*(-(p4(0)/p2(0)+p2(0)/p4(0))+2.d0)) then
c     &                      -(e4max/p2(0)+p2(0)/e4max)+2.d0) then
         ntd2=ntd2+1
         dlips=0.d0
         return
        end if
        q(1:3)=-p3(1:3)-p4(1:3)
        q(0)=sqrt(sqsp**2+q(1)**2+q(2)**2+q(3)**2)
       end if
      else
       if (p4(0) > m4) then
        pp4=sqrt((p4(0)-m4)*(p4(0)+m4))
c        dlips=cps*pp3*pp4/(2.d0-x*(1.d0-ct34))*dx*dct3*dct4
        dlips=cps*pp3*pp4/abs(2.d0-x+(2.d0*p4(0))/ecm*pp3/pp4*ct34)
     &                                                     *dx*dct3*dct4
       else
        dlips=0.d0
       end if
      end if

      contains

      subroutine findy

      real(dp) :: alpha,beta,gamma,delta,sqdel,omx,xmu,mmx

      yp=0.d0
      ym=0.d0

      if (ct34 == 0.d0) then
       iy=1
       yp=(mm-x)/(1.d0-0.5d0*x)
      else
       omx=1.d0-0.5d0*x
c       xmu=sqrt(x**2-4.d0*mu3)

       xmu=sqrt((x-2.d0*mu30)*(x+2.d0*mu30))
       mmx=mm-x
c       alpha=4.d0*omx**2-(ct34*xmu)**2

       alpha=(2.d0*omx-ct34*xmu)*(2.d0*omx+ct34*xmu)
       beta=4.d0*omx*(mm-x)
       gamma=4.d0*(mmx**2+mu4*(ct34*xmu)**2)
       delta=beta**2-alpha*gamma

       if (delta < 0.d0) then
        iy=0
       else if (delta == 0.d0) then
        iy=1
        yp=beta/alpha
       else
        iy=0
        sqdel=sqrt(delta)
        yp=(beta+sqdel)/alpha
        ym=(beta-sqdel)/alpha
        if (abs(2.d0*(mmx-yp*omx)/
     &    (xmu*sqrt((yp-2.d0*mu40)*(yp+2.d0*mu40)))-ct34) < dprec) then
         iy=iy+1
         if (abs(2.d0*(mmx-ym*omx)/
     &    (xmu*sqrt((ym-2.d0*mu40)*(ym+2.d0*mu40)))-ct34) < dprec) then
          iy=iy+1
         end if
        else
         if (abs(2.d0*(mmx-ym*omx)/
     &    (xmu*sqrt((ym-2.d0*mu40)*(ym+2.d0*mu40)))-ct34) < dprec) then
          iy=iy+1
          yp=ym
         end if
        end if
       end if
      end if

      end subroutine findy

      subroutine findy1

      real(dp) :: alpha,mm,beta,gamma,delta,sqdel,omx,xmu,mmx

      yp=0.d0
      ym=0.d0

      mm=mu3+mu4+(1.d0-mu0)*(1.d0+mu0)
      omx=1.d0-0.5d0*x
      xmu=sqrt((x-2.d0*mu30)*(x+2.d0*mu30))
      mmx=mm-x
      delta=1.d0-(2.d0*mu40/mmx)**2*(1.d0-x+mu3)
      if (delta < 0.d0) then
       iy=0
      else if (delta == 0.d0) then
       iy=1
       yp=mmx*omx/(1.d0-x+mu3)
      else
       iy=0
       sqdel=xmu*sqrt(delta)/2.d0
       yp=mmx*(omx+sqdel)/(1.d0-x+mu3)
       ym=mmx*(omx-sqdel)/(1.d0-x+mu3)
       if (abs(2.d0*(mmx-yp*omx)/
     &    (xmu*sqrt((yp-2.d0*mu40)*(yp+2.d0*mu40)))+1.d0) < dprec) then
        iy=iy+1
        if (abs(2.d0*(mmx-ym*omx)/
     &    (xmu*sqrt((ym-2.d0*mu40)*(ym+2.d0*mu40)))+1.d0) < dprec) then
         iy=iy+1
        end if
       else
        if (abs(2.d0*(mmx-ym*omx)/
     &   (xmu*sqrt((ym-2.d0*mu40)*(ym+2.d0*mu40)))+1) < dprec) then
          iy=iy+1
          yp=ym
        end if
       end if
      end if

      end subroutine findy1

      end subroutine kintt3


