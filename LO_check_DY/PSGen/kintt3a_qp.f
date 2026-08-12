      subroutine kintt3a_qp(imom,ecm,ecuta,ctcq,ctc3,m1,m2,m3,m4,
     &                                rand,p1,p2,p3,p4,q,dlips)
c----------------------------------------------------------------------c
c Returns cartezian components of the final state four momenta in the 
c c.m.s. and the phase space normalization of a 2 --> 3 subbrocess 
c with a single t-channel pole ~1/(p1-p3)^2 and a photon four
c momentum q generated according to ~1/[q(0)*(1-beta^2*cos(theta)^2].
c The factor (2*pi)^-5 is not included in dlips.
c 
c Input: imom=1/else => four momenta p3,p4 are/are not calculated,
c        ecm - c.m.s. energy,
c        ecuta - photon energy cut in c.m.s,
c        ctcq,ctc3 - angular cut on cosine of (p1,q) or (p1,p3) angle,
c        m1,m2,m3,m4 - external particle masses,
c        rand(i) - random variables, 0 < rand(i) < 1, i=1,2,...,5
c        p1,p2 - initial state particle four momenta.
c Output: p3,p4,q - final state particle four momenta,
c         dlips - phase space normalization.
c
c Karol Kolodziej, University of Silesia, karol.kolodziej@us.edu.pl
c 21 Jan., 2019; modified on 
c----------------------------------------------------------------------c

      implicit none

      integer, parameter :: dp=kind(1.d0),qp=kind(1.q0)
      integer :: imom,iy,npsp
      real(dp) :: ecuta,ctcq,ctc3,rand(5),dlips
      real(qp) :: ecm,m1,m2,m3,m4,
     &            p1(0:3),p2(0:3),p3(0:3),p4(0:3),q(0:3),
     &                       dprec=1.q-12
      real(dp), parameter :: pi=3.141592653589793238d0,
     &                       dphi=2.d0*pi,cps=dphi**2/8.d0
      real(qp) :: mu30,mu40,mu0,mu3,mu4,mu,mm,xul,yl,yu,x,dx,
     &            ph3,cph3,sph3,phi,cph,sph,
     &            pt3,pp3,db3,beta3,a3,rb3,ct3,dct3,st3,
     &            e3w,pt4,pp4,db4,beta4,a4,da4,rb4,ct,dct,st,ctq3,yp,
     &            ym,y,qtr,beta,beta1,qu

      real :: randy

      if (ecm < m3+m4+ecuta) then
       dlips=0.d0
       return
      end if
c Generate photon energy according to ~1/|q|, q(0)=x*ecm/2:
      qu=(ecm-(m3+m4))*(ecm+(m3+m4))/(2.q0*ecm)
      if (imom == 1) then
       q(0)=ecuta*(qu/ecuta)**rand(1)
      end if
      x=q(0)/(ecm/2.q0)
      dx=q(0)/(ecm/2.q0)*log(qu/ecuta)
c Generate photon polar angle according to ~1/(1-beta^2ct^2):
      beta=sqrt(1.q0-(2.q0*m1/ecm)**2)
c      beta1=((1.q0+beta)/(2.q0*m1/ecm))**2
      beta1=(1.q0+beta*ctcq)/(1.q0-beta*ctcq)
      if (imom == 1) then
c       beta3=beta1**(2.q0*rand(2)-1.q0)
       beta3=beta1**(2.q0*rand(2)-1.q0)
       ct=-(1.q0-beta3)/(1.q0+beta3)/beta
      else
       ct=q(3)/q(0)
      end if
      dct=(1.q0-(beta*ct)**2)/(beta/log(beta1))
      st=sqrt(1.q0-ct**2)
c Generate photon azimuthal angle according to flat distribution:
      if (imom == 1) then
       phi=dphi*rand(3)
       sph=sin(phi)
       cph=cos(phi)
      else
       qtr=q(0)*st
       if (qtr > 0.q0) then
        sph=q(2)/qtr
        cph=q(1)/qtr
       else
        sph=0.q0
        cph=1.q0
       end if
      end if


      mu30=0.q0 !m3/ecm
      mu40=m3/ecm !m4/ecm
      mu0=m4/ecm !sqsp/ecm
      mu3=mu30**2
      mu4=mu40**2
      mu=mu0**2
      mm=mu3+mu4+(1.q0-mu0)*(1.q0+mu0)

      if (p3(0) > m3) then
       pp3=sqrt((p3(0)-m3)*(p3(0)+m3))
      else
       dlips=0.d0
       return
      end if

c Generate phi3 according to flat distribution:
      if (imom == 1) then
       ph3=dphi*rand(4)
       cph3=cos(ph3)
       sph3=sin(ph3)
      else
       pt3=sqrt(p3(1)**2+p3(2)**2)
       if (pt3 > 0.d0) then
        cph3=p3(1)/pt3
        sph3=p3(2)/pt3
       else
        cph3=1.q0
        sph3=0.q0
       end if
      end if

c Generate ct3 in [-ctc3,ctc3]:
      ctq3=1.q0
      call findy
      if (iy > 0) then
       e3w=0.5q0*ecm*yp
       if (e3w > m3) then
c according to ~1/(a3-ct3):
        a3=(-m1**2-m3**2+2.d0*p1(0)*e3w)
     &       /(2.q0*sqrt((p1(0)-m1)*(e3w-m3)*(p1(0)+m1)*(e3w+m3)))
       else
        dlips=0.d0
        return
       end if
       rb3=(a3-ctc3)/(a3+ctc3)
       if (imom == 1) then
        ct3=a3-(a3+ctc3)*rb3**rand(5)
       else
        ct3=p3(3)/pp3
       end if
       dct3=-(a3-ct3)*log(rb3)
      end if

      st3=sqrt((1.q0-ct3)*(1.q0+ct3))

      ctq3=ct*ct3+st*st3*(cph*cph3+sph*sph3)

      if (imom == 1) then
       call findy
       npsp=iy
       if (iy > 0) then
        if (iy == 1) then
         p3(0)=0.5q0*ecm*yp
        else if (iy == 2) then
         call ranlux(randy,1)
         if (randy < 0.5) then
          p3(0)=0.5q0*ecm*yp
         else
          p3(0)=0.5q0*ecm*ym
         end if
        end if
        if (p3(0) > m3) then
         pp3=sqrt((p3(0)-m3)*(p3(0)+m3))
        else
         dlips=0.d0
         return
        end if
        dlips=cps*pp3*q(0)/abs(2.q0-x+(2.q0*p3(0))/ecm*q(0)/pp3*ctq3)
     &                                                    *dx*dct*dct3
       else
        dlips=0.d0
        return
       end if
       if (dlips > 0.q0) then
        p3(1:3)=pp3*(/st3*cph3,st3*sph3,ct3/)
        q(1:3)=q(0)*(/st*cph,st*sph,ct/)
        p4(1:3)=-p3(1:3)-q(1:3)
        p4(0)=sqrt(m4**2+p4(1)**2+p4(2)**2+p4(3)**2)
       end if
      else
       if (p3(0) > m3) then
        pp3=sqrt((p3(0)-m3)*(p3(0)+m3))
        dlips=cps*pp3*q(0)/abs(2.q0-x+(2.q0*p3(0))/ecm*q(0)/pp3*ctq3)
     &                                                    *dx*dct*dct3
       else
        dlips=0.d0
       end if
      end if

      contains

      subroutine findy

      real(qp) :: alpha,bt,gamma,delta,sqdel,omx,xmu,mmx

      yp=0.q0
      ym=0.q0

      if (ctq3 == 0.q0) then
       iy=1
       yp=(mm-x)/(1.q0-0.5q0*x)
      else
       omx=1.q0-0.5q0*x
c       xmu=sqrt(x**2-4.d0*mu3)

       xmu=sqrt((x-2.q0*mu30)*(x+2.q0*mu30))
       mmx=mm-x
c       alpha=4.d0*omx**2-(ctq3*xmu)**2

       alpha=(2.q0*omx-ctq3*xmu)*(2.q0*omx+ctq3*xmu)
       bt=4.q0*omx*(mm-x)
       gamma=4.q0*(mmx**2+mu4*(ctq3*xmu)**2)
       delta=bt**2-alpha*gamma

       if (delta < 0.q0) then
        iy=0
       else if (delta == 0.q0) then
        iy=1
        yp=bt/alpha
       else
        iy=0
        sqdel=sqrt(delta)
        yp=(bt+sqdel)/alpha
        ym=(bt-sqdel)/alpha
        if (abs(2.q0*(mmx-yp*omx)/
     &    (xmu*sqrt((yp-2.q0*mu40)*(yp+2.q0*mu40)))-ctq3) < dprec) then
         iy=iy+1
         if (abs(2.q0*(mmx-ym*omx)/
     &    (xmu*sqrt((ym-2.q0*mu40)*(ym+2.q0*mu40)))-ctq3) < dprec) then
          iy=iy+1
         end if
        else
         if (abs(2.q0*(mmx-ym*omx)/
     &    (xmu*sqrt((ym-2.q0*mu40)*(ym+2.q0*mu40)))-ctq3) < dprec) then
          iy=iy+1
          yp=ym
         end if
        end if
       end if
      end if

      end subroutine findy

      end subroutine kintt3a_qp


