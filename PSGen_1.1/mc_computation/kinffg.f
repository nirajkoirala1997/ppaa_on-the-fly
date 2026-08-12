      subroutine kinffg(imom,ecm,m3,m4,ecut,e3,e7,x,p3,p4,p7,dlips)
c----------------------------------------------------------------------c
c Returns the final state four momenta p3,p4,p7 and the phase space 
c normalization for a 2 --> 3 process in the c.m.s. with p7 being a 
c photon in the final state. Energy of particle 3 is generated according 
c to 1/(p4.p7). The factor (2*pi)^-5 is not included in dlips.
c
c Input: imom=1/(else) => momenta p1 and p2 are/(not) calculated,
c        ecm - the c.m.s. energy,
c        m3,m4 - particle masses,
c        ecut - lower cut on the energy of the photon,
c        x(i) - random variables, 0 < x(i) < 1, i=1,...,5
c
c Karol Kolodziej, University of Silesia, kolodzie@us.edu.pl
c 22 April, 1999; modified on 1 June, 2000, 13 May, 2019
c----------------------------------------------------------------------c

      implicit none
      integer, parameter :: dp=kind(1.d0)

      real(dp) :: ecut,x(5),dlips
      real(dp) :: ecm,m3,m4,e3,e7,p3(0:3),p4(0:3),p7(0:3),px,py,pz,
     &            pp32,pp3,re3,ce3,dele3,e4max,m42,ecm3,e7num

      integer :: ifirst=1,imom

      real(dp) :: e7min,e7max,re7,de7,de3,
     &            ct37,st37,ct3,phi3,phi37,st3,cf3,sf3

      real(dp) :: twopi,dphi3,dphi37,dct3,cps

      save ifirst,dphi3,dphi37,dct3,cps

      if (ifirst == 1 ) then
       ifirst=0
       twopi=8.d0*atan(1.d0)
       dphi3=twopi
       dphi37=twopi
       dct3=2.d0
       cps=dphi3*dphi37/8.d0
      end if

      m42=m4**2

c The energy of particle 3:

      e4max=(ecm+(m4-m3)*(m3+m4)/ecm)/2.d0
      ce3=(ecm+(m3-m4)*(m3+m4)/ecm)/2.d0
      re3=m42/(2.d0*ecm)*(ecut/e4max)/(ce3-m3)
      if (imom == 1) then
       dele3=(ce3-m3)*re3**x(1)
       p3(0)=ce3-dele3
       de3=-dele3*log(re3)
      else
       de3=-(ce3-e3)*log(re3)
      end if

      if (p3(0) < m3) then
       dlips=0.d0
       return
      end if

      if (imom == 1) then
       pp32=(p3(0)-m3)*(p3(0)+m3)
       ecm3=ecm-p3(0)
      else
       pp32=(e3-m3)*(e3+m3)
       ecm3=ecm-e3
      end if

      pp3=sqrt(pp32)

c The photon energy:

      e7num=ecm3**2-pp32-m42

      e7min=max(ecut,e7num/(2.d0*(ecm3+pp3)))

      e7max=e7num/(2.d0*(ecm3-pp3))

      if (e7max <= e7min) then
       dlips=0.d0
       return
      end if

      re7=e7max/e7min

      if (imom == 1) then
       p7(0)=e7min*re7**x(2)
       de7=p7(0)*log(re7)
      else
       de7=e7*log(re7)
      end if

c Phase space normalization:

      dlips=cps*de7*de3*dct3

      if (imom /= 1) return


c energy of particle 4:

      p4(0)=ecm3-p7(0)
     
c cos of the angle between p3 and p7:

      if (pp3 == 0.d0) then
       ct37=1.d0
      else
       ct37=(p4(0)**2-pp32-p7(0)**2-m42)/(2.d0*pp3*p7(0))
      end if

c Angles of particle 3:

      ct3=dct3*x(3)-1.d0

      if (abs(ct37)  > 1.d0) then
       dlips=0.d0
       return
      end if

      if (abs(ct3) < 1.d0) then
       st3=sqrt(1.d0-ct3**2)
      else if (ct3 < -1.d0) then
       ct3=-1.d0
       st3=0.d0
      else if (ct3 > 1.d0) then
       ct3=1.d0
       st3=0.d0
      end if

      phi3=dphi3*x(4)
      cf3=cos(phi3)
      sf3=sin(phi3)  

c The azimuthal angle of particle 7 with respect to particle 3, phi7-phi3:

      phi37=dphi37*x(5)

c Momenta of particles 3, 7 and 4 in the frame where p3(3)=Oz:

      p3(1:)=(/0.d0,0.d0,pp3/)
      st37=sqrt(1.d0-ct37**2)
      p7(1:)=p7(0)*(/st37*cos(phi37),st37*sin(phi37),ct37/)

c Rotation of the frame by angle theta_3 in the xOz plane and by angle
c -phi_3 in the xOy plane:

      p3(1)=pp3*st3*cf3
      p3(2)=pp3*st3*sf3
      p3(3)=pp3*ct3

      px=p7(1)
      py=p7(2)
      pz=p7(3)
      p7(1)=(px*ct3+pz*st3)*cf3-py*sf3
      p7(2)=(px*ct3+pz*st3)*sf3+py*cf3
      p7(3)=-px*st3+pz*ct3

      p4(1:)=-(p3(1:)+p7(1:))

      end



