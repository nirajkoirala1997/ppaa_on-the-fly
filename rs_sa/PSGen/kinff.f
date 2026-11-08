      subroutine kinff(imom,ecm,m1,m2,rand,p1,p2,dlips)
c----------------------------------------------------------------------c
c Returns cartezian components of the final state four momenta and
c the phase space normalization for a 2 --> 2 process in the c.m.s. 
c The factor (2*pi)^-2 is not included in dlips.
c
c Input: imom=1/(else) => momenta p1 and p2 are/(not) calculated
c        ecm - the c.m.s. energy,
c        m1,m2 - particle masses,
c        x(i) - random variables, 0 < x(i) < 1, i=1,2,
c
c Karol Kolodziej, University of Silesia, kolodzie@us.edu.pl
c 24 February, 1999; modified on 23 June, 2005; 29 June, 2007;
c 13 January, 2009;
c 30 January, 2022: formula for e2 and p2(0) changed
c----------------------------------------------------------------------c

      implicit none

      integer, parameter :: dp=kind(1.d0)
      integer :: imom
      real(dp) :: ecm,m1,m2,rand(2),p1(0:3),p2(0:3),pp,dlips
      real(dp), parameter :: pi=3.141592653589793238d0, zero=0.d0,
     &                       dphi=2.d0*pi,dct=2.d0,cps=dct*dphi/4.d0
      real(dp) :: ct,pst,phi,delme,e1,e2

c      if (ecm < m1+m2) then
c       dlips=zero
c      else

c       pp=sqrt((ecm-m1-m2)*(ecm+m1+m2)*(ecm-m1+m2)*(ecm+m1-m2))
c     &                                                /(2.d0*ecm)

c       if (imom == 1) then
c        ct=dct*rand(1)-1.d0
c        st=sqrt(1.d0-ct**2)
c        phi=dphi*rand(2)

c        p1=(/sqrt(m1**2+pp**2),pp*st*cos(phi),pp*st*sin(phi),pp*ct/)
c        p2=(/sqrt(m2**2+pp**2),-p1(1),-p1(2),-p1(3)/)
c       end if

c       dlips=cps*pp/ecm 

c       write(*,*)'dlips1=',dlips,cps,pp,ecm

c      end if

c Modification of 13 January, 2009:

      delme=m1+m2

      if (ecm < delme) then
       dlips=zero
      else
       delme=(m1-m2)*delme/ecm
       if (delme == 0.d0) then
        e1=0.5d0*ecm
        e2=e1
        pp=sqrt((e1-m1)*(e1+m1))
       else
        e1=0.5d0*(ecm+delme)
c        e2=e1-delme
        e2=ecm-e1
        if (e1 < m1 .or. e2 < m2) then 
         dlips=0.d0
         return
        end if
        if (delme > 0.d0) then
         pp=sqrt((e1-m1)*(e1+m1))
        else
         pp=sqrt((e2-m2)*(e2+m2))
        end if
       end if

       if (imom == 1) then
c        ct=dct*rand(1)-1.d0
c        pst=pp*sqrt(1.d0-ct**2)
        ct=dct*rand(1)
        pst=pp*sqrt(ct*(2.d0-ct))
        phi=dphi*rand(2)
c        p1=(/e1,pst*cos(phi),pst*sin(phi),pp*ct/)
        p1=(/e1,pst*cos(phi),pst*sin(phi),pp*(ct-1.d0)/)
c        p2=(/e2,-p1(1),-p1(2),-p1(3)/)
        p2=(/sqrt(m2**2+p1(1)**2+p1(2)**2+p1(3)**2),
     &                            -p1(1),-p1(2),-p1(3)/)
       end if

       dlips=cps*pp/ecm

      end if

      end


