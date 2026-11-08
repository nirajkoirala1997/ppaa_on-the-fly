      subroutine kinin_qp(ecm,m,p,flux)
c----------------------------------------------------------------------c
c Returns four momenta of the initial state particles in the c.m.s.
c and the flux factor for the cscattering process
c
c                    f1(p1) + f2(p2) -> something
c
c The frame is defined such that p1=(E1,0,0,p), p2=(E2,0,0,-p).
c
c Karol Kolodziej, University of Silesia, kolodzie@us.edu.pl
c 15 December, 2005; modiffied on 23 March, 2006;
c                                 23 May, 2007;
c                                 8 Nov. 2018 - quadruple precision
c----------------------------------------------------------------------c

      implicit none

      integer, parameter :: qp=selected_real_kind(32,99)

      real(kind(1.d0)) :: flux
      real(qp) :: ecm,m(2),p(0:3,2),pp,e

c Initial state particle momenta:

      if (m(1) == m(2)) then

       e=ecm/2.d0
       pp=sqrt((e-m(1))*(e+m(1)))

       p(:,1)=(/e,0.q0,0.q0,pp/)
       p(:,2)=(/e,0.q0,0.q0,-pp/)

      else

       pp=sqrt((ecm-m(1)-m(2))*(ecm+m(1)+m(2))
     &        *(ecm-m(1)+m(2))*(ecm+m(1)-m(2)))/(2.q0*ecm)

       p(:,1)=(/sqrt(m(1)**2+pp**2),0.q0,0.q0,pp/)
       p(:,2)=(/sqrt(m(2)**2+pp**2),0.q0,0.q0,-pp/)

      end if

c Flux factor:

      flux=4.d0*ecm*pp

      end
