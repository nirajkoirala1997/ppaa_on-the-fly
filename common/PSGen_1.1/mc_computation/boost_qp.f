      subroutine boost_qp(g,v,p)
c----------------------------------------------------------------------c
c Returns a four vector p boosted back to the frame of velocity v
c
c Input: g - Lorentz gamma factor,
c        v - velocity of the frame,
c        p - four momentum before boost. 
c
c Karol Kolodziej, University of Silesia, kolodzie@us.edu.pl
c 3 July, 2001; modiffied on 7 April, 2005
c 16 Nov. 2018: quadruple precision four momenta implemented
c----------------------------------------------------------------------c

      implicit none

      real(kind(1.q0)) :: g,v(3),p(0:3),p0,vp

      vp=v(1)*p(1)+v(2)*p(2)+v(3)*p(3)

      p0=p(0)
      p(0)=g*(p0+vp)

      p(1:)=p(1:)+g*(g*vp/(g+1.q0)+p0)*v

      end
