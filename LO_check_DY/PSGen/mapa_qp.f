      subroutine mapa_qp(imom,x,q,sqsc,sqsl,sqsu,s,sqs,ds)
c----------------------------------------------------------------------c
c Performs mappings of the s-channel photon and gluon propagators
c
c Karol Kolodziej, University of Silesia, karol.kolodziej@us.edu.pl
c 1 April, 2008; modified on 14 May, 2014
c 16 Nov. 2018: quadrple precision four momenta implemented
c 9 Feb. 2023: ds cznged to quadruple precision
c----------------------------------------------------------------------c

      implicit none
 
      integer, parameter :: dp=kind(1.d0),qp=kind(1.q0)
      integer :: imom
      real(dp) :: x,sqsc,ds
      real(qp) :: q(0:3),sqsl,sqsu,s,sqs,rs

      if (sqsl == 0.d0) sqsl=sqsc

      if (sqsu <= sqsl) then  ! added on 14 May, 2014
       s=0.q0
       sqs=0.q0
       ds=0.d0
       return
      end if

      if (imom == 1) then
       rs=sqsu/sqsl
       s=sqsl*sqsl*rs**(2.d0*x)
       sqs=sqrt(s)
       ds=2.d0*log(rs)*s
      else
       s=(q(0)-q(3))*(q(0)+q(3))-q(1)*q(1)-q(2)*q(2)
       sqs=sqrt(s)
       ds=2.d0*log(sqsu/sqsl)*s
      end if

      end

