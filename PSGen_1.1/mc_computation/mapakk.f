      subroutine mapakk(imom,x,q,sqsc,sqsl,sqsu,s,sqs,ds)
c----------------------------------------------------------------------c
c Performs mappings of the s-channel photon and gluon propagators
c
c Karol Kolodziej, University of Silesia, karol.kolodziej@us.edu.pl
c 1 April, 2008; modified on 14 May, 2014; 11 June, 2019
c----------------------------------------------------------------------c

      implicit none
 
      integer, parameter :: dp=kind(1.d0)
      integer :: imom
      real(dp) :: x,q(0:3),sqsc,sqsl,sqsu,s,sqs,ds
      real(dp) :: rs

      if (sqsl == 0.d0) sqsl=sqsc

      if (sqsu <= sqsl) then  ! added on 14 May, 2014
       s=0.d0
       sqs=0.d0
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
       if (s < sqsl**2) s=sqsl**2 ! 11 June, 2019
       sqs=sqrt(s)
       ds=2.d0*log(sqsu/sqsl)*s
      end if

      end

