      subroutine mapb_qp(imom,x,q,m,wdth,sqsl,sqsu,s,sqs,ds)
c----------------------------------------------------------------------c
c Performs mappings of the s-channel propagators
c
c Karol Kolodziej, University of Silesia, karol.kolodziej@us.edu.pl
c 1 April, 2008; modified
c 16 Nov. 2018: quadruple precision four momenta implemented
c----------------------------------------------------------------------c

      implicit none
 
      integer, parameter :: dp=kind(1.d0),qp=kind(1.q0)
      integer :: imom
      real(dp) :: x,m,wdth,ds
      real(qp) :: q(0:3),sqsl,sqsu,s,sqs
      real(qp) :: wdm,ats,nfct

      if (imom == 1) then
       if (wdth > 0.d0) then
        wdm=wdth*m
        ats=-atan((m-sqsl)*(m+sqsl)/wdm)
        nfct=wdm/(atan((sqsu-m)*(sqsu+m)/wdm)-ats)
        s=wdm*tan(wdm/nfct*x+ats)+m*m
        sqs=sqrt(s)
        ds=(((sqs-m)*(sqs+m))**2+wdm*wdm)/nfct
       else
        ds=(sqsu-sqsl)*(sqsu+sqsl)
        s=ds*x+sqsl*sqsl
        sqs=sqrt(s)
       end if
      else
       s=(q(0)-q(3))*(q(0)+q(3))-q(1)*q(1)-q(2)*q(2)
       sqs=sqrt(s)
       if (wdth > 0.d0) then
        wdm=wdth*m
        ats=-atan((m-sqsl)*(m+sqsl)/wdm)
        nfct=wdm/(atan((sqsu-m)*(sqsu+m)/wdm)-ats)
        ds=(((sqs-m)*(sqs+m))**2+wdm*wdm)/nfct
       else
        ds=(sqsu-sqsl)*(sqsu+sqsl)
       end if
      end if

      end

