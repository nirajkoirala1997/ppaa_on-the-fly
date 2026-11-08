      subroutine mapc_qp(imom,x,q,m,wdth,sqsc,sqsl,sqsu,s,sqs,ds)
c----------------------------------------------------------------------c
c Performs mappings of the s-channel propagators
c
c Karol Kolodziej, University of Silesia, karol.kolodziej@us.edu.pl
c 24 May, 2021, 2008; modified on 
c----------------------------------------------------------------------c

      implicit none
 
      integer, parameter :: dp=kind(1.d0),qp=kind(1.q0)
      integer :: imom
      real(dp) :: x,m,wdth,sqsc,ds
      real(qp) :: wdm,ats,nfct
      real(qp) :: q(0:3),sqsl,sqsu,s,sqs,rs

      if (sqsl == 0.d0) sqsl=sqsc

      if (sqsu <= sqsl) then
       s=0.d0
       sqs=0.d0
       ds=0.d0
       return
      end if

      if (imom == 1) then
       if (m == 0.d0) then
        rs=sqsu/sqsl
        s=sqsl*sqsl*rs**(2.d0*x)
        sqs=sqrt(s)
        ds=2.d0*log(rs)*s
       else
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
       end if
      else
       s=(q(0)-q(3))*(q(0)+q(3))-q(1)*q(1)-q(2)*q(2)
       sqs=sqrt(s)
       if (m == 0.d0) then
        ds=2.d0*log(sqsu/sqsl)*s
       else
        if (wdth > 0.d0) then
         wdm=wdth*m
         ats=-atan((m-sqsl)*(m+sqsl)/wdm)
         nfct=wdm/(atan((sqsu-m)*(sqsu+m)/wdm)-ats)
         ds=(((sqs-m)*(sqs+m))**2+wdm*wdm)/nfct
        else
         ds=(sqsu-sqsl)*(sqsu+sqsl)
        end if
       end if
      end if

      end


