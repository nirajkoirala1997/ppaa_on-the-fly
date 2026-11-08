      subroutine ini_photon(imom,ecm,x,sp,sqsp,dlipsg)

      use fourmom
      use kinparams
      use kincuts 

      implicit none

c Energy and angles of the initial state photon in the c.m.s.:


      integer :: imom

      real(dp), parameter :: dphig=2.d0*pi

      real(dp) :: neg,deg,nctg,ctg,stg,phig,x(3),dlipsg

      real(xp) :: ecm,msum,reg,sp,sqsp,beta,ybta,dctg

      msum=sum(m(3:))
      reg=(ecm-msum)*(ecm+msum)/(2.d0*ecm)/ecut0

      if (reg > 0.d0) then
       neg=1.d0/log(reg)
      else
       dlipsg=0.d0
       return
      end if

      if (imom == 1) then
       p(0,igind)=ecut0*reg**x(3)
      end if

      deg=p(0,igind)/neg

      if (iqprec == 1) then
       beta=sqrt(1.q0-(2.q0*m(1)/ecm)**2)
       ybta=(1.q0+beta*ct0b)/(1.q0-beta*ct0b)
      else
       beta=sqrt(1.d0-(2.d0*m(1)/ecm)**2)
       ybta=(1.d0+beta*ct0b)/(1.d0-beta*ct0b)
      end if

      nctg=beta/log(ybta)

      if (imom == 1) then
       ctg=ybta**(2.d0*x(2)-1.d0)
       ctg=(ctg-1.d0)/(ctg+1.d0)/beta
      else
       ctg=p(3,igind)/p(0,igind)
      end if

      if (ctg > 1.d0) then
       ctg=1.d0
      else if (ctg < -1.d0) then
       ctg=-1.d0
      end if

      if (iqprec == 1) then
       dctg=(1.q0-(beta*ctg)**2)/nctg
      else
       dctg=(1.d0-(beta*ctg)**2)/nctg
      end if

      if (imom == 1) then
       stg=sqrt(1.d0-ctg**2)
       phig=2.d0*pi*x(1)
       p(1,igind)=p(0,igind)*stg*cos(phig)
       p(2,igind)=p(0,igind)*stg*sin(phig)
       p(3,igind)=p(0,igind)*ctg
      end if

      sp=ecm*(ecm-2.d0*p(0,igind))
      sqsp=sqrt(sp)

      dlipsg=0.5d0*p(0,igind)*deg*dctg*dphig

      end
