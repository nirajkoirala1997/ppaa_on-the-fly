      subroutine king2isr(i1,i2,mv1,gv1,imom,ecm,x,ndim,dlips)
 
      use fourmom
      use kinparams
      use engymom
      use kincuts 

      implicit none

      integer :: i1,i2,imom,ndim
      real(dp) :: mv1,gv1,x(ndim),dlips
      real(xp) :: ecm

c Energy and angles of the initial state photon in the c.m.s.:

      real(dp) :: neg,deg,nctg,ctg,stg,phig,dphig=2.d0*pi

      save dphig


      msum=m(i1)+m(i2)
      reg=(ecm-msum)*(ecm+msum)/(2.d0*ecm)/ecut0

      if (reg > 0.d0) then
       neg=1.d0/log(reg)
      else
       dlips=0.d0
       return
      end if

      if (imom == 1) then
       p(0,igind)=ecut0*reg**x(ndim)
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
       ctg=ybta**(2.d0*x(ndim-1)-1.d0)
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
       phig=2.d0*pi*x(ndim-2)
       p(1,igind)=p(0,igind)*stg*cos(phig)
       p(2,igind)=p(0,igind)*stg*sin(phig)
       p(3,igind)=p(0,igind)*ctg
      end if

      sp=ecm*(ecm-2.d0*p(0,igind))
      sqsp=sqrt(sp)

      if (iqprec == 1) then
       call kinff_qp(imom,sqsp,m(i1),m(i2),x(1:2),p(:,i1),p(:,i2),dlips)
      else
       call kinff(imom,sqsp,m(i1),m(i2),x(1:2),p(:,i1),p(:,i2),dlips)
      end if

      dlips=0.5d0*p(0,igind)*deg*dctg*dphig*dlips

      if (dlips <= 0.d0) then
       dlips=0.d0
       return
      end if
      if (imom == 1) then
       v=-p(1:,igind)/sqrt(sp+p(0,igind)**2)
       g=sqrt(sp+p(0,igind)**2)/sqsp
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i1))
        call boost_qp(g,v,p(:,i2))
       else
        call boost(g,v,p(:,i1))
        call boost(g,v,p(:,i2))
       end if
      end if

      end

