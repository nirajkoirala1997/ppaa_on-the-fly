      subroutine king4isr(i1,i2,mv1,gv1,i3,i4,mv2,gv2,
     &                                imom,ecm,x,ndim,dlips)
 
      use fourmom
      use kinparams
      use engymom
      use kincuts 

      implicit none

      integer :: i1,i2,i3,i4,imom,ndim
      real(dp) :: mv1,gv1,mv2,gv2,x(ndim),dlips
      real(xp) :: ecm
      real(dp) :: ds1,ds2
      real(dp) :: dl1,dl2

c Energy and angles of the initial state photon in the c.m.s.:

      real(dp) :: neg,deg,nctg,ctg,stg,phig,dphig=2.d0*pi

      save dphig

      msum=sum(m(3:))
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

      sqsl1=m(i1)+m(i2)
      sqsl2=m(i3)+m(i4)
      sqsu1=sqsp-sqsl2

      if (imom /= 1) then
       q1=p(:,i1)+p(:,i2)
       q2=p(:,i3)+p(:,i4)
      end if


      if (mv1 == 0.d0) then
       if (iqprec == 1) then
        call mapa_qp(imom,x(1),q1,sqscut,sqsl1,sqsu1,s1,sqs1,ds1)
       else
        call mapakk(imom,x(1),q1,sqscut,sqsl1,sqsu1,s1,sqs1,ds1)
       end if
      else
       if (iqprec == 1) then
        call mapb_qp(imom,x(1),q1,mv1,gv1,sqsl1,sqsu1,s1,sqs1,ds1)
       else
        call mapbkk(imom,x(1),q1,mv1,gv1,sqsl1,sqsu1,s1,sqs1,ds1)
       end if
      end if

      sqsu2=sqsp-sqs1

      if (mv2 == 0.d0) then
       if (iqprec == 1) then
        call mapa_qp(imom,x(2),q2,sqscut,sqsl2,sqsu2,s2,sqs2,ds2)
       else
        call mapakk(imom,x(2),q2,sqscut,sqsl2,sqsu2,s2,sqs2,ds2)
       end if
      else
       if (iqprec == 1) then
        call mapb_qp(imom,x(2),q2,mv2,gv2,sqsl2,sqsu2,s2,sqs2,ds2)
       else
        call mapbkk(imom,x(2),q2,mv2,gv2,sqsl2,sqsu2,s2,sqs2,ds2)
       end if
      end if

      if (iqprec == 1) then
       call kinff_qp(imom,sqs1,m(i1),m(i2),x(3:4),p(:,i1),p(:,i2),dl1)
       call kinff_qp(imom,sqs2,m(i3),m(i4),x(5:6),p(:,i3),p(:,i4),dl2)
       call kinff_qp(imom,sqsp,sqs1,sqs2,x(7:8),q1,q2,dlips)
      else
       call kinff(imom,sqs1,m(i1),m(i2),x(3:4),p(:,i1),p(:,i2),dl1)
       call kinff(imom,sqs2,m(i3),m(i4),x(5:6),p(:,i3),p(:,i4),dl2)
       call kinff(imom,sqsp,sqs1,sqs2,x(7:8),q1,q2,dlips)
      end if

      dlips=0.5d0*p(0,igind)*deg*dctg*dphig
     &     *dlips*ds1*ds2*dl1*dl2

      if (dlips <= 0.d0) then
       dlips=0.d0
       return
      end if
      if (imom == 1) then
       v=q1(1:)/q1(0)
       g=q1(0)/sqs1
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i1))
        call boost_qp(g,v,p(:,i2))
       else
        call boost(g,v,p(:,i1))
        call boost(g,v,p(:,i2))
       end if
       if ((p(0,i1)+p(0,i2))**2-(p(1,i1)+p(1,i2))**2-
     &     (p(2,i1)+p(2,i2))**2-(p(3,i1)+p(3,i2))**2 
     &                                         < sqsl1**2) then
        dlips=0.d0
        return
       end if
       v=q2(1:)/q2(0)
       g=q2(0)/sqs2
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i3))
        call boost_qp(g,v,p(:,i4))
       else
        call boost(g,v,p(:,i3))
        call boost(g,v,p(:,i4))
       end if
       if ((p(0,i3)+p(0,i4))**2-(p(1,i3)+p(1,i4))**2-
     &     (p(2,i3)+p(2,i4))**2-(p(3,i3)+p(3,i4))**2 
     &                                         < sqsl2**2) then
        dlips=0.d0
        return
       end if
       v=-p(1:,igind)/sqrt(sp+p(0,igind)**2)
       g=sqrt(sp+p(0,igind)**2)/sqsp
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i1))
        call boost_qp(g,v,p(:,i2))
        call boost_qp(g,v,p(:,i3))
        call boost_qp(g,v,p(:,i4))
       else
        call boost(g,v,p(:,i1))
        call boost(g,v,p(:,i2))
        call boost(g,v,p(:,i3))
        call boost(g,v,p(:,i4))
       end if
      end if

      end

