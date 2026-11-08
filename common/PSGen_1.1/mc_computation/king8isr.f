      subroutine king8isr(i1,i2,mv1,gv1,i3,i4,mv2,gv2,i5,i6,mv3,gv3,
     &                           i7,i8,mv4,gv4,imom,ecm,x,ndim,dlips)
 
      use fourmom
      use kinparams
      use engymom
      use kincuts 

      implicit none

      integer :: i1,i2,i3,i4,i5,i6,i7,i8,imom,ndim
      real(dp) :: mv1,gv1,mv2,gv2,mv3,gv3,mv4,gv4,x(ndim),dlips
      real(xp) :: ecm
      real(dp) :: ds1,ds2,ds3,ds4,ds12,ds34
      real(dp) :: dl1,dl2,dl3,dl4,dl12,dl34

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
      sqsl3=m(i5)+m(i6)
      sqsl4=m(i7)+m(i8)
      sqsl12=sqsl1+sqsl2
      sqsl34=sqsl3+sqsl4
      sqsu12=sqsp-sqsl34

      if (imom /= 1) then
       q1=p(:,i1)+p(:,i2)
       q2=p(:,i3)+p(:,i4)
       q3=p(:,i5)+p(:,i6)
       q4=p(:,i7)+p(:,i8)
       q12=q1+q2
       q34=q3+q4
      end if

      sqsu12=sqsp-sqsl34
      ds12=(sqsu12-sqsl12)*(sqsu12+sqsl12)
      if (imom /= 1) then
       s12=q12(0)**2-q12(1)**2-q12(2)**2-q12(3)**2
      else
       s12=ds12*x(1)+sqsl12**2
      end if
      sqs12=sqrt(s12)

      sqsu34=sqsp-sqs12
      ds34=(sqsu34-sqsl34)*(sqsu34+sqsl34)
      if (imom /= 1) then
       s34=q34(0)**2-q34(1)**2-q34(2)**2-q34(3)**2
      else
       s34=ds34*x(2)+sqsl34**2
      end if
      sqs34=sqrt(s34)

      sqsu1=sqs12-sqsl2
      if (mv1 == 0.d0) then
       if (iqprec == 1) then
        call mapa_qp(imom,x(3),q1,sqscut,sqsl1,sqsu1,s1,sqs1,ds1)
       else
        call mapakk(imom,x(3),q1,sqscut,sqsl1,sqsu1,s1,sqs1,ds1)
       end if
      else
       if (iqprec == 1) then
        call mapb_qp(imom,x(3),q1,mv1,gv1,sqsl1,sqsu1,s1,sqs1,ds1)
       else
        call mapbkk(imom,x(3),q1,mv1,gv1,sqsl1,sqsu1,s1,sqs1,ds1)
       end if
      end if

      sqsu2=sqs12-sqs1
      if (mv2 == 0.d0) then
       if (iqprec == 1) then
        call mapa_qp(imom,x(4),q2,sqscut,sqsl2,sqsu2,s2,sqs2,ds2)
       else
        call mapakk(imom,x(4),q2,sqscut,sqsl2,sqsu2,s2,sqs2,ds2)
       end if
      else
       if (iqprec == 1) then
        call mapb_qp(imom,x(4),q2,mv2,gv2,sqsl2,sqsu2,s2,sqs2,ds2)
       else
        call mapbkk(imom,x(4),q2,mv2,gv2,sqsl2,sqsu2,s2,sqs2,ds2)
       end if
      end if

      sqsu3=sqs34-sqsl4
      if (mv3 == 0.d0) then
       if (iqprec == 1) then
        call mapa_qp(imom,x(5),q3,sqscut,sqsl3,sqsu3,s3,sqs3,ds3)
       else
        call mapakk(imom,x(5),q3,sqscut,sqsl3,sqsu3,s3,sqs3,ds3)
       end if
      else
       if (iqprec == 1) then
        call mapb_qp(imom,x(5),q3,mv3,gv3,sqsl3,sqsu3,s3,sqs3,ds3)
       else
        call mapbkk(imom,x(5),q3,mv3,gv3,sqsl3,sqsu3,s3,sqs3,ds3)
       end if
      end if

      sqsu4=sqs34-sqs3
      if (mv4 == 0.d0) then
       if (iqprec == 1) then
        call mapa_qp(imom,x(6),q4,sqscut,sqsl4,sqsu4,s4,sqs4,ds4)
       else
        call mapakk(imom,x(6),q4,sqscut,sqsl4,sqsu4,s4,sqs4,ds4)
       end if
      else
       if (iqprec == 1) then
        call mapb_qp(imom,x(6),q4,mv4,gv4,sqsl4,sqsu4,s4,sqs4,ds4)
       else
        call mapbkk(imom,x(6),q4,mv4,gv4,sqsl4,sqsu4,s4,sqs4,ds4)
       end if
      end if

      if (iqprec == 1) then
       call kinff_qp(imom,sqs1,m(i1),m(i2),x(7:8),p(:,i1),p(:,i2),dl1)
       call kinff_qp(imom,sqs2,m(i3),m(i4),x(9:10),p(:,i3),p(:,i4),dl2)
       call kinff_qp(imom,sqs3,m(i5),m(i6),x(11:12),p(:,i5),p(:,i6),dl3)
       call kinff_qp(imom,sqs4,m(i7),m(i8),x(13:14),p(:,i7),p(:,i8),dl4)
       call kinff_qp(imom,sqs12,sqs1,sqs2,x(15:16),q1,q2,dl12)
       call kinff_qp(imom,sqs34,sqs3,sqs4,x(17:18),q3,q4,dl34)
       call kinff_qp(imom,sqsp,sqs12,sqs34,x(19:20),q12,q34,dlips)
      else
       call kinff(imom,sqs1,m(i1),m(i2),x(7:8),p(:,i1),p(:,i2),dl1)
       call kinff(imom,sqs2,m(i3),m(i4),x(9:10),p(:,i3),p(:,i4),dl2)
       call kinff(imom,sqs3,m(i5),m(i6),x(11:12),p(:,i5),p(:,i6),dl3)
       call kinff(imom,sqs4,m(i7),m(i8),x(13:14),p(:,i7),p(:,i8),dl4)
       call kinff(imom,sqs12,sqs1,sqs2,x(15:16),q1,q2,dl12)
       call kinff(imom,sqs34,sqs3,sqs4,x(17:18),q3,q4,dl34)
       call kinff(imom,sqsp,sqs12,sqs34,x(19:20),q12,q34,dlips)
      end if

      dlips=0.5d0*p(0,igind)*deg*dctg*dphig
     &     *dlips*ds12*ds34*ds1*ds2*ds3*ds4*dl12*dl34*dl1*dl2*dl3*dl4

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
       v=q3(1:)/q3(0)
       g=q3(0)/sqs3
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i5))
        call boost_qp(g,v,p(:,i6))
       else
        call boost(g,v,p(:,i5))
        call boost(g,v,p(:,i6))
       end if
       if ((p(0,i5)+p(0,i6))**2-(p(1,i5)+p(1,i6))**2-
     &     (p(2,i5)+p(2,i6))**2-(p(3,i5)+p(3,i6))**2 
     &                                         < sqsl3**2) then
        dlips=0.d0
        return
       end if
       v=q4(1:)/q4(0)
       g=q4(0)/sqs4
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i7))
        call boost_qp(g,v,p(:,i8))
       else
        call boost(g,v,p(:,i7))
        call boost(g,v,p(:,i8))
       end if
       if ((p(0,i7)+p(0,i8))**2-(p(1,i7)+p(1,i8))**2-
     &     (p(2,i7)+p(2,i8))**2-(p(3,i7)+p(3,i8))**2 
     &                                         < sqsl4**2) then
        dlips=0.d0
        return
       end if
       v=q12(1:)/q12(0)
       g=q12(0)/sqs12
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
       if ((p(0,i1)+p(0,i2)+p(0,i3)+p(0,i4))**2
     &    -(p(1,i1)+p(1,i2)+p(1,i3)+p(1,i4))**2
     &    -(p(2,i1)+p(2,i2)+p(2,i3)+p(2,i4))**2
     &    -(p(3,i1)+p(3,i2)+p(3,i3)+p(3,i4))**2 < sqsl12**2) then
        dlips=0.d0
        return
       end if
       v=q34(1:)/q34(0)
       g=q34(0)/sqs34
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i5))
        call boost_qp(g,v,p(:,i6))
        call boost_qp(g,v,p(:,i7))
        call boost_qp(g,v,p(:,i8))
       else
        call boost(g,v,p(:,i5))
        call boost(g,v,p(:,i6))
        call boost(g,v,p(:,i7))
        call boost(g,v,p(:,i8))
       end if
       if ((p(0,i5)+p(0,i6)+p(0,i7)+p(0,i8))**2
     &    -(p(1,i5)+p(1,i6)+p(1,i7)+p(1,i8))**2
     &    -(p(2,i5)+p(2,i6)+p(2,i7)+p(2,i8))**2
     &    -(p(3,i5)+p(3,i6)+p(3,i7)+p(3,i8))**2 < sqsl34**2) then
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
        call boost_qp(g,v,p(:,i5))
        call boost_qp(g,v,p(:,i6))
        call boost_qp(g,v,p(:,i7))
        call boost_qp(g,v,p(:,i8))
       else
        call boost(g,v,p(:,i1))
        call boost(g,v,p(:,i2))
        call boost(g,v,p(:,i3))
        call boost(g,v,p(:,i4))
        call boost(g,v,p(:,i5))
        call boost(g,v,p(:,i6))
        call boost(g,v,p(:,i7))
        call boost(g,v,p(:,i8))
       end if
      end if

      end

