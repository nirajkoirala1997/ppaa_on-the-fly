      subroutine kin4_22(i3,i4,mv1,gv1,i5,i6,mv2,gv2,
     &                                  imom,ecm,p,x,ndim,dlips)

      use kinparams

      implicit none
 
      integer :: i3,i4,i5,i6,imom,ndim
      real(dp) :: mv1,gv1,mv2,gv2,x(ndim),dlips,ds1,ds2,dl1,dl2
      real(xp) :: ecm,p(0:3,nexp),q1(0:3),q2(0:3),v(3),g
      real(xp) :: sqsl1,sqsu1,s1,sqs1,sqsl2,sqsu2,s2,sqs2

      sqsl1=m(i3)+m(i4)
      sqsl2=m(i5)+m(i6)
      sqsu1=ecm-sqsl2

      if (imom /= 1) then
       q1=p(:,i3)+p(:,i4)
       q2=p(:,i5)+p(:,i6)
      end if

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim),q1,mv1,gv1,sqsc,sqsl1,sqsu1,s1,sqs1,
     &                                                           ds1)
      else
       call mapckk(imom,x(ndim),q1,mv1,gv1,sqsc,sqsl1,sqsu1,s1,sqs1,ds1)
      end if

      sqsu2=ecm-sqs1

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim-1),q2,mv2,gv2,sqsc,sqsl2,sqsu2,s2,sqs2,
     &                                                              ds2)
      else
       call mapckk(imom,x(ndim-1),q2,mv2,gv2,sqsc,sqsl2,sqsu2,s2,sqs2,
     &                                                             ds2)
      end if

      if (iqprec == 1) then
       call kinff_qp(imom,sqs1,m(i3),m(i4),x(1:2),p(:,i3),p(:,i4),dl1)
       call kinff_qp(imom,sqs2,m(i5),m(i6),x(3:4),p(:,i5),p(:,i6),dl2)
       call kinff_qp(imom,ecm,sqs1,sqs2,x(5:6),q1,q2,dlips)
      else
       call kinff(imom,sqs1,m(i3),m(i4),x(1:2),p(:,i3),p(:,i4),dl1)
       call kinff(imom,sqs2,m(i5),m(i6),x(3:4),p(:,i5),p(:,i6),dl2)
       call kinff(imom,ecm,sqs1,sqs2,x(5:6),q1,q2,dlips)
      end if

      dlips=ds1*ds2*dlips*dl1*dl2

      if (dlips <= 0.d0) then
       dlips=0.d0
       return
      end if

      if (imom == 1) then
       v=q1(1:)/q1(0)
       g=q1(0)/sqs1
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i3))
        call boost_qp(g,v,p(:,i4))
       else
        call boost(g,v,p(:,i3))
        call boost(g,v,p(:,i4))
       end if
       v=q2(1:)/q2(0)
       g=q2(0)/sqs2
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i5))
        call boost_qp(g,v,p(:,i6))
       else
        call boost(g,v,p(:,i5))
        call boost(g,v,p(:,i6))
       end if
      end if

      end
