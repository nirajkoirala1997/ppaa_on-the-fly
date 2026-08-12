      subroutine kin7_222(i3,i4,mv1,gv1,i5,i6,mv2,gv2,i7,i8,mv3,gv3,i9,
     &                                         imom,ecm,p,x,ndim,dlips)

      use kinparams

      implicit none
 
      integer :: i3,i4,i5,i6,i7,i8,i9,imom,ndim
      real(dp) :: mv1,gv1,mv2,gv2,mv3,gv3,
     &            x(ndim),dlips,ds1p,ds1,ds2,ds3,ds4,
     &            dl1p,dl1,dl2,dl3,dl4
      real(xp) :: ecm,p(0:3,nexp),v(3),g
      real(xp), dimension(0:3) :: q1p,q1,q2,q3,q4
      real(xp) :: sqsl1p,sqsu1p,s1p,sqs1p,sqsl1,sqsu1,s1,sqs1,
     &            sqsl2,sqsu2,s2,sqs2,
     &            sqsl3,sqsu3,s3,sqs3,sqsl4,sqsu4,s4,sqs4

      sqsl1p=m(i3)+m(i4)
      sqsl1=sqsl1p+m(i9)
      sqsl2=m(i5)+m(i6)
      sqsl3=m(i7)+m(i8)
      sqsl4=sqsl2+sqsl3

      if (imom /= 1) then
       q1p=p(:,i3)+p(:,i4)
       q1=q1p+p(:,i9)
       q2=p(:,i5)+p(:,i6)
       q3=p(:,i7)+p(:,i8)
       q4=q2+q3
      end if

      sqsu1=ecm-sqsl4

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim),q1,mv1,gv1,sqsc,sqsl1,sqsu1,s1,sqs1,
     &                                                           ds1)
      else
       call mapckk(imom,x(ndim),q1,mv1,gv1,sqsc,sqsl1,sqsu1,s1,sqs1,ds1)
      end if

      sqsu1p=sqs1-m(i9)

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim-1),q1p,mv1,gv1,sqsc,sqsl1p,sqsu1p,
     &                                                  s1p,sqs1p,ds1p)
      else
       call mapckk(imom,x(ndim-1),q1p,mv1,gv1,sqsc,sqsl1p,sqsu1p,
     &                                                  s1p,sqs1p,ds1p)
      end if

      sqsu4=ecm-sqs1

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim-2),q4,mv2+mv3,zero,sqsc,sqsl4,sqsu4,
     &                                                      s4,sqs4,ds4)
      else
       call mapckk(imom,x(ndim-2),q4,mv2+mv3,zero,sqsc,sqsl4,sqsu4,
     &                                                      s4,sqs4,ds4)
      end if

      sqsu2=sqs4-sqsl3

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim-3),q2,mv2,gv2,sqsc,sqsl2,sqsu2,s2,sqs2,
     &                                                             ds2)
      else
       call mapckk(imom,x(ndim-3),q2,mv2,gv2,sqsc,sqsl2,sqsu2,s2,sqs2,
     &                                                            ds2)
      end if

      sqsu3=sqs4-sqs2

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim-4),q3,mv3,gv3,sqsc,sqsl3,sqsu3,s3,sqs3,
     &                                                             ds3)
      else
       call mapckk(imom,x(ndim-4),q3,mv3,gv3,sqsc,sqsl3,sqsu3,s3,sqs3,
     &                                                             ds3)
      end if


      if (iqprec == 1) then
       call kinff_qp(imom,sqs1p,m(i3),m(i4),x(1:2),p(:,i3),p(:,i4),dl1p)
       call kinff_qp(imom,sqs1,sqs1p,m(i9),x(3:4),q1p,p(:,i9),dl1)
       call kinff_qp(imom,sqs2,m(i5),m(i6),x(5:6),p(:,i5),p(:,i6),dl2)
       call kinff_qp(imom,sqs3,m(i7),m(i8),x(7:8),p(:,i7),p(:,i8),dl3)
       call kinff_qp(imom,sqs4,sqs2,sqs3,x(9:10),q2,q3,dl4)
       call kinff_qp(imom,ecm,sqs1,sqs4,x(11:12),q1,q4,dlips)
      else
       call kinff(imom,sqs1p,m(i3),m(i4),x(1:2),p(:,i3),p(:,i4),dl1p)
       call kinff(imom,sqs1,sqs1p,m(i9),x(3:4),q1p,p(:,i9),dl1)
       call kinff(imom,sqs2,m(i5),m(i6),x(5:6),p(:,i5),p(:,i6),dl2)
       call kinff(imom,sqs3,m(i7),m(i8),x(7:8),p(:,i7),p(:,i8),dl3)
       call kinff(imom,sqs4,sqs2,sqs3,x(9:10),q2,q3,dl4)
       call kinff(imom,ecm,sqs1,sqs4,x(11:12),q1,q4,dlips)
      end if

      dlips=ds1p*ds1*ds2*ds3*ds4*dlips*dl1p*dl1*dl2*dl3*dl4

      if (dlips <= 0.d0) then
       dlips=0.d0
       return
      end if

      if (imom == 1) then
       v=q1p(1:)/q1p(0)
       g=q1p(0)/sqs1p
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i3))
        call boost_qp(g,v,p(:,i4))
       else
        call boost(g,v,p(:,i3))
        call boost(g,v,p(:,i4))
       end if
       v=q1(1:)/q1(0)
       g=q1(0)/sqs1
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i3))
        call boost_qp(g,v,p(:,i4))
        call boost_qp(g,v,p(:,i9))
       else
        call boost(g,v,p(:,i3))
        call boost(g,v,p(:,i4))
        call boost(g,v,p(:,i9))
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
       v=q3(1:)/q3(0)
       g=q3(0)/sqs3
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i7))
        call boost_qp(g,v,p(:,i8))
       else
        call boost(g,v,p(:,i7))
        call boost(g,v,p(:,i8))
       end if
       v=q4(1:)/q4(0)
       g=q4(0)/sqs4
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
c       include "test.f"
      end if

      end
