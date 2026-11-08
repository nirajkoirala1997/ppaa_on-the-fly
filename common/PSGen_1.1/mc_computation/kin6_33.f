      subroutine kin6_33(i3,i4,i5,mv1,gv1,mv2,gv2,
     &                   i6,i7,i8,mv3,gv3,mv4,gv4,
     &                                  imom,ecm,p,x,ndim,dlips)

      use kinparams

      implicit none
 
      integer :: i3,i4,i5,i6,i7,i8,imom,ndim
      real(dp) :: mv1,gv1,mv2,gv2,mv3,gv3,mv4,gv4,
     &            x(ndim),dlips,ds1,ds2,ds3,ds4,
     &            dl1,dl2,dl3,dl4
      real(xp) :: ecm,p(0:3,nexp),v(3),g
      real(xp), dimension(0:3) :: q1,q2,q3,q4
      real(xp) :: sqsl1,sqsu1,s1,sqs1,sqsl2,sqsu2,s2,sqs2,
     &            sqsl3,sqsu3,s3,sqs3,sqsl4,sqsu4,s4,sqs4

      sqsl1=m(i4)+m(i5)
      sqsl2=sqsl1+m(i3)
      sqsl3=m(i7)+m(i8)
      sqsl4=sqsl3+m(i6)

      sqsu2=ecm-sqsl4

      if (imom /= 1) then
       q1=p(:,i4)+p(:,i5)
       q2=q1+p(:,i3)
       q3=p(:,i7)+p(:,i8)
       q4=q3+p(:,i6)
      end if

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim),q2,mv2,gv2,sqsc,sqsl2,sqsu2,s2,sqs2,
     &                                                           ds2)
      else
       call mapckk(imom,x(ndim),q2,mv2,gv2,sqsc,sqsl2,sqsu2,s2,sqs2,
     &                                                           ds2)
       end if

      sqsu4=ecm-sqs2

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim-1),q4,mv4,gv4,sqsc,sqsl4,sqsu4,s4,sqs4,
     &                                                              ds4)
      else
       call mapckk(imom,x(ndim-1),q4,mv4,gv4,sqsc,sqsl4,sqsu4,s4,sqs4,
     &                                                              ds4)
      end if

      sqsu1=sqs2-m(i3)

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim-2),q1,mv1,gv1,sqsc,sqsl1,sqsu1,s1,sqs1,
     &                                                              ds1)
      else
       call mapckk(imom,x(ndim-2),q1,mv1,gv1,sqsc,sqsl1,sqsu1,s1,sqs1,
     &                                                              ds1)
      end if

      sqsu3=sqs4-m(i6)

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim-3),q3,mv3,gv3,sqsc,sqsl3,sqsu3,s3,sqs3,
     &                                                              ds3)
      else
       call mapckk(imom,x(ndim-3),q3,mv3,gv3,sqsc,sqsl3,sqsu3,s3,sqs3,
     &                                                              ds3)
      end if

      if (iqprec == 1) then
       call kinff_qp(imom,sqs1,m(i4),m(i5),x(1:2),p(:,i4),p(:,i5),dl1)
       call kinff_qp(imom,sqs2,sqs1,m(i3),x(3:4),q1,p(:,i3),dl2)
       call kinff_qp(imom,sqs3,m(i7),m(i8),x(5:6),p(:,i7),p(:,i8),dl3)
       call kinff_qp(imom,sqs4,sqs3,m(i6),x(7:8),q3,p(:,i6),dl4)
       call kinff_qp(imom,ecm,sqs2,sqs4,x(9:10),q2,q4,dlips)
      else
       call kinff(imom,sqs1,m(i4),m(i5),x(1:2),p(:,i4),p(:,i5),dl1)
       call kinff(imom,sqs2,sqs1,m(i3),x(3:4),q1,p(:,i3),dl2)
       call kinff(imom,sqs3,m(i7),m(i8),x(5:6),p(:,i7),p(:,i8),dl3)
       call kinff(imom,sqs4,sqs3,m(i6),x(7:8),q3,p(:,i6),dl4)
       call kinff(imom,ecm,sqs2,sqs4,x(9:10),q2,q4,dlips)
      end if

      dlips=ds1*ds2*ds3*ds4*dlips*dl1*dl2*dl3*dl4

      if (dlips <= 0.d0) then
       dlips=0.d0
       return
      end if

      if (imom == 1) then
       v=q1(1:)/q1(0)
       g=q1(0)/sqs1
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i4))
        call boost_qp(g,v,p(:,i5))
       else
        call boost(g,v,p(:,i4))
        call boost(g,v,p(:,i5))
       end if
       v=q2(1:)/q2(0)
       g=q2(0)/sqs2
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i3))
        call boost_qp(g,v,p(:,i4))
        call boost_qp(g,v,p(:,i5))
       else
        call boost(g,v,p(:,i3))
        call boost(g,v,p(:,i4))
        call boost(g,v,p(:,i5))
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
        call boost_qp(g,v,p(:,i6))
        call boost_qp(g,v,p(:,i7))
        call boost_qp(g,v,p(:,i8))
       else
        call boost(g,v,p(:,i6))
        call boost(g,v,p(:,i7))
        call boost(g,v,p(:,i8))
       end if
      end if

      end
