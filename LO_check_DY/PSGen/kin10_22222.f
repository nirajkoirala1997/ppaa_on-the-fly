      subroutine kin10_22222(i3,i4,mv1,gv1,i5,i6,mv2,gv2,
     &                     i7,i8,mv3,gv3,i9,i10,mv4,gv4,i11,i12,mv5,gv5,
     &                                    imom,ecm,p,x,ndim,dlips)

      use kinparams

      implicit none
 
      integer :: i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,imom,ndim
      real(dp) :: mv1,gv1,mv2,gv2,mv3,gv3,mv4,gv4,mv5,gv5,
     &            x(ndim),dlips,ds1,ds2,ds3,ds4,ds5,ds6,ds7,ds8,
     &            dl1,dl2,dl3,dl4,dl5,dl6,dl7,dl8
      real(xp) :: ecm,p(0:3,nexp),v(3),g
      real(xp), dimension(0:3) :: q1,q2,q3,q4,q5,q6,q7,q8
      real(xp) :: sqsl1,sqsu1,s1,sqs1,sqsl2,sqsu2,s2,sqs2,
     &            sqsl3,sqsu3,s3,sqs3,sqsl4,sqsu4,s4,sqs4,
     &            sqsl5,sqsu5,s5,sqs5,sqsl6,sqsu6,s6,sqs6,
     &            sqsl7,sqsu7,s7,sqs7,sqsl8,sqsu8,s8,sqs8

      sqsl1=m(i3)+m(i4)
      sqsl2=m(i5)+m(i6)
      sqsl3=sqsl1+sqsl2
      sqsl4=m(i7)+m(i8)
      sqsl5=m(i9)+m(i10)
      sqsl6=sqsl4+sqsl5
      sqsl7=m(i11)+m(i12)
      sqsl8=sqsl6+sqsl7

      if (imom /= 1) then
       q1=p(:,i3)+p(:,i4)
       q2=p(:,i5)+p(:,i6)
       q3=q1+q2
       q4=p(:,i7)+p(:,i8)
       q5=p(:,i9)+p(:,i10)
       q6=q4+q5
       q7=p(:,i11)+p(:,i12)
       q8=q6+q7
      end if

      sqsu3=ecm-sqsl8

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim-7),q3,mv1+mv2,zero,sqsc,sqsl3,sqsu3,
     &                                                   s3,sqs3,ds3)
      else
       call mapckk(imom,x(ndim-7),q3,mv1+mv2,zero,sqsc,sqsl3,sqsu3,
     &                                                   s3,sqs3,ds3)
      end if

      sqsu8=ecm-sqs3

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim),q8,mv1+mv2,zero,sqsc,sqsl8,sqsu8,
     &                                                   s8,sqs8,ds8)
      else
       call mapckk(imom,x(ndim),q8,mv1+mv2,zero,sqsc,sqsl8,sqsu8,
     &                                                   s8,sqs8,ds8)
      end if

      sqsu6=sqs8-sqsl7

c      if (iqprec == 1) then
c       call mapc_qp(imom,x(ndim-1),q6,mv3+mv4,zero,sqsc,sqsl6,sqsu6,
c     &                                                      s6,sqs6,ds6)
c      else
c       call mapckk(imom,x(ndim-1),q6,mv3+mv4,zero,sqsc,sqsl6,sqsu6,
c     &                                                      s6,sqs6,ds6)
c      end if

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim-1),q6,mv3,gv3,sqsc,sqsl6,sqsu6,
     &                                                      s6,sqs6,ds6)
      else
       call mapckk(imom,x(ndim-1),q6,mv3,gv3,sqsc,sqsl6,sqsu6,
     &                                                      s6,sqs6,ds6)
      end if

      sqsu1=sqs3-sqsl2

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim-2),q1,mv1,gv1,sqsc,sqsl1,sqsu1,s1,sqs1,
     &                                                              ds1)
      else
       call mapckk(imom,x(ndim-2),q1,mv1,gv1,sqsc,sqsl1,sqsu1,s1,sqs1,
     &                                                              ds1)
      end if

      sqsu2=sqs3-sqs1

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim-3),q2,mv2,gv2,sqsc,sqsl2,sqsu2,s2,sqs2,
     &                                                              ds2)
      else
       call mapckk(imom,x(ndim-3),q2,mv2,gv2,sqsc,sqsl2,sqsu2,s2,sqs2,
     &                                                              ds2)
      end if

      sqsu4=sqs6-sqsl5

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim-4),q4,mv3,gv3,sqsc,sqsl4,sqsu4,s4,sqs4,
     &                                                              ds4)
      else
       call mapckk(imom,x(ndim-4),q4,mv3,gv3,sqsc,sqsl4,sqsu4,s4,sqs4,
     &                                                              ds4)
      end if

      sqsu5=sqs6-sqs4

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim-5),q5,mv4,gv4,sqsc,sqsl5,sqsu5,s5,sqs5,
     &                                                              ds5)
      else
       call mapckk(imom,x(ndim-5),q5,mv4,gv4,sqsc,sqsl5,sqsu5,s5,sqs5,
     &                                                              ds5)
      end if

      sqsu7=sqs8-sqs6

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim-6),q7,mv5,gv5,sqsc,sqsl7,sqsu7,s7,sqs7,
     &                                                              ds7)
      else
       call mapckk(imom,x(ndim-6),q7,mv5,gv5,sqsc,sqsl7,sqsu7,s7,sqs7,
     &                                                              ds7)
      end if

      if (iqprec == 1) then
       call kinff_qp(imom,sqs1,m(i3),m(i4),x(1:2),p(:,i3),p(:,i4),dl1)
       call kinff_qp(imom,sqs2,m(i5),m(i6),x(3:4),p(:,i5),p(:,i6),dl2)
       call kinff_qp(imom,sqs4,m(i7),m(i8),x(5:6),p(:,i7),p(:,i8),dl3)
       call kinff_qp(imom,sqs5,m(i9),m(i10),x(7:8),p(:,i9),p(:,i10),dl4)
       call kinff_qp(imom,sqs3,sqs1,sqs2,x(9:10),q1,q2,dl5)
       call kinff_qp(imom,sqs6,sqs4,sqs5,x(11:12),q4,q5,dl6)
       call kinff_qp(imom,sqs7,m(i11),m(i12),x(13:14),
     &                                       p(:,i11),p(:,i12),dl7)
       call kinff_qp(imom,sqs8,sqs6,sqs7,x(15:16),q6,q7,dl8)
       call kinff_qp(imom,ecm,sqs3,sqs8,x(17:18),q3,q8,dlips)
      else
       call kinff(imom,sqs1,m(i3),m(i4),x(1:2),p(:,i3),p(:,i4),dl1)
       call kinff(imom,sqs2,m(i5),m(i6),x(3:4),p(:,i5),p(:,i6),dl2)
       call kinff(imom,sqs4,m(i7),m(i8),x(5:6),p(:,i7),p(:,i8),dl3)
       call kinff(imom,sqs5,m(i9),m(i10),x(7:8),p(:,i9),p(:,i10),dl4)
       call kinff(imom,sqs3,sqs1,sqs2,x(9:10),q1,q2,dl5)
       call kinff(imom,sqs6,sqs4,sqs5,x(11:12),q4,q5,dl6)
       call kinff(imom,sqs7,m(i11),m(i12),x(13:14),
     &                                       p(:,i11),p(:,i12),dl7)
       call kinff(imom,sqs8,sqs6,sqs7,x(15:16),q6,q7,dl8)
       call kinff(imom,ecm,sqs3,sqs8,x(17:18),q3,q8,dlips)
      end if

      dlips=ds1*ds2*ds3*ds4*ds5*ds6*ds7*ds8
     &     *dlips*dl1*dl2*dl3*dl4*dl5*dl6*dl7*dl8

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
       v=q3(1:)/q3(0)
       g=q3(0)/sqs3
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i3))
        call boost_qp(g,v,p(:,i4))
        call boost_qp(g,v,p(:,i5))
        call boost_qp(g,v,p(:,i6))
       else
        call boost(g,v,p(:,i3))
        call boost(g,v,p(:,i4))
        call boost(g,v,p(:,i5))
        call boost(g,v,p(:,i6))
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
       v=q5(1:)/q5(0)
       g=q5(0)/sqs5
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i9))
        call boost_qp(g,v,p(:,i10))
       else
        call boost(g,v,p(:,i9))
        call boost(g,v,p(:,i10))
       end if
       v=q6(1:)/q6(0)
       g=q6(0)/sqs6
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i7))
        call boost_qp(g,v,p(:,i8))
        call boost_qp(g,v,p(:,i9))
        call boost_qp(g,v,p(:,i10))
       else
        call boost(g,v,p(:,i7))
        call boost(g,v,p(:,i8))
        call boost(g,v,p(:,i9))
        call boost(g,v,p(:,i10))
       end if
       v=q7(1:)/q7(0)
       g=q7(0)/sqs7
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i11))
        call boost_qp(g,v,p(:,i12))
       else
        call boost(g,v,p(:,i11))
        call boost(g,v,p(:,i12))
       end if
       v=q8(1:)/q8(0)
       g=q8(0)/sqs8
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i7))
        call boost_qp(g,v,p(:,i8))
        call boost_qp(g,v,p(:,i9))
        call boost_qp(g,v,p(:,i10))
        call boost_qp(g,v,p(:,i11))
        call boost_qp(g,v,p(:,i12))
       else
        call boost(g,v,p(:,i7))
        call boost(g,v,p(:,i8))
        call boost(g,v,p(:,i9))
        call boost(g,v,p(:,i10))
        call boost(g,v,p(:,i11))
        call boost(g,v,p(:,i12))
       end if
      end if

      end
