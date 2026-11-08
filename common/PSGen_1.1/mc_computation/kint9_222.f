      subroutine kint9_222(i1,i3,mvt,ctc3,i4,i5,mv1,gv1,i6,i7,mv2,gv2,
     &                    i8,i9,mv3,gv3,i10,i11,imom,ecm,p,x,ndim,dlips)

      use kinparams

      implicit none
 
      integer :: i1,i3,i4,i5,i6,i7,i8,i9,i10,i11,imom,ndim,i
      real(dp) :: mvt,ctc3,mv1,gv1,mv2,gv2,mv3,gv3,x(ndim),dlips,
     &            ds1,ds2,ds4,ds3,ds5,ds6,dsp,dl1,dl2,dl4,dl3,dl5,dl6,
     &            mv12,gv12
      real(dp), parameter :: dphi=2.d0*pi,cps=dphi/4.d0,
     &                       dprec=1.d-6,dprec1=1.d-16
      real(xp) :: ecm,p(0:3,nexp),pp(0:3),q1(0:3),q2(0:3),q4(0:3),
     &            q3(0:3),q5(0:3),q6(0:3),
     &            pp1,pp3,beta3,db3,rb3,act3,ct3,dct3,st3,v(3),g
      real(xp) :: sqslp,sqsup,sp,sqsp,sqsl1,sqsu1,s1,sqs1,
     &            sqsl2,sqsu2,s2,sqs2,sqsl4,sqsu4,s4,sqs4,
     &            sqsl3,sqsu3,s3,sqs3,sqsl5,sqsu5,s5,sqs5,
     &            sqsl6,sqsu6,s6,sqs6
      real(dp) :: phi3

      sqsl1=m(i4)+m(i5)
      sqsl2=m(i6)+m(i7)
      sqsl3=sqsl1+sqsl2
      sqsl4=m(i8)+m(i9)
      sqsl5=m(i10)+m(i11)
      sqsl6=sqsl4+sqsl5
      sqslp=sqsl3+sqsl6
      sqsup=ecm-m(i3)

      if (mv1 >= mv2) then
       mv12=mv1
       gv12=gv1
      else
       mv12=mv2
       gv12=gv2
      end if

      if (imom /= 1) then
       q1=p(:,i4)+p(:,i5)
       q2=p(:,i6)+p(:,i7)
       q3=q1+q2
       q4=p(:,i8)+p(:,i9)
       q5=p(:,i10)+p(:,i11)
       q6=q4+q5
       pp=q3+q6
      end if

      if (iqprec == 1) then
       call mapa_qp(imom,x(1),pp,sqsc,sqslp,sqsup,sp,sqsp,dsp)
      else
       call mapakk(imom,x(1),pp,sqsc,sqslp,sqsup,sp,sqsp,dsp)
      end if
c      if (iqprec == 1) then
c       call mapc_qp(imom,x(1),pp,1.d0,zero,sqsc,sqslp,sqsup,sp,
c     &                                                       sqsp,dsp)
c      else
c       call mapckk(imom,x(1),pp,1.d0,zero,sqsc,sqslp,sqsup,sp,
c     &                                                       sqsp,dsp)
c      end if

      pp1=sqrt(p(1,i1)**2+p(2,i1)**2+p(3,i1)**2)

c Generate ct3 in [-ctc3,ctc3]:
      if (imom == 1) then
       p(0,i3)=(ecm**2+m(i3)**2-sp)/(2.d0*ecm)
      end if
      pp3=sqrt(p(0,i3)**2-m(i3)**2)
      beta3=2.d0*pp1*pp3
     &                /(mvt**2-m(i1)**2-m(i3)**2+2.d0*p(0,i1)*p(0,i3))
      if (iqprec == 1) then
       db3=1.q0-beta3
      else
       db3=1.d0-beta3
      end if
      if (abs(db3) < dprec .and. iqprec /= 1) then
       db3=-((m(i1)+m(i3))**2-mvt**2)/(2.d0*p(0,i1)*p(0,i3))
     &                          +0.5d0*(m(i1)/p(0,i1)+m(i3)/p(0,i3))**2
      end if
      if (iqprec == 1) then
       rb3=(1.q0-(1.q0-db3)*ctc3)/(1.q0+(1.q0-db3)*ctc3)
      else
       rb3=(1.d0-(1.d0-db3)*ctc3)/(1.d0+(1.d0-db3)*ctc3)
      end if
      if (imom == 1) then
       if (iqprec == 1) then
        act3=(1.q0+(1.q0-db3)*ctc3)*rb3**x(2)
        ct3=(1.q0-act3)/beta3
       else
        act3=(1.d0+(1.d0-db3)*ctc3)*rb3**x(2)
        ct3=(1.d0-act3)/beta3
       end if
      else
       ct3=p(3,i3)/pp3
       if (iqprec == 1) then
        act3=1.q0-beta3*ct3
       else
        act3=1.d0-beta3*ct3
       end if
      end if
      if (iqprec == 1) then
       dct3=-(1.q0-(1.q0-db3)*ct3)*log(rb3)/beta3
      else
       dct3=-(1.d0-(1.d0-db3)*ct3)*log(rb3)/beta3
      end if
      if (imom == 1) then
       if (iqprec == 1) then
        st3=sqrt(1.q0-ct3**2)
       else
        if (1.d0-ct3 < dprec1) then 
         st3=sqrt(2.d0*(act3-db3))
        else
         st3=sqrt(1.d0-ct3**2)
        end if
       end if
       phi3=dphi*x(3)
       p(1:3,i3)=(/pp3*st3*cos(phi3),pp3*st3*sin(phi3),pp3*ct3/)
       pp(0)=ecm-p(0,i3)
       pp(1:3)=-p(1:3,i3)
      end if

      sqsu3=sqsp-sqsl6

      if (iqprec == 1) then
       call mapc_qp(imom,x(4),q3,mv12,gv12,sqsc,sqsl3,sqsu3,s3,sqs3,ds3)
       sqsu6=sqsp-sqs3
       call mapc_qp(imom,x(5),q6,1.d0,zero,sqsc,sqsl6,sqsu6,s6,sqs6,ds6)
       sqsu1=sqs3-sqsl2
       call mapc_qp(imom,x(6),q1,mv1,gv1,sqsc,sqsl1,sqsu1,s1,sqs1,ds1)
       sqsu2=sqs3-sqs1
       call mapc_qp(imom,x(7),q2,mv2,gv2,sqsc,sqsl2,sqsu2,s2,sqs2,ds2)
       sqsu4=sqs6-sqsl5
       call mapc_qp(imom,x(8),q4,mv3,gv3,sqsc,sqsl4,sqsu4,s4,sqs4,ds4)
       sqsu5=sqs6-sqs4
       call mapc_qp(imom,x(9),q5,1.d0,zero,sqsc,sqsl5,sqsu5,s5,sqs5,ds5)
      else
       call mapckk(imom,x(4),q3,mv12,gv12,sqsc,sqsl3,sqsu3,s3,sqs3,ds3)
       sqsu6=sqsp-sqs3
       call mapckk(imom,x(5),q6,1.d0,zero,sqsc,sqsl6,sqsu6,s6,sqs6,ds6)
       sqsu1=sqs3-sqsl2
       call mapckk(imom,x(6),q1,mv1,gv1,sqsc,sqsl1,sqsu1,s1,sqs1,ds1)
       sqsu2=sqs3-sqs1
       call mapckk(imom,x(7),q2,mv2,gv2,sqsc,sqsl2,sqsu2,s2,sqs2,ds2)
       sqsu4=sqs6-sqsl5
       call mapckk(imom,x(8),q4,mv3,gv3,sqsc,sqsl4,sqsu4,s4,sqs4,ds4)
       sqsu5=sqs6-sqs4
       call mapckk(imom,x(9),q5,1.d0,zero,sqsc,sqsl5,sqsu5,s5,sqs5,ds5)
      end if

      if (iqprec == 1) then
       call kinff_qp(imom,sqs1,m(i4),m(i5),x(10:11),p(:,i4),p(:,i5),dl1)
       call kinff_qp(imom,sqs2,m(i6),m(i7),x(12:13),p(:,i6),p(:,i7),dl2)
       call kinff_qp(imom,sqs3,sqs1,sqs2,x(14:15),q1,q2,dl3)
       call kinff_qp(imom,sqs4,m(i8),m(i9),x(16:17),p(:,i8),p(:,i9),dl4)
       call kinff_qp(imom,sqs5,m(i10),m(i11),x(18:19),p(:,i10),p(:,i11),
     &                                                              dl5)
       call kinff_qp(imom,sqs6,sqs4,sqs5,x(20:21),q4,q5,dl6)
       call kinff_qp(imom,sqsp,sqs3,sqs6,x(22:23),q3,q6,dlips)
      else
       call kinff(imom,sqs1,m(i4),m(i5),x(10:11),p(:,i4),p(:,i5),dl1)
       call kinff(imom,sqs2,m(i6),m(i7),x(12:13),p(:,i6),p(:,i7),dl2)
       call kinff(imom,sqs3,sqs1,sqs2,x(14:15),q1,q2,dl3)
       call kinff(imom,sqs4,m(i8),m(i9),x(16:17),p(:,i8),p(:,i9),dl4)
       call kinff(imom,sqs5,m(i10),m(i11),x(18:19),p(:,i10),p(:,i11),
     &                                                              dl5)
       call kinff(imom,sqs6,sqs4,sqs5,x(20:21),q4,q5,dl6)
       call kinff(imom,sqsp,sqs3,sqs6,x(22:23),q3,q6,dlips)
      end if

      dlips=dsp*cps*dct3*pp3/ecm*ds1*ds2*ds3*ds4*ds5*ds6
     &                          *dl1*dl2*dl3*dl4*dl5*dl6*dlips

      if (dlips > 0.d0) then
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
         call boost_qp(g,v,p(:,i6))
         call boost_qp(g,v,p(:,i7))
        else
         call boost(g,v,p(:,i6))
         call boost(g,v,p(:,i7))
        end if
        v=q3(1:)/q3(0)
        g=q3(0)/sqs3
        if (iqprec == 1) then
         call boost_qp(g,v,p(:,i4))
         call boost_qp(g,v,p(:,i5))
         call boost_qp(g,v,p(:,i6))
         call boost_qp(g,v,p(:,i7))
        else
         call boost(g,v,p(:,i4))
         call boost(g,v,p(:,i5))
         call boost(g,v,p(:,i6))
         call boost(g,v,p(:,i7))
        end if
        v=q4(1:)/q4(0)
        g=q4(0)/sqs4
        if (iqprec == 1) then
         call boost_qp(g,v,p(:,i8))
         call boost_qp(g,v,p(:,i9))
        else
         call boost(g,v,p(:,i8))
         call boost(g,v,p(:,i9))
        end if
        v=q5(1:)/q5(0)
        g=q5(0)/sqs5
        if (iqprec == 1) then
         call boost_qp(g,v,p(:,i10))
         call boost_qp(g,v,p(:,i11))
        else
         call boost(g,v,p(:,i10))
         call boost(g,v,p(:,i11))
        end if
        v=q6(1:)/q6(0)
        g=q6(0)/sqs6
        if (iqprec == 1) then
         call boost_qp(g,v,p(:,i8))
         call boost_qp(g,v,p(:,i9))
         call boost_qp(g,v,p(:,i10))
         call boost_qp(g,v,p(:,i11))
        else
         call boost(g,v,p(:,i8))
         call boost(g,v,p(:,i9))
         call boost(g,v,p(:,i10))
         call boost(g,v,p(:,i11))
        end if
        v=pp(1:)/pp(0)
        g=pp(0)/sqsp
        if (iqprec == 1) then
         call boost_qp(g,v,p(:,i4))
         call boost_qp(g,v,p(:,i5))
         call boost_qp(g,v,p(:,i6))
         call boost_qp(g,v,p(:,i7))
         call boost_qp(g,v,p(:,i8))
         call boost_qp(g,v,p(:,i9))
         call boost_qp(g,v,p(:,i10))
         call boost_qp(g,v,p(:,i11))
        else
         call boost(g,v,p(:,i4))
         call boost(g,v,p(:,i5))
         call boost(g,v,p(:,i6))
         call boost(g,v,p(:,i7))
         call boost(g,v,p(:,i8))
         call boost(g,v,p(:,i9))
         call boost(g,v,p(:,i10))
         call boost(g,v,p(:,i11))
        end if
c        write(*,*) 'p_f:',p(:,3)+p(:,4)+p(:,5)+p(:,6)+p(:,7)+p(:,8)
c     &                   +p(:,9)
c        do i=3,6
c         write(*,*)'m:',i,
c     &   sqrt(abs(p(0,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2))-m(i),
c     &               m(i)
c        end do
       end if
      else
       dlips=0.d0
      end if

      end
