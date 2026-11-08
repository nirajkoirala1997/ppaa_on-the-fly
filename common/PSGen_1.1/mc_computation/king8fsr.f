      subroutine king8fsr(i1,i2,mv1,gv1,i3,i4,mv2,gv2,i5,i6,mv3,gv3,
     &                           i7,i8,mv4,gv4,imom,ecm,x,ndim,dlips)
 
      use fourmom
      use kinparams
      use engymom
      use kincuts 

      implicit none

      integer :: i1,i2,i3,i4,i5,i6,i7,i8,imom,ndim
      real(dp) :: mv1,gv1,mv2,gv2,mv3,gv3,mv4,gv4,x(ndim),dlips
      real(xp) :: ecm
      real(dp) :: ds1,ds2,ds3,ds4,ds34,ds234
      real(dp) :: dl1,dl2,dl3,dl4,dl34,dl234

      real(dp) :: ecutp

      sqsl1=m(i1)+m(i2)
      sqsl2=m(i3)+m(i4)
      sqsl3=m(i5)+m(i6)
      sqsl4=m(i7)+m(i8)
      sqsl34=sqsl3+sqsl4
      sqsl234=sqsl2+sqsl34

      if (imom /= 1) then
       q1=p(:,i1)+p(:,i2)+p(:,igind)
       q2=p(:,i3)+p(:,i4)
       q3=p(:,i5)+p(:,i6)
       q4=p(:,i7)+p(:,i8)
       q34=q3+q4
       q234=q2+q34
      end if

      sqsu1=ecm-sqsl234
      if (mv1 == 0.d0) then
       if (iqprec == 1) then
        call mapa_qp(imom,x(ndim),q1,sqscut,sqsl1,sqsu1,s1,sqs1,ds1)
       else
        call mapakk(imom,x(ndim),q1,sqscut,sqsl1,sqsu1,s1,sqs1,ds1)
       end if
      else
       if (iqprec == 1) then
        call mapb_qp(imom,x(ndim),q1,mv1,gv1,sqsl1,sqsu1,s1,sqs1,ds1)
       else
        call mapbkk(imom,x(ndim),q1,mv1,gv1,sqsl1,sqsu1,s1,sqs1,ds1)
       end if
      end if

      sqsu234=ecm-sqs1
      ds234=(sqsu234-sqsl234)*(sqsu234+sqsl234)
      if (imom /= 1) then
       s234=q234(0)**2-q234(1)**2-q234(2)**2-q234(3)**2
      else
       s234=ds234*x(ndim-1)+sqsl234**2
      end if
      sqs234=sqrt(s234)

      sqsu2=sqs234-sqsl34
      if (mv2 == 0.d0) then
       if (iqprec == 1) then
        call mapa_qp(imom,x(ndim-2),q2,sqscut,sqsl2,sqsu2,s2,sqs2,ds2)
       else
        call mapakk(imom,x(ndim-2),q2,sqscut,sqsl2,sqsu2,s2,sqs2,ds2)
       end if
      else
       if (iqprec == 1) then
        call mapb_qp(imom,x(ndim-2),q2,mv2,gv2,sqsl2,sqsu2,s2,sqs2,ds2)
       else
        call mapbkk(imom,x(ndim-2),q2,mv2,gv2,sqsl2,sqsu2,s2,sqs2,ds2)
       end if
      end if

      sqsu34=sqs234-sqs2
      ds34=(sqsu34-sqsl34)*(sqsu34+sqsl34)
      if (imom /= 1) then
       s34=q34(0)**2-q34(1)**2-q34(2)**2-q34(3)**2
      else
       s34=ds34*x(ndim-3)+sqsl34**2
      end if
      sqs34=sqrt(s34)

      sqsu3=sqs34-sqsl4
      if (mv3 == 0.d0) then
       if (iqprec == 1) then
        call mapa_qp(imom,x(ndim-4),q3,sqscut,sqsl3,sqsu3,s3,sqs3,ds3)
       else
        call mapakk(imom,x(ndim-4),q3,sqscut,sqsl3,sqsu3,s3,sqs3,ds3)
       end if
      else
       if (iqprec == 1) then
        call mapb_qp(imom,x(ndim-4),q3,mv3,gv3,sqsl3,sqsu3,s3,sqs3,ds3)
       else
        call mapbkk(imom,x(ndim-4),q3,mv3,gv3,sqsl3,sqsu3,s3,sqs3,ds3)
       end if
      end if

      sqsu4=sqs34-sqs3
      if (mv4 == 0.d0) then
       if (iqprec == 1) then
        call mapa_qp(imom,x(ndim-5),q4,sqscut,sqsl4,sqsu4,s4,sqs4,ds4)
       else
        call mapakk(imom,x(ndim-5),q4,sqscut,sqsl4,sqsu4,s4,sqs4,ds4)
       end if
      else
       if (iqprec == 1) then
        call mapb_qp(imom,x(ndim-5),q4,mv4,gv4,sqsl4,sqsu4,s4,sqs4,ds4)
       else
        call mapbkk(imom,x(ndim-5),q4,mv4,gv4,sqsl4,sqsu4,s4,sqs4,ds4)
       end if
      end if

      if (iqprec == 1) then
       call kinff_qp(imom,ecm,sqs1,sqs234,x(ndim-7:ndim-6),q1,q234,
     &                                                         dlips)
      else
       call kinff(imom,ecm,sqs1,sqs234,x(ndim-7:ndim-6),q1,q234,dlips)
      end if

      qq1=sqrt(q1(1)**2+q1(2)**2+q1(3)**2)
      b34g=qq1/q1(0)
      ecutp=ecut0*sqrt((1.d0-b34g)/(1.d0+b34g))
      if (imom /= 1) then
       e34g=sqrt(sqs1**2+qq1**2)
       v34g=q1(1:)/e34g
       g34g=e34g/sqs1
       e3=g34g
     &    *(p(0,i1)-v34g(1)*p(1,i1)-v34g(2)*p(2,i1)-v34g(3)*p(3,i1))
       e7=g34g*(p(0,igind)
     &   -v34g(1)*p(1,igind)-v34g(2)*p(2,igind)-v34g(3)*p(3,igind))
      end if

      if (iqprec == 1) then
       call kinffg_qp(imom,sqs1,m(i1),m(i2),ecutp,e3,e7,
     &                x(ndim-12:ndim-8),p(:,i1),p(:,i2),p(:,igind),dl1)
      else
       call kinffg(imom,sqs1,m(i1),m(i2),ecutp,e3,e7,
     &                x(ndim-12:ndim-8),p(:,i1),p(:,i2),p(:,igind),dl1)
      end if

      if (iqprec == 1) then
       call kinff_qp(imom,sqs2,m(i3),m(i4),x(1:2),p(:,i3),p(:,i4),dl2)
       call kinff_qp(imom,sqs3,m(i5),m(i6),x(3:4),p(:,i5),p(:,i6),dl3)
       call kinff_qp(imom,sqs4,m(i7),m(i8),x(5:6),p(:,i7),p(:,i8),dl4)
       call kinff_qp(imom,sqs234,sqs2,sqs34,x(7:8),q2,q34,dl234)
       call kinff_qp(imom,sqs34,sqs3,sqs4,x(9:10),q3,q4,dl34)
      else
       call kinff(imom,sqs2,m(i3),m(i4),x(1:2),p(:,i3),p(:,i4),dl2)
       call kinff(imom,sqs3,m(i5),m(i6),x(3:4),p(:,i5),p(:,i6),dl3)
       call kinff(imom,sqs4,m(i7),m(i8),x(5:6),p(:,i7),p(:,i8),dl4)
       call kinff(imom,sqs234,sqs2,sqs34,x(7:8),q2,q34,dl234)
       call kinff(imom,sqs34,sqs3,sqs4,x(9:10),q3,q4,dl34)
      end if

      dlips=dlips*ds234*ds34*ds1*ds2*ds3*ds4*dl234*dl34*dl1*dl2*dl3*dl4

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
        call boost_qp(g,v,p(:,igind))
       else
        call boost(g,v,p(:,i1))
        call boost(g,v,p(:,i2))
        call boost(g,v,p(:,igind))
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
       v=q234(1:)/q234(0)
       g=q234(0)/sqs234
       if (iqprec == 1) then
        call boost_qp(g,v,p(:,i3))
        call boost_qp(g,v,p(:,i4))
        call boost_qp(g,v,p(:,i5))
        call boost_qp(g,v,p(:,i6))
        call boost_qp(g,v,p(:,i7))
        call boost_qp(g,v,p(:,i8))
       else
        call boost(g,v,p(:,i3))
        call boost(g,v,p(:,i4))
        call boost(g,v,p(:,i5))
        call boost(g,v,p(:,i6))
        call boost(g,v,p(:,i7))
        call boost(g,v,p(:,i8))
       end if
       if ((p(0,i3)+p(0,i4)+p(0,i5)+p(0,i6)+p(0,i7)+p(0,i8))**2
     &    -(p(1,i3)+p(1,i4)+p(1,i5)+p(1,i6)+p(1,i7)+p(1,i8))**2
     &    -(p(2,i3)+p(2,i4)+p(2,i5)+p(2,i6)+p(2,i7)+p(2,i8))**2
     &    -(p(3,i3)+p(3,i4)+p(3,i5)+p(3,i6)+p(3,i7)+p(3,i8))**2 
     &                                             < sqsl234**2) then
        dlips=0.d0
        return
       end if
      end if

      end

