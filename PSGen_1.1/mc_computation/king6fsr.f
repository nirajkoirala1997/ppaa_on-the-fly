      subroutine king6fsr(i1,i2,mv1,gv1,i3,i4,mv2,gv2,i5,i6,mv3,gv3,
     &                                         imom,ecm,x,ndim,dlips)
 
      use fourmom
      use kinparams
      use engymom
      use kincuts 

      implicit none

      integer :: i1,i2,i3,i4,i5,i6,imom,ndim
      real(dp) :: mv1,gv1,mv2,gv2,mv3,gv3,x(ndim),dlips
      real(xp) :: ecm
      real(dp) :: ds1,ds2,ds3,ds23
      real(dp) :: dl1,dl2,dl3,dl23

      real(dp) :: ecutp

      sqsl1=m(i1)+m(i2)
      sqsl2=m(i3)+m(i4)
      sqsl3=m(i5)+m(i6)
      sqsl23=sqsl2+sqsl3

      if (imom /= 1) then
       q1=p(:,i1)+p(:,i2)+p(:,igind)
       q2=p(:,i3)+p(:,i4)
       q3=p(:,i5)+p(:,i6)
       q23=q2+q3
      end if

      sqsu1=ecm-sqsl23
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

      sqsu23=ecm-sqs1
      ds23=(sqsu23-sqsl23)*(sqsu23+sqsl23)
      if (imom /= 1) then
       s23=q23(0)**2-q23(1)**2-q23(2)**2-q23(3)**2
      else
       s23=ds23*x(ndim-1)+sqsl23**2
      end if
      sqs23=sqrt(s23)

      sqsu2=sqs23-sqsl3
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

      sqsu3=sqs23-sqs2
      if (mv3 == 0.d0) then
       if (iqprec == 1) then
        call mapa_qp(imom,x(ndim-3),q3,sqscut,sqsl3,sqsu3,s3,sqs3,ds3)
       else
        call mapakk(imom,x(ndim-3),q3,sqscut,sqsl3,sqsu3,s3,sqs3,ds3)
       end if
      else
       if (iqprec == 1) then
        call mapb_qp(imom,x(ndim-3),q3,mv3,gv3,sqsl3,sqsu3,s3,sqs3,ds3)
       else
        call mapbkk(imom,x(ndim-3),q3,mv3,gv3,sqsl3,sqsu3,s3,sqs3,ds3)
       end if
      end if

      if (iqprec == 1) then
       call kinff_qp(imom,ecm,sqs1,sqs23,x(ndim-5:ndim-4),q1,q23,dlips)
      else
       call kinff(imom,ecm,sqs1,sqs23,x(ndim-5:ndim-4),q1,q23,dlips)
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
     &                x(ndim-10:ndim-6),p(:,i1),p(:,i2),p(:,igind),dl1)
      else
       call kinffg(imom,sqs1,m(i1),m(i2),ecutp,e3,e7,
     &                x(ndim-10:ndim-6),p(:,i1),p(:,i2),p(:,igind),dl1)
      end if

      if (iqprec == 1) then
       call kinff_qp(imom,sqs2,m(i3),m(i4),x(1:2),p(:,i3),p(:,i4),dl2)
       call kinff_qp(imom,sqs3,m(i5),m(i6),x(3:4),p(:,i5),p(:,i6),dl3)
       call kinff_qp(imom,sqs23,sqs2,sqs3,x(5:6),q2,q3,dl23)
      else
       call kinff(imom,sqs2,m(i3),m(i4),x(1:2),p(:,i3),p(:,i4),dl2)
       call kinff(imom,sqs3,m(i5),m(i6),x(3:4),p(:,i5),p(:,i6),dl3)
       call kinff(imom,sqs23,sqs2,sqs3,x(5:6),q2,q3,dl23)
      end if

      dlips=dlips*ds23*ds1*ds2*ds3*dl23*dl1*dl2*dl3

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
       v=q23(1:)/q23(0)
       g=q23(0)/sqs23
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
       if ((p(0,i3)+p(0,i4)+p(0,i5)+p(0,i6))**2
     &    -(p(1,i3)+p(1,i4)+p(1,i5)+p(1,i6))**2
     &    -(p(2,i3)+p(2,i4)+p(2,i5)+p(2,i6))**2
     &    -(p(3,i3)+p(3,i4)+p(3,i5)+p(3,i6))**2 < sqsl23**2) then
        dlips=0.d0
        return
       end if
      end if

      end

