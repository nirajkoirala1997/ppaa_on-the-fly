      subroutine kin3_2(i3,i4,mv1,gv1,i5,
     &                                  imom,ecm,p,x,ndim,dlips)

      use kinparams

      implicit none
 
      integer :: i3,i4,i5,imom,ndim
      real(dp) :: mv1,gv1,x(ndim),dlips,ds1,dl1
      real(xp) :: ecm,p(0:3,nexp),q1(0:3),v(3),g
      real(xp) :: sqsl1,sqsu1,s1,sqs1

      sqsl1=m(i3)+m(i4)
      sqsu1=ecm-sqsl1

      if (imom /= 1) then
       q1=p(:,i3)+p(:,i4)
      end if

      if (iqprec == 1) then
       call mapc_qp(imom,x(ndim),q1,mv1,gv1,sqsc,sqsl1,sqsu1,s1,sqs1,
     &                                                           ds1)
      else
       call mapckk(imom,x(ndim),q1,mv1,gv1,sqsc,sqsl1,sqsu1,s1,sqs1,ds1)
      end if

      if (iqprec == 1) then
       call kinff_qp(imom,sqs1,m(i3),m(i4),x(1:2),p(:,i3),p(:,i4),dl1)
       call kinff_qp(imom,ecm,sqs1,m(i5),x(3:4),q1,p(:,i5),dlips)
      else
       call kinff(imom,sqs1,m(i3),m(i4),x(1:2),p(:,i3),p(:,i4),dl1)
       call kinff(imom,ecm,sqs1,m(i5),x(3:4),q1,p(:,i5),dlips)
      end if

      dlips=ds1*dlips*dl1

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
      end if

      end
