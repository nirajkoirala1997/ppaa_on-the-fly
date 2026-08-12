      subroutine kintt_rest(imom,sqs,q,p,x,dx,ix,dix,dl)
c----------------------------------------------------------------------c
c Returns cartezian components of the final state four momenta
c p(:,ix(1)),...,p(:,ix(dix)), in the c.m.s. 
c 
c Input: imom=1/else => the four momenta are/are not calculated,
c        sqs - invariant mass of the subsystem that do not form 
c              the t-channel poles,
c        q - four momentum of the particle ix(1),...,ix(dix)) system in
c            c.m.s.
c        x(i) - random variables, 0 < x(i) < 1, i=1,...,dx
c Output: p(:,i), i=ix(1),...,ix(dix) - final state particle four momenta
c         in c.m.s.,
c         dl - phase space normalization, corresponding powers of (2*pi) 
c              are not included
c
c Karol Kolodziej, University of Silesia, karol.kolodziej@us.edu.pl
c 11 July, 2018; modified on 
c----------------------------------------------------------------------c

      use kinparams
c      use mathprms
c      use partcls

      implicit none
      integer :: imom,dx,dix,ix(dix),i
      real(dp) :: sqs,q(0:3),p(0:3,nexp),x(dx),dl,v(3),g
      real(dp) :: sqs1,sqs1l,sqs1u,ds1,s1,q1(0:3),dl1,dl2
      real(dp) :: sqs2,sqs2l,sqs2u,ds2,s2,q2(0:3),dl3
      real(dp) :: sqs3,sqs3l,sqs3u,ds3,s3,q3(0:3),dl4
      real(dp) :: sqs4,sqs4l,sqs4u,ds4,s4,q4(0:3),dl5
      real(dp) :: sqs5,sqs5l,sqs5u,ds5,s5,q5(0:3),dl6
      real(dp) :: sqs6,sqs6l,sqs6u,ds6,s6,q6(0:3),dl7

      if (dix == 2) then
       call kinff(imom,sqs,m(ix(1)),m(ix(2)),x,
     &                                   p(:,ix(1)),p(:,ix(2)),dl) 
       if (dl <= 0.d0) then
        dl=0.d0
        return
       end if
       if (imom == 1) then
        v=q(1:)/q(0)
        g=q(0)/sqs                                                            
        do i=1,dix
         call boost(g,v,p(:,ix(i)))
        end do
       end if
      else if (dix == 3) then
       sqs1l=m(ix(1))
       sqs1u=sqs-(m(ix(2))+m(ix(3)))
       ds1=(sqs1u-sqs1l)*(sqs1u+sqs1l)
       if (ds1 <= 0.d0) then
        dl=0.d0
        return
       end if
       if (imom == 1) then
        s1=ds1*x(1)+sqs1l**2
       else
        q1=p(:,ix(2))+p(:,ix(3))
        s1=(q1(0)-q1(3))*(q1(0)+q1(3))-q1(1)*q1(1)-q1(2)*q1(2)
       end if
       sqs1=sqrt(s1)
       call kinff(imom,sqs,m(ix(1)),sqs1,x(2:3),
     &                                   p(:,ix(1)),q1,dl1) 
       call kinff(imom,sqs1,m(ix(2)),m(ix(3)),x(4:5),
     &                           p(:,ix(2)),p(:,ix(3)),dl2)
       dl=ds1*dl1*dl2
       if (dl <= 0.d0) then
        dl=0.d0
        return
       end if
       if (imom == 1) then
        v=q1(1:)/q1(0)
        g=q1(0)/sqs1
        call boost(g,v,p(:,ix(2)))
        call boost(g,v,p(:,ix(3)))
        v=q(1:)/q(0)
        g=q(0)/sqs                                                            
        do i=1,dix
         call boost(g,v,p(:,ix(i)))
        end do
       end if
      else if (dix == 4) then 
       sqs1l=m(ix(1))+m(ix(2))
       sqs2l=m(ix(3))+m(ix(4))
       sqs1u=sqs-sqs2l
       ds1=(sqs1u-sqs1l)*(sqs1u+sqs1l)
       if (imom == 1) then
        s1=ds1*x(1)+sqs1l**2
       else
        q1=p(:,ix(1))+p(:,ix(2))
        s1=(q1(0)-q1(3))*(q1(0)+q1(3))-q1(1)*q1(1)-q1(2)*q1(2)
       end if
       sqs1=sqrt(s1)
       sqs2u=sqs-sqs1
       ds2=(sqs2u-sqs2l)*(sqs2u+sqs2l)
       if (imom == 1) then
        s2=ds2*x(2)+sqs2l**2
       else
        q2=p(:,ix(3))+p(:,ix(4))
        s2=(q2(0)-q2(3))*(q2(0)+q2(3))-q2(1)*q2(1)-q2(2)*q2(2)
       end if
       sqs2=sqrt(s2)
       call kinff(imom,sqs,sqs1,sqs2,x(3:4),q1,q2,dl1) 
       call kinff(imom,sqs1,m(ix(1)),m(ix(2)),x(5:6),
     &                           p(:,ix(1)),p(:,ix(2)),dl2)
       call kinff(imom,sqs2,m(ix(3)),m(ix(4)),x(7:8),
     &                           p(:,ix(3)),p(:,ix(4)),dl3)
       dl=ds1*ds2*dl1*dl2*dl3
       if (dl <= 0.d0) then
        dl=0.d0
        return
       end if
       if (imom == 1) then
        v=q1(1:)/q1(0)
        g=q1(0)/sqs1
        call boost(g,v,p(:,ix(1)))
        call boost(g,v,p(:,ix(2)))
        v=q2(1:)/q2(0)
        g=q2(0)/sqs2
        call boost(g,v,p(:,ix(3)))
        call boost(g,v,p(:,ix(4)))
        v=q(1:)/q(0)
        g=q(0)/sqs                                                            
        do i=1,dix
         call boost(g,v,p(:,ix(i)))
        end do
       end if
      else if (dix == 5) then
       sqs1l=m(ix(1))+m(ix(2))
       sqs2l=m(ix(3))+m(ix(4))+m(ix(5))
       sqs1u=sqs-sqs2l
       ds1=(sqs1u-sqs1l)*(sqs1u+sqs1l)
       if (imom == 1) then
        s1=ds1*x(1)+sqs1l**2
       else
        q1=p(:,ix(1))+p(:,ix(2))
        s1=(q1(0)-q1(3))*(q1(0)+q1(3))-q1(1)*q1(1)-q1(2)*q1(2)
       end if
       sqs1=sqrt(s1)
       sqs2u=sqs-sqs1
       ds2=(sqs2u-sqs2l)*(sqs2u+sqs2l)
       if (imom == 1) then
        s2=ds2*x(2)+sqs2l**2
       else
        q3=p(:,ix(3))+p(:,ix(4))
        q2=q3+p(:,ix(5))
        s2=(q2(0)-q2(3))*(q2(0)+q2(3))-q2(1)*q2(1)-q2(2)*q2(2)
       end if
       sqs2=sqrt(s2)
       sqs3l=m(ix(3))+m(ix(4))
       sqs3u=sqs2-m(ix(5))
       ds3=(sqs3u-sqs3l)*(sqs3u+sqs3l)
       if (imom == 1) then
        s3=ds3*x(3)+sqs3l**2
       else
        s3=(q3(0)-q3(3))*(q3(0)+q3(3))-q3(1)*q3(1)-q3(2)*q3(2)
       end if
       sqs3=sqrt(s3)
       call kinff(imom,sqs,sqs1,sqs2,x(4:5),q1,q2,dl1) 
       call kinff(imom,sqs1,m(ix(1)),m(ix(2)),x(6:7),
     &                           p(:,ix(1)),p(:,ix(2)),dl2)
       call kinff(imom,sqs2,sqs3,m(ix(5)),x(8:9),
     &                                   q3,p(:,ix(5)),dl3)
       call kinff(imom,sqs3,m(ix(3)),m(ix(4)),x(10:11),
     &                           p(:,ix(3)),p(:,ix(4)),dl4)
       dl=ds1*ds2*ds3*dl1*dl2*dl3*dl4
       if (dl <= 0.d0) then
        dl=0.d0
        return
       end if
       if (imom == 1) then
        v=q1(1:)/q1(0)
        g=q1(0)/sqs1
        call boost(g,v,p(:,ix(1)))
        call boost(g,v,p(:,ix(2)))
        v=q3(1:)/q3(0)
        g=q3(0)/sqs3
        call boost(g,v,p(:,ix(3)))
        call boost(g,v,p(:,ix(4)))
        v=q2(1:)/q2(0)
        g=q2(0)/sqs2
        call boost(g,v,p(:,ix(3)))
        call boost(g,v,p(:,ix(4)))
        call boost(g,v,p(:,ix(5)))
        v=q(1:)/q(0)
        g=q(0)/sqs                                                            
        do i=1,dix
         call boost(g,v,p(:,ix(i)))
        end do
       end if
      else if (dix == 6) then
       sqs1l=m(ix(1))+m(ix(2))+m(ix(3))
       sqs2l=m(ix(4))+m(ix(5))+m(ix(6))
       sqs1u=sqs-sqs2l
       ds1=(sqs1u-sqs1l)*(sqs1u+sqs1l)
       if (imom == 1) then
        s1=ds1*x(1)+sqs1l**2
       else
        q3=p(:,ix(1))+p(:,ix(2))
        q1=q3+p(:,ix(3))
        s1=(q1(0)-q1(3))*(q1(0)+q1(3))-q1(1)*q1(1)-q1(2)*q1(2)
       end if
       sqs1=sqrt(s1)
       sqs2u=sqs-sqs1
       ds2=(sqs2u-sqs2l)*(sqs2u+sqs2l)
       if (imom == 1) then
        s2=ds2*x(2)+sqs2l**2
       else
        q4=p(:,ix(4))+p(:,ix(5))
        q2=q4+p(:,ix(6))
        s2=(q2(0)-q2(3))*(q2(0)+q2(3))-q2(1)*q2(1)-q2(2)*q2(2)
       end if        
       sqs2=sqrt(s2)
       sqs3l=m(ix(1))+m(ix(2))
       sqs3u=sqs1-m(ix(3))
       ds3=(sqs3u-sqs3l)*(sqs3u+sqs3l)
       if (imom == 1) then
        s3=ds3*x(3)+sqs3l**2
       else
        s3=(q3(0)-q3(3))*(q3(0)+q3(3))-q3(1)*q3(1)-q3(2)*q3(2)
       end if
       sqs3=sqrt(s3)
       sqs4l=m(ix(4))+m(ix(5))
       sqs4u=sqs2-m(ix(6))
       ds4=(sqs4u-sqs4l)*(sqs4u+sqs4l)
       if (imom == 1) then
        s4=ds4*x(4)+sqs4l**2
       else
        s4=(q4(0)-q4(3))*(q4(0)+q4(3))-q4(1)*q4(1)-q4(2)*q4(2)
       end if
       sqs4=sqrt(s4)
       call kinff(imom,sqs,sqs1,sqs2,x(5:6),q1,q2,dl1) 
       call kinff(imom,sqs1,sqs3,m(ix(3)),x(7:8),
     &                                   q3,p(:,ix(3)),dl2)
       call kinff(imom,sqs3,m(ix(1)),m(ix(2)),x(9:10),
     &                           p(:,ix(1)),p(:,ix(2)),dl3)
       call kinff(imom,sqs2,sqs4,m(ix(6)),x(11:12),
     &                                   q4,p(:,ix(6)),dl4)
       call kinff(imom,sqs4,m(ix(4)),m(ix(5)),x(13:14),
     &                           p(:,ix(4)),p(:,ix(5)),dl5)
       dl=ds1*ds2*ds3*ds4*dl1*dl2*dl3*dl4*dl5
       if (dl <= 0.d0) then
        dl=0.d0
        return
       end if
       if (imom == 1) then
        v=q3(1:)/q3(0)
        g=q3(0)/sqs3
        call boost(g,v,p(:,ix(1)))
        call boost(g,v,p(:,ix(2)))
        v=q1(1:)/q1(0)
        g=q1(0)/sqs1
        call boost(g,v,p(:,ix(1)))
        call boost(g,v,p(:,ix(2)))
        call boost(g,v,p(:,ix(3)))
        v=q4(1:)/q4(0)
        g=q4(0)/sqs4
        call boost(g,v,p(:,ix(4)))
        call boost(g,v,p(:,ix(5)))
        v=q2(1:)/q2(0)
        g=q2(0)/sqs2
        call boost(g,v,p(:,ix(4)))
        call boost(g,v,p(:,ix(5)))
        call boost(g,v,p(:,ix(6)))
        v=q(1:)/q(0)
        g=q(0)/sqs                                                            
        do i=1,dix
         call boost(g,v,p(:,ix(i)))
        end do
       end if
      else if (dix == 7) then
       sqs1l=m(ix(1))+m(ix(2))+m(ix(3))
       sqs2l=m(ix(4))+m(ix(5))+m(ix(6))+m(ix(7))
       sqs1u=sqs-sqs2l
       ds1=(sqs1u-sqs1l)*(sqs1u+sqs1l)
       if (imom == 1) then
        s1=ds1*x(1)+sqs1l**2
       else
        q3=p(:,ix(1))+p(:,ix(2))
        q1=q3+p(:,ix(3))
        s1=(q1(0)-q1(3))*(q1(0)+q1(3))-q1(1)*q1(1)-q1(2)*q1(2)
       end if
       sqs1=sqrt(s1)
       sqs2u=sqs-sqs1
       ds2=(sqs2u-sqs2l)*(sqs2u+sqs2l)
       if (imom == 1) then
        s2=ds2*x(2)+sqs2l**2
       else
        q4=p(:,ix(4))+p(:,ix(5))
        q5=p(:,ix(6))+p(:,ix(7))
        q2=q4+q5
        s2=(q2(0)-q2(3))*(q2(0)+q2(3))-q2(1)*q2(1)-q2(2)*q2(2)
       end if
       sqs2=sqrt(s2)
       sqs3l=m(ix(1))+m(ix(2))
       sqs3u=sqs1-m(ix(3))
       ds3=(sqs3u-sqs3l)*(sqs3u+sqs3l)
       if (imom == 1) then
        s3=ds3*x(3)+sqs3l**2
       else
        s3=(q3(0)-q3(3))*(q3(0)+q3(3))-q3(1)*q3(1)-q3(2)*q3(2)
       end if
       sqs3=sqrt(s3)
       sqs4l=m(ix(4))+m(ix(5))
       sqs5l=m(ix(6))+m(ix(7))
       sqs4u=sqs2-sqs5l
       ds4=(sqs4u-sqs4l)*(sqs4u+sqs4l)
       if (imom == 1) then
        s4=ds4*x(4)+sqs4l**2
       else
        s4=(q4(0)-q4(3))*(q4(0)+q4(3))-q4(1)*q4(1)-q4(2)*q4(2)
       end if
       sqs4=sqrt(s4)
       sqs5u=sqs2-sqs4
       ds5=(sqs5u-sqs5l)*(sqs5u+sqs5l)
       if (imom == 1) then
        s5=ds5*x(5)+sqs5l**2
       else
        s5=(q5(0)-q5(3))*(q5(0)+q5(3))-q5(1)*q5(1)-q5(2)*q5(2)
       end if
       sqs5=sqrt(s5)
       call kinff(imom,sqs,sqs1,sqs2,x(6:7),q1,q2,dl1) 
       call kinff(imom,sqs1,sqs3,m(ix(3)),x(8:9),
     &                                   q3,p(:,ix(3)),dl2)
       call kinff(imom,sqs3,m(ix(1)),m(ix(2)),x(10:11),
     &                           p(:,ix(1)),p(:,ix(2)),dl3)
       call kinff(imom,sqs2,sqs4,sqs5,x(12:13),
     &                                   q4,q5,dl4)
       call kinff(imom,sqs4,m(ix(4)),m(ix(5)),x(14:15),
     &                           p(:,ix(4)),p(:,ix(5)),dl5)
       call kinff(imom,sqs5,m(ix(6)),m(ix(7)),x(16:17),
     &                           p(:,ix(6)),p(:,ix(7)),dl6)
       dl=ds1*ds2*ds3*ds4*ds5*dl1*dl2*dl3*dl4*dl5*dl6
       if (dl <= 0.d0) then
        dl=0.d0
        return
       end if
       if (imom == 1) then
        v=q3(1:)/q3(0)
        g=q3(0)/sqs3
        call boost(g,v,p(:,ix(1)))
        call boost(g,v,p(:,ix(2)))
        v=q1(1:)/q1(0)
        g=q1(0)/sqs1
        call boost(g,v,p(:,ix(1)))
        call boost(g,v,p(:,ix(2)))
        call boost(g,v,p(:,ix(3)))
        v=q4(1:)/q4(0)
        g=q4(0)/sqs4
        call boost(g,v,p(:,ix(4)))
        call boost(g,v,p(:,ix(5)))
        v=q5(1:)/q5(0)
        g=q5(0)/sqs5
        call boost(g,v,p(:,ix(6)))
        call boost(g,v,p(:,ix(7)))
        v=q2(1:)/q2(0)
        g=q2(0)/sqs2
        call boost(g,v,p(:,ix(4)))
        call boost(g,v,p(:,ix(5)))
        call boost(g,v,p(:,ix(6)))
        call boost(g,v,p(:,ix(7)))
        v=q(1:)/q(0)
        g=q(0)/sqs                                                            
        do i=1,dix
         call boost(g,v,p(:,ix(i)))
        end do
       end if
      else if (dix == 8) then
       sqs1l=m(ix(1))+m(ix(2))+m(ix(3))+m(ix(4))
       sqs2l=m(ix(5))+m(ix(6))+m(ix(7))+m(ix(8))
       sqs1u=sqs-sqs2l
       ds1=(sqs1u-sqs1l)*(sqs1u+sqs1l)
       if (imom == 1) then
        s1=ds1*x(1)+sqs1l**2
       else
        q3=p(:,ix(1))+p(:,ix(2))
        q4=p(:,ix(3))+p(:,ix(4))
        q1=q3+q4
        s1=(q1(0)-q1(3))*(q1(0)+q1(3))-q1(1)*q1(1)-q1(2)*q1(2)
       end if
       sqs1=sqrt(s1)
       sqs2u=sqs-sqs1
       ds2=(sqs2u-sqs2l)*(sqs2u+sqs2l)
       if (imom == 1) then
        s2=ds2*x(2)+sqs2l**2
       else
        q5=p(:,ix(5))+p(:,ix(6))
        q6=p(:,ix(7))+p(:,ix(8))
        q2=q5+q6
        s2=(q2(0)-q2(3))*(q2(0)+q2(3))-q2(1)*q2(1)-q2(2)*q2(2)
       end if
       sqs2=sqrt(s2)
       sqs3l=m(ix(1))+m(ix(2))
       sqs4l=m(ix(3))+m(ix(4))
       sqs3u=sqs1-sqs4l
       ds3=(sqs3u-sqs3l)*(sqs3u+sqs3l)
       if (imom == 1) then
        s3=ds3*x(3)+sqs3l**2
       else
        s3=(q3(0)-q3(3))*(q3(0)+q3(3))-q3(1)*q3(1)-q3(2)*q3(2)
       end if
       sqs3=sqrt(s3)
       sqs4u=sqs1-sqs3
       ds4=(sqs4u-sqs4l)*(sqs4u+sqs4l)
       if (imom == 1) then
        s4=ds4*x(4)+sqs4l**2
       else
        s4=(q4(0)-q4(3))*(q4(0)+q4(3))-q4(1)*q4(1)-q4(2)*q4(2)
       end if
       sqs4=sqrt(s4)
       sqs5l=m(ix(5))+m(ix(6))
       sqs6l=m(ix(7))+m(ix(8))
       sqs5u=sqs2-sqs6l
       ds5=(sqs5u-sqs5l)*(sqs5u+sqs5l)
       if (imom == 1) then
        s5=ds5*x(5)+sqs5l**2
       else
        s5=(q5(0)-q5(3))*(q5(0)+q5(3))-q5(1)*q5(1)-q5(2)*q5(2)
       end if
       sqs5=sqrt(s5)
       sqs6u=sqs2-sqs5
       ds6=(sqs6u-sqs6l)*(sqs6u+sqs6l)
       if (imom == 1) then
        s6=ds6*x(6)+sqs6l**2
       else
        s6=(q6(0)-q6(3))*(q6(0)+q6(3))-q6(1)*q6(1)-q6(2)*q6(2)
       end if
       sqs6=sqrt(s6)
       call kinff(imom,sqs,sqs1,sqs2,x(7:8),q1,q2,dl1) 
       call kinff(imom,sqs1,sqs3,sqs4,x(9:10),q3,q4,dl2)
       call kinff(imom,sqs3,m(ix(1)),m(ix(2)),x(11:12),
     &                           p(:,ix(1)),p(:,ix(2)),dl3)
       call kinff(imom,sqs4,m(ix(3)),m(ix(4)),x(13:14),
     &                           p(:,ix(3)),p(:,ix(4)),dl4)
       call kinff(imom,sqs2,sqs5,sqs6,x(15:16),q5,q6,dl5)
       call kinff(imom,sqs5,m(ix(5)),m(ix(6)),x(17:18),
     &                           p(:,ix(5)),p(:,ix(6)),dl6)
       call kinff(imom,sqs6,m(ix(7)),m(ix(8)),x(19:20),
     &                           p(:,ix(7)),p(:,ix(8)),dl7)
       dl=ds1*ds2*ds3*ds4*ds5*ds6*dl1*dl2*dl3*dl4*dl5*dl6*dl7
       if (dl <= 0.d0) then
        dl=0.d0
        return
       end if
       if (imom == 1) then
        v=q3(1:)/q3(0)
        g=q3(0)/sqs3
        call boost(g,v,p(:,ix(1)))
        call boost(g,v,p(:,ix(2)))
        v=q4(1:)/q4(0)
        g=q4(0)/sqs4
        call boost(g,v,p(:,ix(3)))
        call boost(g,v,p(:,ix(4)))
        v=q1(1:)/q1(0)
        g=q1(0)/sqs1
        call boost(g,v,p(:,ix(1)))
        call boost(g,v,p(:,ix(2)))
        call boost(g,v,p(:,ix(3)))
        call boost(g,v,p(:,ix(4)))
        v=q5(1:)/q5(0)
        g=q5(0)/sqs5
        call boost(g,v,p(:,ix(5)))
        call boost(g,v,p(:,ix(6)))
        v=q6(1:)/q6(0)
        g=q6(0)/sqs6
        call boost(g,v,p(:,ix(7)))
        call boost(g,v,p(:,ix(8)))
        v=q2(1:)/q2(0)
        g=q2(0)/sqs2
        call boost(g,v,p(:,ix(5)))
        call boost(g,v,p(:,ix(6)))
        call boost(g,v,p(:,ix(7)))
        call boost(g,v,p(:,ix(8)))
        v=q(1:)/q(0)
        g=q(0)/sqs                                                            
        do i=1,dix
         call boost(g,v,p(:,ix(i)))
        end do
       end if
      end if

      end
