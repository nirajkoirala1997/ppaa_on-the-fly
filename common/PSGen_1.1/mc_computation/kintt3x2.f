      subroutine kintt3x2(ikin,imom,ecm,p,x,ndim,dlips)
c----------------------------------------------------------------------c
c Returns cartezian components of the final state four momenta in the 
c c.m.s. and the phase space normalization of a 2 --> n process 
c with two t-channel poles ~1/[(p1-p3)^-mt1^2] and ~1/[(p2-p4)^-mt2^2].
c The factor (2*pi)^-(3n-4) is not included in dlips.
c 
c Input: imom=1/else => four momenta p3,p4 are/are not calculated,
c        ecm - c.m.s. energy,
c        x(i) - random variables, 0 < x(i) < 1, i=1,2,...,ndim
c        p(:,i), i=1,2 - initial state particle four momenta.
c Output: p(:,i), i=3,...,nepl - final state particle four momenta,
c         dlips - phase space normalization.
c
c Karol Kolodziej, University of Silesia, karol.kolodziej@us.edu.pl
c 11 July, 2018; modified on 
c----------------------------------------------------------------------c

      use inprms
      use ntnkin
      use kincuts
      use partcls
 
      implicit none
 
      integer :: ikin,imom,ndim
      real(dp) :: ecm,p(0:3,nepl),x(ndim),dlips,v(3),g
      real(dp), dimension(0:3) :: q
      real(dp) :: sqsl1,sqsl2
      real(dp) :: sqsu1,sqsu2
      real(dp) :: sqs1,sqs2
      real(dp) :: ds1,ds2
      real(dp) :: s1,s2
      real(dp) :: dl1,dl2
c      real(dp), dimension(nepl) :: n

      include 'tchcalls.f'

      contains

      subroutine kintt3y2(i1,i3,it1,i2,i4,it2,id)

      integer :: i1,i2,i3,i4,it1,it2,id,i
      real(dp) :: q1(0:3)

      if (id == 1) then
       call kintt3a(imom,ecm,ecuta,ctab,ctc,m(i1),m(i2),m(i3),m(i4),
     &              x,p(:,i1),p(:,i2),p(:,i3),p(:,i4),p(:,ix(1)),dlips)
       return
      end if

      sqsl1=0.d0
      do i=1,id
       sqsl1=sqsl1+m(ix(i))
      end do
      sqsl2=m(i3)+m(i4)
      if (imom /= 1) then
       q1=0.d0
       do i=1,id
        q1=q1+p(:,ix(i))
       end do
      end if
      sqsu1=ecm-sqsl2
      call mapakk(imom,x(1),q1,sqscut,sqsl1,sqsu1,s1,sqs1,ds1)
      call kintt3(imom,ecm,sqs1,
     &             it1,zero,ctc,m(i1),m(i3),it2,zero,ctc,m(i2),m(i4),
     &             x(2:6),p(:,i1),p(:,i2),p(:,i3),p(:,i4),q1,dlips)
      call kintt_rest(imom,sqs1,q1,p,x(7:),ndim-6,ix,id,dl1)
      dlips=dlips*ds1*dl1
      if (dlips <= 0.d0) then
       dlips=0.d0
       return
      end if

      end subroutine kintt3y2
 
      end
