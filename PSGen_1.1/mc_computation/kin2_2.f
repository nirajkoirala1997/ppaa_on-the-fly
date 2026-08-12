      subroutine kin2_2(i3,i4,imom,ecm,p,x,ndim,dlips)

      use kinparams

      implicit none
 
      integer :: i3,i4,imom,ndim
      real(dp) :: x(ndim),dlips
      real(xp) :: ecm,p(0:3,nexp)


      if (iqprec == 1) then
       call kinff_qp(imom,ecm,m(i3),m(i4),x,p(:,i3),p(:,i4),dlips)
      else
       call kinff(imom,ecm,m(i3),m(i4),x,p(:,i3),p(:,i4),dlips)
      end if

      end
