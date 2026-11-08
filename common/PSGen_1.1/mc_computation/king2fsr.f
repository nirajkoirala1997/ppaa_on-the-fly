      subroutine king2fsr(i1,i2,mv1,gv1,imom,ecm,x,ndim,dlips)
 
      use fourmom
      use kinparams
      use engymom
      use kincuts 

      implicit none

      integer :: i1,i2,imom,ndim
      real(dp) :: mv1,gv1,x(ndim),dlips
      real(xp) :: ecm

      if (iqprec == 1) then
       call kinffg_qp(imom,ecm,m(i1),m(i2),ecut0,p(0,i1),p(0,igind),
     &                x,p(:,i1),p(:,i2),p(:,igind),dlips)
      else
       call kinffg(imom,ecm,m(i1),m(i2),ecut0,p(0,i1),p(0,igind),
     &                x,p(:,i1),p(:,i2),p(:,igind),dlips)
      end if

      if (dlips <= 0.d0) then
       dlips=0.d0
       return
      end if

      end

