      subroutine psncls_qp(ikin0,ecm,p,x,ndim,dlinv)
 
      use kinparams
      use psnrmtest
 
      implicit none
 
      integer, parameter :: imom=0
      integer :: ikin0,ndim,ikin
      real(dp) :: x(ndim),dlinv,dlips,dlips0,rlips
      real(qp) :: ecm,p(0:3,nexp)
 
      dlips0=dlinv
      dlinv=0.d0
      dlinvt=0.d0
 
      do ikin=1,nkin
       if (aw(ikin) > 0.d0) then
        call kincls(ikin,imom,ecm,p,x,ndim,dlips)
        if (dlips > 0.d0) then
         dlinvt=dlinvt+aw(ikin)/dlips
         if (ikin == ikin0) then
          rlips=abs(dlips-dlips0)/dlips
          if (rlips > psnac) then
           iwkin(ikin)=iwkin(ikin)+1
           if (rlips > aw_w(ikin,1)) then
            aw_w(ikin,1)=rlips
            aw_w(ikin,2)=dlips0
            aw_w(ikin,3)=dlips
           end if
          end if
          dlips=dlips0
         end if
         dlinv=dlinv+aw(ikin)/dlips
        end if
       end if
      end do
 
      end
