      subroutine psnkk(ikin0,ecm,p,x,ndim,dlinv)
 
      use partcls
      use psnrmtest
 
      implicit none
 
      integer :: ikin0,ndim
      real(8) :: x(ndim),dlinv,dlips0
      real(16) :: ecm,p(0:3,nepl)
 
      dlips0=dlinv
      dlinv=0.d0
      dlinvt=0.d0
 
      call psub1(ikin0,ecm,p,x,ndim,dlinv,dlips0)
 
      end
