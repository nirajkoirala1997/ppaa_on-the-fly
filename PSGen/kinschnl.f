      subroutine kinschnl(ikin,imom,ecm,x,ndim,dlips)
 
      use fourmom
      use kinparams
 
      implicit none
 
      integer :: ikin,imom,ndim
      real(dp) :: x(ndim),dlips
      real(dp) :: ecm
 
      if (ikin == 1) then
       call kin2_2(3,4,
     &                             imom,ecm,p,x,ndim,dlips)
      end if
 
      end
