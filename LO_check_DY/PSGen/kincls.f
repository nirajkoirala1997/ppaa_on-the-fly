      subroutine kincls(ikin,imom,x,ndim,dlips)
 
      use kinparams
      use fourmom
 
      implicit none
 
      integer :: ikin,imom,ndim
      real(8) :: ecm,x(ndim),dlips
 
      ecm=p(0,1)+p(0,2)
 
      if (ikin <= nkins) then
       call kinschnl(ikin,imom,ecm,x,ndim,dlips)
      else
! Add a call to your own made kinematical subroutine here and change 
! the value of parameter nkin in kinparams.f appropriately
      end if
 
      include 'kin_tests.f'
 
      end
