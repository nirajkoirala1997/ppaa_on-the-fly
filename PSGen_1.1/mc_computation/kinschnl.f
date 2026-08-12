      subroutine kinschnl(ikin,imom,ecm,x,ndim,dlips)
 
      use fourmom
      use kinparams
 
      implicit none
 
      integer :: ikin,imom,ndim
      real(dp) :: x(ndim),dlips
      real(dp) :: ecm
 
      if (ikin == 1) then
       call kin7_33(6,7,8,mw,gamw,mt,gamt,3,4,5,mw,gamw,mt,gamt,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 2) then
       call kin7_33(6,7,8,mw,gamw,mt,gamt,3,4,5,mw,zero,mt,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 3) then
       call kin7_33(3,4,5,mw,gamw,mt,gamt,6,7,8,mw,zero,mt,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 4) then
       call kin7_33(6,7,8,mw,zero,mt,zero,3,4,5,mw,zero,mt,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 5) then
       call kin7_222(3,6,mg,zero,4,7,ma,zero,5,8,mz,gamz,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 6) then
       call kin7_222(3,6,mg,zero,4,7,ma,zero,5,8,mz,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 7) then
       call kin7_222(3,6,mg,zero,4,7,mz,gamz,5,8,mz,gamz,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 8) then
       call kin7_222(3,6,mg,zero,4,7,mz,gamz,5,8,mz,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 9) then
       call kin7_222(3,6,mg,zero,5,8,mz,gamz,4,7,mz,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 10) then
       call kin7_222(3,6,mg,zero,4,7,mz,zero,5,8,mz,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 11) then
       call kin7_222(3,6,mg,zero,7,8,mw,gamw,4,5,mw,gamw,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 12) then
       call kin7_222(3,6,mg,zero,7,8,mw,gamw,4,5,mw,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 13) then
       call kin7_222(3,6,mg,zero,4,5,mw,gamw,7,8,mw,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 14) then
       call kin7_222(3,6,mg,zero,7,8,mw,zero,4,5,mw,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 15) then
       call kin7_222(3,6,ma,zero,4,7,ma,zero,5,8,mz,gamz,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 16) then
       call kin7_222(3,6,ma,zero,4,7,ma,zero,5,8,mz,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 17) then
       call kin7_222(3,6,ma,zero,4,7,mz,gamz,5,8,mz,gamz,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 18) then
       call kin7_222(3,6,ma,zero,4,7,mz,gamz,5,8,mz,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 19) then
       call kin7_222(3,6,ma,zero,5,8,mz,gamz,4,7,mz,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 20) then
       call kin7_222(3,6,ma,zero,4,7,mz,zero,5,8,mz,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 21) then
       call kin7_222(3,6,ma,zero,7,8,mw,gamw,4,5,mw,gamw,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 22) then
       call kin7_222(3,6,ma,zero,7,8,mw,gamw,4,5,mw,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 23) then
       call kin7_222(3,6,ma,zero,4,5,mw,gamw,7,8,mw,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 24) then
       call kin7_222(3,6,ma,zero,7,8,mw,zero,4,5,mw,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 25) then
       call kin7_222(4,7,ma,zero,3,6,mz,gamz,5,8,mz,gamz,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 26) then
       call kin7_222(4,7,ma,zero,3,6,mz,gamz,5,8,mz,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 27) then
       call kin7_222(4,7,ma,zero,3,6,mz,zero,5,8,mz,gamz,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 28) then
       call kin7_222(4,7,ma,zero,3,6,mz,zero,5,8,mz,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 29) then
       call kin7_222(4,7,ma,zero,5,8,mz,gamz,3,6,mh,gamh,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 30) then
       call kin7_222(4,7,ma,zero,5,8,mz,gamz,3,6,mh,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 31) then
       call kin7_222(4,7,ma,zero,5,8,mz,zero,3,6,mh,gamh,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 32) then
       call kin7_222(4,7,ma,zero,5,8,mz,zero,3,6,mh,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 33) then
       call kin7_222(3,6,mz,gamz,4,7,mz,gamz,5,8,mz,gamz,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 34) then
       call kin7_222(3,6,mz,gamz,4,7,mz,gamz,5,8,mz,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 35) then
       call kin7_222(3,6,mz,gamz,5,8,mz,gamz,4,7,mz,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 36) then
       call kin7_222(3,6,mz,gamz,4,7,mz,zero,5,8,mz,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 37) then
       call kin7_222(3,6,mz,gamz,7,8,mw,gamw,4,5,mw,gamw,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 38) then
       call kin7_222(3,6,mz,gamz,7,8,mw,gamw,4,5,mw,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 39) then
       call kin7_222(3,6,mz,gamz,4,5,mw,gamw,7,8,mw,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 40) then
       call kin7_222(3,6,mz,gamz,7,8,mw,zero,4,5,mw,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 41) then
       call kin7_222(3,6,mz,zero,4,7,mz,gamz,5,8,mz,gamz,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 42) then
       call kin7_222(3,6,mz,zero,4,7,mz,gamz,5,8,mz,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 43) then
       call kin7_222(3,6,mz,zero,5,8,mz,gamz,4,7,mz,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 44) then
       call kin7_222(3,6,mz,zero,4,7,mz,zero,5,8,mz,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 45) then
       call kin7_222(3,6,mz,zero,7,8,mw,gamw,4,5,mw,gamw,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 46) then
       call kin7_222(3,6,mz,zero,7,8,mw,gamw,4,5,mw,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 47) then
       call kin7_222(3,6,mz,zero,4,5,mw,gamw,7,8,mw,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 48) then
       call kin7_222(3,6,mz,zero,7,8,mw,zero,4,5,mw,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 49) then
       call kin7_222(4,7,mz,gamz,5,8,mz,gamz,3,6,mh,gamh,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 50) then
       call kin7_222(4,7,mz,gamz,5,8,mz,gamz,3,6,mh,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 51) then
       call kin7_222(4,7,mz,gamz,5,8,mz,zero,3,6,mh,gamh,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 52) then
       call kin7_222(4,7,mz,gamz,5,8,mz,zero,3,6,mh,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 53) then
       call kin7_222(5,8,mz,gamz,4,7,mz,zero,3,6,mh,gamh,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 54) then
       call kin7_222(5,8,mz,gamz,4,7,mz,zero,3,6,mh,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 55) then
       call kin7_222(4,7,mz,zero,5,8,mz,zero,3,6,mh,gamh,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 56) then
       call kin7_222(4,7,mz,zero,5,8,mz,zero,3,6,mh,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 57) then
       call kin7_222(7,8,mw,gamw,4,5,mw,gamw,3,6,mh,gamh,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 58) then
       call kin7_222(7,8,mw,gamw,4,5,mw,gamw,3,6,mh,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 59) then
       call kin7_222(7,8,mw,gamw,4,5,mw,zero,3,6,mh,gamh,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 60) then
       call kin7_222(7,8,mw,gamw,4,5,mw,zero,3,6,mh,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 61) then
       call kin7_222(4,5,mw,gamw,7,8,mw,zero,3,6,mh,gamh,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 62) then
       call kin7_222(4,5,mw,gamw,7,8,mw,zero,3,6,mh,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 63) then
       call kin7_222(7,8,mw,zero,4,5,mw,zero,3,6,mh,gamh,9,
     &                             imom,ecm,p,x,ndim,dlips)
      else if (ikin == 64) then
       call kin7_222(7,8,mw,zero,4,5,mw,zero,3,6,mh,zero,9,
     &                             imom,ecm,p,x,ndim,dlips)
      end if
 
      end
