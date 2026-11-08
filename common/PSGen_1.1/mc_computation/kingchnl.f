      subroutine kingchnl(ikin,imom,ecm,x,ndim,dlips)
 
      use fourmom
      use kinparams
 
      implicit none
 
      integer :: ikin,imom,ndim
      real(dp) :: x(ndim),dlips
      real(dp) :: ecm
 
 
      if (ikin == nkinst+1) then
       call king6_33(6,7,8,mw,gamw,mt,gamt,3,4,5,mw,gamw,mt,gamt,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+2) then
       call king6_33(6,7,8,mw,gamw,mt,gamt,3,4,5,mw,zero,mt,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+3) then
       call king6_33(3,4,5,mw,gamw,mt,gamt,6,7,8,mw,zero,mt,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+4) then
       call king6_33(6,7,8,mw,zero,mt,zero,3,4,5,mw,zero,mt,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+5) then
       call king6isr(3,6,mg,zero,4,7,ma,zero,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+6) then
       call king6fsr(3,6,mg,zero,4,7,ma,zero,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+7) then
       call king6fsr(6,3,mg,zero,4,7,ma,zero,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+8) then
       call king6isr(3,6,mg,zero,4,7,ma,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+9) then
       call king6fsr(3,6,mg,zero,4,7,ma,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+10) then
       call king6fsr(6,3,mg,zero,4,7,ma,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+11) then
       call king6isr(3,6,mg,zero,4,7,mz,gamz,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+12) then
       call king6fsr(3,6,mg,zero,4,7,mz,gamz,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+13) then
       call king6fsr(6,3,mg,zero,4,7,mz,gamz,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+14) then
       call king6isr(3,6,mg,zero,4,7,mz,gamz,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+15) then
       call king6fsr(3,6,mg,zero,4,7,mz,gamz,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+16) then
       call king6fsr(6,3,mg,zero,4,7,mz,gamz,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+17) then
       call king6isr(3,6,mg,zero,5,8,mz,gamz,4,7,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+18) then
       call king6fsr(3,6,mg,zero,5,8,mz,gamz,4,7,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+19) then
       call king6fsr(6,3,mg,zero,5,8,mz,gamz,4,7,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+20) then
       call king6isr(3,6,mg,zero,4,7,mz,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+21) then
       call king6fsr(3,6,mg,zero,4,7,mz,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+22) then
       call king6fsr(6,3,mg,zero,4,7,mz,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+23) then
       call king6isr(3,6,mg,zero,7,8,mw,gamw,4,5,mw,gamw,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+24) then
       call king6fsr(3,6,mg,zero,7,8,mw,gamw,4,5,mw,gamw,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+25) then
       call king6fsr(6,3,mg,zero,7,8,mw,gamw,4,5,mw,gamw,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+26) then
       call king6isr(3,6,mg,zero,7,8,mw,gamw,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+27) then
       call king6fsr(3,6,mg,zero,7,8,mw,gamw,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+28) then
       call king6fsr(6,3,mg,zero,7,8,mw,gamw,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+29) then
       call king6isr(3,6,mg,zero,4,5,mw,gamw,7,8,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+30) then
       call king6fsr(3,6,mg,zero,4,5,mw,gamw,7,8,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+31) then
       call king6fsr(6,3,mg,zero,4,5,mw,gamw,7,8,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+32) then
       call king6isr(3,6,mg,zero,7,8,mw,zero,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+33) then
       call king6fsr(3,6,mg,zero,7,8,mw,zero,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+34) then
       call king6fsr(6,3,mg,zero,7,8,mw,zero,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+35) then
       call king6isr(3,6,ma,zero,4,7,ma,zero,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+36) then
       call king6fsr(3,6,ma,zero,4,7,ma,zero,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+37) then
       call king6fsr(6,3,ma,zero,4,7,ma,zero,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+38) then
       call king6isr(3,6,ma,zero,4,7,ma,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+39) then
       call king6fsr(3,6,ma,zero,4,7,ma,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+40) then
       call king6fsr(6,3,ma,zero,4,7,ma,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+41) then
       call king6isr(3,6,ma,zero,4,7,mz,gamz,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+42) then
       call king6fsr(3,6,ma,zero,4,7,mz,gamz,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+43) then
       call king6fsr(6,3,ma,zero,4,7,mz,gamz,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+44) then
       call king6isr(3,6,ma,zero,4,7,mz,gamz,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+45) then
       call king6fsr(3,6,ma,zero,4,7,mz,gamz,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+46) then
       call king6fsr(6,3,ma,zero,4,7,mz,gamz,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+47) then
       call king6isr(3,6,ma,zero,5,8,mz,gamz,4,7,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+48) then
       call king6fsr(3,6,ma,zero,5,8,mz,gamz,4,7,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+49) then
       call king6fsr(6,3,ma,zero,5,8,mz,gamz,4,7,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+50) then
       call king6isr(3,6,ma,zero,4,7,mz,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+51) then
       call king6fsr(3,6,ma,zero,4,7,mz,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+52) then
       call king6fsr(6,3,ma,zero,4,7,mz,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+53) then
       call king6isr(3,6,ma,zero,7,8,mw,gamw,4,5,mw,gamw,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+54) then
       call king6fsr(3,6,ma,zero,7,8,mw,gamw,4,5,mw,gamw,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+55) then
       call king6fsr(6,3,ma,zero,7,8,mw,gamw,4,5,mw,gamw,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+56) then
       call king6isr(3,6,ma,zero,7,8,mw,gamw,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+57) then
       call king6fsr(3,6,ma,zero,7,8,mw,gamw,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+58) then
       call king6fsr(6,3,ma,zero,7,8,mw,gamw,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+59) then
       call king6isr(3,6,ma,zero,4,5,mw,gamw,7,8,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+60) then
       call king6fsr(3,6,ma,zero,4,5,mw,gamw,7,8,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+61) then
       call king6fsr(6,3,ma,zero,4,5,mw,gamw,7,8,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+62) then
       call king6isr(3,6,ma,zero,7,8,mw,zero,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+63) then
       call king6fsr(3,6,ma,zero,7,8,mw,zero,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+64) then
       call king6fsr(6,3,ma,zero,7,8,mw,zero,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+65) then
       call king6isr(4,7,ma,zero,3,6,mz,gamz,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+66) then
       call king6fsr(4,7,ma,zero,3,6,mz,gamz,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+67) then
       call king6fsr(7,4,ma,zero,3,6,mz,gamz,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+68) then
       call king6isr(4,7,ma,zero,3,6,mz,gamz,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+69) then
       call king6fsr(4,7,ma,zero,3,6,mz,gamz,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+70) then
       call king6fsr(7,4,ma,zero,3,6,mz,gamz,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+71) then
       call king6isr(4,7,ma,zero,3,6,mz,zero,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+72) then
       call king6fsr(4,7,ma,zero,3,6,mz,zero,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+73) then
       call king6fsr(7,4,ma,zero,3,6,mz,zero,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+74) then
       call king6isr(4,7,ma,zero,3,6,mz,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+75) then
       call king6fsr(4,7,ma,zero,3,6,mz,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+76) then
       call king6fsr(7,4,ma,zero,3,6,mz,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+77) then
       call king6isr(4,7,ma,zero,5,8,mz,gamz,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+78) then
       call king6fsr(4,7,ma,zero,5,8,mz,gamz,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+79) then
       call king6fsr(7,4,ma,zero,5,8,mz,gamz,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+80) then
       call king6isr(4,7,ma,zero,5,8,mz,gamz,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+81) then
       call king6fsr(4,7,ma,zero,5,8,mz,gamz,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+82) then
       call king6fsr(7,4,ma,zero,5,8,mz,gamz,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+83) then
       call king6isr(4,7,ma,zero,5,8,mz,zero,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+84) then
       call king6fsr(4,7,ma,zero,5,8,mz,zero,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+85) then
       call king6fsr(7,4,ma,zero,5,8,mz,zero,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+86) then
       call king6isr(4,7,ma,zero,5,8,mz,zero,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+87) then
       call king6fsr(4,7,ma,zero,5,8,mz,zero,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+88) then
       call king6fsr(7,4,ma,zero,5,8,mz,zero,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+89) then
       call king6isr(3,6,mz,gamz,4,7,mz,gamz,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+90) then
       call king6fsr(3,6,mz,gamz,4,7,mz,gamz,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+91) then
       call king6fsr(6,3,mz,gamz,4,7,mz,gamz,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+92) then
       call king6isr(3,6,mz,gamz,4,7,mz,gamz,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+93) then
       call king6fsr(3,6,mz,gamz,4,7,mz,gamz,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+94) then
       call king6fsr(6,3,mz,gamz,4,7,mz,gamz,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+95) then
       call king6isr(3,6,mz,gamz,5,8,mz,gamz,4,7,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+96) then
       call king6fsr(3,6,mz,gamz,5,8,mz,gamz,4,7,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+97) then
       call king6fsr(6,3,mz,gamz,5,8,mz,gamz,4,7,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+98) then
       call king6isr(3,6,mz,gamz,4,7,mz,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+99) then
       call king6fsr(3,6,mz,gamz,4,7,mz,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+100) then
       call king6fsr(6,3,mz,gamz,4,7,mz,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+101) then
       call king6isr(3,6,mz,gamz,7,8,mw,gamw,4,5,mw,gamw,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+102) then
       call king6fsr(3,6,mz,gamz,7,8,mw,gamw,4,5,mw,gamw,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+103) then
       call king6fsr(6,3,mz,gamz,7,8,mw,gamw,4,5,mw,gamw,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+104) then
       call king6isr(3,6,mz,gamz,7,8,mw,gamw,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+105) then
       call king6fsr(3,6,mz,gamz,7,8,mw,gamw,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+106) then
       call king6fsr(6,3,mz,gamz,7,8,mw,gamw,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+107) then
       call king6isr(3,6,mz,gamz,4,5,mw,gamw,7,8,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+108) then
       call king6fsr(3,6,mz,gamz,4,5,mw,gamw,7,8,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+109) then
       call king6fsr(6,3,mz,gamz,4,5,mw,gamw,7,8,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+110) then
       call king6isr(3,6,mz,gamz,7,8,mw,zero,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+111) then
       call king6fsr(3,6,mz,gamz,7,8,mw,zero,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+112) then
       call king6fsr(6,3,mz,gamz,7,8,mw,zero,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+113) then
       call king6isr(3,6,mz,zero,4,7,mz,gamz,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+114) then
       call king6fsr(3,6,mz,zero,4,7,mz,gamz,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+115) then
       call king6fsr(6,3,mz,zero,4,7,mz,gamz,5,8,mz,gamz,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+116) then
       call king6isr(3,6,mz,zero,4,7,mz,gamz,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+117) then
       call king6fsr(3,6,mz,zero,4,7,mz,gamz,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+118) then
       call king6fsr(6,3,mz,zero,4,7,mz,gamz,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+119) then
       call king6isr(3,6,mz,zero,5,8,mz,gamz,4,7,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+120) then
       call king6fsr(3,6,mz,zero,5,8,mz,gamz,4,7,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+121) then
       call king6fsr(6,3,mz,zero,5,8,mz,gamz,4,7,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+122) then
       call king6isr(3,6,mz,zero,4,7,mz,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+123) then
       call king6fsr(3,6,mz,zero,4,7,mz,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+124) then
       call king6fsr(6,3,mz,zero,4,7,mz,zero,5,8,mz,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+125) then
       call king6isr(3,6,mz,zero,7,8,mw,gamw,4,5,mw,gamw,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+126) then
       call king6fsr(3,6,mz,zero,7,8,mw,gamw,4,5,mw,gamw,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+127) then
       call king6fsr(6,3,mz,zero,7,8,mw,gamw,4,5,mw,gamw,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+128) then
       call king6isr(3,6,mz,zero,7,8,mw,gamw,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+129) then
       call king6fsr(3,6,mz,zero,7,8,mw,gamw,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+130) then
       call king6fsr(6,3,mz,zero,7,8,mw,gamw,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+131) then
       call king6isr(3,6,mz,zero,4,5,mw,gamw,7,8,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+132) then
       call king6fsr(3,6,mz,zero,4,5,mw,gamw,7,8,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+133) then
       call king6fsr(6,3,mz,zero,4,5,mw,gamw,7,8,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+134) then
       call king6isr(3,6,mz,zero,7,8,mw,zero,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+135) then
       call king6fsr(3,6,mz,zero,7,8,mw,zero,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+136) then
       call king6fsr(6,3,mz,zero,7,8,mw,zero,4,5,mw,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+137) then
       call king6isr(4,7,mz,gamz,5,8,mz,gamz,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+138) then
       call king6fsr(4,7,mz,gamz,5,8,mz,gamz,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+139) then
       call king6fsr(7,4,mz,gamz,5,8,mz,gamz,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+140) then
       call king6isr(4,7,mz,gamz,5,8,mz,gamz,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+141) then
       call king6fsr(4,7,mz,gamz,5,8,mz,gamz,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+142) then
       call king6fsr(7,4,mz,gamz,5,8,mz,gamz,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+143) then
       call king6isr(4,7,mz,gamz,5,8,mz,zero,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+144) then
       call king6fsr(4,7,mz,gamz,5,8,mz,zero,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+145) then
       call king6fsr(7,4,mz,gamz,5,8,mz,zero,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+146) then
       call king6isr(4,7,mz,gamz,5,8,mz,zero,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+147) then
       call king6fsr(4,7,mz,gamz,5,8,mz,zero,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+148) then
       call king6fsr(7,4,mz,gamz,5,8,mz,zero,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+149) then
       call king6isr(5,8,mz,gamz,4,7,mz,zero,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+150) then
       call king6isr(5,8,mz,gamz,4,7,mz,zero,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+151) then
       call king6isr(4,7,mz,zero,5,8,mz,zero,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+152) then
       call king6fsr(4,7,mz,zero,5,8,mz,zero,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+153) then
       call king6fsr(7,4,mz,zero,5,8,mz,zero,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+154) then
       call king6isr(4,7,mz,zero,5,8,mz,zero,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+155) then
       call king6fsr(4,7,mz,zero,5,8,mz,zero,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+156) then
       call king6fsr(7,4,mz,zero,5,8,mz,zero,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+157) then
       call king6isr(7,8,mw,gamw,4,5,mw,gamw,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+158) then
       call king6fsr(7,8,mw,gamw,4,5,mw,gamw,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+159) then
       call king6isr(7,8,mw,gamw,4,5,mw,gamw,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+160) then
       call king6fsr(7,8,mw,gamw,4,5,mw,gamw,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+161) then
       call king6isr(7,8,mw,gamw,4,5,mw,zero,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+162) then
       call king6fsr(7,8,mw,gamw,4,5,mw,zero,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+163) then
       call king6isr(7,8,mw,gamw,4,5,mw,zero,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+164) then
       call king6fsr(7,8,mw,gamw,4,5,mw,zero,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+165) then
       call king6isr(4,5,mw,gamw,7,8,mw,zero,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+166) then
       call king6fsr(4,5,mw,gamw,7,8,mw,zero,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+167) then
       call king6isr(4,5,mw,gamw,7,8,mw,zero,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+168) then
       call king6fsr(4,5,mw,gamw,7,8,mw,zero,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+169) then
       call king6isr(7,8,mw,zero,4,5,mw,zero,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+170) then
       call king6fsr(7,8,mw,zero,4,5,mw,zero,3,6,mh,gamh,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+171) then
       call king6isr(7,8,mw,zero,4,5,mw,zero,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      else if (ikin == nkinst+172) then
       call king6fsr(7,8,mw,zero,4,5,mw,zero,3,6,mh,zero,
     &                             imom,ecm,x,ndim,dlips)
      end if
 
      end
