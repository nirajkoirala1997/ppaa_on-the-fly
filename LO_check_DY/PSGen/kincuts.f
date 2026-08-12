      module kincuts   ! Kinematical cuts

      integer :: icuts,ichan,
     &           ictlb,ictqb,ictgb,ictab,ictpb,ictsb,ictal,ictaq,ictas,
     &           ictqq,ictlq,ictlbq,ictlg,ictll,ictls,
     &           ienl,iemt,ienq,ienb,ieng,iena,
     &           imqq,imjj,imll,imbb,imtth,imbbh,imnb,imbqq,imln,imbln,
     &           impp,
     &           iptl,iptj,iptb,irapl,irapj,irapb,ietal,ietaj,ietab,
     &           idelr
      real(kind(1.d0)) :: ctlb,ctqb,ctgb,ctab,ctpb,ctsb,ctal,ctaq,ctas,
     &                    ctqq,ctlq,ctlbq,ctlg,ctll,ct0b,ctls,
     &                    ecutl,ecutt,ecutq,ecutb,ecutg,ecuta,ecut0,
     &                    mqq,mjj,mll,mtt,mww,mwp1,mwp2,mwm2,mbb,mbbh,
     &                    mbqq,mqq1,mqq2,mblnt1,mblnt2,mln1,mln2,
     &                    mbln1,mbln2,mpp1,mpp2,
     &                    ptl,ptj,ptb,etal,etaj,etab,rapl,rapj,rapb,
     &                    delr
      real(kind(1.d0)) :: ctc,sqscut

      end module kincuts
