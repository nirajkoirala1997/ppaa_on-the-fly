      subroutine cutinfo(iu)
c=============================================================================c
c Writes cut information
c
c Karol Kolodziej, University of Silesia, karol.kolodziej@us.edu.pl
c 8 June, 2010, modiffied on 
c
c=============================================================================c

      use kincuts
      use prtinds

      implicit none

      integer :: iu
      real(dp) :: theta

      write(iu,'(a)') '#-----------------------------------------------'
     &    //'------------------------------#'
      write(iu,'(a)') '# The following cuts have been imposed:'

      if (ictlb == 1) then
       if (nxl > 0) then
        write(iu,'(a)') '# Lepton-beam angle (deg):'
        theta=acos(ctlb)/pi*180.d0
        write(iu,8050) theta,180.d0-theta,ctlb
 8050   format('# ',f11.6,' < Th_lb <',f11.6,'  =>  |cos(Th_lb)|',
     &  ' <',f19.16)
       end if
      end if

      if (ictqb == 1) then
       if (nxq > 0 .or. nxb > 0) then
        write(iu,'(a)') '# Quark-beam angle (deg):'
        theta=acos(ctqb)/pi*180.d0
        write(iu,8052) theta,180.d0-theta,ctqb
 8052   format('# ',f11.6,' < Th_qb <',f11.6,'  =>  |cos(Th_qb)|',
     &  ' <',f19.16)
       end if
      end if

      if (ictgb == 1) then
       if (nxg > 0) then
        write(iu,'(a)') '# Gluon-beam angle (deg):'
        theta=acos(ctgb)/pi*180.d0
        write(iu,8053) theta,180.d0-theta,ctgb
 8053   format('# ',f11.6,' < Th_gb <',f11.6,'  =>  |cos(Th_gb)|',
     &  ' <',f19.16)
       end if
      end if

      if (ictab == 1) then
       if (nxa > 0) then
        write(iu,'(a)') '# Photon-beam angle (deg):'
        theta=acos(ctab)/pi*180.d0
        write(iu,8054) theta,180.d0-theta,ctab
 8054   format('# ',f11.6,' < Th_ab <',f11.6,'  =>  |cos(Th_ab)|',
     &  ' <',f19.16)
       end if
      end if

      if (ictpb == 1) then
       if (nxp > 0) then
        write(iu,'(a)') '# Proton-beam angle (deg):'
        theta=acos(ctpb)/pi*180.d0
        write(iu,8055) theta,180.d0-theta,ctpb
 8055   format('# ',f11.6,' < Th_pb <',f11.6,'  =>  |cos(Th_pb)|',
     &  ' <',f19.16)
       end if
      end if

      if (ictsb == 1) then
       if (nxs > 0) then
        write(iu,'(a)') '# Scalar-beam angle (deg):'
        theta=acos(ctsb)/pi*180.d0
        write(iu,'(a,f10.4,a,f10.4,a,f19.16)') '# ',theta,
     &      ' < Th_sb <',180.d0-theta,'  =>  |cos(Th_sb)| <',ctsb
       end if
      end if

      if (ictal == 1) then
       if (nxa*nxl > 0) then
        write(iu,'(a)') '# Photon-lepton angle (deg):'
        theta=acos(ctal)/pi*180.d0
        write(iu,8056) theta,ctal
 8056   format('# ',f11.6,' < Th_al','  =>  cos(Th_al) <',f19.16)
       end if
      end if

      if (ictaq == 1) then
       if (nxa*nxq > 0 .or. nxa*nxb > 0) then
        write(iu,'(a)') '# Photon-quark angle (deg):'
        theta=acos(ctaq)/pi*180.d0
        write(iu,8058) theta,ctaq
 8058   format('# ',f11.6,' < Th_aq','  =>  cos(Th_aq) <',f19.16)
       end if
      end if

      if (ictas == 1) then
       if (nxa*nxs > 0) then
        write(iu,'(a)') '# Photon-scalar angle (deg):'
        theta=acos(ctas)/pi*180.d0
        write(iu,'(a,f11.6,a,f19.16)') '# ',theta,' < Th_as '//
     &       ' =>  cos(Th_as) <',ctas
       end if
      end if

      if (ictqq == 1) then
       if (nxq+nxb > 1) then
        write(iu,'(a)') '# Quark-quark'' angle (deg):'
        theta=acos(ctqq)/pi*180.d0
        write(iu,8069) theta,ctqq
 8069   format('# ',f11.6,' < Th_qq','  =>  cos(Th_qq) <',f19.16)
       end if
      end if

      if (ictlq == 1) then
       if (nxl*(nxq+nxb) > 0) then
        write(iu,'(a)') '# Lepton-quark angle (deg):'
        theta=acos(ctlq)/pi*180.d0
        write(iu,8057) theta,ctlq
 8057   format('# ',f11.6,' < Th_lq','  =>  cos(Th_lq) <',f19.16)
       end if
      end if

      if (ictlg == 1) then
       if (nxl*nxg > 0) then
        write(iu,'(a)') '# Lepton-gluon angle (deg):'
        theta=acos(ctlg)/pi*180.d0
        write(iu,8158) theta,ctlg
 8158   format('# ',f11.6,' < Th_lg','  =>  cos(Th_lg) <',f19.16)
       end if
      end if

      if (ictls == 1) then
       if (nxl*nxs > 0) then
        write(iu,'(a)') '# Lepton-scalar angle (deg):'
        theta=acos(ctls)/pi*180.d0
        write(iu,8159) theta,ctls
 8159   format('# ',f11.6,' < Th_ls','  =>  cos(Th_ls) <',f19.16)
       end if
      end if

      if (ictll == 1) then
       if (nxl > 1) then
        write(iu,'(a)') '# Lepton-lepton'' angle (deg):'
        theta=acos(ctll)/pi*180.d0
        write(iu,8059) theta,ctll
 8059   format('# ',f11.6,' < Th_ll','  =>  cos(Th_ll) <',f19.16)
       end if
      end if

      if (ienl == 1) then
       if (nxl > 0) then
        write(iu,8230) ecutl
 8230   format('# Lepton energy:    E_l >',f10.4,' GeV')
       end if
      end if

      if (iemt == 1) then
       if (nxn > 0) then
        write(iu,'(a,f10.4,a)')
     &     '# Missing transverse energy:   ET_miss >',ecutt,' GeV'
       end if
      end if

      if (ienq == 1) then 
       if (nxq+nxb > 0) then
        write(iu,8232) ecutq
 8232   format('# Quark energy:    E_q >',f10.4,' GeV')
       end if
      end if

      if (ienb == 1) then 
       if (nxb > 0) then
        write(iu,8332) ecutb
 8332   format('# b-quark energy:    E_b >',f10.4,' GeV')
       end if
      end if

      if (ieng == 1) then 
       if (nxg > 0) then
        write(iu,8233) ecutg
 8233   format('# Gluon energy:    E_g >',f10.4,' GeV')
       end if
      end if

      if (iena == 1) then
       if (nxa > 0) then
        write(iu,8234) ecuta
 8234   format('# Photon energy:    E_a >',f10.4,' GeV')
       end if
      end if

      if (imqq == 1) then
       if (nxq+nxb > 1) then
        write(iu,8400) mqq
 8400   format('# Quark-quark invariant mass:  m(q,q'') >',f10.4,
     &                                                      ' GeV')
       end if
      end if

      if (imbb == 1) then
       if (nxb > 1) then
        write(iu,8402) mbb
 8402   format('# bb'' invariant mass:  m(b,b'') >',f10.4,
     &                                                      ' GeV')
       end if
      end if

      if (impp == 1) then
       if (nxp > 1) then
        write(iu,8404) mpp1,mpp2
 8404   format('# pp~ invariant mass: ',f10.4,' < Q^2(p,p~) <',f10.4,
     &                                                      ' GeV^2')
       end if
      end if

      if (imjj == 1) then
       if (nxq+nxb+nxg > 1) then
        write(iu,'(a,f10.4,a)') '# Jet-jet invariant mass:  m(j,j'') >',
     &       mjj,' GeV'
       end if
      end if

      if (imll == 1) then
       write(iu,8410) mll
 8410  format('# Lepton-lepton invariant mass:  m(l,l'') >',f10.4,
     &                                                      ' GeV')
      end if

      if (iptl == 1) then
       if (nxl > 0) then
        write(iu,'(a,f10.4,a)') '# Lepton transverse momentum:  '
     & //'pT_l  >',ptl,' GeV'
       end if
      end if
      if (iptj == 1) then
       if (nxj > 0) then
        write(iu,'(a,f10.4,a)')'# Light quark or gluon transverse '
     &                        //'momentum:  pT_j  >',ptj,' GeV'
       end if
      end if
      if (iptb == 1) then
       if (nxb > 0) then
        write(iu,'(a,f10.4,a)')'# b-quark transverse '
     &                       //'momentum:  pT_b  >',ptb,' GeV'
       end if
      end if
      if (ietal == 1) then
       if (nxl > 0) then
        write(iu,'(a,f10.4)')'# Lepton pseudorapidity:  '
     &  //'|eta_l|  <',etal
       end if
      end if
      if (ietaj == 1) then
       if (nxj > 0) then
        write(iu,'(a,f10.4)')'# Light quark or gluon pseudorapidity:  '
     &  //'|eta_j|  <',etaj
       end if
      end if
      if (ietab == 1) then
       if (nxb > 0) then
        write(iu,'(a,f10.4)')'# b_quark pseudorapidity:  '
     &  //'|eta_b|  <',etab
       end if
      end if
      if (irapl == 1) then
       if (nxl > 0) then
        write(iu,'(a,f10.4)')'# Lepton rapidity:  '
     &  //'|rap_l|  <',rapl
       end if
      end if
      if (irapj == 1) then
       if (nxj > 0) then
        write(iu,'(a,f10.4)')'# Light quark or gluon rapidity:  '
     &  //'|rap_j|  <',rapj
       end if
      end if
      if (irapb == 1) then
       if (nxb > 0) then
        write(iu,'(a,f10.4)')'# b-quark rapidity:  '
     &  //'|rap_b|  <',rapb
       end if
      end if
      if (idelr == 1) then
       if (nxj+nxb+nxl > 0) then
        write(iu,'(a)')'# Separation in the azimuthal '
     &  //'angle-pseudorapidity plane of two leptons,'
        write(iu,'(a,f10.4)')'# two jets and a lepton and a jet,'
     &  //' Delta R  >',delr
       end if
      end if

      if (imtth == 1) then
       if (imbqq == 1) then
        if (nxb > 1) then
         write(iu,'(a,f10.4,a)') '# bb'' invariant mass:  '//
     &                          '|m(b,b'') - mh| <',mbbh,' GeV'
        end if
       end if
       if (imnb == 1) then
        if (nxq > 1) then
         write(iu,'(a,f10.4,a,f10.4,a)') '# qq'' invariant mass: ',
     &                        mqq1,' GeV < m(q,q'') <',mqq2,' GeV'
         if (imbqq == 1) then
          if (nxb > 0) then
           write(iu,'(a,f10.4,a)') '# bqq'' invariant mass:  '//
     &                 '|m(b,q,q'') - mt| <',mbqq,' GeV'
          end if
         end if
        end if
       end if
       if (imln == 1) then
        if (nxl*nxn > 0) then
         write(iu,'(a,f10.4,a,f10.4,a)') '# ln transverse mass: ',
     &                        mln1,' GeV < mT(l,n) <',mln2,' GeV'
         if (imbln == 1) then
          if (nxb > 0) then
           write(iu,'(a,f8.2,a,f8.2,a)') '# bln transverse mass: mt -'
     &     ,mbln1,' GeV < mT(b,l,n) < mt +', mbln2,' GeV'
          end if
         end if
        end if
       end if
      end if

      write(iu,'(a)') '#-----------------------------------------------'
     &    //'------------------------------#'

      end subroutine cutinfo
