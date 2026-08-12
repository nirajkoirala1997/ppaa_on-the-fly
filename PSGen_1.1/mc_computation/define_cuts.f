      subroutine define_cuts

      use kincuts
      use kinparams

      implicit none

c      real(kind(1.d0)), parameter :: pi=3.141592653589793238d0

      if (icuts == 1) then
c Impose cuts on:
c  transverse momentum of the final state lepton (GeV):
      iptl=0
      ptl=30.d0
c  transverse momentum of the final state jet: light quark or gluon (GeV):
      iptj=0
      ptj=30.d0
c  transverse momentum of the final state b-jet (GeV):
      iptb=0
      ptb=30.d0
c  pseudorapidity of the final state lepton:
      ietal=0
      etal=2.1d0
c  pseudorapidity of the final state jet: light quark or gluon:
      ietaj=0
      etaj=2.4d0
c  pseudorapidity of the final state b-jet:
      ietab=0
      etab=2.4d0
c separation in the azimuthal angle-pseudorapidity plane of two leptons,
c two jets and a lepton and a jet:
      idelr=0
      delr=0.4d0
c  missing transverse energy (GeV):
      iemt=0
      ecutt=20.d0
c  rapidity of the final state lepton:
      irapl=0
      rapl=2.d0
c  rapidity of the final state jet: light quark or gluon:
      irapj=0
      rapj=2.5d0
c  rapidity of the final state b-jet:
      irapb=0
      rapb=2.5d0
c  cos theta(lepton,beam) - ictlb=1/0:
      ictlb=1
c      ctlb=0.985d0   ! 
      ctlb=cos(5.d0/180.d0*pi)
c  cos theta(quark,beam), including b - ictqb=1/0:
      ictqb=0
c      ctqb=0.985d0   ! 
      ctqb=cos(5.d0/180.d0*pi)
c  cos theta(gluon,beam) - ictgb=1/0:
      ictgb=0
c      ctgb=0.985d0   ! 
      ctgb=ctqb
c  cos theta(photon,beam) - ictab=1/0:
      ictab=1
c      ctab=0.985d0
      ctab=cos(5.d0/180.d0*pi)
c  cos theta(proton,beam) - ictpb=1/0:
      ictpb=0
c      ctpb=0.985d0
      ctpb=cos(30.d0/180.d0*pi)
c  cos theta(scalar,beam) - ictsb=1/0:
      ictsb=0
c      ctsb=0.985d0
      ctsb=cos(5.d0/180.d0*pi)
c  cos theta(photon,lepton) - ictal=1/0:
      ictal=1
      ctal=cos(5.d0/180.d0*pi)
c  cos theta(photon,quark), including b - ictaq=1/0:
      ictaq=0
      ctaq=cos(5.d0/180.d0*pi)
c  cos theta(photon,scalar) - ictas=1/0:
      ictas=0
      ctas=cos(5.d0/180.d0*pi)
c  cos theta(quark,quark'), including b - ictqq=1/0:
      ictqq=0
c      ctqq=0.985d0   ! 
      ctqq=cos(10.d0/180.d0*pi)
c  cos theta(lepton,light quark) - ictlq=1/0:
      ictlq=0
c      ctlq=0.985d0   ! 
      ctlq=cos(5.d0/180.d0*pi)
c  cos theta(lepton,b-quark) - ictlbq=1/0:
      ictlbq=0
c      ctlbq=0.985d0   ! 
      ctlbq=cos(1.d0/180.d0*pi)
c  cos theta(lepton,gluon) - ictlg=1/0:
      ictlg=0
c      ctlg=0.985d0   ! 
      ctlg=ctlq
c  cos theta(lepton,scalar) - ictls=1/0:
      ictls=0
c      ctls=0.985d0   ! 
      ctls=cos(5.d0/180.d0*pi)
c  cos theta(lepton,lepton') - ictll=1/0:
      ictll=0
c      ctll=0.985d0   ! 
      ctll=cos(5.d0/180.d0*pi)
c  lepton energy (GeV):
      ienl=1
      ecutl=5.d0
c  light quark energy (GeV):
      ienq=0
      ecutq=10.d0
c  b-quark energy (GeV):
      ienb=0
      ecutb=15.d0
c  gluon energy (GeV):
      ieng=0
      ecutg=1.d0 !ecutq
c  photon energy (GeV):
      iena=1
      ecuta=1.0d0
c  quark-quark' invariant mass, including b (GeV):
      imqq=0
      mqq=10.d0
c  jet-jet' invariant mass, including b (GeV):
      imjj=0
      mjj=10.d0
c  lepton-lepton' invariant mass (GeV):
      imll=0
      mll=10.d0
c  b-b' invarint mass (GeV):
      imbb=0
      mbb=20.d0
c  angular cut for t-channel processes, cos(p_i,p_f)
      ctc=ctlb
      if (ictlb == 0) ctc=1.d0
c invariant mass cut, used in mappings of the photon and gluon propagators 
c that couple to massless fermions relevant if light fermion masses are zero
      sqscut=1.d-8

      call exp_names

      else ! icuts
c no angular cut for t-channel processes, cos(p_i,p_f)
       ctc=1.d0
c invariant mass cut, used in mappings of the photon and gluon propagators 
c that couple to massless fermions relevant if light fermion masses are zero
       sqscut=1.d-8
      end if ! icuts

      if (icuts == 1) then
       if (igind > 0) then
        if (expl(igind) == 'a') then
         if (ictab == 1) then
          ct0b=ctab
         else
          ct0b=1.d0
         end if
         if (iena == 1) ecut0=ecuta
        else if (expl(igind) == 'g') then
         if (ictgb == 1) then
          ct0b=ctgb
         else
          ct0b=1.d0
         end if
         if (ieng == 1) ecut0=ecutg
        end if
       end if
      else
       ecut0=1.d0
      end if



      end
