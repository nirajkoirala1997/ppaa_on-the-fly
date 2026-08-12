      program PSGen
c----------------------------------------------------------------------c
c Generates the multichannel Monte Carlo kinematical routine for 
c the process
c
c                 f_1 f_2 -> f_3 f_4 ... f_nptcl
c
c taking into account poles due to the Feynman propagators of unstable 
c heavy particles which decay into on shell final state particles 
c according to patterns listed in genps.dat, possible t-channel poles 
c and peaks due to the on shell photon/gluon radiation
c
c Karol Kolodziej, University of Silesia, Katowice, Poland
c karol.kolodziej@us.edu.pl, www: https://www.kk.us.edu.pl'
c January 2022 - March 2023; modified 
c----------------------------------------------------------------------c

      use psgen_info
      use proc_name

      implicit none

      character(80) :: process

c Specify a process below. The maximum number of external particles is 11. 
c Standard Model particles include:
c     Quarks:   d u s c b t d~ u~ s~ c~ b~ t~
c     Leptons:  e- mu- ta- e+ mu+ ta+ ve vm vt ve~ vm~ vt~
c     Bosons:   g a z w+ w- h
c Other particles:
c     Nucleons: p, p~, n, n~
c     Mesons: pi-, pi+, pi0, k-, k+, k0, k0~, rho-, rho+, rho0, ome, phi,
c             rho1, rho2
c Define process:

      process='u u~ -> e+ e-' 

c or it will be overwritten, if transfered fro CARLOMAT

      if (len_trim(proc) > 0) process=proc

      write(6,*)'process=',trim(process)

c Generate t-channel kinematics, for processes with identical fermion-antifermion
c pairs in the initial and final state, itchnl=1 (yes, only for double t-channel 
c poles)/2 (yes, both for double and single t-channel poles))/else (no)?

      itchnl=0

c Generate single photon radiation kinematics, for processes with one photon 
c or gluon in the final state, 
c igchnl=1 (yes)/
c       =2 (yes, taking into account poles listed in genps.dat 
c           and information collected in particles.dat)/else (no)?

      igchnl=0

c Use quartic precision for the particle four momenta, masses and widths,
c iquadp=1 (yes)/else (no)?:
 
      iquadp=0

c Short name of the process:

      nicknm='auto'  ! If prcsnm='auto' then it is made automatically.
 
      call read_process(process)
      call genps(nfspt)

      end

