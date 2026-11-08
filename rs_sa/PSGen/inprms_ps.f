      module inprms_ps

      integer, parameter :: dp=kind(1.d0), qp=kind(1.q0)
      real(dp), parameter :: zero=0.d0,
     &                       pi=3.141592653589793238d0,
     &                       twopi=2.d0*pi,
     &                       sqrt2=1.414213562373095d0      

c Conversion costant (pb*GeV^2):

      real(dp), parameter ::   convc=0.3893793656d12

c Intermediate particle masses and widths:'
c
c LOCAL EDIT (mz,gamz,mw,gamw). COMPILE-TIME constants the run card cannot
c reach, set here to ../input.input's values -- same edit and same reasoning
c as ../../PSGen/inprms_ps.f (the born set this directory's own PSGen was
c copied from): mZ is an EXTERNAL ON-SHELL leg (kin3_2 hands m(3) straight to
c kinff), and psgen_init REFUSES TO RUN if this disagrees with the card, since
c an off-shell external leg is a wrong answer, not a slow one. mW is inert
c here (no W appears in this final state) and kept only so this file and
c ../../PSGen/inprms_ps.f/../../kp_sa's card stay comparable.
c
c mt, gamt UNCHANGED FROM UPSTREAM, AND NOT A PHYSICS INPUT HERE. genps.dat's
c new technical `2 z h : mt,zero` entry (see PSGen_1.1/code_generation/genps.dat,
c "z-h TECHNICAL pairing") reuses mt/zero purely as an already-declared
c mass/width pair to drive kin3_2's invariant-mass mapping for the {Z,H}
c sub-system -- gv1=zero is genps.dat's own established convention for a flat
c (non-resonant) mapping (see its repeated ",zero" variants alongside every
c ",gamX" one), not a Breit-Wigner peak. mt=173 GeV sits below the physical
c Z+H threshold (216.4 GeV), which costs VEGAS efficiency, not correctness --
c revisit only if convergence turns out to need it.
      real(dp), parameter ::
     &     ma=zero, gama=zero,                        ! photon
     &     mg=zero, gamg=zero,                        ! gluon
     &     mz=91.188d0,   gamz=1.0d-12,                ! Z
     &     mw=80.3692d0,  gamw=1.0d-12,                ! W
     &     mh=125.2d0,    gamh=0.5829111764397600d-2, ! Higgs
     &     mt=173.0d0,    gamt=1.5d0                   ! top, and this build's technical z-h mass

c Minimum virtuality: the hard lower edge of every s-channel invariant the
c mappings can produce. LOCAL EDIT: made a variable so psgen_interface can set
c it from the run card -- same edit as ../../PSGen/inprms_ps.f, but NOT inert
c here the way it is there: kin3_2 DOES link mapckk (the born build's kin2_2
c never did), so sqsc is the genuine lower edge of the {Z,H} invariant-mass
c window mapckk samples.
      real(dp) :: sqsc=1.d-8

c LOCAL EDIT: the m(Z,gamma) window kin3_2 samples the {3,4} pair in, set by
c psgen_interface from the run card's [cuts] mzgamma_min/max. 0 = no edge
c (the full [m3+m4, ecm-m5] range, as before). Sampling only inside the cut
c window instead of cutting afterwards: the II dipole map leaves m(Z,gamma)
c unchanged, so no real or dipole point inside the window is lost.
      real(dp) :: sqs_win_lo=0.d0, sqs_win_hi=0.d0

c Fermion masses:

      real(dp), parameter :: 
     &       me=0.5109989461d-3, mmu=105.6583745d-3, mtau=1.77686d0,
     &       mu0=5.d-3, md0=10.d-3, ms0=200.d-3, 
     &       mc0=1.275d0, mb0=4.18d0

c Nucleon masses:
c Note that nucleon magnetic moments and vector boson masses are specified
c in nuclff_phok.f in directory carloilb

      real(dp), parameter :: mp=938.272081d-3, mn=939.565413d-3

c Meson masses:

      real(dp), parameter :: mpi=139.57018d-3, mpi0=134.9766d-3,
     &                       mk=493.677d-3,mk0=497.614d-3,
     &                       mrho=775.26d-3,mome=782.65d-3,
     &                       mphi=1019.461d-3,
     &                       mrh1=1465.d-3,mrh2=1720.d-3
c     &                       mrho=774.37d-3,mome=782.65d-3,

c Meson widths:

      real(dp), parameter :: gampi=zero,gmrh=149.1d-3,
     &                       gmom=8.49d-3,
     &                       gmph=4.266d-3,
     &                       gmr1=400.d-3,gmr2=250.d-3
c      real(dp), parameter :: gampi=zero,gmrh=144.21d-3,gmom=8.49d-3

c     &       gampi=1.054571628d-25/2.6033d-8

c Masses for the quadruple precision kinematics, are automatically determined
c from the corresponding masses in double precision defined above:

      real(qp), parameter :: zero_qp=0.q0, ma_qp=0.q0, mg_qp=0.q0
      real(qp) :: mz_qp, mw_qp, mh_qp,
     & me_qp, mmu_qp, mtau_qp,
     & mu0_qp, md0_qp, ms0_qp, mc0_qp, mb0_qp, mt_qp,
     & mp_qp, mn_qp,
     & mpi_qp, mpi0_qp, mk_qp, mk0_qp,
     & mrho_qp, mome_qp, mphi_qp, mrh1_qp, mrh2_qp

c Flags and criteria of the mass and four momenta tests:
      integer :: im_test,ip_test
      real(dp) :: dm_test,dp_test

      real(dp) :: psnac

      end module inprms_ps

