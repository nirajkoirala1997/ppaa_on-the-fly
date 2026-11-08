      module inprms_ps

      integer, parameter :: dp=kind(1.d0), qp=kind(1.q0)
      real(dp), parameter :: zero=0.d0,
     &                       pi=3.141592653589793238d0,
     &                       twopi=2.d0*pi,
     &                       sqrt2=1.414213562373095d0      

c Conversion costant (pb*GeV^2):

      real(dp), parameter ::   convc=0.3893793656d12

c Intermediate particle masses and widths:'
      real(dp), parameter ::  
     &     ma=zero, gama=zero,                        ! photon
     &     mg=zero, gamg=zero,                        ! gluon
     &     mz=91.1876d0,  gamz=2.4952d0,              ! Z 
     &     mw=80.379d0,   gamw=2.085d0,               ! W 
     &     mh=125.18d0,   gamh=0.5829111764397600d-2, ! Higgs
     &     mt=173.0d0,    gamt=1.5d0,                 ! top
     &     sqsc=1.d-8                                 ! minimum virtuality

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

