      module kinparams
 
      use inprms_ps
 
      integer, parameter :: iqprec=0
      integer, parameter :: xp=kind(1.d0)
      integer, parameter :: nexp=4
      integer, parameter :: nfsp=2
      integer, parameter :: nkins=1
      integer, parameter :: nkinst=1
      integer, parameter :: nkinag=1
      integer, parameter :: nkin=1
      integer, parameter :: igind=0
      real(xp) :: m(nexp)
      character(10) :: expc(nexp)
      character(1) :: expl(nexp)
 
      integer :: iwkin(nkin)
      real(kind(1.d0)) :: aw(nkin), xmc(0:nkin)
      real(kind(1.d0)) :: aw_w(nkin,3)
 
      contains
       subroutine exp_masses
       m=(/mu0,mu0,me,me/)
       end subroutine exp_masses
       subroutine exp_names
        expc(1)='u'
        expl(1)='q'
        expc(2)='u~'
        expl(2)='q'
        expc(3)='e+'
        expl(3)='l'
        expc(4)='e-'
        expl(4)='l'
       end subroutine exp_names
 
      end module kinparams
