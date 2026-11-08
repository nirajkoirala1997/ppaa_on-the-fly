      module kinparams
 
      use inprms_ps
 
      integer, parameter :: iqprec=0
      integer, parameter :: xp=kind(1.d0)
      integer, parameter :: nexp=9
      integer, parameter :: nfsp=7
      integer, parameter :: nkins=64
      integer, parameter :: nkinst=64
      integer, parameter :: nkinag=236
      integer, parameter :: nkin=236
      integer, parameter :: igind=9
      real(xp) :: m(nexp)
      character(10) :: expc(nexp)
      character(1) :: expl(nexp)
 
      integer :: iwkin(nkin)
      real(kind(1.d0)) :: aw(nkin), xmc(0:nkin)
      real(kind(1.d0)) :: aw_w(nkin,3)
 
      contains
       subroutine exp_masses
       m=(/me,me,mb0,mmu,zero,mb0,mmu,zero,ma/)
       end subroutine exp_masses
       subroutine exp_names
        expc(1)='e+'
        expl(1)='l'
        expc(2)='e-'
        expl(2)='l'
        expc(3)='b'
        expl(3)='q'
        expc(4)='mu+'
        expl(4)='l'
        expc(5)='vm'
        expl(5)='n'
        expc(6)='b~'
        expl(6)='q'
        expc(7)='mu-'
        expl(7)='l'
        expc(8)='vm~'
        expl(8)='n'
        expc(9)='a'
        expl(9)='a'
       end subroutine exp_names
 
      end module kinparams
