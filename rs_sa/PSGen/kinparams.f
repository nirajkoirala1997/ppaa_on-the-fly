      module kinparams
 
      use inprms_ps
 
      integer, parameter :: iqprec=0
      integer, parameter :: xp=kind(1.d0)
      integer, parameter :: nexp=5
      integer, parameter :: nfsp=3
      include 'nkin_decl.f'
      integer, parameter :: igind=5
      real(xp) :: m(nexp)
      character(10) :: expc(nexp)
      character(1) :: expl(nexp)
 
      integer :: iwkin(nkin)
      real(kind(1.d0)) :: aw(nkin)
      real(kind(1.d0)) :: aw_w(nkin,3)
 
      contains
       subroutine exp_masses
       m=(/mg,mg,mz,mh,mg/)
       end subroutine exp_masses
       subroutine exp_names
        expc(1)='g'
        expl(1)='g'
        expc(2)='g'
        expl(2)='g'
        expc(3)='z'
        expl(3)='z'
        expc(4)='h'
        expl(4)='h'
        expc(5)='g'
        expl(5)='g'
       end subroutine exp_names
 
      end module kinparams
