      subroutine set_indices

      use prtinds

      implicit none

      integer :: i, nlept

c Determine numbers of the final state particles of different kind and
c their indices to be referenced in cut routines
      nxn=0  ! neutrinos
      nxl=0  ! charged leptons
      nxq=0  ! quarks /= b
      nxb=0  ! b-quarks
      nxa=0  ! photons
      nxg=0  ! gluons
      nxs=0  ! pions
      nxj=0  ! light quark or gluon
      nxp=0  ! protons
      ixn=0
      ixl=0
      ixq=0
      ixb=0
      ixa=0
      ixg=0
      ixs=0
      ixj=0
      ixn=0
      ixp=0
      iyl=0
      iyq=0
      iyb=0
      iya=0
      iyg=0
      iys=0
      iyj=0
      iyp=0
      do i=nexp-nfsp+1,nexp
       if (expl(i) == 'n') then
        nxn=nxn+1
        ixn(nxn)=i
        iyn(i)=i
       else if (expl(i) == 'l') then
        nxl=nxl+1
        ixl(nxl)=i
        iyl(i)=i
       else if (expl(i) == 'q' .and. expc(i)(1:1) /= 'b') then
        nxq=nxq+1
        ixq(nxq)=i
        iyq(i)=i
       else if (expc(i)(1:1) == 'b') then
        nxb=nxb+1
        ixb(nxb)=i
        iyb(i)=i
       else if (expl(i) == 'a') then
        nxa=nxa+1
        ixa(nxa)=i
        iya(i)=i
       else if (expl(i) == 'g') then
        nxg=nxg+1
        ixg(nxg)=i
        iyg(i)=i
       else if (expl(i) == 's') then
        nxs=nxs+1
        ixs(nxs)=i
        iys(i)=i
       else if (trim(expc(i)) == 'p' .or. trim(expc(i)) == 'p~') then
        nxp=nxp+1
        ixp(nxp)=i
        iyp(i)=i
       end if
       if ((expl(i) == 'q' .and. expc(i)(1:1) /= 'b') .or.
     &                                    expl(i) == 'g') then
        nxj=nxj+1
        ixj(nxj)=i
        iyj(i)=i
       end if
      end do

      end
