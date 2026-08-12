      subroutine read_process(process)
c----------------------------------------------------------------------c
c Reads particles particles of a process and cheks if they are on a list
c in file 'particles.dat' and whether they couple to the photon or gluon
c
c Karol Kolodziej, University of Silesia, karol.kolodziej@us.edu.pl
c 26 October, 2021 
c modified 
c 1 Feb. 2022: writing of kincls.f improved
c 19 Sept. 2022: writing of kincls.f further improved such that no 
c calls to t-channel kinematics are generated for prcesses without
c t-channej poles, even if itchnl is accidentally set to 1 or 2
c----------------------------------------------------------------------c

      use psgen_info

      implicit none

      character(80) :: process,proc,cline,indp,indq
      character(1), parameter :: space=' '
      character(2), parameter :: separ='->'
      character(2) :: pldat
      character(10) :: ptdat,ind1,ind2,ind3,ind4,mdat,wdat
      integer :: i,j,k,l,n,ptmod,inmod,igvb

      proc=adjustl(process)

      j=0
      l=index(proc,separ)

      if (l == 0) then
       write(*,*)'Wrong separator, or process not entered.'
       write(*,*)'Please, try again!'
       stop
      end if


c Read initial state particles

      do 
       n=index(proc,space)
       j=j+1
       ptcl(j)=proc(1:n-1)
       proc=adjustl(proc(n:l+1))
       if (proc(1:2) == separ) exit
      end do

c No of the initial state particles

      nispt=j

c Read final state 

      k=l+2

      proc=adjustl(process(l+2:))

      do 
       n=index(proc,space)
       if (n > 1) then
        j=j+1
        ptcl(j)=proc(1:n-1)
       else
        exit
       end if
       proc=adjustl(proc(n:))
      end do

c No of particles in the process

      nptcl=j

      write(6,*)'nptcl=',nptcl

c No of particles in the process

      nfspt=nptcl-nispt

c Check if all the particles are on the list in file 'particles.dat' 
c and whether they couple to the photon and gluon

      do i=1,nptcl
       ptmod=1
       inmod=0
       open(unit=66,file='particles.dat',status='old',action='read')
       do
        read(66,'(a10,a10,a10,i2,i2,a2)',end=100) ptdat,mdat,wdat,j,k,
     &                                            pldat
        if (trim(ptcl(i)) == trim(ptdat)) then
         cmass(i)=mdat
         cwdth(i)=wdat
         ipa(i)=j
         ipg(i)=k
         ptclt(i)=adjustl(pldat)
         inmod=1
         exit
        end if
       end do
 100   close(66)
       ptmod=ptmod*inmod
       if (ptmod == 0) then
        write(*,*) 'Warning:',trim(ptcl(i)),
     &             ' is not listed in ''particles.dat'''
       end if
      end do

c Create short name for the process

      if (trim(nicknm) == 'auto') then
       nicknm=''
       do j=1,nispt
        nicknm=trim(adjustl(nicknm))//trim(adjustl(ptcl(j)))
       end do
       nicknm=trim(adjustl(nicknm))//'_'
       do j=nispt+1,nptcl
        nicknm=trim(adjustl(nicknm))//ptcl(j)(1:1)
       end do
      end if

c Check for a photon or gluon in the final state:
      igvb=0
      do i=nispt+1,nptcl,1
       if (ptcl(i) == 'a' .or. ptcl(i) == 'g') then
        igind2=i
        igvb=igvb+1
       end if
      end do
      if (igvb > 0) then
       if (igvb == 1) then
        igchnl=2
       else
        igchnl=0
       end if
      else
       igchnl=0
      end if

      open(unit=29,file='fourmom.f',status='replace',action='write')

      write(29,'(a)') '      module fourmom'
      write(29,'(a)') ' '
      write(ind1,'(i4)') nptcl
      if (iquadp == 1) then
       write(29,'(a)') '      real(kind(1.q0)) :: p(0:3,'//
     &                                    trim(adjustl(ind1))//')'
      else
       write(29,'(a)') '      real(kind(1.d0)) :: p(0:3,'//
     &                                    trim(adjustl(ind1))//')'
      end if
      write(29,'(a)') '      real(kind(1.d0)) :: r_dp(0:3,'//
     &                                    trim(adjustl(ind1))//')'
      write(29,'(a)') '      real(kind(1.q0)) :: r_qp(0:3,'//
     &                                    trim(adjustl(ind1))//')'
      write(29,'(a)') ' '
      write(29,'(a)') '      end module fourmom'
      close(29)

      open(unit=29,file='tchcalls.f',status='replace',action='write')

c Check for particles forming t-channel poles:
      itch1=0
      itch2=0
      ixtch1=0
      ixtch2=0
      do j=3,nptcl
       if (ptcl(j) == ptcl(1)) then
        itch1=itch1+1
        ixtch1(itch1)=j
       else if (ptcl(j) == ptcl(2)) then
        itch2=itch2+1
        ixtch2(itch2)=j
       end if
      end do

      if (itch1 == 0 .and. itch2 == 0) itchnl=0  ! 19 Sept. 2022

      if (itchnl == 1) then
       nkint2=0
c       if (nkint2 > 0) then ! Commented on 21 Jan 2019
       if (nptcl >= 4) then      ! Changed to >= instead of > on 16 Dec. 2022
        if (itch1 > 0 .or. itch2 > 0) then
c        if (itch1 > 0 .and. itch2 > 0) then ! 8 Feb. 2022
         write(ind1,'(i4)') nptcl-4
         write(29,'(a)') '      integer :: ix('//trim(adjustl(ind1))
     &                                                         //')'
         write(29,'(a)') ' '
        end if
        if (itch1 > 0 .and. itch2 > 0) then
         do i=1,itch1
          write(ind1,'(i4)') ixtch1(i)
          do j=1,itch2
           nkint2=nkint2+1
           write(ind3,'(i4)') nkint2
           if (nkint2 == 1) then
            write(29,'(a)') '      if (ikin == nkins+'
     &             //trim(adjustl(ind3))//') then'
           else
            write(29,'(a)') '      else if (ikin == nkins+'
     &            //trim(adjustl(ind3))//') then'
           end if
           write(ind2,'(i4)') ixtch2(j)
           l=0
           do k=3,nptcl
            if (k /= ixtch1(i) .and. k /= ixtch2(j)) then
             write(ind3,'(i4)') k
             l=l+1
             write(ind4,'(i4)') l
             write(29,'(a)') '       ix('//trim(adjustl(ind4))//')='
     &                       //trim(adjustl(ind3))
            end if
           end do
           if (l == 0) ind4='0'
           write(29,'(a)') '       call kintt3y2(1,'
     &          //trim(adjustl(ind1))//',1,2,'//trim(adjustl(ind2))//
     &         ',1,'//trim(adjustl(ind4))//')'
          end do
         end do
        end if
       else
         write(29,'(a)') 'c This file is empty for nptcl <= 4'
       end if
      else
       write(29,'(a)') 'c Assignements below should not play '//
     &                  'a role for itchnl /= 1'
       write(ind1,'(i4)') nptcl-2
       write(29,'(a)') '      integer :: ix('//trim(adjustl(ind1))
     &                                                         //')'
       l=0
       do k=3,nptcl
        l=l+1
        write(ind1,'(i4)') l
        write(ind2,'(i4)') k
        write(29,'(a)') '       ix('//trim(adjustl(ind1))//')='
     &                      //trim(adjustl(ind2))
       end do
       close(29)
      end if

c      

      open(unit=28,file='kincls.f',status='replace',action='write')

      write(28,'(a)') '      subroutine kincls(ikin,imom,x,ndim,'
     &                                                  //'dlips)'

      write(28,'(a)') ' '
c      write(28,'(a)') '      use partcls'
      write(28,'(a)') '      use kinparams'
      write(28,'(a)') '      use fourmom'

      write(28,'(a)') ' '
      write(28,'(a)') '      implicit none'
      write(28,'(a)') ' '
      write(28,'(a)') '      integer :: ikin,imom,ndim'!,ifr=1'

      write(indp,'(i16)') kind(1.d0)
      write(indq,'(i16)') kind(1.q0)
      if (iquadp == 1) then
       write(28,'(a)') '      real('//trim(adjustl(indp))//
     &                ') :: x(ndim),dlips'
       write(28,'(a)') '      real('//trim(adjustl(indq))//
     &                ') :: ecm'
      else
       write(28,'(a)') '      real('//trim(adjustl(indp))//
     &                ') :: ecm,x(ndim),dlips'
      end if
      write(28,'(a)') ' '
      write(28,'(a)') '      ecm=p(0,1)+p(0,2)'
      write(28,'(a)') ' '
c      write(28,'(a)') '      save ifr'
c      write(28,'(a)') ' '
c      write(28,'(a)') '      if (ifr == 1) then'
c      write(28,'(a)') '       call exp_masses'
c      write(28,'(a)') '       ifr=0'
c      write(28,'(a)') '      end if'
c      write(28,'(a)') ' '
      write(28,'(a)') '      if (ikin <= nkins) then'
      write(28,'(a)') '       call kinschnl'//
     &                 '(ikin,imom,ecm,x,ndim,dlips)'
c      write(28,'(a)') '      end if'    ! 1 Feb. 2022


c      if (itchnl == 2 .and. nptcl > 4 .and. nkint2 > 0) then
c      if (itchnl == 2 .and. nptcl > 4) then       ! 8 Feb. 2022
      if (itchnl == 1 .and. nptcl >= 4) then       ! 15 Dec. 2022
c       write(28,'(a)') '      if (ikin > nkins) then'
       write(28,'(a)') '      else if (ikin <= nkinst) then' ! 4 April 2022
       if (iquadp == 1) then
c        write(28,'(a)') '       call kintt3x2_qp'//
        write(28,'(a)') '       call kintchnl'//           ! 8 Feb. 2022
     &            '(ikin,imom,ecm,x,ndim,dlips)'
       else
c        write(28,'(a)') '       call kintt3x2'//
        write(28,'(a)') '       call kintchnl'//           ! 1 Feb. 2022
     &            '(ikin,imom,ecm,x,ndim,dlips)'
       end if
c       write(28,'(a)') '      end if'     ! 1 Feb. 2022
      end if
      if (igchnl == 2 .and. nptcl > 4) then
       write(28,'(a)') '      else if (ikin <= nkinag) then' ! 15 March 2023
       write(28,'(a)') '       call kingchnl'//
     &                 '(ikin,imom,ecm,x,ndim,dlips)'
      end if
      write(28,'(a)') '      else'                          ! 15 March 2023
      write(28,'(a)') '! Add a call to your own made kinematical '
     &  //'subroutine here and change '
      write(28,'(a)') '! the value of parameter nkin in kinparams.f '
     &  //'appropriately'
      write(28,'(a)') '      end if'                        ! 1 Feb. 2022
      write(28,'(a)') ' '
      write(28,'(a)') '      include ''kin_tests.f'''
      write(28,'(a)') ' '
      write(28,'(a)') '      end'

      close(28)

      open(unit=28,file='process_nm.dat',status='replace',
     &                                           action='write')

      write(28,'(a)') '      prcsnm='''//trim(adjustl(process))//''''
      write(28,'(a)') '      nicknm='''//trim(adjustl(nicknm))//''''
      write(28,'(a)') '      prog_vrsn='''//trim(prog_vrsn)//''''

      close(28)

      end
