      subroutine genps(nfsp)
c----------------------------------------------------------------------c
c Generates a multichannel Monte Carlo integration routine for 
c the process:
c
c                 1 + 2 -> 3 + 4 + ... + n
c
c Karol Kolodziej, University of Silesia, karol.kolodziej@us.edu.pl
c
c 11 May, 2021 
c modified October, 2021
c 10 Feb. 2022: generation of mapping patterns in calls to s-channel
c kinematic routines written anew 
c 14-15 Feb. 2022: generation of calls to single tpole kinematics
c written anew
c 18 Feb. 2022: generation of calls to single tpole kinematics improved
c 6 April 2022: a bug in writing variable cinds corrected
c 29 June, 2022: writing of tchcalls.f corrected
c----------------------------------------------------------------------c

      use psgen_info

      implicit none

      integer :: i,j,k,l,m,n,nfsp,ivp(mxks),icolon,nvks,np,
     &           iks,ixp(nfsp),iv,iyp(mxks,nfsp),ip,
     &           icrit,itcrit1,itcrit2,
     &           izp(nfsp),ipl,ipr,ifpx,itp(nfsp),
     &           iup(nfsp),isp(nfsp),ia(mxks),ig(mxks),
     &           iao(mxks),igo(mxks),ia1,ig1,ia2,ig2,ikg,itarg,ikt,
     &           inarg,ivd(mxks),ifpk,ifpl,ifpn,ict1(mxks),ict2(mxks)

      character(78) :: cline
      character(80) :: ckin0,ckin1,ckin2,ckin3,ind,ind1,ind2,ind3,cinds,
     &                 cnp,cnp1,cnp2,cnp3,cindo,ckingi,ckingf,
     &                 ckino1,ckino2,ckino3,tcalli,tcallj,targ0,targ,
     &                 targf,tcallf,indl,indk,indn,ckint,cnt,
     &                 ckint1,cnt1,ckint2,cnt2,ckint3,cnt3,
     &                 targ1,targ11,targ12,targ13,
     &                 targ1f,targ1f1,targ1f2,targ1f3,
     &                 targ2,targ21,targ22,targ23,
     &                 targ2f,targ2f1,targ2f2,targ2f3,ctry(mxks),indg

      character(40), dimension(mxks) :: cvrt,cgam,ck,ci,co,ct10,ct1f,
     &                                  ct1k,ct20,ct2f,ct2k

      character(160) :: cls,ckin,ckino,ckins,cking


      write(ind,'(i4)') nfsp
      ckin0='kin'//trim(adjustl(ind))//'_'
      write(indg,'(i4)') nfsp-1
      ckingi='king'//trim(adjustl(indg))//'isr'//'('
      ckingf='king'//trim(adjustl(indg))//'fsr'//'('
      write(ind,'(i4)') nfsp
      tcalli='       call kint'//trim(adjustl(ind))//'_'
      tcallf='     &                                  '//
     &                                     'imom,ecm,p,x,ndim,dlips)'


      open(unit=7,file='genps.dat',status='old',action='read')
      i=0
      do 
       read(7,'(i1,a78)',end=700) np,cline
       if (np /= 0) then
        i=i+1
        icolon=index(cline,':')
        ivp(i)=np
        cvrt(i)=cline(1:icolon-1)
        cgam(i)=cline(icolon+1:)
       end if
      end do
 700  nvks=i
      close(7)

      iv=0
      iyp=0
      ict1=0
      ict2=0
      ct10=''
      ct1f=''
      ct1k=''
      ct20=''
      ct2f=''
      ct2k=''

      do i=1,nvks
       do k=nispt+1,nptcl
        ifpk=index(cvrt(i),' '//trim(ptcl(k))//' ')
        if (ifpk > 0) then
         do l=k+1,nptcl
          ifpl=index(cvrt(i),' '//trim(ptcl(l))//' ')
          if (ifpl > 0 .and. ifpk /= ifpl) then
           write(indk,'(i4)') k
           write(indl,'(i4)') l
           cinds=trim(adjustl(indk))//','//trim(adjustl(indl))//','
           cindo=trim(adjustl(indl))//','//trim(adjustl(indk))//','
           if (ivp(i) == 2) then
            iv=iv+1
            ivd(iv)=ivp(i)
            write(ind,'(i4)') ivp(i)
            ck(iv)=trim(adjustl(ind))
            ci(iv)=trim(cinds)//trim(adjustl(cgam(i)))//','
            co(iv)=trim(cindo)//trim(adjustl(cgam(i)))//','
            ia(iv)=ipa(k)
            ig(iv)=ipg(k)
            iao(iv)=ipa(l)
            igo(iv)=ipg(l)
            iyp(iv,k-nispt)=k
            iyp(iv,l-nispt)=l
            do m=1,itch1
             if (k == ixtch1(m)) then
              ct10(iv)='(1,'//trim(adjustl(indk))//',zero,ctc,'
              ct1f(iv)=trim(adjustl(indl))//','
              ict1(iv)=1
              exit
             else if (l == ixtch1(m)) then
              ct10(iv)='(1,'//trim(adjustl(indl))//',zero,ctc,'
              ct1f(iv)=trim(adjustl(indk))//','
              ict1(iv)=1
              exit
             end if
            end do
            do m=1,itch2
             if (k == ixtch2(m)) then
              ct20(iv)='(2,'//trim(adjustl(indk))//',zero,ctc,'
              ct2f(iv)=trim(adjustl(indl))//','
              ict2(iv)=1
              exit
             else if (l == ixtch2(m)) then
              ct20(iv)='(2,'//trim(adjustl(indl))//',zero,ctc,'
              ct2f(iv)=trim(adjustl(indk))//','
              ict2(iv)=1
              exit
             end if
            end do
           else
            do n=l+1,nptcl
             ifpn=index(cvrt(i),' '//trim(ptcl(n))//' ')
             if (ifpn > 0 .and. ifpn /= ifpk .and. ifpn /= ifpl) then
              write(indn,'(i4)') n
c              cinds=trim(cinds)//','//trim(adjustl(indn))//','
              cinds=trim(cinds)//trim(adjustl(indn))//','   ! 6 April 2022
              cindo=trim(cindo)//trim(adjustl(indn))//','
              if (ivp(i) == 3) then
               iv=iv+1
               ivd(iv)=ivp(i)
               write(ind,'(i4)') ivp(i)
               ck(iv)=trim(adjustl(ind))
               ci(iv)=trim(cinds)//trim(adjustl(cgam(i)))//','
               co(iv)=trim(cindo)//trim(adjustl(cgam(i)))//','
               ia(iv)=ipa(k)
               ig(iv)=ipg(k)
               iao(iv)=ipa(l)
               igo(iv)=ipg(l)
               iyp(iv,k-nispt)=k
               iyp(iv,l-nispt)=l
               iyp(iv,n-nispt)=n
               do m=1,itch1
                if (k == ixtch1(m)) then
                 ct10(iv)='(1,'//trim(adjustl(indk))//',zero,ctc,'
                 ct1f(iv)=trim(adjustl(indl))//','//trim(adjustl(indn))
     &                  //','!mw,zero'
                 ct1k(iv)='2'
                 ict1(iv)=1
                 exit
                else if (l == ixtch1(m)) then
                 ct10(iv)='(1,'//trim(adjustl(indl))//',zero,ctc,'
                 ct1f(iv)=trim(adjustl(indk))//','//trim(adjustl(indn))
     &                  //','!mw,zero'
                 ct1k(iv)='2'
                 ict1(iv)=1
                 exit
                else if (n == ixtch1(m)) then
                 ct10(iv)='(1,'//trim(adjustl(indn))//',zero,ctc,'
                 ct1f(iv)=trim(adjustl(indk))//','//trim(adjustl(indl))
     &                  //','!mw,zero'
                 ct1k(iv)='2'
                 ict1(iv)=1
                 exit
                end if
               end do
               do m=1,itch2
                if (k == ixtch2(m)) then
                 ct20(iv)='(2,'//trim(adjustl(indk))//',zero,ctc,'
                 ct2f(iv)=trim(adjustl(indl))//','//trim(adjustl(indn))
     &                  //','!mw,zero'
                 ct2k(iv)='2'
                 ict2(iv)=1
                 exit
                else if (l == ixtch2(m)) then
                 ct20(iv)='(2,'//trim(adjustl(indl))//',zero,ctc,'
                 ct2f(iv)=trim(adjustl(indk))//','//trim(adjustl(indn))
     &                  //','!mw,zero'
                 ct2k(iv)='2'
                 ict2(iv)=1
                 exit
                else if (n == ixtch2(m)) then
                 ct20(iv)='(2,'//trim(adjustl(indn))//',zero,ctc,'
                 ct2f(iv)=trim(adjustl(indk))//','//trim(adjustl(indl))
     &                  //','!mw,zero'
                 ct2k(iv)='2'
                 ict2(iv)=1
                 exit
                end if
               end do
              end if
             end if
            end do
           end if
          end if
         end do
        end if
       end do
      end do

c      write(*,*) 'ixtch1:', ixtch1(:itch1)
c      write(*,*) 'ixtch2:', ixtch2(:itch2)

c      do i=1,iv
c       write(*,*) i, iyp(i,:)
c       write(*,*) i,ia(i),ci(i)
c       write(*,*) i,iao(i),co(i)
c       write(*,*) i, ck(i)
c       if (ict1(i) == 1) then
c        write(*,*) 't:'
c        write(*,*) i, ct10(i)
c        write(*,*) i, ct1f(i)
c        write(*,*) i, ct1k(i)
c       end if
c       if (ict2(i) == 1) then
c        write(*,*) 't:'
c        write(*,*) i, ct20(i)
c        write(*,*) i, ct2f(i)
c        write(*,*) i, ct2k(i)
c       end if
c      end do




      open(unit=19,file='schcalls.txt',status='replace',action='write')
      open(unit=20,file='schcalls2.txt',status='replace',action='write')
      open(unit=191,file='gchcalls.txt',status='replace',action='write')
      open(unit=201,file='gchcalls2.txt',status='replace',
     &                                                   action='write')

      iks=0
      ikg=0
      ikt=0
      ckin=''
      ckino=''
      cnp=''
      ckint='     &'
      cnt=''
      targ1=''
      targ1f=''
      targ2=''
      targ2f=''
      ctry=''



      if (nfsp == 2) then
       iks=iks+1
       ckins=trim(ckin0)//'2(3,4,'
       write(19,'(a)') ckins(:80)
       write(20,'(a)') ckins(81:)
      else if (nfsp == 3 .and. igind2 > nispt) then
       iks=iks+1
       write(ind,'(i4)') igind2
       do i=1,iv
        ckins=trim(ckin0)//trim(ck(i))//'('//trim(adjustl(ci(i)))
     &        //trim(adjustl(ind))//','
        write(19,'(a)') ckins(:80)
        write(20,'(a)') ckins(81:)
        ckin=trim(adjustl(ci(i)))
        ckino=trim(adjustl(co(i)))
        ipl=0
        ipr=0
        do ip=1,nfsp
         if (isp(ip) > 0) then
          ipl=ipl+1
         else if (isp(ip) == 0) then
          ipr=nispt+ip
         end if
        end do
        call photon_radiation
       end do
      else
       do i=1,iv
        do j=i+1,iv
         icrit=0
         itcrit1=0
         itcrit2=0
         do ip=1,nfsp
          icrit=icrit+iyp(i,ip)*iyp(j,ip)
         end do
         if (icrit == 0) then
          izp=iyp(i,:)+iyp(j,:)
          ckin=trim(ckin)//trim(adjustl(ci(i)))//trim(adjustl(ci(j)))
          ckino=trim(ckino)//trim(adjustl(co(i)))//trim(adjustl(ci(j)))
          cnp=trim(cnp)//trim(adjustl(ck(i)))//trim(adjustl(ck(j)))
          if (ict1(i) == 0 .and. ict2(i) == 0 .or. itcrit1 == 1
     &                                        .or. itcrit2 == 1) then
           ckint=trim(ckint)//trim(adjustl(ci(i)))
           cnt=trim(cnt)//trim(adjustl(ck(i)))
          else
           if (ict1(i) > 0) then
            targ1=trim(targ1)//ct10(i)
            targ1f=trim(targ1f)//trim(ct1f(i))
            itcrit1=1
           end if
           if (ict2(i) > 0) then
            targ2=trim(targ2)//ct20(i)
            targ2f=trim(targ2f)//trim(ct2f(i))
            itcrit2=1
           end if
          end if
          if (ict1(j) == 0 .and. ict2(j) == 0 .or. itcrit1 == 1
     &                                        .or. itcrit2 == 1) then
           ckint=trim(ckint)//trim(adjustl(ci(j)))
           cnt=trim(cnt)//trim(adjustl(ck(j)))
          else 
           if (ict1(j) > 0) then
            targ1=trim(targ1)//ct10(j)
            targ1f=trim(targ1f)//trim(ct1f(j))
            itcrit1=1
           end if
           if (ict2(j) > 0) then
            targ2=trim(targ2)//ct20(j)
            targ2f=trim(targ2f)//trim(ct2f(j))
            itcrit2=1
           end if
          end if
          ipl=0
          ipr=0
          do ip=1,nfsp
           if (izp(ip) > 0) then
            ipl=ipl+1
           else if (izp(ip) == 0) then
            ipr=nispt+ip
           end if
          end do
          if (nfsp-ipl <= 1) then
           if (ipr > 0) then
            call photon_radiation
            write(ind,'(i4)') ipr
            ckin=trim(ckin)//trim(adjustl(ind))//','
            ckint=trim(ckint)//trim(adjustl(ind))//','
           end if
           iks=iks+1
           ckins=trim(ckin0)//trim(cnp)//'('//trim(ckin)
           write(19,'(a)') ckins(:80)
           write(20,'(a)') ckins(81:)
           ckin=''
           ckino=''
           cnp=''
           call single_tpole
           ckint='     &'
           cnt=''
           targ1=''
           targ1f=''
           targ2=''
           targ2f=''
          else
           ckin1=ckin
           ckino1=ckino
           cnp1=cnp
           ckint1=ckint
           cnt1=cnt
           targ11=targ1
           targ1f1=targ1f
           targ21=targ2
           targ2f1=targ2f
           do k=j+1,iv
            icrit=0
            do ip=1,nfsp
             icrit=icrit+izp(ip)*iyp(k,ip)
            end do
            if (icrit == 0) then
             itp=izp(:)+iyp(k,:)
             ckin=trim(ckin)//trim(adjustl(ci(k)))
             ckino=trim(ckino)//trim(adjustl(ci(k)))  ! 16 Feb. 2023
             cnp=trim(cnp)//trim(adjustl(ck(k)))
             if (ict1(k) == 0 .and. ict2(k) == 0 .or. itcrit1 == 1
     &                                           .or. itcrit2 == 1) then
              ckint=trim(ckint)//trim(adjustl(ci(k)))
              cnt=trim(cnt)//trim(adjustl(ck(k)))
             else 
              if (ict1(k) > 0) then
               targ1=trim(targ1)//ct10(k)
               targ1f=trim(targ1f)//trim(ct1f(k))
               itcrit1=1
              end if
              if (ict2(k) > 0) then
               targ2=trim(targ2)//ct20(k)
               targ2f=trim(targ2f)//trim(ct2f(k))
               itcrit2=1
              end if
             end if
             ipl=0
             ipr=0
             do ip=1,nfsp
              if (itp(ip) > 0) then
               ipl=ipl+1
              else if (itp(ip) == 0) then
               ipr=nispt+ip
              end if
             end do
             if (nfsp-ipl <= 1) then
              iks=iks+1
              ckin=trim(ckin0)//trim(cnp)//'('//trim(ckin)
              if (ipr > 0) then
               call photon_radiation
               write(ind,'(i4)') ipr
               ckin=trim(ckin)//trim(adjustl(ind))//','
               ckint=trim(ckint)//trim(adjustl(ind))//','
              end if
              write(19,'(a)') ckin(:80)
              write(20,'(a)') ckin(81:)
              ckin=trim(ckin1)
              ckino=trim(ckino1)
              cnp=trim(cnp1)
              call single_tpole
              ckint=trim(ckint1)
              cnt=trim(cnt1)
              targ1=trim(targ11)
              targ1f=trim(targ1f1)
              targ2=trim(targ21)
              targ2f=trim(targ2f1)
             else
              ckin2=ckin
              ckino2=ckino
              cnp2=cnp
              ckint2=ckint
              cnt2=cnt
              targ12=targ1
              targ1f2=targ1f
              targ22=targ2
              targ2f2=targ2f
              do l=k+1,iv
               icrit=0
               do ip=1,nfsp
                icrit=icrit+itp(ip)*iyp(l,ip)
               end do
               if (icrit == 0) then
                iup=itp(:)+iyp(l,:)
                ckin=trim(ckin)//trim(adjustl(ci(l)))
                ckino=trim(ckino)//trim(adjustl(ci(l)))   ! 16 Feb. 2023
                cnp=trim(cnp)//trim(adjustl(ck(l)))
                if (ict1(l) == 0 .and. ict2(l) == 0 .or. itcrit1 == 1
     &                                          .or. itcrit2 == 1) then
                 ckint=trim(ckint)//trim(adjustl(ci(l)))
                 cnt=trim(cnt)//trim(adjustl(ck(l)))
                else
                 if (ict1(l) > 0) then
                  targ1=trim(targ1)//ct10(l)
                  targ1f=trim(targ1f)//trim(ct1f(l))
                  itcrit1=1
                 end if
                 if (ict2(l) > 0) then
                  targ2=trim(targ2)//ct20(l)
                  targ2f=trim(targ2f)//trim(ct2f(l))
                  itcrit2=1
                 end if
                end if
                ipl=0
                ipr=0
                do ip=1,nfsp
                 if (iup(ip) > 0) then
                  ipl=ipl+1
                 else if (iup(ip) == 0) then
                  ipr=nispt+ip
                 end if
                end do
                if (nfsp-ipl <= 1) then
                 iks=iks+1
                 ckin=trim(ckin0)//trim(cnp)//'('//trim(ckin)
                 if (ipr > 0) then
                  call photon_radiation
                  write(ind,'(i4)') ipr
                  ckin=trim(ckin)//trim(adjustl(ind))//','
                  ckint=trim(ckint)//trim(adjustl(ind))//','
                 end if
                 write(19,'(a)') ckin(:80)
                 write(20,'(a)') ckin(81:)
                 ckin=trim(ckin2)
                 ckino=trim(ckino2)
                 cnp=trim(cnp2)
                 call single_tpole
                 ckint=trim(ckint2)
                 cnt=trim(cnt2)
                 targ1=trim(targ12)
                 targ1f=trim(targ1f2)
                 targ2=trim(targ22)
                 targ2f=trim(targ2f2)
                else
                 ckin3=ckin
                 ckino3=ckino
                 cnp3=cnp
                 ckint3=ckint
                 cnt3=cnt
                 targ13=targ1
                 targ1f3=targ1f
                 targ23=targ2
                 targ2f3=targ2f
                 do m=l+1,iv
                  icrit=0
                  do ip=1,nfsp
                   icrit=icrit+iup(ip)*iyp(m,ip)
                  end do
                  if (icrit == 0) then
                   isp=iup(:)+iyp(m,:)
                   ckin=trim(ckin)//trim(adjustl(ci(m)))
                   ckino=trim(ckino)//trim(adjustl(ci(m)))  ! 16 Feb. 2023
                   cnp=trim(cnp)//trim(adjustl(ck(m)))
                   if (ict1(m) == 0 .and. ict2(m) == 0 .or. itcrit1 == 1
     &                                           .or. itcrit2 == 1) then
                    ckint=trim(ckint)//trim(adjustl(ci(m)))
                    cnt=trim(cnt)//trim(adjustl(ck(m)))
                   else 
                    if (ict1(m) > 0) then
                     targ1=trim(targ1)//ct10(m)
                     targ1f=trim(targ1f)//trim(ct1f(m))
                     itcrit1=1
                    end if
                    if (ict2(m) > 0) then
                     targ2=trim(targ2)//ct20(m)
                     targ2f=trim(targ2f)//trim(ct2f(m))
                     itcrit2=1
                    end if
                   end if
                   ipl=0
                   ipr=0
                   do ip=1,nfsp
                    if (isp(ip) > 0) then
                     ipl=ipl+1
                    else if (isp(ip) == 0) then
                     ipr=nispt+ip
                    end if
                   end do
                   if (nfsp-ipl <= 1) then
                    iks=iks+1
                    ckin=trim(ckin0)//trim(cnp)//'('//trim(ckin)
                    if (ipr > 0) then
                     call photon_radiation
                     write(ind,'(i4)') ipr
                     ckin=trim(ckin)//trim(adjustl(ind))//','
                     ckint=trim(ckint)//trim(adjustl(ind))//','
                    end if
                    write(19,'(a)') ckin(:80)
                    write(20,'(a)') ckin(81:)
                    ckin=trim(ckin3)
                    ckino=trim(ckino3)
                    cnp=trim(cnp3)
                    call single_tpole
                    ckint=trim(ckint3)
                    cnt=trim(cnt3)
                    targ1=trim(targ13)
                    targ1f=trim(targ1f3)
                    targ2=trim(targ23)
                    targ2f=trim(targ2f3)
                   end if
                  end if
                 end do
                 ckin=trim(ckin2)
                 ckino=trim(ckino2)
                 cnp=trim(cnp2)
                 ckint=trim(ckint2)
                 cnt=trim(cnt2)
                 targ1=trim(targ12)
                 targ1f=trim(targ1f2)
                 targ2=trim(targ22)
                 targ2f=trim(targ2f2)
                end if
               end if
              end do
              ckin=trim(ckin1)
              ckino=trim(ckino1)
              cnp=trim(cnp1)
              ckint=trim(ckint1)
              cnt=trim(cnt1)
              targ1=trim(targ11)
              targ1f=trim(targ1f1)
              targ2=trim(targ21)
              targ2f=trim(targ2f1)
             end if
            end if
           end do
           ckin=''
           ckino=''
           cnp=''
           ckint='     &'
           cnt=''
           targ1=''
           targ1f=''
           targ2=''
           targ2f=''
          end if
         else
          ckin=''
          ckino=''
          cnp=''
          ckint='     &'
          cnt=''
          targ1=''
          targ1f=''
          targ2=''
          targ2f=''
         end if
        end do
       end do
      end if

      close(19)
      close(20)
      close(191)
      close(201)

      if (itch1 > 0 .or. itch2 > 0) then
       write(29,'(a)') '      end if'
      end if


      open(unit=19,file='schcalls.txt',status='old',action='read')
      open(unit=18,file='schcalls2.txt',status='old',action='read')
      open(unit=20,file='kinschnl.f',status='replace',action='write')
      open(unit=191,file='gchcalls.txt',status='old',action='read')
      open(unit=181,file='gchcalls2.txt',status='old',action='read')
      open(unit=201,file='kingchnl.f',status='replace',action='write')

      write(20,'(a)') '      subroutine kinschnl(ikin,imom,ecm,'
     &                                                //'x,ndim,dlips)'
      write(20,'(a)') ' '
      write(20,'(a)') '      use fourmom'
      write(20,'(a)') '      use kinparams'
      write(20,'(a)') ' '
      write(20,'(a)') '      implicit none'
      write(20,'(a)') ' '
      write(20,'(a)') '      integer :: ikin,imom,ndim'
      write(20,'(a)') '      real(dp) :: x(ndim),dlips'
      if (iquadp == 1) then
       write(20,'(a)') '      real(qp) :: ecm' ! 20 Dec. 2022
c       write(20,'(a)') '      real(qp) :: ecm,p(0:3,nexp)'
      else
       write(20,'(a)') '      real(dp) :: ecm' ! 20 Dec. 2022
c       write(20,'(a)') '      real(dp) :: ecm,p(0:3,nexp)'
      end if
      write(20,'(a)') ' '

      if (igchnl == 2) then
       write(201,'(a)') '      subroutine kingchnl(ikin,imom,ecm,'
     &                                                //'x,ndim,dlips)'
       write(201,'(a)') ' '
       write(201,'(a)') '      use fourmom'
       write(201,'(a)') '      use kinparams'
c       write(201,'(a)') '      use engymom'
       write(201,'(a)') ' '
       write(201,'(a)') '      implicit none'
       write(201,'(a)') ' '
       write(201,'(a)') '      integer :: ikin,imom,ndim'
       write(201,'(a)') '      real(dp) :: x(ndim),dlips'
       if (iquadp == 1) then
        write(201,'(a)') '      real(qp) :: ecm' ! 20 Dec. 2022
c        write(201,'(a)') '      real(qp) :: ecm,p(0:3,nexp)'
       else
        write(201,'(a)') '      real(dp) :: ecm' ! 20 Dec. 2022
c        write(201,'(a)') '      real(dp) :: ecm,p(0:3,nexp)'
       end if
       write(201,'(a)') ' '
c       write(201,'(a)') '      ecm=ee'
c       write(201,'(a)') '      p=pp'
       write(201,'(a)') ' '
      else
       write(201,'(a)') '      subroutine kingchnl'
       write(201,'(a)') 'c This subroutine is empty'
       write(201,'(a)') '      end'
      end if
       iks=0

      do 
       read(19,'(a)',end=190) ckin(:80)
       read(18,'(a)',end=190) ckin(81:)
       iks=iks+1
       write(ind,'(i10)') iks
       if (iks == 1) then
        write(20,'(a)') '      if (ikin == '
     &                             //trim(adjustl(ind))//') then'
       else
        write(20,'(a)') '      else if (ikin == '
     &                             //trim(adjustl(ind))//') then'
       end if
       cls='      call '//trim(ckin)
       call wrt72(20,' '//trim(cls))
       write(20,'(a)') '     &                             '//
     &                                    'imom,ecm,p,x,ndim,dlips)'
      end do

 190  close(19)
      close(18)

      nkins2=iks

      write(20,'(a)') '      end if'
      write(20,'(a)') ' '
      write(20,'(a)') '      end'
      
      close(20)

      if (igchnl == 2) then
       ikg=0
       do 
        read(191,'(a)',end=195) ckin(:80)
        read(181,'(a)',end=195) ckin(81:)
        ikg=ikg+1
        write(ind,'(i10)') ikg
        if (ikg == 1) then
         write(201,'(a)') '      if (ikin == nkinst+'
     &                              //trim(adjustl(ind))//') then'
        else
         write(201,'(a)') '      else if (ikin == nkinst+'
     &                              //trim(adjustl(ind))//') then'
        end if
        cls='      call '//trim(ckin)
        call wrt72(201,' '//trim(cls))
        write(201,'(a)') '     &                             '//
     &                                    'imom,ecm,x,ndim,dlips)'
       end do
 
 195   close(191)
       close(181)

       nking2=ikg

       write(201,'(a)') '      end if'
       write(201,'(a)') ' '
c       write(201,'(a)') '      if (imom == 1) pp=p'
c       write(201,'(a)') ' '
       write(201,'(a)') '      end'
      
       close(201)
      end if

c      if (ikg > 0) then
       open(unit=287,file='engymom.f',status='replace',action='write')
       write(287,'(a)') '      module engymom'
       write(287,'(a)') ' '
       write(287,'(a)') '      use inprms_ps'
c       write(287,'(a)') '      use partcls'
       write(287,'(a)') ' '
       if (iquadp == 1) then
        write(287,'(a)') '      real(qp) :: '    ! 21 Dec. 2022
c        write(287,'(a)') '      real(qp) :: ecm,p(0:3,nepl),'    ! 4 April 2022
c        write(287,'(a)') '      real(qp) :: ecm,p(0:3,0:nepl),' ! 8 Feb. 2022
       else
        write(287,'(a)') '      real(dp) :: '    ! 21 Dec. 2022
c        write(287,'(a)') '      real(dp) :: ecm,p(0:3,nepl),'    ! 4 April 2022
c        write(287,'(a)') '      real(dp) :: ecm,p(0:3,0:nepl),' ! 8 Feb. 2022
       end if
       write(287,'(a)') '     &            msum,reg,sp,sqsp,g,v(3),'
       write(287,'(a)') '     &            beta,ybta,dctg,'
       write(287,'(a)') '     &            q1(0:3),q2(0:3),q3(0:3),'//
     &                                   'q4(0:3),'
       write(287,'(a)') '     &            q23(0:3),q12(0:3),q34(0:3),'
     &                                 //'q234(0:3),'
       write(287,'(a)') '     &            sqsl1,sqsl2,sqsl3,sqsl4,'
       write(287,'(a)') '     &            sqsl23,sqsl12,sqsl34,'//
     &                                    'sqsl234,'
       write(287,'(a)') '     &            sqsu1,sqsu2,sqsu3,sqsu4,'
       write(287,'(a)') '     &            sqsu23,sqsu12,sqsu34,'//
     &                                    'sqsu234,'
       write(287,'(a)') '     &            s1,s2,s3,s4,'
       write(287,'(a)') '     &            s23,s12,s34,s234,'
       write(287,'(a)') '     &            sqs1,sqs2,sqs3,sqs4,'
       write(287,'(a)') '     &            sqs23,sqs12,sqs34,sqs234,'
       write(287,'(a)') '     &            qq1,b34g,e34g,v34g(3),g34g,'
       write(287,'(a)') '     &            e3,e7'
       write(287,'(a)') ' '
       write(287,'(a)') '      end module engymom'
       close(287)      
c      end if

      call wrkinprms(nfsp)

      write(*,*) 'Code has been generated. Go to ../mc_computation '
     &           //'to test it.'

      open(unit=19,file='kinkk1.f',status='replace',action='write')
      write(19,'(a)')'      subroutine kinkk1'
      write(19,'(a)')'c This subroutine is empty'
      write(19,'(a)')'      end subroutine kinkk1'
      close(19)

      contains
       subroutine photon_radiation
        integer :: in
        character(160) :: ckinga

        in=index(ckin,'(')                       ! 16 Feb. 2023>
        if (in > 0) then
         ckinga=trim(ckin(in+1:))
        else
         ckinga=trim(ckin)
        end if                                   ! 16 Feb. 2023<
        if (ipa(1) == 1 .and. ipa(2) == 1) then
         ikg=ikg+1
c         cking=trim(ckingi)//trim(ckin)
         in=index(cnp,'3')                      ! 20 Feb. 2023>
         if (in > 0) then
          cking='king'//trim(adjustl(indg))//'_'//trim(cnp)//'('
     &          //trim(ckinga)
          write(191,'(a)') trim(cking(:80))
          write(201,'(a)') cking(81:)
          return
         else                                    ! 20 Feb. 2023
          cking=trim(ckingi)//trim(ckinga)   ! 16 Feb. 2023
          write(191,'(a)') trim(cking(:80))
          write(201,'(a)') cking(81:)
         end if
        end if
        if (trim(ptcl(ipr)) == 'a') then
         if (ia(i) == 1) then
          ikg=ikg+1
c          cking=trim(ckingf)//trim(ckin)
          cking=trim(ckingf)//trim(ckinga)   ! 16 Feb. 2023
          write(191,'(a)') trim(cking(:80))
          write(201,'(a)') cking(81:)
         end if
         if (iao(i) == 1) then
          ikg=ikg+1
          cking=trim(ckingf)//trim(ckino)
          write(191,'(a)') trim(cking(:80))
          write(201,'(a)') cking(81:)
         end if
        end if
       end subroutine photon_radiation

      subroutine single_tpole

       integer :: it,itwrt1,itwrt2,itnowr
       character(80) :: ctarg1,ctarg2

       itwrt1=0
       itwrt2=0
       itnowr=0

       if (nfsp >= 8 .and. (trim(cnt) == '2' .or. trim(cnt) == '22')
     &     .or. (len_trim(targ1) == 0 .and. len_trim(targ2) == 0)
     &     .or. (nfsp == 6 .or. nfsp == 7) .and. trim(cnt) == '2'
     &     .or. nfsp == 9 .and. trim(cnt) == '32')
     &                               then  ! 29 June, 2022; 15 Sept, 2022
        itnowr=1
        return
       end if

       if (itcrit1 == 1) then 
        ctarg1=trim(ckint)//trim(targ1f)
        do it=1,ikt
         if (trim(ctarg1) == trim(ctry(it))) then
          itwrt1=1
          exit
         end if
        end do
        if (itwrt1 == 0) then
         ikt=ikt+1
         ctry(ikt)=ctarg1
        end if
       end if
       if (itcrit2 == 1) then
        ctarg2=trim(ckint)//trim(targ2f)
        do it=1,ikt
         if (trim(ctarg2) == trim(ctry(it))) then
          itwrt2=1
          exit
         end if
        end do
        if (itwrt2 == 0) then
         ikt=ikt+1
         ctry(ikt)=ctarg2
        end if
       end if

       if (itcrit1 == 1 .and. itwrt1 == 0) then 
        nkint2=nkint2+1
        write(ind3,'(i4)') nkint2
        if (nkint2 == 1) then
         write(29,'(a)') '      if (ikin == nkins+'
     &               //trim(adjustl(ind3))//') then'
        else
         write(29,'(a)') '      else if (ikin == nkins+'
     &               //trim(adjustl(ind3))//') then'
        end if
        write(29,'(a)') trim(tcalli)//trim(cnt)//trim(targ1)
        write(29,'(a)') trim(ctarg1)
        write(29,'(a)') trim(tcallf)
       end if
       if (itcrit2 == 1 .and. itwrt2 == 0) then
        nkint2=nkint2+1
        write(ind3,'(i4)') nkint2
        if (nkint2 == 1) then
         write(29,'(a)') '      if (ikin == nkins+'
     &               //trim(adjustl(ind3))//') then'
        else
         write(29,'(a)') '      else if (ikin == nkins+'
     &               //trim(adjustl(ind3))//') then'
        end if
        write(29,'(a)') trim(tcalli)//trim(cnt)//trim(targ2)
        write(29,'(a)') trim(ctarg2)
        write(29,'(a)') trim(tcallf)
       end if

       itcrit1=0
       itcrit2=0

      end subroutine single_tpole

      end
