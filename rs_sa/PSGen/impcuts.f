      subroutine impcuts(p,istf)

      use inprms_ps
      use kincuts
      use prtinds

      implicit none

      integer :: istf,i,j,ifrst=1
      real(dp) :: p(0:3,nexp),pij(0:3),emisst,q2pp
      real(dp), dimension(nexp-nfsp+1:nexp) :: pst,pp,phi,eta,rap

      save ifrst

c      if (ifrst == 1) then
c       call define_cuts
c       call set_indices
cc       call cutinfo(6)
c       ifrst=0
c      end if

      istf=1
      if (ictlb == 1) then
       do i=1,nxl
        if (abs(p(3,ixl(i))/sqrt(p(1,ixl(i))**2+p(2,ixl(i))**2
     &                                 +p(3,ixl(i))**2)) > ctlb) then
         istf=0
         return
        end if
       end do
      end if
      if (ictqb == 1) then
       do i=nexp-nfsp+1,nexp
        if (iyq(i) == i .or. iyb(i) == i) then
         if (abs(p(3,i)/sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2)) > ctqb) 
     &                                                           then
          istf=0
          return
         end if
        end if
       end do
      end if
      if (ictgb == 1) then
       do i=1,nxg
        if (abs(p(3,ixg(i))/sqrt(p(1,ixg(i))**2+p(2,ixg(i))**2
     &                               +p(3,ixg(i))**2)) > ctgb) then
         istf=0
         return
        end if
       end do
      end if
      if (ictab == 1) then
       do i=1,nxa
        if (abs(p(3,ixa(i))/sqrt(p(1,ixa(i))**2+p(2,ixa(i))**2
     &                               +p(3,ixa(i))**2)) > ctab) then
         istf=0
         return
        end if
       end do
      end if
      if (ictsb == 1) then
       do i=1,nxs
        if (abs(p(3,ixs(i))/sqrt(p(1,ixs(i))**2+p(2,ixs(i))**2
     &                               +p(3,ixs(i))**2)) > ctsb) then
         istf=0
         return
        end if
       end do
      end if
      if (ictpb == 1) then
       do i=1,nxp
        if (abs(p(3,ixp(i))/sqrt(p(1,ixp(i))**2+p(2,ixp(i))**2
     &                               +p(3,ixp(i))**2)) > ctpb) then
         istf=0
         return
        end if
       end do
      end if
      if (ictal == 1) then
       do i=1,nxa
        do j=1,nxl
         if ((p(1,ixa(i))*p(1,ixl(j))+p(2,ixa(i))*p(2,ixl(j))
     &       +p(3,ixa(i))*p(3,ixl(j)))/sqrt((p(1,ixa(i))**2
     &       +p(2,ixa(i))**2+p(3,ixa(i))**2)*(p(1,ixl(j))**2
     &       +p(2,ixl(j))**2+p(3,ixl(j))**2)) > ctal) then
          istf=0
          return
         end if
        end do
       end do
      end if
      if (ictaq == 1) then
       do i=1,nxa
        do j=nexp-nfsp+1,nexp
         if (iyq(j) == j .or. iyb(j) == j) then
          if ((p(1,ixa(i))*p(1,j)+p(2,ixa(i))*p(2,j)+p(3,ixa(i))*p(3,j))
     &        /sqrt((p(1,ixa(i))**2+p(2,ixa(i))**2+p(3,ixa(i))**2)
     &        *(p(1,j)**2+p(2,j)**2+p(3,j)**2)) > ctaq) then
           istf=0
           return
          end if
         end if
        end do
       end do
      end if
      if (ictas == 1) then
       do i=1,nxa
        do j=1,nxs
         if ((p(1,ixa(i))*p(1,ixs(j))+p(2,ixa(i))*p(2,ixs(j))
     &       +p(3,ixa(i))*p(3,ixs(j)))/sqrt((p(1,ixa(i))*p(1,ixa(i))
     &       +p(2,ixa(i))*p(2,ixa(i))+p(3,ixa(i))*p(3,ixa(i)))
     &      *(p(1,ixs(j))*p(1,ixs(j))+p(2,ixs(j))*p(2,ixs(j))
     &       +p(3,ixs(j))*p(3,ixs(j)))) > ctas) then
          istf=0
          return
         end if
        end do
       end do
      end if
      if (ictqq == 1) then
       do i=nexp-nfsp+1,nexp-1
        if (iyq(i) == i .or. iyb(i) == i) then
         do j=i+1,nexp
          if (iyq(j) == j .or. iyb(j) == j) then
           if ((p(1,i)*p(1,j)+p(2,i)*p(2,j)+p(3,i)*p(3,j))/sqrt(
     &         (p(1,i)*p(1,i)+p(2,i)*p(2,i)+p(3,i)*p(3,i))
     &        *(p(1,j)*p(1,j)+p(2,j)*p(2,j)+p(3,j)*p(3,j))) > ctqq)
     &                                                           then
            istf=0
            return
           end if
          end if
         end do
        end if
       end do
      end if
      if (ictlq == 1) then
       do i=1,nxl
        do j=1,nxq
         if ((p(1,ixl(i))*p(1,ixq(j))+p(2,ixl(i))*p(2,ixq(j))
     &       +p(3,ixl(i))*p(3,ixq(j)))/sqrt((p(1,ixl(i))*p(1,ixl(i))
     &       +p(2,ixl(i))*p(2,ixl(i))+p(3,ixl(i))*p(3,ixl(i)))
     &      *(p(1,ixq(j))*p(1,ixq(j))+p(2,ixq(j))*p(2,ixq(j))
     &       +p(3,ixq(j))*p(3,ixq(j)))) > ctlq) then
          istf=0
          return
         end if
        end do
       end do
      end if
      if (ictlbq == 1) then
       do i=1,nxl
        do j=1,nxb
         if ((p(1,ixl(i))*p(1,ixb(j))+p(2,ixl(i))*p(2,ixb(j))
     &       +p(3,ixl(i))*p(3,ixb(j)))/sqrt((p(1,ixl(i))*p(1,ixl(i))
     &       +p(2,ixl(i))*p(2,ixl(i))+p(3,ixl(i))*p(3,ixl(i)))
     &      *(p(1,ixb(j))*p(1,ixb(j))+p(2,ixb(j))*p(2,ixb(j))
     &       +p(3,ixb(j))*p(3,ixb(j)))) > ctlbq) then
          istf=0
          return
         end if
        end do
       end do
      end if
      if (ictlg == 1) then
       do i=1,nxl
        do j=1,nxg
         if ((p(1,ixl(i))*p(1,ixg(j))+p(2,ixl(i))*p(2,ixg(j))
     &       +p(3,ixl(i))*p(3,ixg(j)))/sqrt((p(1,ixl(i))*p(1,ixl(i))
     &       +p(2,ixl(i))*p(2,ixl(i))+p(3,ixl(i))*p(3,ixl(i)))
     &      *(p(1,ixg(j))*p(1,ixg(j))+p(2,ixg(j))*p(2,ixg(j))
     &       +p(3,ixg(j))*p(3,ixg(j)))) > ctlg) then
          istf=0
          return
         end if
        end do
       end do
      end if
      if (ictls == 1) then
       do i=1,nxl
        do j=1,nxs
         if ((p(1,ixl(i))*p(1,ixs(j))+p(2,ixl(i))*p(2,ixs(j))
     &       +p(3,ixl(i))*p(3,ixs(j)))/sqrt((p(1,ixl(i))*p(1,ixl(i))
     &       +p(2,ixl(i))*p(2,ixl(i))+p(3,ixl(i))*p(3,ixl(i)))
     &      *(p(1,ixs(j))*p(1,ixs(j))+p(2,ixs(j))*p(2,ixs(j))
     &       +p(3,ixs(j))*p(3,ixs(j)))) > ctls) then
          istf=0
          return
         end if
        end do
       end do
      end if
      if (ictll == 1) then
       do i=1,nxl-1
        do j=i+1,nxl
         if ((p(1,ixl(i))*p(1,ixl(j))+p(2,ixl(i))*p(2,ixl(j))
     &       +p(3,ixl(i))*p(3,ixl(j)))/sqrt((p(1,ixl(i))*p(1,ixl(i))
     &       +p(2,ixl(i))*p(2,ixl(i))+p(3,ixl(i))*p(3,ixl(i)))
     &       *(p(1,ixl(j))*p(1,ixl(j))+p(2,ixl(j))*p(2,ixl(j))
     &        +p(3,ixl(j))*p(3,ixl(j)))) > ctll) then
          istf=0
          return
         end if
        end do
       end do
      end if
      if (ienl == 1) then
       do i=1,nxl
        if (p(0,ixl(i)) < ecutl) then
         istf=0
         return
        end if
       end do
      end if
      if (iemt == 1) then
       emisst=0.d0
       do i=1,nxn
        emisst=emisst+sqrt(p(1,ixn(i))**2+p(2,ixn(i))**2)
        if (emisst < ecutt) then
         istf=0
         return
        end if
       end do
      end if
      if (ienq == 1) then
       do i=1,nxq
        if (p(0,ixq(i)) < ecutq) then
         istf=0
         return
        end if
       end do
      end if
      if (ienb == 1) then
       do i=1,nxb
        if (p(0,ixb(i)) < ecutb) then
         istf=0
         return
        end if
       end do
      end if  ! 31 July, 2009
      if (ieng == 1) then
       do i=1,nxg
        if (p(0,ixg(i)) < ecutg) then
         istf=0
         return
        end if
       end do
      end if
      if (iena == 1) then
       do i=1,nxa
        if (p(0,ixa(i)) < ecuta) then
         istf=0
         return
        end if
       end do
      end if
      if (imqq == 1) then
       do i=nexp-nfsp+1,nexp-1
        if (iyq(i) == i .or. iyb(i) == i) then
         do j=i+1,nexp
          if (iyq(j) == j .or. iyb(j) == j) then
           pij=p(:,i)+p(:,j)
           if ((pij(0)-pij(3))*(pij(0)+pij(3))-pij(1)*pij(1)
     &         -pij(2)*pij(2) < mqq*mqq) then
            istf=0
            return
           end if
          end if
         end do
        end if
       end do
      end if
      if (imbb == 1) then
       do i=nexp-nfsp+1,nexp-1
        if (iyb(i) == i) then
         do j=i+1,nexp
          if (iyb(j) == j) then
           pij=p(:,i)+p(:,j)
           if ((pij(0)-pij(3))*(pij(0)+pij(3))-pij(1)*pij(1)
     &         -pij(2)*pij(2) < mbb*mbb) then
            istf=0
            return
           end if
          end if
         end do
        end if
       end do
      end if
      if (impp == 1) then
       do i=nexp-nfsp+1,nexp-1
        if (iyp(i) == i) then
         do j=i+1,nexp
          if (iyp(j) == j) then
           pij=p(:,i)+p(:,j)
           q2pp=(pij(0)-pij(3))*(pij(0)+pij(3))-pij(1)*pij(1)
     &         -pij(2)*pij(2) 
           if (q2pp < mpp1 .or. mpp2 < q2pp) then
            istf=0
            return
           end if
          end if
         end do
        end if
       end do
      end if
      if (imjj == 1) then
       do i=nexp-nfsp+1,nexp-1
        if (iyj(i) == i .or. iyb(i) == i) then
         do j=i+1,nexp
          if (iyj(j) == j .or. iyb(j) == j) then
           pij=p(:,i)+p(:,j)
           if ((pij(0)-pij(3))*(pij(0)+pij(3))-pij(1)*pij(1)
     &         -pij(2)*pij(2) < mjj*mjj) then
            istf=0
            return
           end if
          end if
         end do
        end if
       end do
      end if
      if (imll == 1) then
       do i=nexp-nfsp+1,nexp
        if (iyl(i) == i) then
         do j=i+1,nexp
          if (iyl(j) == j) then
           pij=p(:,i)+p(:,j)
           if ((pij(0)-pij(3))*(pij(0)+pij(3))-pij(1)*pij(1)
     &         -pij(2)*pij(2) < mll*mll) then
            istf=0
            return
           end if
          end if
         end do
        end if
       end do
      end if
      if (iptl == 1 .or. iptj == 1 .or. irapl == 1 .or. irapj == 1 
     &   .or. ietal == 1 .or. ietaj == 1 .or. idelr == 1
     &   .or. iptb == 1 .or. irapb == 1 .or. ietab == 1) then
c Rapidities, pseudorapidities, azimuthal angles:
       do i=nexp-nfsp+1,nexp
        pst(i)=sqrt(p(1,i)**2+p(2,i)**2)
        pp(i)=sqrt(pst(i)**2+p(3,i)**2)
        if (pst(i) > 0.d0) then
         eta(i)=0.5d0*log((pp(i)+p(3,i))/(pp(i)-p(3,i)))
         rap(i)=0.5d0*log((p(0,i)+p(3,i))/(p(0,i)-p(3,i)))
         if (p(2,i)/pst(i) >= 0.d0) then
          phi(i)=acos(p(1,i)/pst(i))
         else
          phi(i)=twopi-acos(p(1,i)/pst(i))
         end if
        else
         if (p(3,i) < 0.d0) then
          phi(i)=pi
          eta(i)=-0.5d0*log(huge(1.d0))
         else
          phi(i)=0.d0
          eta(i)=0.5d0*log(huge(1.d0))
         end if
         rap(i)=eta(i)
        end if
        if (iptl == 1 .and. iyl(i) == i) then
         if (abs(pst(i)) < ptl) then
          istf=0
          return
         end if
        end if
        if (iptj == 1 .and. iyj(i) == i) then
         if (abs(pst(i)) < ptj) then
          istf=0
          return
         end if
        end if
        if (iptb == 1 .and. iyb(i) == i) then
         if (abs(pst(i)) < ptb) then
          istf=0
          return
         end if
        end if
        if (ietal == 1 .and. iyl(i) == i) then
         if (abs(eta(i)) > etal) then
          istf=0
          return
         end if
        end if
        if (ietaj == 1 .and. iyj(i) == i) then
         if (abs(eta(i)) > etaj) then
          istf=0
          return
         end if
        end if
        if (ietab == 1 .and. iyb(i) == i) then
         if (abs(eta(i)) > etab) then
          istf=0
          return
         end if
        end if
        if (irapl == 1 .and. iyl(i) == i) then
         if (abs(rap(i)) > rapl) then
          istf=0
          return
         end if
        end if
        if (irapj == 1 .and. iyj(i) == i) then
         if (abs(rap(i)) > rapj) then
          istf=0
          return
         end if
        end if
        if (irapb == 1 .and. iyb(i) == i) then
         if (abs(rap(i)) > rapb) then
          istf=0
          return
         end if
        end if
       end do
      end if
      if (idelr == 1) then
       do i=nexp-nfsp+1,nexp
        if (iyl(i) == i .or. iyj(i) == i .or. iyb(i) == i) then
         do j=i+1,nexp
          if ((iyl(j) == j .or. iyj(j) == j .or. iyb(j) == j) .and.
     &          sqrt((eta(i)-eta(j))**2+(phi(i)-phi(j))**2) < delr) then
           istf=0
           return
          end if
         end do
        end if
       end do
      end if

      end subroutine impcuts
