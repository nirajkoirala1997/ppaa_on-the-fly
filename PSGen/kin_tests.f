      if (dlips > 0) then
       if (imom == 1) then
        if (im_test == 1) then
         call mass_test
        end if
        if (ip_test == 1) then
         call momenta_test
        end if
       end if
      end if

      contains
       subroutine mass_test

       integer :: i
       real(kind(1.d0)) :: dm

       do i=1,nexp
        dm=sqrt(p(0,i)**2-p(1,i)**2-p(2,i)**2-p(3,i)**2)
        if (abs(dm-m(i)) > dm_test) then
   	 write(*,'(a,i6)') 'Mass test for kinematics:',ikin
	 write(*,'(a,i3,a,e24.16,e24.15)') 
     &              ' ptcl:',i,':',dm,m(i)
        end if
       end do

       end subroutine mass_test

       subroutine momenta_test

       integer :: i
       real(kind(1.d0)), dimension(0:3) :: dpi,dpf

       dpi=0.d0
       dpf=0.d0

       do i=1,nexp-nfsp
        dpi=dpi+p(:,i)
       end do
       do i=nexp-nfsp+1,nexp
        dpf=dpf+p(:,i)
       end do
       if (abs(dpf(0)-dpi(0)) > dp_test) then
	write(*,'(a,i6)') 'Energy test for kinematics:',ikin
        write(*,'(3(e23.15))') dpi(0),dpf(0),dpi(0)-dpf(0)
       end if
       if (abs(dpf(1)-dpi(1))+abs(dpf(2)-dpi(2))+abs(dpf(3)-dpi(3))
     &                                                   > dp_test) then
	write(*,'(a,i6)') 'Momentum test for kinematics:',ikin
	write(*,'(a,3(e20.12))') 'mom_in:',dpi(1:3)
	write(*,'(a,3(e20.12))') 'mom_f: ',dpf(1:3)
       end if

       end subroutine momenta_test
