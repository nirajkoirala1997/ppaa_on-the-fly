      subroutine psgen(ikin,ecm,x1,x2,x,ndim,flux,dlips)
c----------------------------------------------------------------------c
c Calculates the particle four momenta in the centre of mass system, 
c delivered in module fourmom, initial flux factor 'flux' and phase 
c space normalization 'dlips' of the process
c
c                 f_1 f_2 -> f_3 f_4 ... f_nexp,
c
c where
c ikin is kinematics selected out of 'nkin' kinematical routines 
c        generated,
c ecm is the centre of mass energy, which is reduced to sqrt(x1*x2)*ecm if 
c x1,x2, the initial particle momentum fractions are smaller then 1, 
c        i.e. 0 < x1,x2 < 1,
c x(ndim), with ndim=3*(nexp-2)-4, are random numbers in (0,1)
c
c Karol Kolodziej, University of Silesia, karol.kolodziej@us.edu.pl
c 20 October, 2022; modified 
c 8 March, 2023: call to param_trnas shifted to main program, cps 
c defined at the statge of compilation
c----------------------------------------------------------------------c

      use kinparams
      use kincuts
      use fourmom
      use psnrmtest

      implicit none
 
      integer :: ikin,ndim
      real(kind(1.d0)) :: ecm,x1,x2,x(ndim),flux,dlips

      integer :: imom,i,istf
      real(kind(1.d0)) :: ecm0=0.d0,cps=(twopi)**(4-3*nfsp),x1x2,ecmp,
     &                    dlinv,vz,gz,r(0:3,nexp),p_dp(0:3,nexp)
      real(kind(1.q0)) :: ecm_qp

      save ecm0,cps

      x1x2=x1*x2

      if (x1x2 < 1.d0) then
       ecmp=sqrt(x1x2)*ecm
      else
       ecmp=ecm
      end if

      if (ecmp /= ecm0) then
       if (iqprec == 1) then
        ecm_qp=ecmp
        call kinin_qp(ecm_qp,m(1:2),p(:,1:2),flux)
       else
        call kinin(ecmp,m(1:2),p(:,1:2),flux)
       end if
       ecm0=ecmp
      end if

      call kincls(ikin,1,x,ndim,dlinv)

      if (dlinv > zero) then
       if (icuts == 1) then
        if (x1x2 < 1.d0) then
c Boost momenta along the z axis back to the lab system to impose cuts
         vz=(x1-x2)/(x1+x2)
         gz=0.5d0*(x1+x2)/(ecmp/ecm)
         do i=1,nexp
          if (iqprec == 1) then
           r_qp(:,i)=(/gz*(p(0,i)+vz*p(3,i)),p(1,i),p(2,i),
     &                               gz*(p(3,i)+vz*p(0,i))/)
          else
           r_dp(:,i)=(/gz*(p(0,i)+vz*p(3,i)),p(1,i),p(2,i),
     &                             gz*(p(3,i)+vz*p(0,i))/)
          end if
         end do
         if (iqprec == 1) then
          r_dp=r_qp
         end if
         call impcuts(r_dp,istf)
        else
         p_dp=p
         call impcuts(p_dp,istf)
        end if
        if (istf == 0) then
         dlips=0.d0
         return
        end if
       end if

       call psncls(ikin,x,ndim,dlinv)

       if (dlinv > zero) then
        dlips=cps/dlinv
       else
        dlips=zero
       end if
       if (dlinvt > zero) then
        dlipst=cps/dlinvt
       else
        dlipst=zero
       end if
      end if




      end
