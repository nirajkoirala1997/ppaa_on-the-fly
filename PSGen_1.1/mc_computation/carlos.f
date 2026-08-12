      subroutine carlos(f,ndim,nkin,xmc,po,xmx,nmx,ncall,int,var,sigma)
c----------------------------------------------------------------------c
c The Monte Carlo integration routine
c
c Input: f - integrand,
c        ndim - No. of dimensions,
c        nkin - No. of kinematical channels,
c        xmc - weights for the channel selection,
c        po - a part of the original weights xmc that is to be taken
c             for the next integration, po should be in [0.d0,0.5.d0],
c        ncall - max. number of calls.
c Output:
c        int - integral,
c        var -variance,
c        sigma = sqrt(variance/ncall) - standard deviation.
c
c Karol Kolodziej, University of Silesia, kolodzie@us.edu.pl
c 21 Feb, 2002; modiffied on 
c 6 March, 2002; 20 March, 2002;
c 8 May, 2008: addaptaion of integration weights implemented and
c              error estimate changed;
c 10 May, 2008: old error estimate restored;
c 4 Sept. 2008;
c 24 March 2011: error estimation for distribution bins added 
c 19 Nov., 2020: a criterion introduced which allows to discard a call 
c to the integrand if its value exceeds the maximum found for a given 
c kinematical channel by factor 'mxfctr' specified in carlocom.f
c----------------------------------------------------------------------c

      use drivekk
      use distribs

      implicit none

      integer, parameter :: dp=kind(1.d0)
      integer :: ndim,nkin,ncall,nmx(nkin)
      real(dp) :: f,xmc(0:nkin),xi(nkin),int,sigma,xm0(0:nkin),po,
     &            xmx(nkin)
      integer :: i,j
      real(dp) :: var,a,a2
      real :: x(ndim),rand(1)

      int=0.d0
      var=0.d0
      if (iscan == 0) then
       xi=0.d0
       xm0=xmc
      end if

      do i=1,ncall

c Choose the kinematics:

       if (iscan == 0) then
        if (irand == 1) then
         call random_number(rand)
        else if (irand == 2) then
         call ranlux(rand,1)
        end if
c        call random_number(rand)
c        call ranlux(rand,1)

        do j=1,nkin
c         if (dble(rand) < xmc(j)) then
         if (dble(rand(1)) < xmc(j)) then
          ikin=j
          exit
         end if
        end do
       end if

c Generate random argument:

       if (irand == 1) then
        call random_number(x)
       else if (irand == 2) then
        call ranlux(x,ndim)
       end if
c       call random_number(x)
c       call ranlux(x,ndim)

       a=f(dble(x),ndim)/ncall

       if (iscan == 1 .and. a > xmx(ikin)) then 
        xmx(ikin)=a
       end if
       if (iscan == 0) then ! 4 Sept. 2008
        if (a > mxfctr*xmx(j)) then
         nmx(j)=nmx(j)+1
         a=0.d0
        else if (a > xmx(j)) then
         xmx(j)=a
        end if
        xi(j)=xi(j)+a
       end if                ! 4 Sept. 2008

       a2=a*a
       int=int+a
       var=var+a2

c Distributions:

       if (idis == 1 .and. iscan == 0) then
        do j=1,ndist
         if (nbs(j) > 0) then
          if (ibin(j) > 0) then
           di(ibin(j),j)=di(ibin(j),j)+a
           vdi(ibin(j),j)=vdi(ibin(j),j)+a2
          end if
         else
          exit
         end if
        end do
       end if

      end do


      if (iscan == 0) then
       xi=xi/int
      end if
ckk>
c      int=int/ncall
c      var=var/ncall
ckk<
      var=ncall*var-int*int
c      sigma=sqrt(abs(var)/ncall)
ckk>
      var=var/ncall
      sigma=sqrt(abs(var))
ckk<

c Redefine weights for the next integration:


      if (iscan == 0) then
       xmc(0)=0.d0
       do j=1,nkin
        xmc(j)=xmc(j-1)+xi(j)
       end do
       xmc=po*xm0+(1.d0-po)*xmc
      end if

      end

