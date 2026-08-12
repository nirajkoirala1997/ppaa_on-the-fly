!c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~[Misc functions]
      double precision function PK(PKtype,Partontype,fxntype,x)
      implicit none

      character*10 PKtype,Partontype,fxntype
      double precision one_minus_x,x,temp,log1mxbyx,log1mx
      double precision  CA,CF,NF,TR,PI

      parameter(PI=3.141592653589793238462643D0)
      parameter(CA=3d0)
      parameter(CF = 4.0d0/3.0d0)
      parameter(NF = 5d0)
      parameter(TR = 0.5d0)

      one_minus_x = (1.0d0-x)
      log1mxbyx = dLog((1d0-x)/x)
      log1mx = dLog(1d0 - x)

!*************************************************
      if (trim(Partontype) .eq. 'gg') then  !*****
!*************************************************

!...........................................c
      if (trim(PKtype) .eq. 'P') then !.....c
!...........................................c

      if (trim(fxntype) .eq. 'Plus') then

         PK = 2d0*CA*1d0/one_minus_x

      elseif (trim(fxntype) .eq. 'Regular') then

         PK = 2d0*CA*(one_minus_x/x -1d0 + x*one_minus_x)

      elseif (trim(fxntype) .eq. 'Delta') then

         PK = 11d0/6d0*CA - 2d0/3d0*NF*TR 
      else
      print*,"Invalid function type P,Kb or Kt"
      stop
      endif ! function type for Pab ends here

!.................................................c
        elseif (trim(PKtype) .eq. 'Kb') then!.....c
!.................................................c

      if (trim(fxntype) .eq. 'Plus') then

         PK = 2d0*CA*1d0/one_minus_x*log1mxbyx

      elseif (trim(fxntype) .eq. 'Regular') then

         PK = 2d0*CA*(one_minus_x/x -1d0 + x*one_minus_x)*log1mxbyx

      elseif (trim(fxntype) .eq. 'Delta') then

         PK = -((50d0/9d0-Pi**2)*CA - 16d0/9d0*TR*NF)

      else
      print*,"Invalid function type P,Kb or Kt"
      stop
      endif ! function type for Kbar_ab ends here

!.................................................c
        elseif (trim(PKtype) .eq. 'Kt') then!.....c
!.................................................c

      if (trim(fxntype) .eq. 'Plus') then

         PK = CA*2d0/one_minus_x*log1mx

      elseif (trim(fxntype) .eq. 'Regular') then

         PK = 2d0*CA*(one_minus_x/x -1d0 + x*one_minus_x)*log1mx 

      elseif (trim(fxntype) .eq. 'Delta') then

         PK = -CA*PI*PI/3d0

      else
      print*,"Invalid function type P,Kb or Kt"
      stop
      endif ! function type for Ktilda_ab ends here

        else 
      Print*,"Invalid PKtype P,Kbar or Ktilde"
      stop
      endif  !Here PKtype ends

!****************************************************
      elseif (trim(Partontype) .eq. 'qq') then !*****
!****************************************************
!...........................................c
      if (trim(PKtype) .eq. 'P') then
!...........................................c

      if (trim(fxntype) .eq. 'Plus') then

         PK = CF*(1.0d0+x*x)/(1.0d0-x)

      elseif (trim(fxntype) .eq. 'Regular') then

         PK = -CF*(1.0d0+x)

      elseif (trim(fxntype) .eq. 'Delta') then

         PK = 0d0 
      else
      print*,"Invalid function type P,Kb or Kt"
      stop
      endif ! function type for Pab ends here

!...........................................c
        elseif (trim(PKtype) .eq. 'Kb') then
!...........................................c

      if (trim(fxntype) .eq. 'Plus') then

        PK= CF*2.0d0/one_minus_x*log1mxbyx

      elseif (trim(fxntype) .eq. 'Regular') then

         PK = CF*(-(1.0d0+x)*log1mxbyx + one_minus_x)

      elseif (trim(fxntype) .eq. 'Delta') then

         PK = -CF*(5.0d0-pi*pi)

      else
      print*,"Invalid function type P,Kb or Kt"
      stop
      endif ! function type for Kbar_ab ends here

!...........................................c
        elseif (trim(PKtype) .eq. 'Kt') then
!...........................................c

      if (trim(fxntype) .eq. 'Plus') then

         PK = CF*2.0d0/one_minus_x*log1mx

      elseif (trim(fxntype) .eq. 'Regular') then

         PK = -CF*(1.0d0+x)*log1mx

      elseif (trim(fxntype) .eq. 'Delta') then

         PK =  -CF*pi*pi/3.0d0

      else
      print*,"Invalid function type P,Kb or Kt"
      stop
      endif ! function type for Ktilda_ab ends here

        else 
      Print*,"Invalid PKtype P,Kbar or Ktilde"
      stop
      endif  !Here PKtype ends

!****************************************************
      elseif (trim(Partontype) .eq. 'qg') then !*****
!****************************************************
!...........................................c
      if (trim(PKtype) .eq. 'P') then
!...........................................c

      if (trim(fxntype) .eq. 'Plus') then

         PK = 0d0 

      elseif (trim(fxntype) .eq. 'Regular') then

         PK = CF*(1d0+ one_minus_x**2)/x 

      elseif (trim(fxntype) .eq. 'Delta') then

         PK = 0d0
      else
      print*,"Invalid function type P,Kb or Kt"
      stop
      endif ! function type for Pab ends here

!...........................................c
        elseif (trim(PKtype) .eq. 'Kb') then
!...........................................c

      if (trim(fxntype) .eq. 'Plus') then

        PK= 0d0 

      elseif (trim(fxntype) .eq. 'Regular') then

         PK = CF*(1d0+ one_minus_x**2)*log1mxbyx/x + CF*x

      elseif (trim(fxntype) .eq. 'Delta') then

         PK = 0d0 

      else
      print*,"Invalid function type P,Kb or Kt"
      stop
      endif ! function type for Kbar_ab ends here

!...........................................c
        elseif (trim(PKtype) .eq. 'Kt') then
!...........................................c

      if (trim(fxntype) .eq. 'Plus') then

         PK = 0d0 

      elseif (trim(fxntype) .eq. 'Regular') then

         PK = CF*(1d0+ one_minus_x**2)*log1mx/x

      elseif (trim(fxntype) .eq. 'Delta') then

         PK = 0d0 

      else
      print*,"Invalid function type P,Kb or Kt"
      stop
      endif ! function type for Ktilda_ab ends here

        else 
      Print*,"Invalid PKtype P,Kbar or Ktilde"
      stop
      endif  !Here PKtype ends

!****************************************************
      elseif (trim(Partontype) .eq. 'gq') then !*****
!****************************************************
!...........................................c
      if (trim(PKtype) .eq. 'P') then
!...........................................c

      if (trim(fxntype) .eq. 'Plus') then

         PK = 0d0 

      elseif (trim(fxntype) .eq. 'Regular') then

         PK = TR*(x**2 + one_minus_x**2) 

      elseif (trim(fxntype) .eq. 'Delta') then

         PK = 0d0
      else
      print*,"Invalid function type P,Kb or Kt"
      stop
      endif ! function type for Pab ends here

!...........................................c
        elseif (trim(PKtype) .eq. 'Kb') then
!...........................................c

      if (trim(fxntype) .eq. 'Plus') then

        PK= 0d0 

      elseif (trim(fxntype) .eq. 'Regular') then

        PK = TR*(x**2 + one_minus_x**2)*log1mxbyx + TR*2d0*x*one_minus_x

      elseif (trim(fxntype) .eq. 'Delta') then

         PK = 0d0

      else
      print*,"Invalid function type P,Kb or Kt"
      stop
      endif ! function type for Kbar_ab ends here

!...........................................c
        elseif (trim(PKtype) .eq. 'Kt') then
!...........................................c

      if (trim(fxntype) .eq. 'Plus') then

         PK = 0d0 

      elseif (trim(fxntype) .eq. 'Regular') then

         PK = TR*(x**2 + one_minus_x**2)*log1mx

      elseif (trim(fxntype) .eq. 'Delta') then

         PK = 0d0 

      else
      print*,"Invalid function type P,Kb or Kt"
      stop
      endif ! function type for Ktilda_ab ends here

        else 
      Print*,"Invalid PKtype P,Kbar or Ktilde"
      stop
      endif  !Here PKtype ends

      else
      print*,"Invalid Partontype gg/qq/qg/gq"
      stop
      endif ! main channel condition ends here

      return
      end
!_________________________________________________________________c 
      double precision function xint_h(x,PKtype)
      implicit none
      double precision x, dilog
      double precision CA,CF,PI,hbarc2
      external dilog
      character*10 PKtype
      parameter (PI=3.14159265358979d0)
      parameter(CA = 3d0)
      parameter(CF = 4d0/3d0)
!... here only plus function is taken, the colour factor needs to be 
!... considered carefully in the integrand.f

!... For quark quark splitting
      if ( trim(PKtype) .eq. 'Pqq') then
      xint_h =CF*(3d0-2d0*x-x**2-4d0*dLog(1d0-x))/2

      elseif (trim(PKtype) .eq. 'Kbarqq') then
      xint_h = -CF*2d0*(-0.5d0*(dLog(-1d0 + x)*(-2d0*dLog(-1d0 + 1d0/x)
     .       +dLog(-1d0+x)-2d0*dLog(x)))+dilog(1d0-x))

      elseif (trim(PKtype) .eq. 'Ktqq') then
      xint_h = -CF*dLog(1d0-x)**2

!... For gluon gluon splitting
      elseif (trim(PKtype) .eq. 'Pgg') then
      xint_h = -2d0*CA*dLog(1d0 - x) 

      elseif (trim(PKtype) .eq. 'Kbargg') then
      xint_h = -0.5d0*PI**2 -0.5d0*dLog(1d0-x)**2-dilog(1-x)  
      xint_h = 2d0*CA*xint_h

      elseif (trim(PKtype) .eq. 'Ktgg') then
      xint_h = -dLog(1d0-x)**2 
      xint_h = CA*xint_h

      else
      print*,'Invalid xintPK type'
      stop
      endif      
      end
!############################################################# 
       double precision function dilog(x)
       implicit double precision (a-z)
       parameter (pi6=1.644934066848226d+00)
       parameter (een=1.d+00)
       parameter (vier=0.25d+00)
       parameter (b2=+2.7777777777777778D-02)
       parameter (b3=-2.7777777777777778D-04)
       parameter (b4=+4.7241118669690098D-06)
       parameter (b5=-9.1857730746619641D-08)
       parameter (b6=+1.8978869988971001D-09)
       parameter (b7=-4.0647616451442256D-11)
       parameter (b8=+8.9216910204564523D-13)
1      if(x.lt.0.d0)go to 3
       if(x.gt.0.5d0)go to 4
       z=-dlog(1.d0-x)
7      z2=z*z
       dilog=z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b8+b7)+b6)
     1 +b5)+b4)+b3)+b2)+een)-z2*vier
       if(x.gt.een)dilog=-dilog-.5d0*u*u+2.d0*pi6
       return
3      if(x.gt.-een)go to 5
       y=een/(een-x)
       z=-dlog(een-y)
       z2=z*z
       u=dlog(y)
       dilog=z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b8+b7)+b6)
     1 +b5)+b4)+b3)+b2)+een)-z2*vier-u*(z+.5d0*u)-pi6
       return
4      if(x.ge.een)go to 10
       y=een-x
       z=-dlog(x)
6      u=dlog(y)
       z2=z*z
       dilog=-z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b8+b7)+b6)
     1 +b5)+b4)+b3)+b2)+een-u)+z2*vier+pi6
       if(x.gt.een)dilog=-dilog-.5d0*z*z+pi6*2.d0
       return
5      y=een/(een-x)
       z=-dlog(y)
       z2=z*z
       dilog=-z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b8+b7)+b6)
     1 +b5)+b4)+b3)+b2)+een)-z2*vier
       return
10     if(x.eq.een)go to 20
       xx=1.d0/x
       if(x.gt.2.d0)go to 11
       z=dlog(x)
       y=1.d0-xx
       go to 6
11     u=dlog(x)
       z=-dlog(1.d0-xx)
       go to 7
20     dilog=pi6
       return
       end
