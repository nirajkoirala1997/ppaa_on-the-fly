!c---------------------------------------------------------------------
!c     Initial state dipole for the case of Drell- Yan.
!c---------------------------------------------------------------------
module dipoles
   use iso_fortran_env, only: real64
   use globals_mod
   use misc_utils_mod
   use constants_mod
   use openloops
   implicit none
   public :: dipole_initial_initial 
contains

   subroutine dipole_initial_initial(AL, p_ex, ptilde, dip)
      use iso_fortran_env, only: real64
      implicit none
      real(real64), intent(in) :: AL
      real(real64), intent(in) :: p_ex(0:3,5)
      real(real64), intent(in) :: ptilde(0:3,4,2)
      real(real64), intent(out) :: dip(5,2)
   
      real(real64) :: x_i_ab, m2_ew
      real(real64) :: pa(0:3), pb(0:3), p_i(0:3)
      real(real64) :: born(5)
      real(real64) :: Born_cc_tmp(4,4)
      real(real64) :: Born_cc_ai_b(4,4,5), Born_cc_bi_a(4,4,5)
      real(real64) :: p_ex_leg1(0:3,4), p_ex_leg2(0:3,4)
      integer :: i
   
      pa = p_ex(:,1)
      pb = p_ex(:,2)
      p_i = p_ex(:,5)
   
      p_ex_leg1 = ptilde(:,:,1)
      p_ex_leg2 = ptilde(:,:,2)
   
      x_i_ab = (dot(pa,pb)-dot(p_i,pa)-dot(p_i,pb))/dot(pa,pb) 
   
      do i = 1, 5
         call evaluate_ccmatrix(id_LO(i), p_ex_leg1, born(i), Born_cc_tmp, m2_ew)
         Born_cc_ai_b(:,:,i) = Born_cc_tmp
         dip(i,1) = -1.0 / (2.0 * x_i_ab * dot(pa, p_i)) * (8.0 * PI * AL * CF) &
                    * (2.0 / (1.0 - x_i_ab) - (1.0 + x_i_ab)) * Born_cc_ai_b(1,2,i) / CF
      enddo
   
      do i = 1, 5
         call evaluate_ccmatrix(id_LO(i), p_ex_leg2, born(i), Born_cc_tmp, m2_ew)
         Born_cc_bi_a(:,:,i) = Born_cc_tmp

         dip(i,2) = -1.0 / (2.0 * x_i_ab * dot(pb, p_i)) * (8.0 * PI * AL * CF) &
                    * (2.0 / (1.0 - x_i_ab) - (1.0 + x_i_ab)) * Born_cc_bi_a(1,2,i) / CF
      enddo
   
   end subroutine dipole_initial_initial
   end module dipoles 

!c---------------------------------------------------------------------
!c     Initial state dipole for the case of Drell- Yan gq channel
!
!       function dipole_gq_q(k,p)
!      implicit double precision (a-h,o-z)
!      parameter(PI=3.141592653589793238D0)
!      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),
!     .           p6(0:3),p7(0:3),p8(0:3),p9(0:3),p(0:3,1:5)
!      common/usedalpha/AL,ge
!      call p2dtop1d_5(p,p1,p2,p3,p4,p5)
!
!      s12=2.d0*dot(p1,p2) 
!      s13=2.d0*dot(p1,p3)
!      s14=2.d0*dot(p1,p4)
!      s15=2.d0*dot(p1,p5)
!      s23=2.d0*dot(p2,p3)
!      s24=2.d0*dot(p2,p4)
!      s25=2.d0*dot(p2,p5)
!      s34=2.d0*dot(p3,p4)
!      s35=2.d0*dot(p3,p5)
!      s45=2.d0*dot(p4,p5)
!
!      Tr = 0.5d0
!
!      if ( k .eq. 1 ) then ! Dipole leg 1 
!        call reducemomenta2(1,p1,p2,p3,p4,p5,p6,p7,p8,p9)
!        Born= born_uu2ee(1,p6,p7,p8,p9) 
!
!        dipole_gq_q=
!     -   (-8.*Al*Born*Pi*(s12**2 - 2.*s12*s15 + 2.*s15**2 - 2.*s12*s25 +
!     -      4.*s15*s25 + 2.*s25**2)*Tr)/(s12*s15*(s12 - s15 - s25))
!
!      else if ( k .eq. 2 ) then ! Diople leg 2
!        call reducemomenta2(2,p1,p2,p3,p4,p5,p6,p7,p8,p9)
!        Born= born_uu2ee(2,p6,p7,p8,p9) 
!
!        dipole_gq_q=
!     -   (-8.*Al*Born*Pi*(s12**2 - 2.*s12*s15 + 2.*s15**2 - 2.*s12*s25 +
!     -      4.*s15*s25 + 2.*s25**2)*Tr)/(s12*(s12 - s15 - s25)*s25)
!
!      endif  
!      return
!      end 
!c---------------------------------------------------------------------
!
!
!c---------------------------------------------------------------------c
!c       Initial State splitting with initial spectator                c
!c            REDUCED MOMENTA                                          c
!c---------------------------------------------------------------------c
!       subroutine reducemomenta2(k,p1,p2,p3,p4,p5,p15,p2til,p3til,p4til)
!       implicit double precision (a-h,o-z)
!
!       dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),
!     .   p15(0:3),p2til(0:3),p3til(0:3),p4til(0:3),ak(0:3),aktil(0:3)
!     .   ,akadd(0:3),p25(0:3),p1til(0:3),diff(0:3)
!      s12=2.d0*dot(p1,p2)
!      s13=2.d0*dot(p1,p3)
!      s14=2.d0*dot(p1,p4)
!      s15=2.d0*dot(p1,p5)
!      s23=2.d0*dot(p2,p3)
!      s24=2.d0*dot(p2,p4)
!      s25=2.d0*dot(p2,p5)
!      s34=2.d0*dot(p3,p4)
!      s35=2.d0*dot(p3,p5)
!      s45=2.d0*dot(p4,p5)
!
!        x512=(s12-s15-s25)/s12
!        if (k .eq .1) then   ! this choice is for leg1
!
!        do i=0,3
!
!        p15(i)=x512*p1(i)
!        p2til(i)=p2(i)
!        ak(i)=p1(i)+p2(i)-p5(i)
!        aktil(i)= p15(i)+p2(i)
!        akadd(i)=ak(i)+aktil(i)
!        enddo
!        do i=0,3
!        p3til(i)= p3(i)-2d0*dot(p3,akadd)*akadd(i)/dot(akadd,akadd)
!     .            + 2d0*dot(p3,ak)*aktil(i)/dot(ak,ak)
!        p4til(i)= p4(i)-2d0*dot(p4,akadd)*akadd(i)/dot(akadd,akadd)
!     .            + 2d0*dot(p4,ak)*aktil(i)/dot(ak,ak)
!        enddo
!        else if (k .eq. 2) then  ! this choice is for leg2
!
!        do i=0,3
!
!        p25(i)=x512*p2(i)
!        p1til(i)=p1(i)
!        ak(i)=p1(i)+p2(i)-p5(i)
!        aktil(i)= p25(i)+p1(i)
!        akadd(i)=ak(i)+aktil(i)
!
!        p15(i) = p1til(i) ! in output p6,p7,p8,p9 numbered accordingly
!        p2til(i)= p25(i)
!
!        enddo
!        do i=0,3
!        p3til(i)= p3(i)-2d0*dot(p3,akadd)*akadd(i)/dot(akadd,akadd)
!     .            + 2d0*dot(p3,ak)*aktil(i)/dot(ak,ak)
!        p4til(i)= p4(i)-2d0*dot(p4,akadd)*akadd(i)/dot(akadd,akadd)
!     .            + 2d0*dot(p4,ak)*aktil(i)/dot(ak,ak)
!        enddo
!        endif
!
!        end
!c---------------------------------------------------------------------
