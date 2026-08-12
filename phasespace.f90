!##################################################################################!
module kinematics_mod
  use iso_fortran_env, only: real64
  use globals_mod, only : s
  use misc_utils_mod
  use constants_mod
  use misc_utils_mod
  implicit none
  public :: kinvar3,kinvar2,kinvar2_PK,reduceps
contains

!##################################################################################!
  subroutine kinvar2(xx, xxinvmass, p1, p2, p3, p4) !##############################!
!##################################################################################!
    implicit none
    real(real64), intent(in)  :: xx(10)
    real(real64), intent(out) :: xxinvmass
    real(real64), intent(out) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3)
    real(real64) :: m1,m2,m3,m4

    real(real64) :: xa, xb, v, omv, srs2,s34
    m1 = mass(1)
    m2 = mass(1)
    m3 = mass(1)
    m4 = mass(1)

    ! Extract variables
    xa = xx(1)
    xb = xx(2)
    v  = xx(3)
    omv = 1.0_real64 - v

    srs2 = 0.5_real64 * sqrt(s)

    ! Incoming parton 4-vectors
    p1(0) = srs2 * xa
    p1(1) = 0.0_real64
    p1(2) = 0.0_real64
    p1(3) = p1(0)

    p2(0) = srs2 * xb
    p2(1) = 0.0_real64
    p2(2) = 0.0_real64
    p2(3) = -p2(0)

    ! Outgoing parton 4-vectors
    p3(0) = srs2 * (xa * v + xb * omv)
    p3(1) = sqrt(s * xa * xb * v * omv)
    p3(2) = 0.0_real64
    p3(3) = srs2 * (xa * v - xb * omv)

     p4 = p1 + p2 - p3

    ! Invariant mass of p3 + p4 system
    s34 = 2.0_real64 * dot(p3, p4)
    xxinvmass = sqrt(s34)

  end subroutine kinvar2
!====================================================!
   subroutine kinvar2_PK(xa,xb,xc,Qmass,p1,p2,p3,p4) !
!====================================================!
  implicit none 
  real(real64), intent(in)  :: xa, xb, xc 
  real(real64), intent(out) :: Qmass 
  real(real64), intent(out) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3)
  real(real64) :: m1,m2,m3,m4

  real(real64) :: v, omv, srs2
  m1 = mass(1)
  m2 = mass(1)
  m3 = mass(1)
  m4 = mass(1)

   v = xc
   omv=1d0-v 
   srs2=0.5*dsqrt(s)

!incoming parton 4-vectors
   p1(0)=srs2*xa
   p1(1)=0d0
   p1(2)=0d0
   p1(3)=p1(0)

   p2(0)=srs2*xb
   p2(1)=0d0
   p2(2)=0d0
   p2(3)=-p2(0)

!outgoing parton 4-vectors
   p3(0)=srs2*(xa*v+xb*omv)
   p3(1)=dsqrt(s*xa*xb*v*omv)
   p3(2)=0d0
   p3(3)=srs2*(xa*v-xb*omv)

    p4 = p1 + p2 - p3 

    Qmass = dsqrt(2.0d0*dot(p3,p4))

   end subroutine kinvar2_PK

!##################################################################################!
  subroutine kinvar3(xx, xxjac, xinvmass, p1, p2, p3, p4, p5, unphy) !#############!
!##################################################################################!
    implicit none

    real(real64), intent(in)  :: xx(10)
    real(real64), intent(out) :: xxjac, xinvmass
    real(real64), intent(out) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3)
    integer,      intent(out) :: unphy

    ! Local variables
    real(real64) :: xa, xb, v, w, srs2, s34
    real(real64) :: sp, rsp, xjac
    real(real64) :: ct, st, phi, xjac3, xjac4
    real(real64) :: cphi, sphi, e5min, e5max, xjac5, e5, p5m
    real(real64) :: p5t, p5x, p5y, p5z
    real(real64) :: sigma, tau, amp, amm, a1, b1, a2
    real(real64) :: e3min, e3max, xjac6, e3, p3m
    real(real64) :: e4, p4m, czeta, zeta, szeta
    real(real64) :: p3t, p3x, p3y, p3z
    real(real64) :: p4t, p4x, p4y, p4z
    real(real64) :: beta, gamma
    real(real64) :: m1,m2,m3,m4,m5 

    m1 = mass(1)
    m2 = mass(1)
    m3 = mass(1)
    m4 = mass(1)
    m5 = mass(1)

    ! Extract input
    xa = xx(1)
    xb = xx(2)
    xjac = 1d0
    sp = xa * xb * s
    rsp = dsqrt(sp)
    srs2 = 0.5d0 * sqrt(s)

    ! Incoming parton momenta
    p1(0) = xa * srs2
    p1(1) = 0d0
    p1(2) = 0d0
    p1(3) = p1(0)

    p2(0) = xb * srs2
    p2(1) = 0d0
    p2(2) = 0d0
    p2(3) = -p2(0)

    ! Random values
    v = xx(3)
    w = xx(4)

    unphy = 0
    if (sp < (m3 + m4)**2 - m5**2) then
      unphy = 1
      return
    end if

    ! Angular variables
    ct = -1d0 + 2d0*v
    xjac3 = 2d0
    st = sqrt(1d0 - ct*ct)

    phi = 2d0*pi*w
    xjac4 = 2d0*pi
    cphi = cos(phi)
    sphi = sin(phi)

    ! Final state particle p5
    e5min = m5
    e5max = 0.5d0*(rsp + (m5**2 - (m3+m4)**2)/rsp)
    xjac5 = e5max - e5min
    e5 = xjac5*xx(5) + e5min
    p5m = sqrt(abs(e5**2 - m5**2))

    p5t = e5
    p5x = p5m * st
    p5y = 0d0
    p5z = p5m * ct

    sigma = sqrt(sp) - e5
    tau = sigma**2 - p5m**2
    amp = m3 + m4
    amm = m3 - m4
    a1 = sigma * (tau + amp*amm)
    b1 = (tau - amp**2)*(tau - amm**2)
    a2 = p5m * sqrt(b1)
    e3min = 0.5d0 / tau * (a1 - a2)
    e3max = 0.5d0 / tau * (a1 + a2)

    xjac6 = e3max - e3min
    e3 = xjac6 * xx(6) + e3min
    p3m = sqrt(abs(e3**2 - m3**2))

    e4 = rsp - e3 - e5
    p4m = sqrt(abs(e4**2 - m4**2))
    czeta = (p4m**2 - p3m**2 - p5m**2)/(2d0 * p3m * p5m)

    if (czeta >= 1.0d0) then
      unphy = unphy + 1
      return
    end if

    zeta = acos(czeta)
    szeta = sin(zeta)

    ! Final state p3
    p3t = e3
    p3x = p3m * (ct * cphi * szeta + st * czeta)
    p3y = p3m * (sphi * szeta)
    p3z = p3m * (ct * czeta - st * cphi * szeta)

    ! Final state p4
    p4t = e4
    p4x = -p3x - p5x
    p4y = -p3y - p5y
    p4z = -p3z - p5z

    ! Boost back to lab frame
    beta = (xa - xb) / (xa + xb)
    gamma = 1d0 / sqrt(1d0 - beta*beta)

    p3(0) = gamma * (p3t + beta * p3z)
    p3(1) = p3x
    p3(2) = p3y
    p3(3) = gamma * (p3z + beta * p3t)

    p5(0) = gamma * (p5t + beta * p5z)
    p5(1) = p5x
    p5(2) = p5y
    p5(3) = gamma * (p5z + beta * p5t)

    p4(0) = p1(0) + p2(0) - p3(0) - p5(0)
    p4(1) = p1(1) + p2(1) - p3(1) - p5(1)
    p4(2) = p1(2) + p2(2) - p3(2) - p5(2)
    p4(3) = p1(3) + p2(3) - p3(3) - p5(3)

    ! Jacobian
    xxjac = xjac3 * xjac4 * xjac5 * xjac6

    ! Invariant mass
    s34 = 2d0 * dot(p3, p4)
    xinvmass = sqrt(s34)
    return
  end subroutine kinvar3

  !==================================================!

  subroutine reduceps(p1, p2, p3, p4, p5, ptilde)
    implicit none
    real(real64), intent(in)  :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3)
    real(real64), intent(out) :: ptilde(0:3, 1:4, 1:2)

    ! Local variables
    real(real64) :: s12, s15, s25, x512
    real(real64) :: K(0:3), Ktil(0:3), K_plus_Kt(0:3)
    integer :: i
    integer, parameter :: iemitter = 2

    s12 = 2.0_real64 * dot(p1, p2)
    s15 = 2.0_real64 * dot(p1, p5)
    s25 = 2.0_real64 * dot(p2, p5)
    x512 = (s12 - s15 - s25) / s12

    do i = 1, iemitter
      select case (i)
      case (1)  ! emitter is p1
        ptilde(:,1,i) = x512 * p1
        ptilde(:,2,i) = p2

        K    = p1 + p2 - p5
        Ktil = ptilde(:,1,i) + p2
        K_plus_Kt = K + Ktil

        ptilde(:,3,i) = p3 - 2.0_real64 * dot(p3, K_plus_Kt) * K_plus_Kt / dot(K_plus_Kt, K_plus_Kt) &
                            + 2.0_real64 * dot(p3, K) * Ktil / dot(K, K)
        ptilde(:,4,i) = p4 - 2.0_real64 * dot(p4, K_plus_Kt) * K_plus_Kt / dot(K_plus_Kt, K_plus_Kt) &
                            + 2.0_real64 * dot(p4, K) * Ktil / dot(K, K)

      case (2)  ! emitter is p2
        ptilde(:,1,i) = p1
        ptilde(:,2,i) = x512 * p2

        K    = p1 + p2 - p5
        Ktil = ptilde(:,2,i) + p1
        K_plus_Kt = K + Ktil

        ptilde(:,3,i) = p3 - 2.0_real64 * dot(p3, K_plus_Kt) * K_plus_Kt / dot(K_plus_Kt, K_plus_Kt) &
                            + 2.0_real64 * dot(p3, K) * Ktil / dot(K, K)
        ptilde(:,4,i) = p4 - 2.0_real64 * dot(p4, K_plus_Kt) * K_plus_Kt / dot(K_plus_Kt, K_plus_Kt) &
                            + 2.0_real64 * dot(p4, K) * Ktil / dot(K, K)
      end select
    end do
  end subroutine reduceps

end module kinematics_mod





!##################################################################################!

!c---------------------------------------------------------------------
!        subroutine resetmomenta(p1,p2,p3,p4,p5)
!        implicit none
!        double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
!        integer i
!        do i=0,3
!        p1(i)=0d0
!        p2(i)=0d0
!        p3(i)=0d0
!        p4(i)=0d0
!        p5(i)=0d0
!        enddo
!        return
!        end
!c---------------------------------------------------------------------
         subroutine p1d_to_p2d_5(p1,p2,p3,p4,p5,p)
         implicit double precision (a-h,o-z)

         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p(0:3,1:5)
        do i=0,3
         p(i,1)=p1(i)
         p(i,2)=p2(i)
         p(i,3)=p3(i)
         p(i,4)=p4(i)
         p(i,5)=p5(i)
        enddo
         end
!c---------------------------------------------------------------------
!
!c---------------------------------------------------------------------
!         subroutine p2dtop1d_5(p,p1,p2,p3,p4,p5)
!         implicit double precision (a-h,o-z)
!
!         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p(0:3,1:5)
!        do i=0,3
!         p1(i)=p(i,1)
!         p2(i)=p(i,2)
!         p3(i)=p(i,3)
!         p4(i)=p(i,4)
!         p5(i)=p(i,5)
!        enddo
!         end
!c---------------------------------------------------------------------
!c---------------------------------------------------------------------
!         subroutine printmomenta(p1,p2,p3,p4)
!                 implicit double precision (a-h,o-z)
!                 dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
!                 write(*,*)"p1= ",p1
!                 write(*,*)"p2= ",p2
!                 write(*,*)"p3= ",p3
!                 write(*,*)"p4= ",p4
!         end
subroutine print_momenta(p)
  implicit none
  real(kind(1.d0)), intent(in) :: p(0:3,4)
  integer :: i

  do i = 1, 4
    print *, 'p', i, ' =',p(0,i), ',', p(1,i), ',', p(2,i), ',', p(3,i)
  end do
end subroutine print_momenta

!c---------------------------------------------------------------------
!c---------------------------------------------------------------------
!         subroutine p2d_to_p1d_4(p,p1,p2,p3,p4)
!         implicit double precision (a-h,o-z)
!
!         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
!        do i=0,3
!         p1(i)=p(i,1)
!         p2(i)=p(i,2)
!         p3(i)=p(i,3)
!         p4(i)=p(i,4)
!        enddo
!         end
!c---------------------------------------------------------------------
!c---------------------------------------------------------------------
         subroutine p1d_to_p2d_4(p1,p2,p3,p4,p)
         implicit double precision (a-h,o-z)

         dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p(0:3,1:4)
        do i=0,3
         p(i,1)=p1(i)
         p(i,2)=p2(i)
         p(i,3)=p3(i)
         p(i,4)=p4(i)
        enddo
         end
!c---------------------------------------------------------------------
!
