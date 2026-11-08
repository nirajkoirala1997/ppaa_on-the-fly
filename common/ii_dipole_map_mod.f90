!#################################################################################!
! Catani-Seymour initial-initial momentum mapping -- reduced (Born-level)
! kinematics from the real kinematics, for one emitter leg.
!
! WHAT THIS IS. Pure Lorentz algebra on five momenta: given the real 5-point
! kinematics (two massless initial legs, three final legs, the last of which is
! the parton being subtracted), it builds the two-momentum-conserving,
! on-shell-preserving REDUCED kinematics `ptilde` -- one set per choice of
! which initial leg is the emitter. Nothing here evaluates an amplitude or
! knows what the three final legs physically are.
!
! WHY IT IS IN common/ AND NOT IN A STUDY DIRECTORY. Same argument
! dipole_kernels_mod.f90 and kp_kernels_mod.f90 both make: this mapping is
! process-independent Catani-Seymour machinery, not physics specific to
! g g -> Z gamma. It works identically for any process with two massless
! initial legs and exactly two "spectator" final legs plus one radiated
! parton -- which is also why its shape is fixed at exactly two spectator
! final legs (p3, p4) rather than a general list: that IS this process's
! topology (and the g g -> Z H stand-in's), and generalising it to an
! arbitrary final-state count is not a problem this build has.
!
! WHERE THE NUMBERS COME FROM. Transcribed unchanged from `reduceps` in the
! production reference's phasespace.f90
! (/mnt/cluster_storage/FO_all/mcpc2/pp2VH_mg5_gg/phasespace.f90) -- algebra
! untouched, only the name and the local `dot` made explicit rather than
! inherited from a shared misc_utils_mod. Validated indirectly through
! checks/check_dipole_ii.f90 (the mapped kinematics feed evaluate_born_cc/sc,
! and the resulting dipole values are compared against the live reference).
!
! NOT IN $(COMMON_OBJS), same reasoning as dipole_kernels_mod.f90 and
! kp_kernels_mod.f90: named explicitly in the OBJECTS list of whichever
! program actually links it.
!#################################################################################!

module ii_dipole_map_mod

  use iso_fortran_env, only : dp => real64

  implicit none
  private

  public :: ii_dipole_map

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! p1, p2   the two initial (massless) legs
! p3, p4   the two spectator final legs (Z, H for this build)
! p5       the parton being subtracted (the radiated gluon)
!
! ptilde(0:3, 1:4, 1:2)   the reduced kinematics, one set per emitter choice:
!                         ptilde(:,:,1) has p1 as the rescaled emitter, p2 the
!                         spectator; ptilde(:,:,2) has them swapped. Legs 1,2
!                         of the third index are the two initial legs of that
!                         reduced kinematics, legs 3,4 are p3, p4 boosted/
!                         reflected into it.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  subroutine ii_dipole_map(p1, p2, p3, p4, p5, ptilde)
  real(dp), intent(in)  :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3)
  real(dp), intent(out) :: ptilde(0:3, 1:4, 1:2)

  real(dp) :: s12, s15, s25, x512
  real(dp) :: K(0:3), Ktil(0:3), K_plus_Kt(0:3)
  integer  :: i
  integer, parameter :: iemitter = 2

  s12  = 2.0_dp * dot(p1, p2)
  s15  = 2.0_dp * dot(p1, p5)
  s25  = 2.0_dp * dot(p2, p5)
  x512 = (s12 - s15 - s25) / s12

  do i = 1, iemitter
     select case (i)
     case (1)  ! emitter is p1
        ptilde(:,1,i) = x512 * p1
        ptilde(:,2,i) = p2

        K    = p1 + p2 - p5
        Ktil = ptilde(:,1,i) + p2
        K_plus_Kt = K + Ktil

        ptilde(:,3,i) = p3 - 2.0_dp * dot(p3, K_plus_Kt) * K_plus_Kt / dot(K_plus_Kt, K_plus_Kt) &
                           + 2.0_dp * dot(p3, K) * Ktil / dot(K, K)
        ptilde(:,4,i) = p4 - 2.0_dp * dot(p4, K_plus_Kt) * K_plus_Kt / dot(K_plus_Kt, K_plus_Kt) &
                           + 2.0_dp * dot(p4, K) * Ktil / dot(K, K)

     case (2)  ! emitter is p2
        ptilde(:,1,i) = p1
        ptilde(:,2,i) = x512 * p2

        K    = p1 + p2 - p5
        Ktil = ptilde(:,2,i) + p1
        K_plus_Kt = K + Ktil

        ptilde(:,3,i) = p3 - 2.0_dp * dot(p3, K_plus_Kt) * K_plus_Kt / dot(K_plus_Kt, K_plus_Kt) &
                           + 2.0_dp * dot(p3, K) * Ktil / dot(K, K)
        ptilde(:,4,i) = p4 - 2.0_dp * dot(p4, K_plus_Kt) * K_plus_Kt / dot(K_plus_Kt, K_plus_Kt) &
                           + 2.0_dp * dot(p4, K) * Ktil / dot(K, K)
     end select
  end do

  end subroutine ii_dipole_map

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! Lorentz scalar product, +--- metric -- same convention every other module in
! this build uses. Kept local, same reasoning as dipole_kernels_mod.f90.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  real(dp) function dot(p, q) result(d)
  real(dp), intent(in) :: p(0:3), q(0:3)
  d = p(0)*q(0) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
  end function dot

end module ii_dipole_map_mod
!#################################################################################!
