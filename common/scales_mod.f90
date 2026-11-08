!#####################################################################
! Renormalization/factorization scale variation grid.
!
! Single source of truth for the (muR,muF) variation set. Every
! consumer -- the weight calculation, the histogram columns, and the
! result printout -- indexes the same table, so there is exactly one
! place to change if the grid is ever extended.
!
! Ported near-verbatim from pp2VH_mg5_gg/src/core/scales_mod.f90 (branch
! nlo, HEAD a240659): this module has no process-specific content, so
! there is nothing here for onshell/decay to differ on. Do not reorder
! the grid -- tools/ggzh_analysis-v5.1's MUR_FACTOR/MUF_FACTOR tables and
! the row order of pp2VH_mg5_gg/ref_data/LO_total.txt both assume this
! exact order and this exact icentral.
!#####################################################################

module scales_mod

  use iso_fortran_env, only: real64
  implicit none
  private
  public :: nscale, icentral
  public :: mur_factors, muf_factors
  public :: scale_label, scale_pair_label

  ! Standard 7-point variation: both scales varied by {1/2,1,2} with the
  ! two combinations differing by a factor 4 dropped.
  integer, parameter :: nscale = 7

  real(real64), parameter :: mur_factors(nscale) = &
       (/ 0.5d0, 0.5d0, 1.0d0, 1.0d0, 1.0d0, 2.0d0, 2.0d0 /)
  real(real64), parameter :: muf_factors(nscale) = &
       (/ 0.5d0, 1.0d0, 0.5d0, 1.0d0, 2.0d0, 1.0d0, 2.0d0 /)

  ! Index of the central (1,1) choice.
  integer, parameter :: icentral = 4

  ! Human/machine readable summary of the grid, used in output headers.
  character(len=*), parameter :: scale_label = &
       '(0.5,0.5) (0.5,1) (1,0.5) (1,1) (1,2) (2,1) (2,2)'

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function scale_pair_label(is) result(lbl)
  integer, intent(in) :: is
  character(len=16) :: lbl
  write(lbl,'(A,F4.2,A,F4.2,A)') '(',mur_factors(is),',',muf_factors(is),')'
  end function scale_pair_label

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module scales_mod

!#####################################################################
