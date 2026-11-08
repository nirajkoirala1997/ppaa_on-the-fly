!#################################################################################!
! Real emission minus dipoles for p p -> e+ e-.
!
!   yy(1), yy(2)  beam momentum fractions, flat
!   yy(3:7)       PSGen final-state variables
!   yy(8)         PSGen channel
!
! Real and dipoles are each alpha_s^1, so scale variations replay one amplitude
! call with AL/AL_c. m(ee) is unchanged by the dipole map, but each term is still
! judged at its own kinematics, and each fills histograms at its own kinematics.
!
! Channels and which dipole subtracts them:
!   q q~ (both beam orderings): q -> q g on either leg
!   g q   : gluon on beam 1, Born q~ q        q g   : gluon on beam 2, Born q q~
!   g q~  : gluon on beam 1, Born q q~        q~ g  : gluon on beam 2, Born q~ q
!#################################################################################!
  double precision function frs(yy, vwgt)

  use run_state_mod,     only : ecm_had, xmur, xmuf, xcentral_factor, scale_dynamic
  use run_card_mod,      only : hbarc2_in, del_s, del_c
  use cuts_mod,          only : event_cut
  use olp_mod,           only : olp_real, olp_set_alphas, NFLAVOUR, NREAL_CHANNEL
  use dipole_mod,        only : dipoles_evaluate
  use ii_dipole_map_mod, only : ii_dipole_map
  use psgen_interface,   only : psgen_event, psgen_ndim, psgen_nexp, psgen_nkin
  use observables_mod,   only : obs_mass
  use analysis_mod,      only : analysis_fill
  use scales_mod,        only : nscale, icentral, mur_factors, muf_factors
  use pdfvar_mod,        only : nvar, pdfvar_central, pdfvar_select, pdfvar_restore

  implicit none

  real(8), intent(in) :: yy(16)
  real(8), intent(in) :: vwgt

  real(8) :: xa, xb, xr(psgen_ndim), plab(0:3,psgen_nexp), flux, dlips
  real(8) :: ptilde(0:3,4,2)
  real(8) :: mur_central, muf_central, mur, muf, factor_mur, factor_muf, scale_mass
  real(8) :: alphas_central, alphas_ratio, alphasPDF
  real(8) :: f1(-6:6), f2(-6:6), rsp, dot12, dot15, dot25
  real(8) :: amp(NREAL_CHANNEL)
  real(8) :: dipole_qq(NFLAVOUR,2), dipole_qqbar(NFLAVOUR,2), dipole_qbarq(NFLAVOUR,2)
  real(8) :: real_term, dipole_term(2), lum_qqbar
  real(8) :: weight(nvar), w_real(nvar), w_dip1(nvar), w_dip2(nvar)
  integer :: channel, flavour, ivar, member, member_active
  logical :: keep_real, keep_dip(2)

  external alphasPDF

    frs = 0d0

    xa = yy(1)
    xb = yy(2)
    if (xa <= 0d0 .or. xb <= 0d0 .or. xa >= 1d0 .or. xb >= 1d0) return

    xr(1:psgen_ndim) = yy(3:2+psgen_ndim)
    channel = min(1 + int(yy(3+psgen_ndim)*psgen_nkin), psgen_nkin)

    call psgen_event(ecm_had, xa, xb, channel, xr, plab, flux, dlips)
    if (dlips <= 0d0) return

! Technical cut against the numerically delicate soft and collinear corners.
    rsp   = sqrt(xa*xb)*ecm_had
    dot12 = dot4(plab(:,1), plab(:,2))
    dot15 = dot4(plab(:,1), plab(:,5))
    dot25 = dot4(plab(:,2), plab(:,5))
    if (plab(0,5) < del_s*rsp) return
    if (dot15 < del_c*dot12 .or. dot25 < del_c*dot12) return

    call ii_dipole_map(plab(:,1), plab(:,2), plab(:,3), plab(:,4), plab(:,5), ptilde)

    keep_real   = event_cut(plab)
    keep_dip(1) = event_cut(ptilde(:,:,1))
    keep_dip(2) = event_cut(ptilde(:,:,2))
    if (.not. keep_real .and. .not. any(keep_dip)) return

    if (scale_dynamic) then
       scale_mass  = obs_mass(plab(:,3) + plab(:,4))
       mur_central = xcentral_factor*scale_mass
       muf_central = mur_central
    else
       mur_central = xcentral_factor*xmur
       muf_central = xcentral_factor*xmuf
    endif

    alphas_central = alphasPDF(mur_central)
    call olp_set_alphas(alphas_central, mur_central)

    amp = 0d0
    if (keep_real) then
! With massless quarks s, b equal d and c equals u: only channels of the two
! distinct flavours in each group of five are evaluated.
       do channel = 1, NREAL_CHANNEL
          flavour = mod(channel-1, NFLAVOUR) + 1
          if (flavour <= 2) then
             call olp_real(channel, plab, amp(channel))
          else
             amp(channel) = amp(channel - 2)
          endif
       end do
    endif
    call dipoles_evaluate(alphas_central, plab, ptilde, dipole_qq,            &
                          dipole_qqbar, dipole_qbarq)

    member_active = pdfvar_central()
    do ivar = 1, nvar
       if (ivar <= nscale) then
          factor_mur = mur_factors(ivar)
          factor_muf = muf_factors(ivar)
          member     = pdfvar_central()
       else
          factor_mur = 1d0
          factor_muf = 1d0
          member     = ivar - nscale - 1
       endif
       if (member /= member_active) then
          call pdfvar_select(member)
          member_active = member
       endif

       muf = factor_muf*muf_central
       mur = factor_mur*mur_central
       call pdf(xa, muf, f1)
       call pdf(xb, muf, f2)
       alphas_ratio = alphasPDF(mur)/alphas_central

       real_term   = 0d0
       dipole_term = 0d0
       do flavour = 1, NFLAVOUR
          lum_qqbar = f1(-flavour)*f2(flavour) + f1(flavour)*f2(-flavour)
          real_term = real_term + lum_qqbar*amp(flavour)                      &
               + f1(0)*f2( flavour)*amp( 5+flavour)                           &
               + f1(flavour)*f2(0) *amp(10+flavour)                           &
               + f1(0)*f2(-flavour)*amp(15+flavour)                           &
               + f1(-flavour)*f2(0)*amp(20+flavour)

          dipole_term(1) = dipole_term(1) + lum_qqbar*dipole_qq(flavour,1)    &
               + f1(0)*f2( flavour)*dipole_qbarq(flavour,1)                   &
               + f1(0)*f2(-flavour)*dipole_qqbar(flavour,1)
          dipole_term(2) = dipole_term(2) + lum_qqbar*dipole_qq(flavour,2)    &
               + f1(flavour)*f2(0) *dipole_qqbar(flavour,2)                   &
               + f1(-flavour)*f2(0)*dipole_qbarq(flavour,2)
       end do

       w_real(ivar) = 0d0
       w_dip1(ivar) = 0d0
       w_dip2(ivar) = 0d0
       if (keep_real)   w_real(ivar) =  alphas_ratio*real_term
       if (keep_dip(1)) w_dip1(ivar) = -alphas_ratio*dipole_term(1)
       if (keep_dip(2)) w_dip2(ivar) = -alphas_ratio*dipole_term(2)
    end do
    if (member_active /= pdfvar_central()) call pdfvar_restore()

    w_real = hbarc2_in/flux*w_real*dlips
    w_dip1 = hbarc2_in/flux*w_dip1*dlips
    w_dip2 = hbarc2_in/flux*w_dip2*dlips
    weight = w_real + w_dip1 + w_dip2

    if (any(weight /= weight)) return
    frs = weight(icentral)

    if (keep_real)   call analysis_fill(plab, psgen_nexp, w_real*vwgt)
    if (keep_dip(1)) call analysis_fill(ptilde(:,:,1), 4, w_dip1*vwgt)
    if (keep_dip(2)) call analysis_fill(ptilde(:,:,2), 4, w_dip2*vwgt)

  contains

    real(8) function dot4(p, q) result(d)
    real(8), intent(in) :: p(0:3), q(0:3)
    d = p(0)*q(0) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
    end function dot4

  end function frs
!#################################################################################!
