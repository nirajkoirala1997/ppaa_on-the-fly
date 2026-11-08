!#################################################################################!
! Virtual + Catani-Seymour I-operator for p p -> e+ e-, on Born kinematics.
!
! In units of n = alpha_s(mu)/(2 pi) CF Born, with L = log(mu^2/shat), in Recola's
! normalisation (the COLI one; OpenLoops' with polenorm = 1):
!
!   virtual = -2/eps^2 - 3/eps + (A - L^2 - 3 L),   A = -8 + 4 pi^2/3 at L = 0
!   I       = +2/eps^2 + 3/eps + (10 - 4 pi^2/3 + L^2 + 3 L)
!
! so their finite sum is independent of mu. The poles are asserted on every event.
! Scale variations replay the one central loop call through those formulas.
!#################################################################################!
  double precision function fvirt(yy, vwgt)

  use run_state_mod,   only : ecm_had, xmur, xmuf, xcentral_factor, scale_dynamic
  use run_card_mod,    only : hbarc2_in, virtual_terms
  use cuts_mod,        only : event_cut
  use olp_mod,         only : olp_virtual, olp_set_alphas, NFLAVOUR
  use psgen_interface, only : psgen_event, psgen_ndim, psgen_nexp, psgen_nkin
  use observables_mod, only : obs_mass
  use analysis_mod,    only : analysis_fill
  use scales_mod,      only : nscale, icentral, mur_factors, muf_factors
  use pdfvar_mod,      only : nvar, pdfvar_central, pdfvar_select, pdfvar_restore

  implicit none

  real(8), intent(in) :: yy(16)
  real(8), intent(in) :: vwgt

  real(8), parameter :: PI_VALUE = 3.141592653589793238462643d0, CF = 4d0/3d0
  real(8), parameter :: POLE_TOLERANCE = 1d-6
  real(8), parameter :: I_CONSTANT = 10d0 - 4d0*PI_VALUE**2/3d0

  real(8) :: xa, xb, xr(psgen_ndim), plab(0:3,psgen_nexp), flux, dlips, shat
  real(8) :: mur_central, muf_central, mur, muf, factor_mur, factor_muf, scale_mass
  real(8) :: alphas_central, alphas_var, alphasPDF, log_central, log_var
  real(8) :: f1(-6:6), f2(-6:6)
  real(8) :: born(NFLAVOUR), finite(NFLAVOUR), pole1, pole2, normalisation
  real(8) :: virt_constant(NFLAVOUR), coefficient(NFLAVOUR), luminosity_sum
  real(8) :: weight(nvar)
  integer :: channel, flavour, ivar, member, member_active

  external alphasPDF

    fvirt = 0d0

    xa = yy(1)
    xb = yy(2)
    if (xa <= 0d0 .or. xb <= 0d0 .or. xa >= 1d0 .or. xb >= 1d0) return

    xr(1:psgen_ndim) = yy(3:2+psgen_ndim)
    channel = min(1 + int(yy(3+psgen_ndim)*psgen_nkin), psgen_nkin)

    call psgen_event(ecm_had, xa, xb, channel, xr, plab, flux, dlips)
    if (dlips <= 0d0) return
    if (.not. event_cut(plab)) return

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

    shat = xa*xb*ecm_had*ecm_had
    log_central = log(mur_central**2/shat)

! With massless quarks s, b equal d and c equals u.
    do flavour = 1, NFLAVOUR
       if (flavour <= 2) then
          call olp_virtual(flavour, plab, born(flavour), finite(flavour), pole1, pole2)
          normalisation = alphas_central/(2d0*PI_VALUE)*CF*born(flavour)
          if (abs(pole2/normalisation + 2d0) > POLE_TOLERANCE .or.            &
              abs(pole1/normalisation + 3d0) > POLE_TOLERANCE) then
             print '(A,2ES16.8)', 'ERROR fvirt: poles do not cancel the I-operator: ', &
                                  pole2/normalisation, pole1/normalisation
             stop 1
          endif
          virt_constant(flavour) = finite(flavour)/normalisation              &
                                 + log_central**2 + 3d0*log_central
       else
          born(flavour)          = born(flavour-2)
          virt_constant(flavour) = virt_constant(flavour-2)
       endif
    end do

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
       alphas_var = alphasPDF(mur)
       log_var    = log(mur**2/shat)

       coefficient = 0d0
       if (trim(virtual_terms) /= 'iop')  coefficient = coefficient + virt_constant - log_var**2 - 3d0*log_var
       if (trim(virtual_terms) /= 'virt') coefficient = coefficient + I_CONSTANT + log_var**2 + 3d0*log_var

       luminosity_sum = 0d0
       do flavour = 1, NFLAVOUR
          luminosity_sum = luminosity_sum                                     &
               + (f1(-flavour)*f2(flavour) + f1(flavour)*f2(-flavour))        &
                 *coefficient(flavour)*born(flavour)
       end do
       weight(ivar) = hbarc2_in/flux*alphas_var/(2d0*PI_VALUE)*CF*luminosity_sum*dlips
    end do
    if (member_active /= pdfvar_central()) call pdfvar_restore()

    if (any(weight /= weight)) return
    fvirt = weight(icentral)

    call analysis_fill(plab, psgen_nexp, weight*vwgt)

  end function fvirt
!#################################################################################!
