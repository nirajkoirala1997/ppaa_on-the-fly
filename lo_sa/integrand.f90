!#################################################################################!
! The leading-order integrand for p p -> gamma gamma.
!
!   yy(1), yy(2)   beam momentum fractions, flat
!   yy(3), yy(4)   PSGen final-state variables
!   yy(5)          PSGen channel (one channel at this multiplicity)
!
! The Born has no power of alpha_s, so scale variation enters only through muF.
! One q q~ amplitude per flavour serves both beam orderings; see cuts_mod.
!#################################################################################!
  double precision function flo_LO(yy, vwgt)

  use run_state_mod,   only : ecm_had, xmuf, xcentral_factor, scale_dynamic
  use run_card_mod,    only : hbarc2_in
  use cuts_mod,        only : event_cut
  use olp_mod,         only : olp_born, NFLAVOUR
  use psgen_interface, only : psgen_event, psgen_ndim, psgen_nexp, psgen_nkin
  use observables_mod, only : obs_mass
  use analysis_mod,    only : analysis_fill
  use scales_mod,      only : nscale, icentral, muf_factors
  use pdfvar_mod,      only : nvar, pdfvar_central, pdfvar_select, pdfvar_restore

  implicit none

! Sized to vsup.f's MXDIM: VEGAS passes its whole X array.
  real(8), intent(in) :: yy(16)
  real(8), intent(in) :: vwgt

  real(8) :: xa, xb, xr(psgen_ndim), plab(0:3,psgen_nexp), flux, dlips
  real(8) :: muf_central, muf, factor_muf
  real(8) :: f1(-6:6), f2(-6:6), born(NFLAVOUR), luminosity_sum
  real(8) :: weight(nvar)
  integer :: channel, flavour, ivar, member, member_active

    flo_LO = 0d0

    xa = yy(1)
    xb = yy(2)
    if (xa <= 0d0 .or. xb <= 0d0 .or. xa >= 1d0 .or. xb >= 1d0) return

    xr(1:psgen_ndim) = yy(3:2+psgen_ndim)
    channel = min(1 + int(yy(3+psgen_ndim)*psgen_nkin), psgen_nkin)

    call psgen_event(ecm_had, xa, xb, channel, xr, plab, flux, dlips)
    if (dlips <= 0d0) return
    if (.not. event_cut(plab)) return

    if (scale_dynamic) then
       muf_central = xcentral_factor*obs_mass(plab(:,3) + plab(:,4))
    else
       muf_central = xcentral_factor*xmuf
    endif

    do flavour = 1, NFLAVOUR
       call olp_born(flavour, plab, born(flavour))
    end do

    member_active = pdfvar_central()
    do ivar = 1, nvar
       if (ivar <= nscale) then
          factor_muf = muf_factors(ivar)
          member     = pdfvar_central()
       else
          factor_muf = 1d0
          member     = ivar - nscale - 1
       endif
       if (member /= member_active) then
          call pdfvar_select(member)
          member_active = member
       endif

       muf = factor_muf*muf_central
       call pdf(xa, muf, f1)
       call pdf(xb, muf, f2)

       luminosity_sum = 0d0
       do flavour = 1, NFLAVOUR
          luminosity_sum = luminosity_sum                                     &
               + (f1(-flavour)*f2(flavour) + f1(flavour)*f2(-flavour))*born(flavour)
       end do

! flux = 2*shat from PSGen; dlips already carries its (2*pi) factors.
       weight(ivar) = hbarc2_in/flux*luminosity_sum*dlips
    end do
    if (member_active /= pdfvar_central()) call pdfvar_restore()

    if (any(weight /= weight)) return
    flo_LO = weight(icentral)

    call analysis_fill(plab, psgen_nexp, weight*vwgt)

  end function flo_LO
!#################################################################################!
