!#################################################################################!
! The Catani-Seymour KP term for p p -> e+ e-, on Born kinematics.
!
!   yy(1), yy(2)   Born momentum fractions xa, xb, flat
!   yy(3), yy(4)   PSGen final-state variables
!   yy(5)          PSGen channel
!   yy(6)          the convolution variable
!
! For each beam, with Born fraction xB and h(z) = pdf(xB/z)/z,
!
!   Int_xB^1 dz [ regular(z) h(z) + plus(z) (h(z) - h(1)) ]
!     - h(1) Int_0^xB dz plus(z)  + delta h(1)
!
! the middle term being the part of the plus prescription below xB, where the PDF
! vanishes. L = log(muF^2/shat) is independent of z in these variables.
!
! A gluon on one beam enters with the quark and antiquark of the other beam and
! the q q~ Born of that flavour, per the one-amplitude shortcut of cuts_mod.
!#################################################################################!
  double precision function fkp(yy, vwgt)

  use run_state_mod,     only : ecm_had, xmur, xmuf, xcentral_factor, scale_dynamic
  use run_card_mod,      only : hbarc2_in
  use cuts_mod,          only : event_cut
  use olp_mod,           only : olp_born, NFLAVOUR
  use psgen_interface,   only : psgen_event, psgen_ndim, psgen_nexp, psgen_nkin
  use observables_mod,   only : obs_mass
  use analysis_mod,      only : analysis_fill
  use dy_kp_kernels_mod, only : qq_plus, qq_regular, qq_delta, gq_regular
  use scales_mod,        only : nscale, icentral, mur_factors, muf_factors
  use pdfvar_mod,        only : nvar, pdfvar_central, pdfvar_select, pdfvar_restore

  implicit none

  real(8), intent(in) :: yy(16)
  real(8), intent(in) :: vwgt

  real(8), parameter :: TWO_PI = 6.283185307179586476925287d0

  real(8) :: xa, xb, xr(psgen_ndim), plab(0:3,psgen_nexp), flux, dlips, shat
  real(8) :: mur_central, muf_central, mur, muf, factor_mur, factor_muf, scale_mass
  real(8) :: alphasPDF, scale_log, alphas_var
  real(8) :: fraction(2), z_convolution(2), jac_convolution(2), z_endpoint(2), jac_endpoint(2)
  real(8) :: f1(-6:6), f2(-6:6), f1_z(-6:6), f2_z(-6:6)
  real(8) :: born(NFLAVOUR), lum_born, lum_z(2), lum_gluon(2), beam_sum
  real(8) :: plus0(2), plus1(2), reg0(2), reg1(2), glu0(2), glu1(2), end0(2), end1(2)
  real(8) :: plus_value, endpoint_value, regular_value, gluon_value
  real(8) :: flavour_sum, weight(nvar)
  integer :: channel, flavour, ivar, member, member_active, beam

  external alphasPDF

    fkp = 0d0

    xa = yy(1)
    xb = yy(2)
    if (xa <= 0d0 .or. xb <= 0d0 .or. xa >= 1d0 .or. xb >= 1d0) return
    if (yy(6) <= 0d0 .or. yy(6) >= 1d0) return

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

    shat = xa*xb*ecm_had*ecm_had
    fraction = [xa, xb]
    do beam = 1, 2
       z_convolution(beam)   = fraction(beam) + (1d0 - fraction(beam))*yy(6)
       jac_convolution(beam) = 1d0 - fraction(beam)
       z_endpoint(beam)      = fraction(beam)*yy(6)
       jac_endpoint(beam)    = fraction(beam)
       if (z_convolution(beam) >= 1d0) return
       call qq_plus   (z_convolution(beam), plus0(beam), plus1(beam))
       call qq_regular(z_convolution(beam), reg0(beam),  reg1(beam))
       call gq_regular(z_convolution(beam), glu0(beam),  glu1(beam))
       call qq_plus   (z_endpoint(beam),    end0(beam),  end1(beam))
    end do

! With massless quarks s, b equal d and c equals u.
    do flavour = 1, NFLAVOUR
       if (flavour <= 2) then
          call olp_born(flavour, plab, born(flavour))
       else
          born(flavour) = born(flavour-2)
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
       alphas_var = alphasPDF(mur)
       scale_log  = log(muf*muf/shat)

       call pdf(xa, muf, f1)
       call pdf(xb, muf, f2)
       call pdf(xa/z_convolution(1), muf, f1_z)
       call pdf(xb/z_convolution(2), muf, f2_z)
       f1_z = f1_z/z_convolution(1)
       f2_z = f2_z/z_convolution(2)

       flavour_sum = 0d0
       do flavour = 1, NFLAVOUR
          lum_born     = f1(-flavour)*f2(flavour)   + f1(flavour)*f2(-flavour)
          lum_z(1)     = f1_z(-flavour)*f2(flavour) + f1_z(flavour)*f2(-flavour)
          lum_z(2)     = f1(-flavour)*f2_z(flavour) + f1(flavour)*f2_z(-flavour)
          lum_gluon(1) = f1_z(0)*(f2(flavour) + f2(-flavour))
          lum_gluon(2) = f2_z(0)*(f1(flavour) + f1(-flavour))

          beam_sum = 0d0
          do beam = 1, 2
             plus_value     = plus0(beam) + scale_log*plus1(beam)
             endpoint_value = end0(beam)  + scale_log*end1(beam)
             regular_value  = reg0(beam)  + scale_log*reg1(beam)
             gluon_value    = glu0(beam)  + scale_log*glu1(beam)
             beam_sum = beam_sum                                              &
                  + jac_convolution(beam)*( regular_value*lum_z(beam)         &
                                          + plus_value*(lum_z(beam) - lum_born)&
                                          + gluon_value*lum_gluon(beam) )     &
                  - jac_endpoint(beam)*endpoint_value*lum_born                &
                  + qq_delta()*lum_born
          end do
          flavour_sum = flavour_sum + born(flavour)*beam_sum
       end do

       weight(ivar) = hbarc2_in/flux*dlips*alphas_var/TWO_PI*flavour_sum
    end do
    if (member_active /= pdfvar_central()) call pdfvar_restore()

    if (any(weight /= weight)) return
    fkp = weight(icentral)

    call analysis_fill(plab, psgen_nexp, weight*vwgt)

  end function fkp
!#################################################################################!
