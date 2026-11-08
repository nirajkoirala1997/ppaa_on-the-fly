! R / sum(D) in the soft and collinear limits, on hand-built momenta. Each channel
! must approach 1 in its own singular limit.
program check_dipole_limits
  use iso_fortran_env,   only : real64
  use run_card_mod,      only : read_input_file
  use olp_mod,           only : olp_init, olp_real, NFLAVOUR, NREAL_CHANNEL
  use dipole_mod,        only : dipoles_evaluate
  use ii_dipole_map_mod, only : ii_dipole_map
  implicit none
  real(real64), parameter :: ALPHAS = 0.118d0, SQRTS = 500d0
  real(real64) :: lambda, momenta(0:3,5), ptilde(0:3,4,2), amp(NREAL_CHANNEL)
  real(real64) :: d_qq(NFLAVOUR,2), d_qqbar(NFLAVOUR,2), d_qbarq(NFLAVOUR,2)
  integer :: step, channel

  call read_input_file('../../input.input', 'rs')
  call olp_init(ALPHAS, 100d0, .false., .true.)

  print '(A)', 'SOFT: q q~  R/(D1+D2)'
  do step = 1, 5
     lambda = 10d0**(-step)
     call build(lambda*SQRTS/2d0, 0.6d0, momenta)
     call evaluate()
     print '(ES10.1,5F16.10)', lambda, (amp(channel)/(d_qq(channel,1)+d_qq(channel,2)), channel=1,5)
  end do

  print '(A)', 'COLLINEAR TO BEAM 1, z = 0.3: q q~ R/D1; g q R/D1; g q~ R/D1'
  do step = 1, 5
     lambda = 10d0**(-step)
     call build(0.3d0*SQRTS/2d0, cos(lambda), momenta)
     call evaluate()
     print '(ES10.1,6F16.10)', lambda, amp(2)/d_qq(2,1), amp(1)/d_qq(1,1), amp(7)/d_qbarq(2,1), &
          amp(6)/d_qbarq(1,1), amp(17)/d_qqbar(2,1), amp(16)/d_qqbar(1,1)
  end do

  print '(A)', 'COLLINEAR TO BEAM 2, z = 0.3: q q~ R/D2; q g R/D2; q~ g R/D2'
  do step = 1, 5
     lambda = 10d0**(-step)
     call build(0.3d0*SQRTS/2d0, -cos(lambda), momenta)
     call evaluate()
     print '(ES10.1,6F16.10)', lambda, amp(2)/d_qq(2,2), amp(1)/d_qq(1,2), amp(12)/d_qqbar(2,2), &
          amp(11)/d_qqbar(1,2), amp(22)/d_qbarq(2,2), amp(21)/d_qbarq(1,2)
  end do

contains

  subroutine evaluate()
  integer :: ichannel
  call ii_dipole_map(momenta(:,1), momenta(:,2), momenta(:,3), momenta(:,4), momenta(:,5), ptilde)
  do ichannel = 1, NREAL_CHANNEL
     call olp_real(ichannel, momenta, amp(ichannel))
  end do
  call dipoles_evaluate(ALPHAS, momenta, ptilde, d_qq, d_qqbar, d_qbarq)
  end subroutine evaluate

! Partonic CM frame; the e+ e- pair decays isotropically-fixed in its rest frame.
  subroutine build(energy5, cos_theta5, p)
  real(real64), intent(in)  :: energy5, cos_theta5
  real(real64), intent(out) :: p(0:3,5)
  real(real64) :: sin_theta5, pair(0:3), mass_pair, lepton(0:3), dot_vec, factor
  real(real64), parameter :: COS_LEPTON = 0.3d0, PHI_LEPTON = 1.1d0
  real(real64) :: sin_lepton

  sin_theta5 = sqrt(max(0d0, 1d0 - cos_theta5**2))
  if (abs(cos_theta5) > 0.9d0) sin_theta5 = sin(acos(abs(cos_theta5)))
  p(:,1) = [SQRTS/2d0, 0d0, 0d0,  SQRTS/2d0]
  p(:,2) = [SQRTS/2d0, 0d0, 0d0, -SQRTS/2d0]
  p(:,5) = energy5*[1d0, sin_theta5*0.8d0, sin_theta5*0.6d0, cos_theta5]
  pair = p(:,1) + p(:,2) - p(:,5)
  mass_pair = sqrt(pair(0)**2 - sum(pair(1:3)**2))
  sin_lepton = sqrt(1d0 - COS_LEPTON**2)
  lepton = 0.5d0*mass_pair*[1d0, sin_lepton*cos(PHI_LEPTON), sin_lepton*sin(PHI_LEPTON), COS_LEPTON]
  dot_vec = sum(pair(1:3)*lepton(1:3))
  factor  = (lepton(0) + dot_vec/(pair(0) + mass_pair))/mass_pair
  p(0,3)   = (pair(0)*lepton(0) + dot_vec)/mass_pair
  p(1:3,3) = lepton(1:3) + pair(1:3)*factor
  p(:,4)   = pair - p(:,3)
  end subroutine build

end program check_dipole_limits
