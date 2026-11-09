! Pointwise checks of the Recola q q~ -> gamma gamma amplitudes at random points:
!   - Born / analytic (2/3) e^4 Q^4 (t/u + u/t), which has no identical-photon 1/2;
!   - beam swap, q(p1) q~(p2) against q(p2) q~(p1), for Born, loop and real;
!   - exchange of the two photons;
!   - flavour equalities s = b = d and c = u;
!   - virtual poles in units of alpha_s/(2 pi) CF Born at mu^2 = shat (expect -2, -3).
!
!     make clean; make PROG=check_diphoton_amplitudes; ./check_diphoton_amplitudes
program check_diphoton_amplitudes
  use iso_fortran_env, only: real64
  use run_card_mod,    only: read_input_file
  use olp_mod,         only: olp_init, olp_set_alphas, olp_born, olp_virtual, olp_real, &
                             olp_alpha_ew
  implicit none

  real(real64), parameter :: ALPHAS_FIXED = 0.118d0, ROOT_S = 500d0
  real(real64), parameter :: PI_VALUE = 3.141592653589793238462643d0, CF = 4d0/3d0
  integer, parameter :: NPOINTS = 5, DOWN = 1, UP = 2, STRANGE = 3, CHARM = 4, BOTTOM = 5
  integer, parameter :: REAL_QQBAR_DOWN = 1, REAL_QQBAR_UP = 2

  real(real64) :: born_momenta(0:3,4), real_momenta(0:3,5), swapped_born(0:3,4), swapped_real(0:3,5)
  real(real64) :: charge(5), coupling, mandelstam_t, mandelstam_u, analytic
  real(real64) :: value_a, value_b, born, finite, pole1, pole2, finite_b, born_b, pole1_b, pole2_b
  integer :: ipoint, flavour, channel

  charge = [-1d0/3d0, 2d0/3d0, -1d0/3d0, 2d0/3d0, -1d0/3d0]

  call read_input_file('../input.input', 'lo')
  call olp_init(ALPHAS_FIXED, ROOT_S, .true., .true.)
  call olp_set_alphas(ALPHAS_FIXED, ROOT_S)
  coupling = 4d0*PI_VALUE*olp_alpha_ew()

  do ipoint = 1, NPOINTS
     call random_born(born_momenta)
     call random_real(real_momenta)
     print '(A,I0)', 'POINT ', ipoint

     mandelstam_t = -2d0*dot4(born_momenta(:,1), born_momenta(:,3))
     mandelstam_u = -2d0*dot4(born_momenta(:,1), born_momenta(:,4))
     do flavour = DOWN, UP
        analytic = 2d0/3d0*coupling**2*charge(flavour)**4 &
                   *(mandelstam_t/mandelstam_u + mandelstam_u/mandelstam_t)
        call olp_born(flavour, born_momenta, value_a)
        print '(A,I2,F20.15)', '  born/analytic  flavour', flavour, value_a/analytic
     end do

     swapped_born = born_momenta(:,[2,1,3,4])
     swapped_real = real_momenta(:,[2,1,3,4,5])
     do flavour = DOWN, UP
        call olp_virtual(flavour, born_momenta, born, finite, pole1, pole2)
        call olp_virtual(flavour, swapped_born, born_b, finite_b, pole1_b, pole2_b)
        print '(A,I2,2F20.15)', '  beam swap born, virt ratio   ', flavour, born_b/born, finite_b/finite
        print '(A,I2,2F20.12)', '  poles eps^-2, eps^-1         ', flavour, &
             pole2/(ALPHAS_FIXED/(2d0*PI_VALUE)*CF*born), pole1/(ALPHAS_FIXED/(2d0*PI_VALUE)*CF*born)

        call olp_virtual(flavour, born_momenta(:,[1,2,4,3]), born_b, finite_b, pole1_b, pole2_b)
        print '(A,I2,2F20.15)', '  photon swap born, virt ratio ', flavour, born_b/born, finite_b/finite
     end do

     do channel = REAL_QQBAR_DOWN, REAL_QQBAR_UP
        call olp_real(channel, real_momenta, value_a)
        call olp_real(channel, swapped_real, value_b)
        print '(A,I2,F20.15)', '  beam swap real ratio         ', channel, value_b/value_a
        call olp_real(channel, real_momenta(:,[1,2,4,3,5]), value_b)
        print '(A,I2,F20.15)', '  photon swap real ratio       ', channel, value_b/value_a
     end do

     call compare_flavours(STRANGE, DOWN)
     call compare_flavours(BOTTOM, DOWN)
     call compare_flavours(CHARM, UP)
  end do

contains

  subroutine compare_flavours(flavour_heavy, flavour_light)
  integer, intent(in) :: flavour_heavy, flavour_light
  real(real64) :: born_h, finite_h, born_l, finite_l, ignored_1, ignored_2, real_h, real_l
  call olp_virtual(flavour_heavy, born_momenta, born_h, finite_h, ignored_1, ignored_2)
  call olp_virtual(flavour_light, born_momenta, born_l, finite_l, ignored_1, ignored_2)
  call olp_real(flavour_heavy, real_momenta, real_h)
  call olp_real(flavour_light, real_momenta, real_l)
  print '(A,2I2,3F20.15)', '  flavour ratio born virt real ', flavour_heavy, flavour_light, &
       born_h/born_l, finite_h/finite_l, real_h/real_l
  end subroutine compare_flavours

  real(real64) function dot4(a, b) result(product_value)
  real(real64), intent(in) :: a(0:3), b(0:3)
  product_value = a(0)*b(0) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)
  end function dot4

  subroutine beams(momenta)
  real(real64), intent(inout) :: momenta(0:,:)
  momenta(:,1) = 0.5d0*ROOT_S*[1d0, 0d0, 0d0,  1d0]
  momenta(:,2) = 0.5d0*ROOT_S*[1d0, 0d0, 0d0, -1d0]
  end subroutine beams

  function random_direction() result(direction)
  real(real64) :: direction(3), random_pair(2), cos_theta, sin_theta, phi
  call random_number(random_pair)
  cos_theta = 2d0*random_pair(1) - 1d0
  sin_theta = sqrt(1d0 - cos_theta**2)
  phi = 2d0*PI_VALUE*random_pair(2)
  direction = [sin_theta*cos(phi), sin_theta*sin(phi), cos_theta]
  end function random_direction

  subroutine random_born(momenta)
  real(real64), intent(out) :: momenta(0:3,4)
  real(real64) :: direction(3)
  call beams(momenta)
  direction = random_direction()
  momenta(:,3) = 0.5d0*ROOT_S*[1d0,  direction]
  momenta(:,4) = 0.5d0*ROOT_S*[1d0, -direction]
  end subroutine random_born

! The parton is hard and isotropic; the photons decay isotropically in the
! recoiling system's rest frame.
  subroutine random_real(momenta)
  real(real64), intent(out) :: momenta(0:3,5)
  real(real64) :: direction(3), system(0:3), photon_rest(0:3), velocity(3), random_fraction
  real(real64) :: parton_energy, system_mass, gamma_factor, velocity_dot
  call beams(momenta)
  call random_number(random_fraction)
  parton_energy = ROOT_S*(0.1d0 + 0.2d0*random_fraction)
  direction = random_direction()
  momenta(:,5) = parton_energy*[1d0, direction]
  system = momenta(:,1) + momenta(:,2) - momenta(:,5)
  system_mass = sqrt(dot4(system, system))
  direction = random_direction()
  photon_rest = 0.5d0*system_mass*[1d0, direction]
  velocity = system(1:3)/system(0)
  gamma_factor = system(0)/system_mass
  velocity_dot = dot_product(velocity, photon_rest(1:3))
  momenta(0,3) = gamma_factor*(photon_rest(0) + velocity_dot)
  momenta(1:3,3) = photon_rest(1:3) &
       + velocity*(gamma_factor*photon_rest(0) + gamma_factor**2/(gamma_factor + 1d0)*velocity_dot)
  momenta(:,4) = system - momenta(:,3)
  end subroutine random_real

end program check_diphoton_amplitudes
