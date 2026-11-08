!#####################################################################
! Recola 1 amplitudes for p p -> e+ e- at NLO QCD.
!
! One module for every program: lo_sa, kp_sa, virt_sa and rs_sa only differ in
! which process groups they ask olp_init to generate.
!
! Flavour index f = 1..5 is the PDG code: d, u, s, c, b. Legs 3 and 4 are always
! e+ and e-; leg 5 is the radiated parton.
!
! Real-emission channels use the same numbering as DYN's ol_interface:
!    1..5    q     q~  -> e+ e- g
!    6..10   g     q   -> e+ e- q
!   11..15   q     g   -> e+ e- q
!   16..20   g     q~  -> e+ e- q~
!   21..25   q~    g   -> e+ e- q~
!
! Every returned value is summed over final-state and averaged over
! initial-state colour and helicity by Recola.
!
! The loop interference needs muIR to follow the event, which Recola allows only
! after set_dynamic_settings_rcl(1).
!#####################################################################

module olp_mod

  use iso_fortran_env, only: real64
  use recola
  use run_card_mod,    only: mZ_in, wZ_in, mW_in, wW_in, mtop_in, wtop_in, gfermi
  implicit none
  private
  public :: olp_init, olp_set_alphas, olp_born, olp_virtual, olp_real
  public :: olp_set_collier_dir
  public :: NFLAVOUR, NREAL_CHANNEL

  integer, parameter :: NFLAVOUR = 5, NREAL_CHANNEL = 25
  character(len=1), parameter :: quark_name(NFLAVOUR) = ['d','u','s','c','b']

  ! Recola process numbers: Born and loop are 1..5, real channels follow.
  integer, parameter :: REAL_ID_OFFSET = NFLAVOUR

  logical, save :: generated = .false., born_is_loop = .false., have_real = .false.
  character(len=256), save :: collier_dir = ''

contains

  subroutine olp_set_collier_dir(dir)
  character(len=*), intent(in) :: dir
  collier_dir = dir
  end subroutine olp_set_collier_dir

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Set parameters and generate the requested groups. with_born is always
! implied: Born is needed by every program. with_loop generates Born and
! loop together; with_real adds the 25 real-emission channels.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine olp_init(alphas_ref, mu_ref, with_loop, with_real)
  real(real64), intent(in) :: alphas_ref, mu_ref
  logical,      intent(in) :: with_loop, with_real

  integer :: flavour, channel
  character(len=3) :: born_order

  call set_output_file_rcl('*')
  if (len_trim(collier_dir) > 0) then
     call execute_command_line('mkdir -p '//trim(collier_dir))
     call set_collier_output_dir_rcl(trim(collier_dir))
  endif
  call set_print_level_amplitude_rcl(0)
  call set_print_level_squared_amplitude_rcl(0)
  call set_print_level_correlations_rcl(0)
  call set_print_level_RAM_rcl(0)

  call set_complex_mass_scheme_rcl
  call set_pole_mass_z_rcl  (mZ_in,   wZ_in)
  call set_pole_mass_w_rcl  (mW_in,   wW_in)
  call set_pole_mass_top_rcl(mtop_in, wtop_in)
  call set_pole_mass_up_rcl     (0d0)
  call set_pole_mass_down_rcl   (0d0)
  call set_pole_mass_charm_rcl  (0d0, 0d0)
  call set_pole_mass_strange_rcl(0d0)
  call set_pole_mass_bottom_rcl (0d0, 0d0)
  call use_gfermi_scheme_rcl(a=alpha_gmu_complex())
  call set_alphas_rcl(alphas_ref, mu_ref, -1)

  if (with_loop) then
     call set_dynamic_settings_rcl(1)
     call set_mu_ir_rcl(mu_ref)
     call set_mu_uv_rcl(mu_ref)
     call set_delta_ir_rcl(0d0, 0d0)
     call set_compute_ir_poles_rcl(1)
  endif

  born_order = 'LO'
  if (with_loop) born_order = 'NLO'
  do flavour = 1, NFLAVOUR
     call define_process_rcl(flavour,                                         &
          quark_name(flavour)//' '//quark_name(flavour)//'~ -> e+ e-',        &
          trim(born_order))
     if (with_loop) then
        call unselect_all_gs_powers_LoopAmpl_rcl(flavour)
        call select_gs_power_LoopAmpl_rcl(flavour, 2)
     endif
  end do

  if (with_real) then
     do channel = 1, NREAL_CHANNEL
        call define_process_rcl(REAL_ID_OFFSET + channel,                     &
                                trim(real_process_string(channel)), 'LO')
     end do
  endif

  call generate_processes_rcl

  born_is_loop = with_loop
  have_real    = with_real
  generated    = .true.
  end subroutine olp_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! G_mu alpha built from the COMPLEX W and Z masses, as OpenLoops does under
! ew_scheme = 1 and as DYN therefore used. Recola's own G_mu option uses real
! masses, which gives an alpha smaller by a factor 1.000484 at these inputs.
  real(real64) function alpha_gmu_complex() result(alpha)
  real(real64), parameter :: PI_VALUE = 3.141592653589793238462643d0
  complex(real64) :: mass2_w, mass2_z

  mass2_w = cmplx(mW_in**2, -mW_in*wW_in, kind=real64)
  mass2_z = cmplx(mZ_in**2, -mZ_in*wZ_in, kind=real64)
  alpha = sqrt(2d0)/PI_VALUE*gfermi*abs(mass2_w*(1d0 - mass2_w/mass2_z))
  end function alpha_gmu_complex

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function real_process_string(channel) result(process)
  integer, intent(in) :: channel
  character(len=32) :: process
  integer :: flavour
  character(len=1) :: q

  flavour = mod(channel-1, NFLAVOUR) + 1
  q = quark_name(flavour)
  select case ((channel-1)/NFLAVOUR)
    case (0) ; process = q//' '//q//'~ -> e+ e- g'
    case (1) ; process = 'g '//q//' -> e+ e- '//q
    case (2) ; process = q//' g -> e+ e- '//q
    case (3) ; process = 'g '//q//'~ -> e+ e- '//q//'~'
    case (4) ; process = q//'~ g -> e+ e- '//q//'~'
    case default ; process = ''
  end select
  end function real_process_string

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Per-event coupling and scales. muIR follows muR for the loop.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine olp_set_alphas(alphas, mu)
  real(real64), intent(in) :: alphas, mu
  call require_generated('olp_set_alphas')
  call set_alphas_rcl(alphas, mu, -1)
  if (born_is_loop) then
     call set_mu_ir_rcl(mu)
     call set_mu_uv_rcl(mu)
  endif
  end subroutine olp_set_alphas

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine olp_born(flavour, momenta, answer)
  integer,      intent(in)  :: flavour
  real(real64), intent(in)  :: momenta(0:3,4)
  real(real64), intent(out) :: answer

  call require_generated('olp_born')
  call compute_process_rcl(flavour, momenta, 'LO')
  call get_squared_amplitude_rcl(flavour, 0, 'LO', answer)
  if (answer /= answer) answer = 0d0
  end subroutine olp_born

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! The renormalised one-loop QCD interference 2 Re(M0* M1): finite part and
! the coefficients of the IR poles, in Recola's normalisation.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine olp_virtual(flavour, momenta, born, finite, pole1, pole2)
  integer,      intent(in)  :: flavour
  real(real64), intent(in)  :: momenta(0:3,4)
  real(real64), intent(out) :: born, finite, pole1, pole2

  call require_generated('olp_virtual')
  if (.not. born_is_loop) then
     print '(A)','ERROR olp_virtual: olp_init was called without with_loop'
     stop 1
  endif
  call compute_process_rcl(flavour, momenta, 'NLO')
  call get_squared_amplitude_rcl(flavour, 0, 'LO', born)
  call get_squared_amplitude_rcl(flavour, 1, 'NLO', finite)
  call get_squared_amplitude_rcl(flavour, 1, 'NLO-IR1', pole1)
  call get_squared_amplitude_rcl(flavour, 1, 'NLO-IR2', pole2)
  end subroutine olp_virtual

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine olp_real(channel, momenta, answer)
  integer,      intent(in)  :: channel
  real(real64), intent(in)  :: momenta(0:3,5)
  real(real64), intent(out) :: answer

  call require_generated('olp_real')
  if (.not. have_real) then
     print '(A)','ERROR olp_real: olp_init was called without with_real'
     stop 1
  endif
  call compute_process_rcl(REAL_ID_OFFSET + channel, momenta, 'LO')
  call get_squared_amplitude_rcl(REAL_ID_OFFSET + channel, 1, 'LO', answer)
  if (answer /= answer) answer = 0d0
  end subroutine olp_real

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine require_generated(caller)
  character(len=*), intent(in) :: caller
  if (.not. generated) then
     print '(A)','ERROR '//caller//': called before olp_init'
     stop 1
  endif
  end subroutine require_generated

end module olp_mod
