!#####################################################################
! Small shared utilities, and the run card.
!
! WHY THE CARD IS KEYED. The previous input.input was POSITIONAL: values
! were read in a fixed order with no names, so inserting, removing or
! reordering a line silently shifted every subsequent value into the
! wrong variable -- and a run that reads mee_min into muf_fixed does not
! crash, it just produces a wrong number. The card now uses the same
! [section] / key = value format as the KALEU Drell-Yan code, parsed by
! input_card_mod, which makes an unknown key, a missing required key and
! a duplicate key all hard errors.
!
! Every key this module reads is REQUIRED unless the comment says
! otherwise. Defaults are worth having only where there is one obviously
! right answer; for anything that moves the cross section a default is
! just a silent assumption, so the card has to say it.
!#####################################################################

module misc_utils_mod

  use iso_fortran_env, only: real64
  use input_card_mod,  only: card_load,                                       &
                             card_get_int, card_get_real, card_get_logical,   &
                             card_require_int, card_require_real, card_require_str
  implicit none
  private
  public :: dot, read_input_file
  public :: cuts_label, scales_label
!------------------[ Run card ]------------------
  public :: seed, n_adapt, it_adapt, n_prod, it_prod
  public :: out_dir, out_file
  public :: ecm
  public :: pdf_set, pdf_member
  public :: scale_dynamic, mur_fixed, muf_fixed
  public :: mZ_in, wZ_in, mW_in, wW_in, mH_in, mtop_in
  public :: gfermi
  public :: mee_min, mee_max

!------------------[ Variable declarations ]------------------
  integer      :: seed, n_adapt, it_adapt, n_prod, it_prod
  character(len=128) :: out_dir, out_file
  real(real64) :: ecm
  character(len=128) :: pdf_set
  integer      :: pdf_member
  logical      :: scale_dynamic
  real(real64) :: mur_fixed, muf_fixed
  real(real64) :: mZ_in, wZ_in, mW_in, wW_in, mH_in, mtop_in
  real(real64) :: gfermi
  real(real64) :: mee_min, mee_max

contains

  !==========================================================
  real(real64) function dot(p, q) result(d)
    ! Lorentz scalar product: dot(p, q) = p0*q0 - p1*q1 - ...
    real(real64), intent(in) :: p(0:3), q(0:3)
    d = p(0)*q(0) - p(1)*q(1) - p(2)*q(2) - p(3)*q(3)
  end function dot
  !==========================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Read and validate the card.
!
! card_report_unused is NOT called here: the [histograms] rows are read
! later, by analysis_mod, and anything still unclaimed at that point is
! a typo. main.f90 calls it once everything has had its turn.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine read_input_file(filename)
  character(len=*), intent(in) :: filename

  call card_load(filename)

!--- [run] -----------------------------------------------------------
! n_adapt/it_adapt build the VEGAS grid and are thrown away; n_prod and
! it_prod produce the quoted number AND every histogram. it_prod is not
! cosmetic: one production iteration is the sample unit behind every
! error bar in the report, so it wants to be large enough (>~10) for the
! spread across iterations to mean something.
  seed     = card_require_int('run','seed')
  n_adapt  = card_require_int('run','n_adapt')
  it_adapt = card_require_int('run','it_adapt')
  n_prod   = card_require_int('run','n_prod')
  it_prod  = card_require_int('run','it_prod')

  if (seed     <= 0) call fail('[run] seed must be > 0')
  if (n_adapt  <= 0) call fail('[run] n_adapt must be > 0')
  if (it_adapt <= 0) call fail('[run] it_adapt must be > 0')
  if (n_prod   <= 0) call fail('[run] n_prod must be > 0')
  if (it_prod  <  2) call fail('[run] it_prod must be >= 2 -- ' //             &
                               'the error is the spread ACROSS production ' // &
                               'iterations and needs at least two of them')

!--- [output] --------------------------------------------------------
  out_dir  = card_require_str('output','dir')
  out_file = card_require_str('output','file')

!--- [beam] ----------------------------------------------------------
  ecm = card_require_real('beam','ecm')
  if (ecm <= 0d0) call fail('[beam] ecm must be > 0')

!--- [pdf] -----------------------------------------------------------
! member defaults to 0: the central member is the only one this build
! uses, and naming it in every card would be noise.
  pdf_set    = card_require_str('pdf','set')
  pdf_member = card_get_int('pdf','member',0)

!--- [scales] --------------------------------------------------------
! dynamic = yes -> muR = muF = m(e+e-) event by event.
!
! muR IS INERT AT THIS ORDER. The tree amplitude carries no power of gs,
! so nothing in the integrand may read it; it is carried only so a run
! can be labelled and compared against a code that does vary it.
  scale_dynamic = card_get_logical('scales','dynamic',.false.)
  mur_fixed     = card_require_real('scales','mur_fixed')
  muf_fixed     = card_require_real('scales','muf_fixed')
  if (mur_fixed <= 0d0 .or. muf_fixed <= 0d0) &
     call fail('[scales] mur_fixed and muf_fixed must be > 0')

!--- [masses] --------------------------------------------------------
! These reach OpenLoops through constants_mod, and wZ is physics rather
! than bookkeeping: the Z propagator is resonant in the middle of the
! range this run generates.
  mZ_in   = card_require_real('masses','mZ')
  wZ_in   = card_require_real('masses','wZ')
  mW_in   = card_require_real('masses','mW')
  wW_in   = card_require_real('masses','wW')
  mH_in   = card_get_real    ('masses','mH',   125.20d0)
  mtop_in = card_get_real    ('masses','mtop', 172.57d0)

  if (mZ_in <= 0d0 .or. mW_in <= 0d0) call fail('[masses] mZ and mW must be > 0')
  if (wZ_in <= 0d0) call fail('[masses] wZ = 0 is rejected: the Z propagator ' // &
                              'is resonant inside the generated range')

!--- [ew] ------------------------------------------------------------
  gfermi = card_require_real('ew','gfermi')
  if (gfermi <= 0d0) call fail('[ew] gfermi must be > 0')

!--- [cuts] ----------------------------------------------------------
! The only selection this build applies. mee_min is NOT optional: the
! leptons are massless and the amplitude contains gamma* -> e+ e-, so
! the integrand diverges as 1/m(e+e-)^2 and a zero here does not give a
! large cross section, it gives a meaningless one.
!
! mee_max = 0 means "no upper edge" and resolves to ecm.
  mee_min = card_require_real('cuts','mee_min')
  mee_max = card_get_real    ('cuts','mee_max', 0d0)
  if (mee_min <= 0d0) call fail('[cuts] mee_min must be > 0')
  if (mee_max <= 0d0) mee_max = ecm
  if (mee_max <= mee_min) call fail('[cuts] mee_max must exceed mee_min')
  if (mee_max > ecm)      mee_max = ecm

end subroutine read_input_file

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! One-line descriptions for the report header. Two runs whose labels
! differ describe different quantities and must not be merged, which is
! the whole reason these travel with the data.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function cuts_label() result(lab)
  character(len=128) :: lab
  write(lab,'(A,F0.3,A,F0.3)') 'mee in [', mee_min, ',', mee_max
  lab = trim(lab)//']'
  end function cuts_label

  function scales_label() result(lab)
  character(len=128) :: lab
  if (scale_dynamic) then
     lab = 'dynamic: muR = muF = m(e+e-)'
  else
     write(lab,'(A,F0.4,A,F0.4)') 'fixed: muR = ', mur_fixed, ', muF = ', muf_fixed
  endif
  end function scales_label

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine fail(msg)
  character(len=*), intent(in) :: msg
  print '(A)','ERROR input card: '//trim(msg)
  stop 1
  end subroutine fail

end module misc_utils_mod

!#####################################################################
