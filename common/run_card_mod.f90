!#####################################################################
! The run card: what THIS run's input.input means.
!
! The layer above card_parser_mod. The parser knows the file format and
! no section names; this module knows which keys exist, what type each
! one is, which are required, and what counts as a legal value. The two
! other consumers of the card own their own sections the same way --
! cuts_mod reads [cuts], analysis_mod reads the [histograms] rows.
!
! WHY THE CARD IS KEYED. The previous input.input was POSITIONAL: values
! were read in a fixed order with no names, so inserting, removing or
! reordering a line silently shifted every subsequent value into the
! wrong variable -- and a run that reads mee_min into muf_fixed does not
! crash, it just produces a wrong number. The card now uses the same
! [section] / key = value format as the KALEU build next door, parsed by
! card_parser_mod, which makes an unknown key, a missing required key and
! a duplicate key all hard errors.
!
! Every key this module reads is REQUIRED unless the comment says
! otherwise. Defaults are worth having only where there is one obviously
! right answer; for anything that moves the cross section a default is
! just a silent assumption, so the card has to say it.
!#####################################################################

module run_card_mod

  use iso_fortran_env, only: real64
  use card_parser_mod, only: card_load,                                       &
                             card_get_int, card_get_real, card_get_logical,   &
                             card_get_str,                                    &
                             card_require_int, card_require_real, card_require_str
  implicit none
  private
  public :: read_input_file
  public :: scales_label, psgen_label
!------------------[ Run card ]------------------
  public :: seed, n_adapt, it_adapt, n_prod, it_prod, epsrel, resume_run
  public :: out_dir, out_file
  public :: ecm
  public :: pdf_set, pdf_member, pdf_variations
  public :: scale_dynamic, mur_fixed, muf_fixed, central_factor
  public :: mZ_in, wZ_in, mW_in, wW_in, mH_in, wH_in, mtop_in, wtop_in
  public :: gfermi, aewm1_in
  public :: pi_in, hbarc2_in
  public :: sqsc_in
  public :: olp_provider, recola_process
  public :: del_s, del_c
  public :: real_generator_in
  public :: virtual_terms
  public :: iso_mode, iso_n, iso_epsilon, iso_et_max, iso_delta_0, iso_m_hard

!------------------[ Variable declarations ]------------------
  integer      :: seed, n_adapt, it_adapt, n_prod, it_prod
! Relative-accuracy target for the production pass, one per [stats_*] block.
! 0 -- the default when the key is absent -- means no early exit: run exactly
! it_prod iterations, which is what every card did before this key existed.
! Non-zero turns it_prod into a CEILING.
  real(real64) :: epsrel
! Every [stats_*] block, all read regardless of which program is running --
! see read_input_file's header. n_adapt/it_adapt/n_prod/it_prod above are
! the SELECTED block, copied over by program name.
  integer :: n_adapt_lo,   it_adapt_lo,   n_prod_lo,   it_prod_lo
  integer :: n_adapt_kp,   it_adapt_kp,   n_prod_kp,   it_prod_kp
  integer :: n_adapt_virt, it_adapt_virt, n_prod_virt, it_prod_virt
  integer :: n_adapt_rs,   it_adapt_rs,   n_prod_rs,   it_prod_rs
  real(real64) :: epsrel_lo, epsrel_kp, epsrel_virt, epsrel_rs
! Opt-in, and only opt-in: without this key a run never looks for a
! checkpoint, so a stray file left by a killed job can never silently
! change what a later run does.
  logical      :: resume_run
! 256 rather than 128 for the directory: a batch job's collection point is
! routinely a long absolute path, and main REFUSES to truncate one (a report
! written somewhere other than where the launcher looks for it reads as a failed
! job), so a short buffer turns into a hard stop at the worst moment.
  character(len=256) :: out_dir
  character(len=128) :: out_file
  real(real64) :: ecm
  character(len=128) :: pdf_set
  integer      :: pdf_member
! Evaluate every member of the set at the central scale, alongside the
! 7-point scale grid. Requires pdf_member = 0 -- see pdfvar_mod's
! pdfvar_init, which enforces that: members 1.. are error directions,
! not a prediction, and centring the band on one would measure the
! uncertainty around the wrong central value.
  logical      :: pdf_variations
  logical      :: scale_dynamic
  real(real64) :: mur_fixed, muf_fixed
! Multiplies the central scale (dynamic or fixed) before the 7-point
! grid factors are applied, so the band can be re-centred without
! moving the grid's own shape. Default 1.0 -- unchanged from before
! this key existed.
  real(real64) :: central_factor
  real(real64) :: mZ_in, wZ_in, mW_in, wW_in, mH_in, wH_in, mtop_in, wtop_in
  real(real64) :: gfermi
! 1/alpha_em, the analytic chain's (ggza_analytic/) independent EW input --
! see ggza_phys_onshell_mod.f90. Must match mg5_input/param_card.inc's own
! AEWM1 for the analytic provider to agree with Recola/MadGraph, the same
! reason gfermi above is read from the card rather than hardcoded.
  real(real64) :: aewm1_in
!------------------[ Constants ]------------------
! These are in the card so that EVERY number the run uses is in one file
! and travels with the report. They are not physics choices, though, and
! that is the whole difficulty: a mistyped pi is invisible. It reaches
! alpha through alpha_gfermi, and four electroweak vertices turn a wrong
! digit into a flat offset on every distribution at once -- the same
! signature as the b b~ flavour bug and the alpha gap between the two
! providers this build used to carry -- both of which took an external
! comparison to catch, and neither of which is catchable any more now
! that only one provider computes the full amplitude. See olp_mod.f90.
!
! So they are checked against the references below rather than trusted.
  real(real64) :: pi_in, hbarc2_in

! Reference values, used ONLY to validate the card -- never as inputs.
! Both are exact: PI_REF to well past double precision, and HBARC2_REF
! because hbar and c have been DEFINED rather than measured since the 2019
! SI redefinition, so (hbar c)^2 = 0.3893793721 GeV^2 mbarn carries no
! experimental error to round away.
  real(real64), parameter :: PI_REF     = 3.141592653589793238462643d0
  real(real64), parameter :: HBARC2_REF = 0.3893793721d12
!------------------[ Phase-space sampling ]------------------
! The beam momentum fractions are not configurable: x1 and x2 are thrown
! flat in the integrand, which is where they live -- PSGen takes them as
! arguments and never generates them.
  real(real64) :: sqsc_in
!------------------[ Amplitude provider ]------------------
! Which OLP registers and evaluates the amplitudes. Both are linked; this
! only decides which registration block in main runs. Validated at read
! time rather than at registration, so a typo fails before the PDF set is
! loaded instead of a second later.
  character(len=32) :: olp_provider
  character(len=32) :: recola_process

!------------------[ rs_sa-only, read here so every piece can share one card ]---
! [subtraction], [olp] real_generator and [virtual] terms are each read
! by exactly one program (rs_sa, rs_sa, virt_sa respectively) -- but
! this module is linked into all four, and card_report_unused hard-errors on
! any key its own binary never reads. Reading them HERE, unconditionally,
! with the same defaults their private readers used to carry, marks them used
! for every program without inventing a per-section ignore mechanism: a
! program that has no use for del_s simply never looks at it again.
  real(real64) :: del_s, del_c
  character(len=32) :: real_generator_in
  character(len=32) :: virtual_terms
! [isolation] -- applied by rs_sa only (photon_isolation_mod.f90 takes its
! values from here), but read here for the same reason as the block above: all
! four programs share one card, and a Born-kinematics piece has no parton to
! isolate against, so it would otherwise leave these keys unread.
  integer      :: iso_mode
  real(real64) :: iso_n, iso_epsilon, iso_et_max, iso_delta_0, iso_m_hard

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Read and validate the card.
!
! card_report_unused is NOT called here: the [histograms] rows are read
! later, by analysis_mod, and anything still unclaimed at that point is
! a typo. main.f90 calls it once everything has had its turn.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine read_input_file(filename, program)
  character(len=*), intent(in) :: filename
! Which of lo|kp|virt|rs this binary is -- selects which [stats_*] block
! becomes n_adapt/it_adapt/n_prod/it_prod below. See the [stats_*] block
! further down for why there are four of them in one file rather than one.
  character(len=*), intent(in) :: program

  call card_load(filename)

!--- [run] -----------------------------------------------------------
  seed = card_require_int('run','seed')
  if (seed <= 0) call fail('[run] seed must be > 0')

! Resume the PRODUCTION pass from its last completed iteration if a
! checkpoint for this seed is sitting in [output] dir. Defaults to no,
! and a run that says no behaves exactly as it did before checkpointing
! existed. A crash during adaptation is not covered -- checkpointing
! only starts once production does -- so such a run simply begins again.
  resume_run = card_get_logical('run','resume',.false.)

!--- [stats_lo] [stats_kp] [stats_virt] [stats_rs] --------------------
! ONE FILE, FOUR STATISTICS BLOCKS. Every piece shares one seed, one set
! of cuts, one set of physics values -- but lo_sa, kp_sa, virt_sa and
! rs_sa have very different per-call cost and per-call variance, so the
! call counts that make each converge are not the same number and forcing
! them to be would either waste hours on lo_sa or under-sample rs_sa.
! ALL FOUR BLOCKS ARE READ HERE, UNCONDITIONALLY, regardless of which
! `program` this binary is -- exactly like [subtraction]/[olp] real_generator
! above: every key in the file has to be claimed by SOMEONE or
! card_report_unused hard-errors, and this is the one module linked into
! all four binaries. n_adapt/it_adapt build the VEGAS grid and are thrown
! away; n_prod and it_prod produce the quoted number AND every histogram.
! it_prod is not cosmetic: one production iteration is the sample unit
! behind every error bar in the report, so it wants to be large enough
! (>~10) for the spread across iterations to mean something.
  call read_stats_block('stats_lo',   n_adapt_lo,   it_adapt_lo,             &
                                       n_prod_lo,    it_prod_lo,   epsrel_lo)
  call read_stats_block('stats_kp',   n_adapt_kp,   it_adapt_kp,             &
                                       n_prod_kp,    it_prod_kp,   epsrel_kp)
  call read_stats_block('stats_virt', n_adapt_virt, it_adapt_virt,           &
                                       n_prod_virt,  it_prod_virt, epsrel_virt)
  call read_stats_block('stats_rs',   n_adapt_rs,   it_adapt_rs,             &
                                       n_prod_rs,    it_prod_rs,   epsrel_rs)

  select case (trim(program))
    case ('lo')
      n_adapt = n_adapt_lo   ; it_adapt = it_adapt_lo
      n_prod  = n_prod_lo    ; it_prod  = it_prod_lo
      epsrel  = epsrel_lo
    case ('kp')
      n_adapt = n_adapt_kp   ; it_adapt = it_adapt_kp
      n_prod  = n_prod_kp    ; it_prod  = it_prod_kp
      epsrel  = epsrel_kp
    case ('virt')
      n_adapt = n_adapt_virt ; it_adapt = it_adapt_virt
      n_prod  = n_prod_virt  ; it_prod  = it_prod_virt
      epsrel  = epsrel_virt
    case ('rs')
      n_adapt = n_adapt_rs   ; it_adapt = it_adapt_rs
      n_prod  = n_prod_rs    ; it_prod  = it_prod_rs
      epsrel  = epsrel_rs
    case default
      call fail('read_input_file: program must be lo, kp, virt or rs, not "'// &
                trim(program)//'"')
  end select

!--- [output] --------------------------------------------------------
! Just the run directory. Which file each piece writes inside it is
! fixed in that program's own main.f90, not a card key -- see the
! [output] section's own comment in onshell/input.input.
  out_dir  = card_require_str('output','dir')

!--- [beam] ----------------------------------------------------------
  ecm = card_require_real('beam','ecm')
  if (ecm <= 0d0) call fail('[beam] ecm must be > 0')

!--- [pdf] -----------------------------------------------------------
! member defaults to 0: the central member is the only one this build
! uses, and naming it in every card would be noise.
  pdf_set    = card_require_str('pdf','set')
  pdf_member = card_get_int('pdf','member',0)
! Opt-in: without this key a run behaves exactly as it did before the
! PDF-member axis existed -- one member, no [pdf] variations key needed
! at all in an old card.
  pdf_variations = card_get_logical('pdf','variations',.false.)
  if (pdf_variations .and. pdf_member /= 0) &
     call fail('[pdf] variations = yes requires [pdf] member = 0 -- ' //      &
               'members 1.. are error directions, not the central PDF')

!--- [scales] --------------------------------------------------------
! dynamic = yes -> muR = muF = the invariant mass of the final state with
! partons excluded, event by event. Every study in this build computes
! that same quantity: the 2 -> 2 ones sum legs 3..nexp (there is no parton
! to drop), and the real-emission one sums legs 3 and 4 explicitly,
! skipping the radiated parton. It matches the sibling reference's
! scale_invmass, which sums legs 3..n skipping partons -- the scale
! definition has to agree with it or the comparison moves the target with
! no other symptom.
!
! muR IS LIVE, and this comment used to say the opposite. The q q~ tree
! skeleton this card reader came from had four electroweak vertices and no
! power of gs, so muR could not enter and muR-independence was a free
! correctness check. Every study here is now loop-induced and carries
! alpha_s(muR)^2, so a total that does NOT move when the scale moves means
! a frozen alpha_s -- a bug, and the exact inversion of the old rule. See
! run_state_mod.f90's header and olp_set_alphas in olp_mod.f90.
  scale_dynamic = card_get_logical('scales','dynamic',.false.)
  mur_fixed     = card_require_real('scales','mur_fixed')
  muf_fixed     = card_require_real('scales','muf_fixed')
  if (mur_fixed <= 0d0 .or. muf_fixed <= 0d0) &
     call fail('[scales] mur_fixed and muf_fixed must be > 0')

! Re-centres the 7-point grid (scales_mod::mur_factors/muf_factors) without
! moving its shape: with central_factor = 0.5 the grid runs from Q/4 to Q
! rather than from Q/2 to 2Q. Optional -- an old card with no such key gets
! 1.0, which is the grid this build always had.
  central_factor = card_get_real('scales','central_factor',1d0)
  if (central_factor <= 0d0) call fail('[scales] central_factor must be > 0')

!--- [masses] --------------------------------------------------------
! These go straight to whichever provider is in use -- olp_mod reads them
! from here -- and the widths are physics rather than bookkeeping: the
! provider runs in the complex-mass scheme built from exactly these
! numbers, and both the Z and the W are resonant inside the generated
! range. They are also compiled into PSGen's Breit-Wigner mappings -- see
! the mass check in psgen_init.
!
! Every mass has its width beside it. wtop = 0 is a real choice, not an
! omission: the top is only ever an internal t-channel propagator in the
! W+W- topologies, never an external leg, so it is treated as stable.
  mZ_in   = card_require_real('masses','mZ')
  wZ_in   = card_require_real('masses','wZ')
  mW_in   = card_require_real('masses','mW')
  wW_in   = card_require_real('masses','wW')
  mH_in   = card_get_real    ('masses','mH',   125.20d0)
  wH_in   = card_get_real    ('masses','wH',   6.382339d-3)
  mtop_in = card_get_real    ('masses','mtop', 172.57d0)
  wtop_in = card_get_real    ('masses','wtop', 0d0)

  if (mZ_in <= 0d0 .or. mW_in <= 0d0) call fail('[masses] mZ and mW must be > 0')
  if (mH_in <= 0d0 .or. mtop_in <= 0d0) call fail('[masses] mH and mtop must be > 0')
  if (wH_in <  0d0 .or. wtop_in <  0d0) call fail('[masses] wH and wtop must be >= 0')
  if (wZ_in <= 0d0) call fail('[masses] wZ = 0 is rejected: the Z propagator ' // &
                              'is resonant inside the generated range')
! Same reasoning, and it is the addition over the Drell-Yan build: this
! final state is same-flavour, so it resonates in the W too.
  if (wW_in <= 0d0) call fail('[masses] wW = 0 is rejected: this same-flavour ' // &
                              'final state is W-resonant')

!--- [ew] ------------------------------------------------------------
  gfermi = card_require_real('ew','gfermi')
  if (gfermi <= 0d0) call fail('[ew] gfermi must be > 0')

  aewm1_in = card_require_real('ew','aewm1')
  if (aewm1_in <= 0d0) call fail('[ew] aewm1 must be > 0')

!--- [constants] -----------------------------------------------------
! In the card so that every number the run used is in one file and is
! echoed into every report. CHECKED, NOT TRUSTED: neither of these is a
! quantity anyone should be varying, and a wrong digit in either one
! rescales the cross section smoothly with nothing to give it away.
!
! pi enters through alpha (see alpha_gfermi in olp_mod) and therefore
! four times over. The tolerance is loose enough for a legitimately
! short spelling -- 3.14159265 sits at 1.1e-9 -- and tight enough that
! any single mistyped digit in the first eight places is fatal.
  pi_in = card_get_real('constants','pi', PI_REF)
  if (pi_in <= 0d0) call fail('[constants] pi must be > 0')
  if (abs(pi_in - PI_REF) > 1d-8*PI_REF)                                      &
     call fail('[constants] pi is not pi -- this is a typo, not a physics ' // &
               'choice. It enters alpha four times over and would shift '  // &
               'every distribution by a flat factor.')

! hbarc2 converts GeV^-2 to fb and multiplies the answer directly, so an
! error here is a pure rescaling. The 1e-3 window catches the realistic
! mistake, which is a value in the wrong unit -- pb and nb are away by
! factors of 1000.
!
! NOTE the exact value replaced a rounded 0.3894e12 on 2026-08-15. The old
! one was 5.3e-5 (0.0053 %) high, so everything this build produced before
! that date carries the same flat excess -- under the statistical error,
! but do not merge results from either side of it without allowing for it.
  hbarc2_in = card_get_real('constants','hbarc2', HBARC2_REF)
  if (hbarc2_in <= 0d0) call fail('[constants] hbarc2 must be > 0')
  if (abs(hbarc2_in - HBARC2_REF) > 1d-3*HBARC2_REF)                          &
     call fail('[constants] hbarc2 is far from 0.3893793721e12 GeV^2 fb -- '// &
               'check the unit. It multiplies the cross section directly.')

!--- [psgen] ---------------------------------------------------------
! x1 and x2 are thrown flat on (0,1), jacobian 1, and handed to PSGen as
! they are. There is no key for it: the mapping is not a choice this build
! offers, and at 13.6 TeV VEGAS is left to find the small-x region on its
! own grid.
!
! sqsc is PSGen's minimum virtuality: the hard lower edge of EVERY
! s-channel invariant its mappings can produce. It is NOT KALEU's `ecut`,
! which bounds widthless propagators only. Here it bounds the
! Breit-Wigner channels too, so raising it above the smallest invariant
! the cuts do NOT bound -- m(ve ve~) and m(e ve), both unbounded below in
! this run -- deletes phase space that no channel generates back. The
! upstream default leaves every physical region reachable; it is exposed
! only so the efficiency cost of that can be measured.
  sqsc_in = card_get_real('psgen','sqsc', 1d-8)
  if (sqsc_in <= 0d0) call fail('[psgen] sqsc must be > 0')

! Which OLP supplies the amplitudes. Defaults to recola: this process has
! no tree, so Recola is registered with 'NLO' to generate the loop diagrams
! and needs nothing installed. madgraph loads an already-built library,
! nothing to install either.
! `gosam` IS ACCEPTED HERE BUT ONLY ../onshell LINKS IT. This file is
! shared by both studies and cannot tell which one is running; a study
! whose build does not carry a provider rejects it in its own olp_mod
! (add_process_olp: 'unknown OLP') rather than here, which is a clear
! stop at startup and not a wrong number.
  olp_provider = trim(adjustl(card_get_str('olp','provider','recola')))
  if (olp_provider /= 'recola' .and. olp_provider /= 'madgraph' .and. &
      olp_provider /= 'gosam' .and. olp_provider /= 'analytic') &
     call fail('[olp] provider must be recola, madgraph, gosam or analytic, not "'// &
               trim(olp_provider)//'"')

! WHICH AMPLITUDE RECOLA IS ASKED FOR. This is a PHYSICS choice, not a
! code path: 'full' and 'ww' are different quantities and produce
! different cross sections from the same card.
!
!   full  every diagram. The gauge-invariant answer, and the only one
!         that is an observable.
!   ww    only the diagrams carrying an s-channel W+ W- pair, selected
!         through Recola's own bracket syntax. The same subset
!         MadGraph's linked library computes, so this is the setting
!         that makes a runLO comparison against [olp] provider =
!         madgraph a comparison of like with like. It is NOT gauge
!         invariant -- see main.f90 where the string is registered.
!
! Read only when provider = recola. MadGraph's process is fixed inside
! its own prebuilt library and no card key can move it.
  recola_process = trim(adjustl(card_get_str('olp','recola_process','full')))
  if (recola_process /= 'full' .and. recola_process /= 'ww') &
     call fail('[olp] recola_process must be full or ww, not "'// &
               trim(recola_process)//'"')

!--- [subtraction] -----------------------------------------------------------
! rs_sa's technical cut -- see rs_sa/integrand.f90's header for del_s vs
! del_c. Defaults match rs_sa/subtraction_mod.f90's former own, from Stage
! 7's measurement.
  del_s = card_get_real('subtraction','del_s', 1.0d-5)
  del_c = card_get_real('subtraction','del_c', 1.0d-7)
  if (del_s <= 0d0) call fail('[subtraction] del_s must be > 0')
  if (del_c <= 0d0) call fail('[subtraction] del_c must be > 0')

!--- [olp] real_generator (rs_sa only) ----------------------------------------
! There is no ext_Z key any more: rs_sa's diagram selection is fixed -- see
! rs_sa/olp_mod.f90's DIAGRAM SELECTION header. A card still carrying it is
! rejected by card_report_unused, like any other unknown key.
! rs_sa's own main.f90 does the recola/madgraph select-case (and the
! madgraph-forces-recola warning) itself; this only pins the type and the
! recognised values so a typo here is caught the same way every other [olp]
! key is.
  real_generator_in = trim(adjustl(card_get_str('olp','real_generator','recola')))
  if (real_generator_in /= 'recola' .and. real_generator_in /= 'madgraph') &
     call fail('[olp] real_generator must be recola or madgraph, not "'// &
               trim(real_generator_in)//'"')

!--- [virtual] terms (virt_sa only) -------------------------------------------
! Which pieces flo_virt assembles -- virt | iop | virt_iop. Validated here;
! virt_sa's own virt_terms_init does the term-switch wiring, same as before
! when this lived at [run] contribution.
  virtual_terms = trim(adjustl(card_get_str('virtual','terms','virt_iop')))
  if (virtual_terms /= 'virt' .and. virtual_terms /= 'iop' .and. &
      virtual_terms /= 'virt_iop') &
     call fail('[virtual] terms must be virt, iop or virt_iop, not "'// &
               trim(virtual_terms)//'"')

!--- [isolation] (applied by rs_sa only) -------------------------------------
! Frixione smooth cone. Defaults are MATRIX's ppeexa03 ones, unchanged from
! when photon_isolation_mod.f90 read these keys itself -- see that module's
! header for what each parameter means.
  iso_mode    = card_get_int ('isolation','mode',    1)
  iso_n       = card_get_real('isolation','n',       1.0d0)
  iso_epsilon = card_get_real('isolation','epsilon', 0.5d0)
  iso_et_max  = card_get_real('isolation','et_max',  5.0d0)
  iso_delta_0 = card_get_real('isolation','delta_0', 0.4d0)
  iso_m_hard  = card_get_real('isolation','m_hard',  91.188d0)
  if (iso_mode < 1 .or. iso_mode > 3) call fail('[isolation] mode must be 1, 2 or 3')
  if (iso_n       <= 0d0) call fail('[isolation] n must be > 0')
  if (iso_epsilon <= 0d0) call fail('[isolation] epsilon must be > 0')
  if (iso_delta_0 <= 0d0) call fail('[isolation] delta_0 must be > 0')

end subroutine read_input_file

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! One [stats_*] section, same four required keys and the same validation
! [run] used to carry directly.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine read_stats_block(section, n_adapt_out, it_adapt_out,            &
                                        n_prod_out,  it_prod_out, epsrel_out)
  character(len=*), intent(in)  :: section
  integer,          intent(out) :: n_adapt_out, it_adapt_out
  integer,          intent(out) :: n_prod_out,  it_prod_out
  real(real64),     intent(out) :: epsrel_out

  n_adapt_out  = card_require_int(section,'n_adapt')
  it_adapt_out = card_require_int(section,'it_adapt')
  n_prod_out   = card_require_int(section,'n_prod')
  it_prod_out  = card_require_int(section,'it_prod')

! card_get_real, NOT card_require_real: a card written before this key existed
! has no epsrel line, gets 0, and reproduces its old numbers exactly.
  epsrel_out   = card_get_real(section,'epsrel', 0d0)

  if (n_adapt_out  <= 0) call fail('['//trim(section)//'] n_adapt must be > 0')
  if (it_adapt_out <= 0) call fail('['//trim(section)//'] it_adapt must be > 0')
  if (n_prod_out   <= 0) call fail('['//trim(section)//'] n_prod must be > 0')
  if (it_prod_out  <  2) call fail('['//trim(section)//'] it_prod must be >= 2 -- ' // &
                               'the error is the spread ACROSS production ' //         &
                               'iterations and needs at least two of them')
  if (epsrel_out <  0d0) call fail('['//trim(section)//'] epsrel must be >= 0 ' // &
                               '-- 0 means no early exit')
  if (epsrel_out >= 1d0) call fail('['//trim(section)//'] epsrel must be < 1 -- ' // &
                               'a 100% target is not a convergence criterion')

  end subroutine read_stats_block

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! One-line descriptions for the report header. Two runs whose labels
! differ describe different quantities and must not be merged, which is
! the whole reason these travel with the data.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function psgen_label() result(lab)
  character(len=128) :: lab
  lab = 'PSGen + VEGAS, x12=flat'
  end function psgen_label

! Two runs whose labels differ describe different quantities and must
! never be merged (tools/merge_runs.py's CONSISTENCY check), so the
! 7-point grid and any non-default central_factor have to change this
! string exactly when they change the quantity computed.
!
! THE FINAL STATE IS NAMED GENERICALLY, AND HAS TO BE. This module is
! shared by every study, and they have different final states -- Z gamma
! upstairs, Z H in the validation stand-ins, Z H + parton in the
! real-emission one. This string said 'm(e+ e- ve ve~)' until 2026-09-14,
! the final state of the decay-channel build this card reader was written
! for, so every report any of these studies has ever written names a
! process none of them computes. What they DO all compute is the line
! below: the final-state invariant mass with partons excluded. Which
! process produced a report is already recorded, and merge-checked, by
! '# mode =' and '# amplitude =' -- this field only has to describe the
! scale CHOICE, and now does.
  function scales_label() result(lab)
  character(len=128) :: lab
  character(len=32)  :: cf
! F12.4 + adjustl rather than F0.4: gfortran's minimal-width F0.4 drops the
! leading zero, so a re-centred grid would be labelled '.5000' -- which reads
! as a typo in a field that is the merge key. The separating spaces go in the
! concatenation below, not in cf, because trim() would eat a trailing one.
  cf = ''
  if (abs(central_factor - 1d0) > 1d-12) then
     write(cf,'(F12.4)') central_factor
     cf = adjustl(cf)
  endif
  if (scale_dynamic) then
     if (len_trim(cf) > 0) then
        lab = 'dynamic: muR = muF = '//trim(cf)// &
              ' x m(final state, partons excluded), 7-point grid'
     else
        lab = 'dynamic: muR = muF = m(final state, partons excluded), 7-point grid'
     endif
  else
! THE EFFECTIVE SCALES, central_factor ALREADY FOLDED IN -- not the raw card
! values. central_factor multiplies the fixed scales exactly as it multiplies
! the dynamic one (see any integrand's scale block: mur_c = xcentral_factor
! *xmur), so a label quoting mur_fixed alone would be identical for two runs
! that computed different numbers, and merge_runs.py would pool them without
! complaint. This field is the guard against that, so it has to move whenever
! the scale actually does.
     write(lab,'(A,F0.4,A,F0.4,A)') 'fixed: muR = ', central_factor*mur_fixed, &
                                    ', muF = ', central_factor*muf_fixed,      &
                                    ', 7-point grid'
  endif
  end function scales_label

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine fail(msg)
  character(len=*), intent(in) :: msg
  print '(A)','ERROR input card: '//trim(msg)
  stop 1
  end subroutine fail

end module run_card_mod

!#####################################################################
