!#####################################################################
! Checkpoint and resume for the VEGAS production pass.
!
! A production campaign here is 240M events. Before this module a job
! killed partway through one had exactly one option: start again at
! iteration zero. This writes enough state at every production iteration
! boundary that the next launch can pick up at the next iteration and
! land on a report identical to the uninterrupted one, bit for bit --
! which is the only standard worth holding it to, since anything less
! means the RNG stream or the grid came back merely plausible.
!
! WHAT IS SAVED, and why it is exactly this:
!
!   the VEGAS grid + running sums   VEGAS_GET_STATE (vsup.f)
!   the RNG position                BRM48O          (vsup.f)
!   the histogram accumulators      analysis_write_state
!
! Those three are the entire live state at an iteration boundary. Note
! what is NOT here: the per-iteration scratch inside vegas_m, and hcur
! inside each histogram, both provably zero at the boundary; and
! avgi/sd/chi2, rebuilt in full from the running sums by the next
! iteration. See the checkpoint-hook comment in vsup.f.
!
! SCOPE: THE PRODUCTION PASS ONLY. A crash during adaptation still costs
! a full restart. Adaptation is the cheaper pass, and resuming it would
! need a per-iteration hook inside vegas_m's own ISTAT=3 loop that does
! not exist. The consequence is benign and worth stating: a crash before
! production begins leaves no checkpoint at all, so the next launch
! simply runs normally.
!
! TRIGGERED ONLY BY `[run] resume = yes`. There is deliberately no
! auto-detection of a stray checkpoint file: without the key, a run
! behaves precisely as it did before this module existed.
!#####################################################################

module checkpoint_mod

  use iso_fortran_env, only: real64
  use card_parser_mod, only: card_dump
  use analysis_mod,    only: analysis_write_state, analysis_read_state

  implicit none
  private
  public :: checkpoint_configure, checkpoint_exists
  public :: checkpoint_resume, checkpoint_delete, checkpoint_snapshot

  ! Must match XI's declared shape in vegas_m (vsup.f). Not derived from
  ! anything: vegas_m's arrays are fixed, and a mismatch here would
  ! overrun them through the ENTRY point rather than be caught.
  integer, parameter :: NDMX = 50, MXDIM = 16

  integer, parameter :: IU = 79

  logical, save :: enabled  = .false.
  logical, save :: want_res = .false.
  character(len=512), save :: cpfile = ''

  ! Identity of the run this checkpoint belongs to, checked on resume.
  integer, save :: my_seed = 0, my_ndim = 0, my_nprod = 0, my_itprod = 0

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Called once from main, AFTER analysis_book_from_card: checkpoint_resume
! restores into already-booked histograms and cannot create them.
!
! The file is seed-tagged for the same reason reports are (LO_<seed>.dat):
! a campaign runs many seeds into one out_dir, and an untagged checkpoint
! would have every one of them overwriting the others' state.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine checkpoint_configure(resume_requested, seed, out_dir, ndim, &
                                  n_prod, it_prod)
  logical,          intent(in) :: resume_requested
  integer,          intent(in) :: seed, ndim, n_prod, it_prod
  character(len=*), intent(in) :: out_dir

  character(len=32) :: tag

  want_res  = resume_requested
  my_seed   = seed
  my_ndim   = ndim
  my_nprod  = n_prod
  my_itprod = it_prod

  write(tag,'(A,I0,A)') 'checkpoint_', seed, '.dat'
  cpfile  = trim(out_dir)//'/'//trim(tag)
  enabled = .true.

  end subroutine checkpoint_configure

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! A checkpoint is only ever LOOKED FOR when the card asked for one. The
! want_res test first is what structurally rules out ever resuming from
! a file nobody mentioned.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  logical function checkpoint_exists() result(yes)
  yes = .false.
  if (.not. want_res) return
  inquire(file=trim(cpfile), exist=yes)
  end function checkpoint_exists

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Take a snapshot. Called at every production iteration boundary, via
! the external wrapper at the bottom of this file.
!
! status='replace' every time: one rolling file per seed, so disk use is
! bounded no matter how long the campaign runs, and there is never a
! question of which of several files is the live one.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine checkpoint_snapshot()
  integer :: nd, ndo, ng, mds, npg, ndm, it
  integer :: ijkl, ntot, ntot2
  integer :: ios, j
  real(real64) :: calls, dxg, dv2g, xnd, xjac, si, si2, swgt, schi
  real(real64) :: xi(NDMX,MXDIM)

  if (.not. enabled) return

  call vegas_get_state(nd, ndo, ng, mds, npg, ndm, it,           &
                       calls, dxg, dv2g, xnd, xjac,              &
                       si, si2, swgt, schi, xi)
  call brm48o(ijkl, ntot, ntot2)

  open(unit=IU, file=trim(cpfile), status='replace', action='write', &
       iostat=ios)
  if (ios /= 0) then
     ! Not fatal on purpose. A checkpoint that cannot be written is a
     ! lost safety net, not a wrong answer -- killing a run that is
     ! otherwise healthy would be the worse outcome by far.
     print '(A)','WARNING checkpoint: cannot write '//trim(cpfile)
     return
  endif

  write(IU,'(A)')   '# VEGAS production state. Written once per iteration.'
  write(IU,'(A)')   '# ==== CHECKPOINT ===='
  write(IU,'(A,I0)')'# seed = ',            my_seed
  write(IU,'(A,I0)')'# ndim = ',            my_ndim
  write(IU,'(A,I0)')'# n_prod = ',          my_nprod
  write(IU,'(A,I0)')'# it_prod = ',         my_itprod
  write(IU,'(A,I0)')'# iterations_done = ', it

  call card_dump(IU)

  write(IU,'(A)')     '# ==== RNG ===='
  write(IU,'(3I12)')  ijkl, ntot, ntot2
  write(IU,'(A)')     '# ==== END RNG ===='

  write(IU,'(A)')      '# ==== VEGAS ===='
  write(IU,'(7I12)')   nd, ndo, ng, mds, npg, ndm, it
  write(IU,'(E26.17)') calls, dxg, dv2g, xnd, xjac, si, si2, swgt, schi
  do j = 1, MXDIM
     write(IU,'(E26.17)') xi(:,j)
  end do
  write(IU,'(A)')      '# ==== END VEGAS ===='

  write(IU,'(A)') '# ==== HISTOGRAMS ===='
  call analysis_write_state(IU)
  write(IU,'(A)') '# ==== END HISTOGRAMS ===='

  write(IU,'(A)') '# ==== END CHECKPOINT ===='

  close(IU)

  end subroutine checkpoint_snapshot

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Restore. Called before brm48i and before vsup.
!
! EVERY mismatch is fatal, never a silent fall back to a fresh run. The
! whole value of a resume is that its answer is indistinguishable from
! an uninterrupted one; a resume that quietly patched over a changed
! card would produce a number built from two different setups and no
! way to tell from the report.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine checkpoint_resume(start_it)
  integer, intent(out) :: start_it

  integer :: nd, ndo, ng, mds, npg, ndm, it
  integer :: ijkl, ntot, ntot2
  integer :: ios, j, it_cap
  real(real64) :: calls, dxg, dv2g, xnd, xjac, si, si2, swgt, schi
  real(real64) :: xi(NDMX,MXDIM)

  open(unit=IU, file=trim(cpfile), status='old', action='read', iostat=ios)
  if (ios /= 0) then
     print '(A)','ERROR checkpoint_resume: cannot read '//trim(cpfile)
     stop 1
  endif

  call skip_to(IU, '# ==== CHECKPOINT ====')

  call expect_int(IU, '# seed = ',            my_seed,   'seed')
  call expect_int(IU, '# ndim = ',            my_ndim,   'ndim')
  call expect_int(IU, '# n_prod = ',          my_nprod,  'n_prod')
  ! it_prod is a CEILING, not a property of the data on disk. With an epsrel
  ! target a run stops early and is meant to be continued by RAISING it, so a
  ! ceiling larger than the one the checkpoint was written under is fine.
  ! Only a smaller one is an error: it would discard iterations already
  ! accumulated. n_prod above stays a strict equality -- that one IS a
  ! property of the accumulated data.
  call read_int(IU, '# it_prod = ', it_cap)
  if (it_cap > my_itprod) then
     print '(A,I0,A,I0,A)','ERROR checkpoint_resume: it_prod is ',       &
           my_itprod, ' for this run but the checkpoint was written '//  &
           'under a ceiling of ', it_cap,                                &
           ' -- lowering it would discard iterations already on disk'
     stop 1
  endif
  call read_int  (IU, '# iterations_done = ', start_it)

  if (start_it >= my_itprod) then
     print '(A,I0,A,I0,A)','ERROR checkpoint_resume: the checkpoint '//  &
           'already holds ', start_it, ' of ', my_itprod,                &
           ' production iterations -- the ceiling is exhausted. Raise '//&
           'it_prod in this program''s [stats_*] block and rerun to '// &
           'add iterations; lowering epsrel alone cannot buy any.'
     stop 1
  endif
  if (start_it < 1) then
     print '(A,I0)','ERROR checkpoint_resume: nonsensical '// &
           'iterations_done = ', start_it
     stop 1
  endif

  call check_card(IU)

  call skip_to(IU, '# ==== RNG ====')
  read(IU,'(3I12)') ijkl, ntot, ntot2

  call skip_to(IU, '# ==== VEGAS ====')
  read(IU,'(7I12)')   nd, ndo, ng, mds, npg, ndm, it
  read(IU,'(E26.17)') calls, dxg, dv2g, xnd, xjac, si, si2, swgt, schi
  do j = 1, MXDIM
     read(IU,'(E26.17)') xi(:,j)
  end do

  call skip_to(IU, '# ==== HISTOGRAMS ====')
  call analysis_read_state(IU)

  close(IU)

  ! Once, with the real values. brm48i replays the stream by generating
  ! and discarding ntot2*1e9 + ntot numbers, so calling it twice -- once
  ! to seed, once to position -- would pay that cost for nothing.
  call brm48i(ijkl, ntot, ntot2)

  call vegas_put_state(nd, ndo, ng, mds, npg, ndm, it,           &
                       calls, dxg, dv2g, xnd, xjac,              &
                       si, si2, swgt, schi, xi)

  write(*,'(A,I0,A,I0)') '  checkpoint restored from '//trim(cpfile)// &
                         ' -- resuming after iteration ', start_it,    &
                         ' of ', my_itprod

  end subroutine checkpoint_resume

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Called after the report is safely written, not before: a crash during
! the report write still leaves something to resume from.
!
! Beyond tidiness, this used to make "is this seed still in flight"
! answerable from the presence of the file alone. With an epsrel target
! that is no longer true: a run that CONVERGED also keeps its checkpoint,
! deliberately, so a rerun with a tighter epsrel can continue from it.
! The test is now "checkpoint present AND report absent".
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine checkpoint_delete()
  integer :: ios
  logical :: there

  if (.not. enabled) return
  inquire(file=trim(cpfile), exist=there)
  if (.not. there) return

  open(unit=IU, file=trim(cpfile), status='old', iostat=ios)
  if (ios /= 0) return
  close(IU, status='delete')

  end subroutine checkpoint_delete

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! The card check. card_dump writes to a unit, so a fresh dump goes to a
! scratch unit and the two are compared line by line.
!
! This catches everything the four scalar fields above do not: cuts,
! masses, widths, PDF, scales, histogram binning. A resume across an
! edited card would otherwise blend two setups into one number.
! [stats_*] lines are skipped: this program's own n_prod/it_prod are
! checked above, and other programs' call counts are not its business.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine check_card(iunit)
  integer, intent(in) :: iunit

  integer, parameter :: ISCR = 80
  character(len=512) :: a, b
  integer :: iosa, iosb

  open(unit=ISCR, status='scratch', action='readwrite')
  call card_dump(ISCR)
  rewind(ISCR)

  call skip_to(iunit, '# ==== INPUT CARD ====')
  call skip_to(ISCR,  '# ==== INPUT CARD ====')

  do
     do
        read(iunit,'(A)',iostat=iosa) a
        if (iosa /= 0 .or. index(a, '# [stats_') /= 1) exit
     end do
     do
        read(ISCR, '(A)',iostat=iosb) b
        if (iosb /= 0 .or. index(b, '# [stats_') /= 1) exit
     end do
     if (iosa /= 0 .or. iosb /= 0) exit
     if (trim(a) /= trim(b)) then
        print '(A)','ERROR checkpoint_resume: the run card has changed '// &
                    'since the checkpoint was written.'
        print '(A)','  checkpoint: '//trim(a)
        print '(A)','  now:        '//trim(b)
        stop 1
     endif
     if (trim(a) == '# ==== END INPUT CARD ====') exit
  end do

  if (iosa /= 0 .or. iosb /= 0) then
     print '(A)','ERROR checkpoint_resume: truncated input-card block'
     stop 1
  endif

  close(ISCR)

  end subroutine check_card

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Small readers. The format is `# key = value` on its own line, in a
! fixed order, so these stay this simple deliberately -- a checkpoint is
! written and read by this file alone and never hand-edited.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine skip_to(iunit, marker)
  integer,          intent(in) :: iunit
  character(len=*), intent(in) :: marker
  character(len=512) :: line
  integer :: ios
  do
     read(iunit,'(A)',iostat=ios) line
     if (ios /= 0) then
        print '(A)','ERROR checkpoint: "'//trim(marker)//'" not found'
        stop 1
     endif
     if (trim(line) == marker) return
  end do
  end subroutine skip_to

  subroutine read_int(iunit, prefix, val)
  integer,          intent(in)  :: iunit
  character(len=*), intent(in)  :: prefix
  integer,          intent(out) :: val
  character(len=512) :: line
  integer :: ios
  read(iunit,'(A)',iostat=ios) line
  if (ios /= 0 .or. index(line, prefix) /= 1) then
     print '(A)','ERROR checkpoint: expected a line beginning "'// &
                 trim(prefix)//'", got "'//trim(line)//'"'
     stop 1
  endif
  read(line(len(prefix)+1:),*,iostat=ios) val
  if (ios /= 0) then
     print '(A)','ERROR checkpoint: unreadable value in "'//trim(line)//'"'
     stop 1
  endif
  end subroutine read_int

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! The seed is checked HERE, against the resolved runtime value, and not
! left to the card comparison above. `./runLO <seed>` overrides the
! card, so [run] seed is typically one fixed placeholder across a whole
! campaign -- meaning the card text is identical between seeds and would
! happily accept seed 1002 resuming from seed 1001's checkpoint sitting
! in the same out_dir. Only this test separates them.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine expect_int(iunit, prefix, want, what)
  integer,          intent(in) :: iunit, want
  character(len=*), intent(in) :: prefix, what
  integer :: got
  call read_int(iunit, prefix, got)
  if (got /= want) then
     print '(A,I0,A,I0)','ERROR checkpoint_resume: '//trim(what)//     &
           ' is ', want, ' for this run but the checkpoint holds ', got
     stop 1
  endif
  end subroutine expect_int

end module checkpoint_mod

!#####################################################################
! F77-callable hook for vsup.f, same pattern as analysis_mod's
! histo_end_iteration: the fixed-form integrator must not have to know
! how this compiler mangles a module procedure name.
!#####################################################################

subroutine checkpoint_write_state()
  use checkpoint_mod, only : checkpoint_snapshot
  implicit none
  call checkpoint_snapshot()
end subroutine checkpoint_write_state
