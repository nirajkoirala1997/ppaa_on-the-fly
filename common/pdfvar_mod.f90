!#####################################################################
! PDF member variation grid.
!
! The second uncertainty axis, alongside scales_mod. That module fixes
! the (muR,muF) grid at compile time; this one cannot, because the
! number of members only becomes known once LHAPDF has loaded the set.
! So the two are combined here into a single flat variation index that
! every consumer -- the weight vector, the histogram storage, the
! output columns -- runs over:
!
!   1 .. nscale             the 7-point scale grid, at the card's member
!   nscale+1 .. nvar        the central scale, at PDF member 0 .. npdf
!
! Member 0 is deliberately kept in the second block even though it
! duplicates the central scale entry: with `member = 0` the two must
! agree to the last digit, which is a free end-to-end check that the
! two blocks really do share a central.
!
! WHY THIS IS CHEAP. The squared matrix element sees the PDF set only
! through alpha_s (this build's own coupling-ratio replay in each
! integrand). So a member whose alpha_s matches the central one reuses
! the central matrix element exactly, and costs nothing but two PDF
! evaluations. Only genuine alpha_s members (the trailing pair of a
! `_pdfas` set) change alpha_s and are covered by the same replay --
! no fresh amplitude call is needed for them either, since the Born's
! only muR/alpha_s dependence is the overall alpha_s^igs power.
!
! Ported near-verbatim from pp2VH_mg5_gg/src/core/pdfvar_mod.f90 (branch
! nlo, HEAD a240659).
!#####################################################################

module pdfvar_mod

  use iso_fortran_env, only: real64
  use scales_mod, only: nscale
  implicit none
  private
  public :: pdfvar_init, pdfvar_select, pdfvar_restore
  public :: pdfvar_active, pdfvar_nmem, pdfvar_central, pdfvar_set
  public :: nvar, MAXVAR, pdfvar_index, pdfvar_label

  ! Upper bound on the flat variation count, in the style of the other
  ! MAX* limits in this code: a static ceiling keeps the per-event weight
  ! array and the run accumulators on the stack, while the histograms
  ! (which dominate the memory) allocate the exact nvar they need.
  integer, parameter :: MAXVAR = 128

  ! Set once by pdfvar_init and read-only thereafter. `nvar` is exported
  ! because every array in the chain is dimensioned by it, but protected
  ! so only this module can set it.
  integer, protected :: nvar = nscale  ! total variations; nscale when PDF variation is off
  integer :: npdf = 0             ! largest member index, i.e. NumMembers-1
  integer :: icentral_mem = 0     ! the member the scale grid is evaluated at
  logical :: active = .false.
  character(len=128) :: setname = ''

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Load the set, select the central member, and size the variation grid.
! Replaces a bare initpdfsetbyname/initPDF pair -- every run goes
! through here so that `nvar` is always consistent with what LHAPDF
! actually has loaded.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine pdfvar_init(set_name, central_member, on)
  character(len=*), intent(in) :: set_name
  integer,          intent(in) :: central_member
  logical,          intent(in) :: on

  integer :: nmem

  setname      = set_name
  icentral_mem = central_member
  active       = on

  ! Without this every member prints a banner the first time it is
  ! touched, so the first event of a 43-member run would bury its log
  ! under 43 blocks of LHAPDF boilerplate -- once per parallel seed.
  if (active) call setlhaparm('SILENT')

  call initpdfsetbyname(set_name)

  ! LHAPDF5 convention, kept by LHAPDF6's compatibility layer: this is
  ! the largest member INDEX, one less than NumMembers. Members are
  ! therefore 0..npdf inclusive, and a 1..npdf loop would quietly drop
  ! the last one -- which for a _pdfas set is an alpha_s variation.
  call numberPDF(nmem)
  npdf = nmem

  if (central_member < 0 .or. central_member > npdf) then
     print '(A,I0,A,A,A,I0,A)','ERROR pdfvar_init: member ',central_member, &
           ' is outside set ',trim(set_name),' (0..',npdf,')'
     stop 1
  endif

  if (active) then
     ! With variations on, the scale block and the PDF block have to share
     ! one central value, or the two bands are quoted around different
     ! numbers and cannot be combined.
     if (central_member /= 0) then
        print '(A,I0,A)','ERROR pdfvar_init: [pdf] variations = yes requires member = 0, got ', &
              central_member,'.'
        print '(A)','       Members 1.. are error directions, not the central PDF, so the'
        print '(A)','       uncertainty would be measured around the wrong central value.'
        stop 1
     endif
     nvar = nscale + npdf + 1
     if (nvar > MAXVAR) then
        print '(A,I0,A,I0,A)','ERROR pdfvar_init: nvar = ',nvar, &
              ' exceeds MAXVAR = ',MAXVAR,' -- raise MAXVAR in pdfvar_mod.f90'
        stop 1
     endif
  else
     nvar = nscale
  endif

  call initPDF(icentral_mem)

  end subroutine pdfvar_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Flat index of PDF member `imem`. Kept here so the layout is defined
! in exactly one place; the output header advertises the same mapping.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer function pdfvar_index(imem) result(iv)
  integer, intent(in) :: imem
  iv = nscale + 1 + imem
  end function pdfvar_index

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine pdfvar_select(imem)
  integer, intent(in) :: imem
  call initPDF(imem)
  end subroutine pdfvar_select

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Put the central member back. Every routine that walks the members
! must end with this: everything outside the variation loop (and the
! next event's scale block) assumes the central member is the active one.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine pdfvar_restore()
  call initPDF(icentral_mem)
  end subroutine pdfvar_restore

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  logical function pdfvar_active() result(val)
  val = active
  end function pdfvar_active

  integer function pdfvar_nmem() result(val)
  val = npdf
  end function pdfvar_nmem

  integer function pdfvar_central() result(val)
  val = icentral_mem
  end function pdfvar_central

  function pdfvar_set() result(val)
  character(len=128) :: val
  val = setname
  end function pdfvar_set

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Machine-readable description of the flat index, for the report header.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function pdfvar_label() result(lbl)
  character(len=64) :: lbl
  if (active) then
     write(lbl,'(A,I0,A,I0)') 'scale:1-',nscale,' pdfmem:0-',npdf
  else
     write(lbl,'(A,I0)') 'scale:1-',nscale
  endif
  end function pdfvar_label

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module pdfvar_mod

!#####################################################################
