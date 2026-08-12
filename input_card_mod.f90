!#####################################################################
! Run card parser.
!
! Reads a sectioned, keyed input file:
!
!     [section]
!     key = value          ! trailing comment
!
! and, in the [histograms] section only, repeated definition lines:
!
!     name  obs = pt   leg = 3   nbins = 30   xmin = 0.0   xmax = 600.0
!
! WHY KEYED. The previous format was positional: values were read in a
! fixed order with no names, so inserting, removing or reordering a line
! silently shifted every subsequent value into the wrong variable. Here a
! value is bound to its key, order does not matter, and three classes of
! mistake are caught explicitly rather than producing wrong numbers:
!
!   * an unknown key            -> reported by card_report_unused
!   * a missing required key    -> reported by the card_require_* getters
!   * a duplicate key           -> reported at parse time
!
! Comments start with '!' or '#'. Blank lines are ignored. Values may be
! quoted if they contain spaces.
!#####################################################################

module input_card_mod

  use iso_fortran_env, only: real64
  implicit none
  private
  public :: card_load, card_report_unused, card_dump
  public :: card_get_int, card_get_real, card_get_str, card_get_logical
  public :: card_require_int, card_require_real, card_require_str
  public :: card_nhist, card_hist_name
  public :: card_hist_int, card_hist_real, card_hist_str, card_hist_intlist

  integer, parameter :: MAXENT  = 256
  integer, parameter :: LSEC    = 32
  integer, parameter :: LKEY    = 32
  integer, parameter :: LVAL    = 128

  type :: entry_type
     character(len=LSEC) :: section = ''
     character(len=LKEY) :: name    = ''   ! non-empty only for [histograms] rows
     character(len=LKEY) :: key     = ''
     character(len=LVAL) :: value   = ''
     integer :: iline = 0
     logical :: used  = .false.
  end type entry_type

  type(entry_type), save :: ents(MAXENT)
  integer, save :: nents = 0

  ! Distinct histogram names, in file order.
  character(len=LKEY), save :: hnames(MAXENT)
  integer, save :: nhist = 0

  character(len=256), save :: card_file = ''

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine card_load(filename)
  character(len=*), intent(in) :: filename

  character(len=512) :: line, raw
  character(len=LSEC) :: section
  integer :: iu, ios, iline, ic, ieq
  character(len=LKEY) :: hname

  card_file = filename
  nents = 0
  nhist = 0
  section = ''

  iu = 91
  open(unit=iu, file=filename, status='old', action='read', iostat=ios)
  if (ios /= 0) then
     print*,'ERROR input card: cannot open ',trim(filename)
     stop 1
  endif

  iline = 0
  do
     read(iu,'(A)',iostat=ios) raw
     if (ios /= 0) exit
     iline = iline + 1

     line = raw
     ! strip comments
     ic = index(line,'!') ; if (ic > 0) line = line(1:ic-1)
     ic = index(line,'#') ; if (ic > 0) line = line(1:ic-1)
     line = adjustl(line)
     if (len_trim(line) == 0) cycle

     ! section header
     if (line(1:1) == '[') then
        ic = index(line,']')
        if (ic < 3) then
           print '(A,I0,A)','ERROR input card line ',iline,': malformed section header'
           stop 1
        endif
        section = lower(trim(adjustl(line(2:ic-1))))
        cycle
     endif

     if (len_trim(section) == 0) then
        print '(A,I0,A)','ERROR input card line ',iline,': entry before any [section]'
        stop 1
     endif

     if (trim(section) == 'histograms') then
        ! <name> key=value key=value ...
        call split_first_token(line, hname)
        if (len_trim(hname) == 0) then
           print '(A,I0,A)','ERROR input card line ',iline,': histogram row has no name'
           stop 1
        endif
        call register_hist(hname)
        call parse_kv_list(line, section, hname, iline)
     else
        ieq = index(line,'=')
        if (ieq <= 1) then
           print '(A,I0,A)','ERROR input card line ',iline, &
                 ': expected "key = value" in section ['//trim(section)//']'
           print '(A)','   got: '//trim(raw)
           stop 1
        endif
        call add_entry(section, '', lower(trim(adjustl(line(1:ieq-1)))), &
                       dequote(trim(adjustl(line(ieq+1:)))), iline)
     endif
  end do
  close(iu)

  if (nents == 0) then
     print*,'ERROR input card: ',trim(filename),' contains no entries'
     stop 1
  endif

  end subroutine card_load

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Parse "key = value key = value ..." from the remainder of a line.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine parse_kv_list(rest, section, hname, iline)
  character(len=*), intent(in) :: rest, section, hname
  integer,          intent(in) :: iline

  character(len=512) :: buf
  character(len=LKEY) :: k
  character(len=LVAL) :: v
  integer :: ieq, i, n

  buf = adjustl(rest)
  do
     if (len_trim(buf) == 0) exit
     ieq = index(buf,'=')
     if (ieq <= 1) then
        print '(A,I0,A)','ERROR input card line ',iline, &
              ': expected "key = value" pairs after histogram name'
        stop 1
     endif
     k = lower(trim(adjustl(buf(1:ieq-1))))
     buf = adjustl(buf(ieq+1:))
     ! value runs to the next whitespace
     n = len_trim(buf)
     i = index(buf(1:max(n,1)),' ')
     if (i == 0) then
        v = dequote(trim(buf))
        buf = ''
     else
        v = dequote(trim(buf(1:i-1)))
        buf = adjustl(buf(i:))
     endif
     call add_entry(section, hname, k, v, iline)
  end do

  end subroutine parse_kv_list

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine add_entry(section, name, key, value, iline)
  character(len=*), intent(in) :: section, name, key, value
  integer,          intent(in) :: iline
  integer :: i

  if (len_trim(key) == 0) then
     print '(A,I0,A)','ERROR input card line ',iline,': empty key'
     stop 1
  endif

  do i = 1, nents
     if (trim(ents(i)%section) == trim(section) .and. &
         trim(ents(i)%name)    == trim(name)    .and. &
         trim(ents(i)%key)     == trim(key)) then
        print '(A,I0,A,I0)','ERROR input card line ',iline, &
              ': duplicate key "'//trim(key)//'" in ['//trim(section)// &
              '], first defined on line ',ents(i)%iline
        stop 1
     endif
  end do

  if (nents >= MAXENT) then
     print*,'ERROR input card: too many entries, increase MAXENT'
     stop 1
  endif

  nents = nents + 1
  ents(nents)%section = section
  ents(nents)%name    = name
  ents(nents)%key     = key
  ents(nents)%value   = value
  ents(nents)%iline   = iline
  ents(nents)%used    = .false.

  end subroutine add_entry

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine register_hist(hname)
  character(len=*), intent(in) :: hname
  integer :: i
  do i = 1, nhist
     if (trim(hnames(i)) == trim(hname)) return
  end do
  nhist = nhist + 1
  hnames(nhist) = hname
  end subroutine register_hist

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Lookup. Returns the index of a matching entry, or 0, and marks it used.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer function find(section, name, key) result(idx)
  character(len=*), intent(in) :: section, name, key
  integer :: i
  character(len=LSEC) :: s
  character(len=LKEY) :: k
  ! Keys and sections are stored lowercased, so fold the query too --
  ! otherwise a caller asking for 'mZ' would silently miss the stored
  ! 'mz' and fall back to its default.
  s = lower(trim(section))
  k = lower(trim(key))
  idx = 0
  do i = 1, nents
     if (trim(ents(i)%section) == trim(s)    .and. &
         trim(ents(i)%name)    == trim(name) .and. &
         trim(ents(i)%key)     == trim(k)) then
        ents(i)%used = .true.
        idx = i
        return
     endif
  end do
  end function find

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Typed getters. The card_get_* forms take a default; the card_require_*
! forms abort when the key is absent.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer function card_get_int(section, key, default) result(val)
  character(len=*), intent(in) :: section, key
  integer,          intent(in) :: default
  integer :: i, ios
  val = default
  i = find(section,'',key)
  if (i == 0) return
  read(ents(i)%value,*,iostat=ios) val
  if (ios /= 0) call bad_value(i,'integer')
  end function card_get_int

  real(real64) function card_get_real(section, key, default) result(val)
  character(len=*), intent(in) :: section, key
  real(real64),     intent(in) :: default
  integer :: i, ios
  val = default
  i = find(section,'',key)
  if (i == 0) return
  read(ents(i)%value,*,iostat=ios) val
  if (ios /= 0) call bad_value(i,'real')
  end function card_get_real

  function card_get_str(section, key, default) result(val)
  character(len=*), intent(in) :: section, key, default
  character(len=LVAL) :: val
  integer :: i
  val = default
  i = find(section,'',key)
  if (i == 0) return
  val = ents(i)%value
  end function card_get_str

  logical function card_get_logical(section, key, default) result(val)
  character(len=*), intent(in) :: section, key
  logical,          intent(in) :: default
  integer :: i
  character(len=LVAL) :: v
  val = default
  i = find(section,'',key)
  if (i == 0) return
  v = lower(trim(ents(i)%value))
  select case (trim(v))
    case ('yes','true','t','1','on')  ; val = .true.
    case ('no','false','f','0','off') ; val = .false.
    case default ; call bad_value(i,'logical (yes/no)')
  end select
  end function card_get_logical

  integer function card_require_int(section, key) result(val)
  character(len=*), intent(in) :: section, key
  integer :: i, ios
  i = find(section,'',key)
  if (i == 0) call missing(section,key)
  read(ents(i)%value,*,iostat=ios) val
  if (ios /= 0) call bad_value(i,'integer')
  end function card_require_int

  real(real64) function card_require_real(section, key) result(val)
  character(len=*), intent(in) :: section, key
  integer :: i, ios
  i = find(section,'',key)
  if (i == 0) call missing(section,key)
  read(ents(i)%value,*,iostat=ios) val
  if (ios /= 0) call bad_value(i,'real')
  end function card_require_real

  function card_require_str(section, key) result(val)
  character(len=*), intent(in) :: section, key
  character(len=LVAL) :: val
  integer :: i
  i = find(section,'',key)
  if (i == 0) call missing(section,key)
  val = ents(i)%value
  end function card_require_str

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Histogram accessors.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer function card_nhist() result(n)
  n = nhist
  end function card_nhist

  function card_hist_name(ih) result(nm)
  integer, intent(in) :: ih
  character(len=LKEY) :: nm
  nm = hnames(ih)
  end function card_hist_name

  integer function card_hist_int(ih, key, default) result(val)
  integer,          intent(in) :: ih, default
  character(len=*), intent(in) :: key
  integer :: i, ios
  val = default
  i = find('histograms', hnames(ih), key)
  if (i == 0) return
  read(ents(i)%value,*,iostat=ios) val
  if (ios /= 0) call bad_value(i,'integer')
  end function card_hist_int

  real(real64) function card_hist_real(ih, key, default) result(val)
  integer,          intent(in) :: ih
  character(len=*), intent(in) :: key
  real(real64),     intent(in) :: default
  integer :: i, ios
  val = default
  i = find('histograms', hnames(ih), key)
  if (i == 0) return
  read(ents(i)%value,*,iostat=ios) val
  if (ios /= 0) call bad_value(i,'real')
  end function card_hist_real

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! A comma-separated list of integers, e.g.  legs = 3,4,5,6
!
! Returns nvals = 0 when the key is absent, so a caller can tell an
! omitted list from an empty one. A value runs to the next whitespace
! (parse_kv_list), so the list must not contain spaces -- which is also
! what keeps it distinguishable from the next key on the row.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine card_hist_intlist(ih, key, vals, nvals)
  integer,          intent(in)  :: ih
  character(len=*), intent(in)  :: key
  integer,          intent(out) :: vals(:)
  integer,          intent(out) :: nvals

  character(len=LVAL) :: buf
  integer :: i, ic, ios

  nvals = 0
  i = find('histograms', hnames(ih), key)
  if (i == 0) return

  buf = adjustl(ents(i)%value)
  do
     if (len_trim(buf) == 0) exit
     if (nvals >= size(vals)) then
        print '(A,I0,A)','ERROR input card line ',ents(i)%iline, &
              ': too many entries in "'//trim(key)//'"'
        stop 1
     endif
     ic = index(buf,',')
     nvals = nvals + 1
     if (ic == 0) then
        read(buf,*,iostat=ios) vals(nvals)
        buf = ''
     else
        read(buf(1:ic-1),*,iostat=ios) vals(nvals)
        buf = adjustl(buf(ic+1:))
     endif
     if (ios /= 0) call bad_value(i,'comma-separated integer list')
  end do

  end subroutine card_hist_intlist

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function card_hist_str(ih, key, default) result(val)
  integer,          intent(in) :: ih
  character(len=*), intent(in) :: key, default
  character(len=LVAL) :: val
  integer :: i
  val = default
  i = find('histograms', hnames(ih), key)
  if (i == 0) return
  val = ents(i)%value
  end function card_hist_str

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Call after every getter has run. Any entry nobody asked for is almost
! certainly a typo or a stale field, so it is an error rather than a
! warning -- silently ignoring it is exactly the failure mode the keyed
! format exists to prevent.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine card_report_unused()
  integer :: i, nbad
  nbad = 0
  do i = 1, nents
     if (ents(i)%used) cycle
     nbad = nbad + 1
     if (nbad == 1) print '(A)','ERROR input card: unrecognised entries in '//trim(card_file)//':'
     if (len_trim(ents(i)%name) > 0) then
        print '(A,I0,A)','   line ',ents(i)%iline,'  ['//trim(ents(i)%section)//'] '// &
              trim(ents(i)%name)//'  '//trim(ents(i)%key)//' = '//trim(ents(i)%value)
     else
        print '(A,I0,A)','   line ',ents(i)%iline,'  ['//trim(ents(i)%section)//'] '// &
              trim(ents(i)%key)//' = '//trim(ents(i)%value)
     endif
  end do
  if (nbad > 0) then
     print '(A)','   (check spelling; every key must be one the code reads)'
     stop 1
  endif
  end subroutine card_report_unused

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine card_dump(iunit)
  integer, intent(in) :: iunit
  integer :: i
  write(iunit,'(A)') '# ==== INPUT CARD ===='
  do i = 1, nents
     if (len_trim(ents(i)%name) > 0) then
        write(iunit,'(A)') '# ['//trim(ents(i)%section)//'] '//trim(ents(i)%name)//' '// &
              trim(ents(i)%key)//' = '//trim(ents(i)%value)
     else
        write(iunit,'(A)') '# ['//trim(ents(i)%section)//'] '// &
              trim(ents(i)%key)//' = '//trim(ents(i)%value)
     endif
  end do
  write(iunit,'(A)') '# ==== END INPUT CARD ===='
  write(iunit,*)
  end subroutine card_dump

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! helpers
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine missing(section, key)
  character(len=*), intent(in) :: section, key
  print '(A)','ERROR input card '//trim(card_file)//': required key "'//trim(key)// &
        '" is missing from section ['//trim(section)//']'
  stop 1
  end subroutine missing

  subroutine bad_value(i, want)
  integer,          intent(in) :: i
  character(len=*), intent(in) :: want
  print '(A,I0,A)','ERROR input card line ',ents(i)%iline,': "'//trim(ents(i)%key)// &
        ' = '//trim(ents(i)%value)//'" is not a valid '//trim(want)
  stop 1
  end subroutine bad_value

  subroutine split_first_token(line, tok)
  character(len=*), intent(inout) :: line
  character(len=*), intent(out)   :: tok
  integer :: i
  line = adjustl(line)
  i = index(trim(line),' ')
  if (i == 0) then
     tok  = trim(line)
     line = ''
  else
     tok  = line(1:i-1)
     line = adjustl(line(i:))
  endif
  end subroutine split_first_token

  function dequote(s) result(o)
  character(len=*), intent(in) :: s
  character(len=LVAL) :: o
  integer :: n
  o = s
  n = len_trim(o)
  if (n >= 2) then
     if ((o(1:1) == '"'  .and. o(n:n) == '"') .or. &
         (o(1:1) == "'"  .and. o(n:n) == "'")) o = o(2:n-1)
  endif
  end function dequote

  function lower(s) result(o)
  character(len=*), intent(in) :: s
  character(len=len(s)) :: o
  integer :: i, ic
  o = s
  do i = 1, len_trim(o)
     ic = iachar(o(i:i))
     if (ic >= 65 .and. ic <= 90) o(i:i) = achar(ic+32)
  end do
  end function lower

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module input_card_mod

!#####################################################################
