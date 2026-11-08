! PSGen's own 3-body measure against the massless closed form, with and without
! the m(ee) window. Int dlips over uniform xr must equal
!   Phi_3 = Int dm2/(2 pi) (s-m2)/(8 pi s) * 1/(8 pi)  over the m(ee)^2 range.
program check_ps_volume
  use iso_fortran_env, only : real64
  use psgen_interface, only : psgen_init, psgen_event, psgen_ndim, psgen_nexp
  implicit none
  real(real64), parameter :: PI = 3.141592653589793238462643d0, ECM = 13600d0
  integer, parameter :: NPOINT = 2000000
  real(real64) :: xbeam, shat, window_low, window_high, xr(psgen_ndim)
  real(real64) :: plab(0:3,psgen_nexp), flux, dlips, sum_w, sum_w2, mean, error, exact
  integer :: i, iarg
  character(len=32) :: arg

  xbeam = 0.05d0
  window_low = 0d0 ; window_high = 0d0
  do iarg = 1, command_argument_count()
     call get_command_argument(iarg, arg)
     if (iarg == 1) read(arg,*) window_low
     if (iarg == 2) read(arg,*) window_high
  end do
  call psgen_init(6, 1d-8, window_low, window_high)
  shat = xbeam*xbeam*ECM*ECM
  sum_w = 0d0 ; sum_w2 = 0d0
  do i = 1, NPOINT
     call random_number(xr)
     call psgen_event(ECM, xbeam, xbeam, 1, xr, plab, flux, dlips)
     sum_w = sum_w + dlips ; sum_w2 = sum_w2 + dlips*dlips
  end do
  mean  = sum_w/NPOINT
  error = sqrt((sum_w2/NPOINT - mean*mean)/NPOINT)
  exact = primitive(merge(shat, min(window_high**2,shat), window_high <= 0d0)) &
        - primitive(window_low**2)
  write(*,'(A,F10.2,2F10.2)') 'RESULT sqrt(shat), window: ', sqrt(shat), window_low, window_high
  write(*,'(A,ES18.10,A,ES10.3,A,ES18.10,A,F12.8,A,F8.2)') 'RESULT MC ', mean, ' +- ', error, &
       ' exact ', exact, ' ratio ', mean/exact, ' pull ', (mean-exact)/error
contains
  real(real64) function primitive(m2) result(v)
  real(real64), intent(in) :: m2
  v = (m2 - m2*m2/(2d0*shat))/(2d0*PI)/(8d0*PI)/(8d0*PI)
  end function primitive
end program check_ps_volume
