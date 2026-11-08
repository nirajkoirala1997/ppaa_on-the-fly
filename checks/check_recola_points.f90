! Prints Recola's Born, one-loop and real values at the fixed phase-space points
! in the file named on the command line. The file holds, per point: muR, then
! the 4 Born momenta, then the 5 real-emission momenta, each as E px py pz.
!
!     ./check_recola_points points.dat > recola.out
program check_recola_points
  use iso_fortran_env, only: real64
  use run_card_mod,    only: read_input_file
  use olp_mod,         only: olp_init, olp_set_alphas, olp_born, olp_virtual, &
                             olp_real, NFLAVOUR, NREAL_CHANNEL
  implicit none

  real(real64), parameter :: ALPHAS_FIXED = 0.118d0
  real(real64) :: renorm_scale, born_momenta(0:3,4), real_momenta(0:3,5)
  real(real64) :: born, born_tree, finite, pole1, pole2, real_value
  integer :: npoints, ipoint, flavour, channel, iunit, ileg
  character(len=256) :: points_file

  call get_command_argument(1, points_file)
  call read_input_file('../input.input', 'lo')
  call olp_init(ALPHAS_FIXED, 100d0, .true., .true.)

  open(newunit=iunit, file=trim(points_file), status='old')
  read(iunit,*) npoints
  do ipoint = 1, npoints
     read(iunit,*) renorm_scale
     do ileg = 1, 4
        read(iunit,*) born_momenta(:,ileg)
     end do
     do ileg = 1, 5
        read(iunit,*) real_momenta(:,ileg)
     end do

     call olp_set_alphas(ALPHAS_FIXED, renorm_scale)
     do flavour = 1, NFLAVOUR
        call olp_born(flavour, born_momenta, born_tree)
        call olp_virtual(flavour, born_momenta, born, finite, pole1, pole2)
        print '(A,2I4,ES25.16)', 'RESULT born ', ipoint, flavour, born_tree
        print '(A,2I4,ES25.16)', 'RESULT virt ', ipoint, flavour, finite
        print '(A,2I4,ES25.16)', 'RESULT pole1', ipoint, flavour, pole1
        print '(A,2I4,ES25.16)', 'RESULT pole2', ipoint, flavour, pole2
     end do
     do channel = 1, NREAL_CHANNEL
        call olp_real(channel, real_momenta, real_value)
        print '(A,2I4,ES25.16)', 'RESULT real ', ipoint, channel, real_value
     end do
  end do
  close(iunit)
end program check_recola_points
