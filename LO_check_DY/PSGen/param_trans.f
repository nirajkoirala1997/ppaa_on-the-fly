      subroutine param_trans(unit)

      use kinparams

      implicit none

      integer :: unit

c      aw=1.d0/nkin   ! 4 May 2023
c      iwkin=0
c      aw_w=0.d0

      if (iqprec == 1) then
       call qp_masses
      end if
      call exp_masses
      call exp_names
      call define_cuts
      call set_indices
      call logo_ps(unit)
      call cutinfo(unit)
      call print_masses(unit)

      end
