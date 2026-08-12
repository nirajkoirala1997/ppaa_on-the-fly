      subroutine qp_masses

      use inprms_ps

      implicit none

c Masses in quadruple precision:
      me_qp=qpmass(me)
      mmu_qp=qpmass(mmu)
      mtau_qp=qpmass(mtau)
      mu0_qp=qpmass(mu0)
      md0_qp=qpmass(md0)
      ms0_qp=qpmass(ms0)
      mc0_qp=qpmass(mc0) 
      mb0_qp=qpmass(mb0) 
      mt_qp=qpmass(mt)
      mz_qp=qpmass(mz)
      mw_qp=qpmass(mw)
      mh_qp=qpmass(mh)
      mp_qp=qpmass(mp)
      mn_qp=qpmass(mn)
      mpi_qp=qpmass(mpi)
      mpi0_qp=qpmass(mpi0)
      mk_qp=qpmass(mk)
      mk0_qp=qpmass(mk0)
      mrho_qp=qpmass(mrho)
      mome_qp=qpmass(mome)
      mphi_qp=qpmass(mphi)
      mrh1_qp=qpmass(mrh1)
      mrh2_qp=qpmass(mrh2)

      contains

       function qpmass(mdp)
       real(dp) :: mdp
       real(qp) :: qpmass
       character(80) :: ind

       write(ind,'(e23.16)') mdp
       open(unit=12,status='scratch',action='readwrite')
       write(12,'(a)') trim(ind)
       rewind(12)
       read(12,'(e23.16)') qpmass
       close(12)

       end function qpmass

      end
