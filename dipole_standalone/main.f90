program Drell_yan_dipoleSubtraction 
  use iso_fortran_env
  use globals_mod
  use misc_utils_mod, only : pt_real,its_real,it_max,real_file,ecm,pdf_NLO,ecm
  use misc_utils_mod, only : invmass_base,step_size,iprint,run_tag,alew,read_input_file
  use printframes
  use constants_mod

! LOCAL PARAMETERS
    real(real64)       :: ans, sd, chi2
    character(len=100) :: mode
    real(real64), allocatable :: ai_nlo3(:), err(:)
! EXTERNAL FUNCTIONS
    external fnlo3

   call read_input_file('../input.input')

      npt_real= pt_real

      s=ecm*ecm
      ge = 1d0/alew
      if(iprint == 1) call output(run_tag,real_file)
! ~~~~~~~~~~~~~~~~[--------------------------]~~~~~~~~~~~~~~~~~~~c        

        call init_mass
        call initpdfsetbyname(pdf_NLO)
        Call initPDF(0)
        call ol_LO_init
        call ol_NLO_1real_init

          print*,"  ----------------------------------"
          print*,"  |Initializing Dipole Subtraction  |"
          print*,"  ----------------------------------"
          print*," "
          print*," "
          print*," "
          call printframe1(pt_real,its_real)
          allocate(ai_nlo3(it_max), err(it_max)) 


          xq = invmass_base
          do j = 1,it_max
            
          call printframe2(xq)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call brm48i(40,0,0) 
            call vsup(6,npt_real,its_real,fnlo3,ans,sd,chi2)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ai_nlo3(j) = ans
               err(j) = sd

            mode  = "[real-dipole]"
            call printframe3(mode,ans,sd,chi2)
            xq = xq + step_size 
          enddo

 ! Print results to stdout
  xq = invmass_base
  call printframe0(mode)
  do j = 1, it_max
     write(*,'(i7,3e27.15)') int(xq), ai_nlo3(j), err(j)
     xq = xq + step_size
  end do

  ! Write to file
  if (iprint == 1) then
     xq = invmass_base
     open(unit=20, file='../summary/'//trim(run_tag)//'/'//trim(real_file), status='unknown')
     do j = 1, it_max
        write(20,*) xq, ai_nlo3(j), err(j)
        xq = xq + step_size
     end do
     close(20)
  end if
end program Drell_yan_dipoleSubtraction
