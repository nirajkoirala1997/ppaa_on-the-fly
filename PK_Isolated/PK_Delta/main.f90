program PK_delta 
  use iso_fortran_env
  use globals_mod
  use printframes
  use constants_mod

! Local variables
  character(len=50)  ::  headline,headline1,headline2
  character(len=100) ::  command,run_tag,dir_path,filename,filename1
  integer            ::  iprint
  real(real64) :: PKPlus,err_plus,chi2,ai_lo2,sd
  real(real64), allocatable :: PKDel(:), err_Del(:)
  external fnlo2_PKDel 

      !input data card
      open(unit=10,file='../../run.vegas.dat',status='unknown')
      do i=1,6     
        read (10,*)
      enddo

      read (10,*) pt_delta,pt_delta           ! vegas points     
      read (10,*) its_delta,its_delta          ! vegas iterations 
      npt_delta = pt_delta
      close(10)
      
      open(unit=10,file='../../param_card.dat',status='unknown')    
      read (10,*) ge          ! [ 1/Alpha_ew ]
      close(10)

      open(unit=15,file='../../run.machine.dat',status='unknown')
      read (15,*) mid           ! machine id Tevatron:0 LHC:1
      read (15,*) ecm           ! ecm
      read (15,*) pdf_name_NLO          ! lhapdf set
      read (15,*) it_max        ! it_max no of q for distribution
      read (15,*) xq            ! initialise xq value
      read (15,*) xincr         ! increment in Gev from xq 
      read (15,*) run_tag
      read (15,*) iprint        ! to print data in file
      close(15)

! ~~~~~~~~~~~~~~~~[Writing in a file to store]~~~~~~~~~~~~~~~~~~~c        

      open(unit=20,file='../../output_files.dat',status='unknown')
      do i = 1,7
      read (20,*)
      enddo
      read (20,*) filename
      close(20)
      if (iprint == 1) call output(run_tag,filename)            
      dir_path = '../summary/'//trim(run_tag)//'/'//trim(filename)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c        
      delta = delta_cut 

      call init_mass
      call initpdfsetbyname(pdf_name_NLO)
      Call initPDF(0)
      call ol_LO_init
      allocate(PKDel(it_max), err_Del(it_max))
      s=ecm*ecm
      xq_initial = xq
      

!      headline = "  PK delta     "
!      call printframe0(headline)

      call printframe1(pt_delta,its_delta)   ! Prints Vegas points
      
      do j=1,it_max
       call printframe2(xq)

       call brm48i(40,0,0) 
       call vsup(3,npt_delta,its_delta,fnlo2_PKDel,ai_lo2,sd,chi2)

       PKDel(j) = ai_lo2
       err_Del(j) = sd

       headline = "Delta"
       call printframe3(headline,ai_lo2,sd,chi2)   
     call printframe0(headline)
     write(*,'(i7,e27.15,e27.15)')int(xq),ai_lo2,sd
       xq = xq + xincr
      enddo


! ~~~~~~~~~~~Writing in a file to compare~~~~~~~~~~~~c        
       if (iprint == 1 ) then
         open(unit=21,file=dir_path,status='unknown')
           xq = xq_initial
           do j=1, it_max
             write(21,*)int(xq),PKDel(j),err_Del(j)
             xq = xq + xincr
           enddo
         close(21)
       endif
end program PK_delta
! ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
