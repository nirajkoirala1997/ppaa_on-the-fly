program PK_regular 
  use iso_fortran_env
  use globals_mod
  use printframes
  use constants_mod

! Local variables
  character(len=50)  ::  headline,headline1,headline2
  character(len=100) ::  command,run_tag,dir_path,filename,filename1
  integer            ::  iprint
  real(real64) :: PKPlus,err_plus,chi2,ai_lo2,sd
  real(real64), allocatable :: PKReg(:), err_Reg(:)
  external fnlo2_PKReg 

!      external flo2_Plus,flo2_PKDel,flo2_PKReg,xint_PlusA,xint_Plus_h
!      external flo2_PlusA,flo2_PlusB

      !input data card
      open(unit=10,file='../../run.vegas.dat',status='unknown')
      do i=1,6     
        read (10,*)
      enddo
      read (10,*) pt_regular           ! vegas points     
      read (10,*) its_regular          ! vegas iterations 
      npt_regular = pt_regular
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
      do i = 1,6
      read (20,*)
      enddo
      read (20,*) filename
      close(20)
      if (iprint == 1) call output(run_tag,filename)            
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c        
      dir_path = '../summary/'//trim(run_tag)//'/'//trim(filename)
      delta = delta_cut 

      call init_mass
      call initpdfsetbyname(pdf_name_NLO)
      Call initPDF(0)
      call ol_LO_init
       s=ecm*ecm
      allocate(PKReg(it_max), err_Reg(it_max))
      xq_initial = xq

!      headline = "  PK regular     "
!      call printframe0(headline)

      call printframe1(pt_regular,its_regular)   ! Prints Vegas points

    do j=1,it_max
      call printframe2(xq)

      call brm48i(40,0,0) 
      call vsup(4,npt_regular,its_regular,fnlo2_PKReg,ai_lo2,sd,chi2)

      PKReg(j) = ai_lo2
      err_Reg(j) = sd

      headline = "Regular"
      call printframe3(headline,ai_lo2,sd,chi2)   
      xq = xq + xincr
     enddo

     xq = xq_initial
     do j = 1,it_max
        write(*,*)int(xq),PKReg(j),err_Reg(j)
        xq = xq + xincr
     enddo


! ~~~~~~~~~~~Writing in a file to compare~~~~~~~~~~~~c        
       if (iprint == 1 ) then
         open(unit=21,file=dir_path,status='unknown')
           xq = xq_initial
           do j = 1,it_max
              write(21,*)int(xq),PKReg(j),err_Reg(j)
              xq = xq + xincr
           enddo
         close(21)
       endif
end program PK_regular
! ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
