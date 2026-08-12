program intPK_Plus
  use iso_fortran_env
  use globals_mod
  use printframes
  use constants_mod

! Local variables
  character(len=50)  ::  headline,headline1,headline2
  character(len=100) ::  run_tag, dir_path, filename, filename1, run_tag1
  integer            ::  iprint
  real(real64) :: chi2,ai_lo2_a,ai_lo2_h,sd
  real(real64) :: ai_lo2A,ai_lo2B,sdA,sdB
  real(real64), allocatable :: PKPlus(:),err_Plus(:)
  external xint_PlusA,xint_Plus_h,flo2_PlusA,flo2_PlusB
  
      !input data card
      open(unit=10,file='../../run.vegas.dat',status='unknown')
      do i=1,12     
        read (10,*)
      enddo
      read (10,*) pt_plus, pt_plus2           ! vegas points     
      read (10,*) its_plus, its_plus2          ! vegas iterations 
      npt_plus = pt_plus
      npt_plus2 = pt_plus2
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
      do i = 1,5
      read (20,*)
      enddo
      read (20,*) filename
      close(20)
      if (iprint .eq. 1) call output(run_tag,filename)            
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~c        
      dir_path = '../summary/'//trim(run_tag)//'/'//trim(filename)

      call init_mass
      call initpdfsetbyname(pdf_name_NLO)
      Call initPDF(0)
      call ol_LO_init
       s=ecm*ecm

     allocate(PKPlus(it_max),err_Plus(it_max))
      xq_initial = xq 
      
      do j=1,it_max
      call printframe2(xq)

!      headline = "     [+] distribution     "
!      call printframe0(headline)
!... ...... Plus implementation type

      if (iselect_integrand == 1 ) then

!      headline1 = "[+] with integrand level subtraction"
!      call printframe1(pt1,its1)   ! Prints Vegas points
!      call printframe0(headline1)
!
!!      -------------------------------------------------
!         call brm48i(40,0,0) 
!         call vsup(2,npt1,its1,flo2_Plus,ai_lo2,sd,chi2)
!!      -------------------------------------------------
!
!         PKPlus   = ai_lo2
!         err_plus = sd

     continue


      elseif (iselect_modified == 1 ) then
      headline1 = "[+] with intXh(x) "

        !==================================================!
         call printframe4('[+] with integrated components')
         call printframe1(pt_plus,its_plus)   ! Prints Vegas points
         call brm48i(40,0,0) 
         call vsup(4,npt_plus,its_plus,xint_PlusA,ai_lo2_a,sdA,chi2)
        !==================================================!
        !==================================================!
         call printframe1(pt_plus2,its_plus2)   ! Prints Vegas points
         call brm48i(40,0,0) 
         call vsup(3,npt_plus2,its_plus2,xint_Plus_h,ai_lo2_h,sdB,chi2)
        !==================================================!

         PKPlus(j)   = ai_lo2_a + ai_lo2_h 
         err_plus(j) = dsqrt(sdA**2 + sdB**2)

      headline = "Plus"
      call printframe3(headline,PKPlus(j),err_plus(j),chi2)   

       elseif (iselect_integral == 1 ) then 

        continue

         headline1 = "PlusA distribution"
         headline2 = "PlusB distribution"
         call printframe1(pt_plus,its_plus)   ! Prints Vegas points

        !==================================================!
         call printframe4(headline1)
         call brm48i(40,0,0)
         call vsup(4,npt_plus,its_plus,flo2_PlusA,ai_lo2A,sdA,chi2)
        !==================================================!
        !==================================================!
         call printframe4(headline2)
         call brm48i(40,0,0)
         call vsup(4,npt_plus,its_plus,flo2_PlusB,ai_lo2B,sdB,chi2)
        !==================================================!

         PKPlus(j) = (ai_lo2A - ai_lo2B)
         err_plus(j) = dsqrt(sdA**2 + sdB**2)

      endif
      
           xq = xq + xincr
      enddo
      xq = xq_initial
      do i=1,it_max
           write(*,*)int(xq),PKPlus(i),err_Plus(i)
           xq = xq + xincr
      enddo


! ~~~~~~~~~~~Writing in a file to compare~~~~~~~~~~~~c        
       if (iprint == 1 ) then
         open(unit=21,file=dir_path,status='unknown')
           xq = xq_initial
           do i=1,it_max
               write(21,*)int(xq),PKPlus(i),err_Plus(i)
               xq = xq + xincr
           enddo
         close(21)
       endif
end program intPK_Plus  
! ~~~~~~~~~~~----------------------------~~~~~~~~~~~~c        
