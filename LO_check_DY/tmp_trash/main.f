      program uU2eE_LO 
      use openloops
      implicit double precision (a-h,o-z)
      dimension x(10),y(10),ai_lo2(1:50),err(0:50)
      parameter (pi=3.14159265358979d0)

      common/energy/s
      common/amass/am1,am2,am3,am4,am5
      common/usedalpha/AL,ge   
      common/distribution/xq
      common/prc_id/id_LO,id_NLO_1r,id_NLO_1loop
      character*50 name,mode
      character*100 run_tag,filename
      external flo2_LO

      !input data card
      open(unit=10,file='../run.vegas.dat',status='unknown')    
      do i=1,9
      read (10,*)
      enddo
      read (10,*) pt1          ! vegas points     LO 2 body
      read (10,*) its1          ! vegas iterations LO 2 body
      npt1 = pt1
      close(10)

      open(unit=10,file='../param_card.dat',status='unknown')    
      read (10,*) ge       ! [ 1/Alpha_ew ]
      close(10)
      ge = 1d0/128d0

      open(unit=15,file='../run.machine.dat',status='unknown')
      read (15,*) mid           ! machine id Tevatron:0 LHC:1
      read (15,*) ecm           ! ecm
      read (15,*) name          !lhapdf set
      read (15,*) it_max        !lhapdf set
      read (15,*) xq_initial
      read (15,*) step_size         !step
      read (15,*) run_tag           !run_tag saves the output 
      read (15,*) iprint            !save data in output file ../summary
      close(15)


      open(unit=20,file='../output_files.dat',status='unknown')
      read (20,*)
      read (20,*)
      read (20,*)
      read (20,*) filename
      close(20)

c ~~~~~~~~~~~~~~~~[Writing in a file to store]~~~~~~~~~~~~~~~~~~~c        
      open(unit=20,file='../output_files.dat',status='unknown')
      read (20,*) 
      read (20,*) 
      read (20,*) 
      read (20,*) filename
      close(20)
      if(iprint .eq. 1) call output(run_tag,filename)
c ~~~~~~~~~~~~~~~~[--------------------------]~~~~~~~~~~~~~~~~~~~c        


        call initpdfsetbyname(name)
        Call initPDF(0)
        call ol_LO_init(id_LO)
        call ol_NLO_1loop_init(id_NLO_1loop)
      
      am1=0.0d0
      am2=0.0d0
      am3=0d0
      am4=0d0
      am5=0d0
      s=ecm*ecm

        print*," "
        print*," "
        print*,"____________________________________"
        Print*," Calculating LO_DY"
        print*,"____________________________________"
        print*,"````````````````````````````````````"
        print*," "
        print*," "

      call printframe1(pt1,its1)   ! vegas points print 

        xq = xq_initial
        do j=1,it_max

        call printframe2(xq)

         call brm48i(40,0,0) ! initialize random number generator
         call vsup(3,npt1,its1,flo2_LO,ans,sd,chi2)

            ai_lo2(j) = ans
              err(j)  = sd

            mode  = "Leading Order"
            call printframe3(mode,ans,sd,chi2)
         xq = xq + step_size
        enddo

         xq = xq_initial

         call printframe4(mode)

        do j=1,it_max
          write(*,'(i7,3e27.15)')int(xq),ai_lo2(j),err(j)
          xq = xq + step_size
        enddo

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (iprint .eq. 0) goto 123
        open(unit=20,file='../summary/'//trim(run_tag)//'/'
     .  //trim(filename),status='unknown')
         xq = xq_initial
         do i=1,it_max
          write(20,*)xq,ai_lo2(i),err(i)
          xq = xq + step_size
         enddo
         close(20)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       end
