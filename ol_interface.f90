!c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Masses, widths and the electroweak scheme. One routine, called by every
! registration path below, because Born, real and virtual must run at the same
! electroweak input or their sum is meaningless.
!
! ORDER MATTERS HERE. ew_scheme must be set BEFORE alpha: ew_scheme = 0 is what
! makes OpenLoops accept an alpha at all, and setting alpha first has it
! silently discarded. ew_scheme = 1 (let OpenLoops derive alpha from Gmu itself)
! gives a cross section 0.0969% higher -- see constants_mod::alpha_gfermi.
       subroutine ol_ew_setup
        use constants_mod, only : mass, decayW, alpha_gfermi
        use openloops

        call set_parameter("use_cms", 1)          ! complex-mass scheme
        call set_parameter("ew_scheme", 0)        ! take alpha as given ...
        call set_parameter("alpha", alpha_gfermi())   ! ... and this is it

        call set_parameter("mass(23)",  mass(23))
        call set_parameter("width(23)", decayW(23))
        call set_parameter("mass(24)",  mass(24))
        call set_parameter("width(24)", decayW(24))
        call set_parameter("mass(25)",  mass(25))
        call set_parameter("width(25)", decayW(25))
        call set_parameter("mass(6)",   mass(6))
        call set_parameter("width(6)",  decayW(6))
        call set_parameter("mass(5)",   mass(5))  ! massless b: 5-flavour, PDF only

       return
       end subroutine ol_ew_setup
!c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       subroutine ol_LO_init
        use globals_mod, only : id_LO
        use constants_mod
        use openloops

      
        call set_parameter("order_ew", 2)
        call set_parameter("order_qcd", 0)

        call ol_ew_setup

        call set_parameter("verbose",-1)
        call set_parameter("model","sm")
                                               ! OPENLOOPS MANUAL
!Page 16 :http://link.springer.com/10.1140/epjc/s10052-019-7306-2


       id_LO(2) = register_process("2 -2 -> 11 -11", 1)  ! u ū
       id_LO(1) = register_process("1 -1 -> 11 -11", 1)  ! d d̄

       id_LO(4) = register_process("4 -4 -> 11 -11", 1)  ! c c̄
       id_LO(3) = register_process("3 -3 -> 11 -11", 1)  ! s s̄

       id_LO(5) = register_process("5 -5 -> 11 -11", 1)  ! b b̄

      return
      end subroutine ol_LO_init
   !=============================
       subroutine ol_NLO_1real_init
        use globals_mod, only : id_NLO_1r
        use constants_mod
        use openloops

      
        call set_parameter("order_ew", 2)
        call set_parameter("order_qcd", 1)

        call ol_ew_setup

        call set_parameter("verbose",-1)
        call set_parameter("model","sm") 

       id_NLO_1r(2) = register_process("2 -2 -> 11 -11 21", 1)  ! u ū
       id_NLO_1r(1) = register_process("1 -1 -> 11 -11 21", 1)  ! d d̄
       id_NLO_1r(4) = register_process("4 -4 -> 11 -11 21", 1)  ! c c̄
       id_NLO_1r(3) = register_process("3 -3 -> 11 -11 21", 1)  ! s s̄
       id_NLO_1r(5) = register_process("5 -5 -> 11 -11 21", 1)  ! b b̄

      return
      end subroutine ol_NLO_1real_init
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       subroutine ol_NLO_init_1loop
        use globals_mod, only : id_NLO_1loop
        use constants_mod
        use openloops

      
!        call set_parameter("order_ew", 2)
!        call set_parameter("order_qcd", 1)

        call ol_ew_setup
        call set_parameter("ckmorder", 0)

        call set_parameter("ct_on", 0)
        call set_parameter("iop_on",1)
!        call set_parameter("r2_on", 1)


        call set_parameter("verbose",-1)
        call set_parameter("model","sm") 

       id_NLO_1loop(2) = register_process("2 -2 -> 11 -11", 11)  ! u ū
       id_NLO_1loop(1) = register_process("1 -1 -> 11 -11", 11)  ! d d̄
                                                              
       id_NLO_1loop(4) = register_process("4 -4 -> 11 -11", 11)  ! c c̄
       id_NLO_1loop(3) = register_process("3 -3 -> 11 -11", 11)  ! s s̄
                                                              
       id_NLO_1loop(5) = register_process("5 -5 -> 11 -11", 11)  ! b b̄

      return
      end subroutine ol_NLO_init_1loop
