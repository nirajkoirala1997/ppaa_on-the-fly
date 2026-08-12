!c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       subroutine ol_LO_init
        use globals_mod, only : id_LO
        use constants_mod
        use openloops

      
        call set_parameter("order_ew", 2)
        call set_parameter("order_qcd", 0)

        call set_parameter("mass(23)", mass(23))
        call set_parameter("mass(24)", mass(24))
        call set_parameter("width(23)", decayW(23))
        call set_parameter("width(24)", decayW(24))
        call set_parameter("ew_scheme",2)

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


        call set_parameter("mass(23)", mass(23))
        call set_parameter("mass(24)", mass(24))
        call set_parameter("width(23)", decayW(23))
        call set_parameter("width(24)", decayW(24))
        call set_parameter("ew_scheme",2)

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

        call set_parameter("mass(23)", mass(23))
        call set_parameter("mass(24)", mass(24))
        call set_parameter("width(23)", decayW(23))
        call set_parameter("width(24)", decayW(24))
        call set_parameter("ckmorder", 0)
        call set_parameter("ew_scheme",2)

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
