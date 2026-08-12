module matrix_elements_mod
  use openloops
  use globals_mod, only : ge, id_LO, id_NLO_1r, id_NLO_1loop
  use misc_utils_mod
  use constants_mod
  implicit none
  private

  public :: evaluate_born,evaluate_1real,evaluate_1loop

contains

  !========================================================!
  subroutine evaluate_born(p1, p2, p3, p4,answer)
    real(8), intent(out) :: answer(5) 
    real(8), intent(in) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3)

    ! Internal variables
    real(8) :: s12, s13, s23
    real(8) ::  e
    real(8) :: p_ex(0:3,4)
!    real(8) :: qu2,gew, qe, qu, cw, sw, cw2, sw2, zprop
    integer :: i
    real(8) :: qe, qu, qd, qu2, qd2    

    ! Initialize coupling
    e = dsqrt(ge * 4.d0 * PI)

    ! Mandelstam variables
    s13 = 2.d0 * dot(p1, p3) ! t
    s23 = 2.d0 * dot(p1, p4) ! u
    s12 = 2.d0 * dot(p1, p2) ! s
    qe = e
    qu = 2d0/3d0 
    qd = -1d0/3d0
    qu2 = qu**2 
    qd2 = qd**2 
    answer(2) = (2*e**4*qu2*(-2*s13*s23 + s12*(s13 + s23)))/(3d0*s12**2)
    answer(1) = (2*e**4*qd2*(-2*s13*s23 + s12*(s13 + s23)))/(3d0*s12**2)
    answer(3) = answer(1)
    answer(4) = answer(2)
    answer(5) = answer(1)
    answer = 0d0

! Set OpenLoops parameters and evaluate
!    call SET_PARAMETER('alpha_s', AL)
!    call SET_PARAMETER('mu', scale)

    call p1d_to_p2d_4(p1, p2, p3, p4, p_ex)
    do i = 1,5
    call evaluate_tree(id_LO(i), p_ex, answer(i))
    enddo

  end subroutine evaluate_born
  !========================================================!

  !========================================================!
  subroutine evaluate_1real(AL, p1, p2, p3, p4, p5,answer)
    real(8), intent(out) :: answer(5) 
    real(8), intent(in) :: AL, p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3)

    ! Internal variables
    real(8) ::  e
    real(8) :: p_ex(0:3,5)
    integer :: i

    ! Initialize coupling
    e = dsqrt(ge * 4.d0 * PI)

! Set OpenLoops parameters and evaluate
    call SET_PARAMETER('alpha_s', AL)
!    call SET_PARAMETER('mu', scale)
    call p1d_to_p2d_5(p1, p2, p3, p4 ,p5, p_ex)
    do i = 1,5
    call evaluate_tree(id_NLO_1r(i), p_ex, answer(i))
    enddo
  end subroutine evaluate_1real
  !========================================================!


  !========================================================!
  subroutine evaluate_1loop(AL,mur,p1, p2, p3, p4, answer)
    real(8), intent(out) :: answer(5) 
    real(8), intent(in) :: AL, p1(0:3), p2(0:3), p3(0:3), p4(0:3), mur

    ! Internal variables
    real(8) ::  e
    real(8) :: p_ex(0:3,5)
    integer :: i
    real(8) :: am210 , am211(0:2), acc

    ! Initialize coupling
    e = dsqrt(ge * 4.d0 * PI)

      call SET_PARAMETER('alpha_s', AL)
      call SET_PARAMETER('mu', mur)

      call set_parameter("check_poles",1)
      call set_parameter("ckmorder",6)
      call set_parameter("polenorm",1)
      call set_parameter("ct_on", 1)
      call set_parameter("iop_on",1)
      call set_parameter("r2_on", 1)

!      call evaluate_tree(id_LO,p_ex,born)
!      call evaluate_cc(id_LO,p_ex,amp_210,amp2cc,amp2ewcc)

!      call evaluate_loopct(id_NLO_1loop,p_ex,aborn,am2Rct)

! Set OpenLoops parameters and evaluate
     answer = 0.0
     call p1d_to_p2d_4(p1, p2, p3, p4 , p_ex)
    do i = 1,5
     call evaluate_loop(id_NLO_1loop(i), p_ex, am210, am211, acc)
     answer(i) = am211(0)
    enddo


  end subroutine evaluate_1loop
  !========================================================!



end module matrix_elements_mod
