module PK_functions_mod
  use constants_mod
  implicit none
  private

  public :: PK, dilog, xint_h

contains

  function PK(PKtype, Partontype, fxntype, x) result(res)
    implicit none
    character(len=*), intent(in) :: PKtype, Partontype, fxntype
    double precision, intent(in) :: x
    double precision :: res
    double precision :: one_minus_x, log1mxbyx, log1mx

    one_minus_x = 1.0d0 - x
    log1mxbyx = dlog((1d0 - x) / x)
    log1mx = dlog(1d0 - x)

    select case (trim(Partontype))
  !=============!
    case ('gg')
  !=============!
      select case (trim(PKtype))
      case ('P')
        select case (trim(fxntype))
        case ('Plus')
          res = 2d0 * CA / one_minus_x
        case ('Regular')
          res = 2d0 * CA * (one_minus_x / x - 1d0 + x * one_minus_x)
        case ('Delta')
          res = 11d0/6d0 * CA - 2d0/3d0 * NF * TR
        case default
          call invalid_fxntype()
        end select

      case ('Kb')
        select case (trim(fxntype))
        case ('Plus')
          res = 2d0 * CA / one_minus_x * log1mxbyx
        case ('Regular')
          res = 2d0 * CA * (one_minus_x / x - 1d0 + x * one_minus_x) * log1mxbyx
        case ('Delta')
          res = -((50d0/9d0 - PI**2) * CA - 16d0/9d0 * TR * NF)
        case default
          call invalid_fxntype()
        end select

      case ('Kt')
        select case (trim(fxntype))
        case ('Plus')
          res = 2d0 * CA / one_minus_x * log1mx
        case ('Regular')
          res = 2d0 * CA * (one_minus_x / x - 1d0 + x * one_minus_x) * log1mx
        case ('Delta')
          res = -CA * PI**2 / 3d0
        case default
          call invalid_fxntype()
        end select

      case default
        call invalid_PKtype()

      end select

  !=============!
    case ('qq')
  !=============!
      select case (trim(PKtype))
        case ('P')
          select case (trim(fxntype))
             case ('Plus')
               res = CF * (1d0 + x*x) / one_minus_x
             case ('Regular' , 'Delta')
               res = 0d0 
             case default
               call invalid_fxntype()
          end select

        case ('Kb')
          select case (trim(fxntype))
             case ('Plus')
               res = CF * 2d0 / one_minus_x * log1mxbyx
             case ('Regular')
               res = CF * (-(1d0 + x) * log1mxbyx + one_minus_x)
             case ('Delta')
               res = -CF * (5d0 - PI*PI)
          case default
            call invalid_fxntype()
          end select

        case ('Kt')
          select case (trim(fxntype))
             case ('Plus')
               res = CF * 2d0 / one_minus_x * log1mx
             case ('Regular')
               res = -CF * (1d0 + x) * log1mx
             case ('Delta')
               res = -CF * PI*PI / 3d0
          case default
            call invalid_fxntype()
          end select

        case default
          call invalid_PKtype()
      end select

  !=============!
    case ('qg')
  !=============!
      select case (trim(PKtype))
      case ('P')
        select case (trim(fxntype))
        case ('Plus', 'Delta')
          res = 0d0
        case ('Regular')
          res = CF * (1d0 + one_minus_x**2) / x
        case default
          call invalid_fxntype()
        end select

      case ('Kb')
        select case (trim(fxntype))
        case ('Plus', 'Delta')
          res = 0d0
        case ('Regular')
          res = CF * (1d0 + one_minus_x**2) * log1mxbyx / x + CF * x
        case default
          call invalid_fxntype()
        end select

      case ('Kt')
        select case (trim(fxntype))
        case ('Plus', 'Delta')
          res = 0d0
        case ('Regular')
          res = CF * (1d0 + one_minus_x**2) * log1mx / x
        case default
          call invalid_fxntype()
        end select

      case default
        call invalid_PKtype()
      end select

  !=============!
    case ('gq')
  !=============!
      select case (trim(PKtype))
      case ('P')
        select case (trim(fxntype))
        case ('Plus', 'Delta')
          res = 0d0
        case ('Regular')
          res = TR * (x**2 + one_minus_x**2)
        case default
          call invalid_fxntype()
        end select

      case ('Kb')
        select case (trim(fxntype))
        case ('Plus', 'Delta')
          res = 0d0
        case ('Regular')
          res = TR * (x**2 + one_minus_x**2) * log1mxbyx + 2d0 * TR * x * one_minus_x
        case default
          call invalid_fxntype()
        end select

      case ('Kt')
        select case (trim(fxntype))
        case ('Plus', 'Delta')
          res = 0d0
        case ('Regular')
          res = TR * (x**2 + one_minus_x**2) * log1mx
        case default
          call invalid_fxntype()
        end select

      case default
        call invalid_PKtype()
      end select

    case default
      print *, "Invalid Partontype. Expected: gg, qq, qg, gq"
      stop
    end select

  end function PK

  subroutine invalid_fxntype()
    print *, "Invalid function type. Use: Plus, Regular, or Delta"
    stop
  end subroutine

  subroutine invalid_PKtype()
    print *, "Invalid PKtype. Use: P, Kb, or Kt"
    stop
  end subroutine

 function xint_h(x, PKtype) result(val)
  implicit none
  ! Arguments
  double precision, intent(in) :: x
  character(len=*), intent(in) :: PKtype
  double precision :: val


  ! Main logic
  select case (trim(PKtype))

  case ('Pqq')
    val = CF * (3d0 - 2d0*x - x**2 - 4d0*dlog(1d0 - x)) / 2d0

  case ('Kbarqq')
    val = -CF* (PI**2 + dlog(1.0-x)**2 + 2.0*dilog(1.0-x) )

  case ('Ktqq')
    val = -CF * log(1d0 - x)**2

  case ('Pgg')
    val = -2d0 * CA * log(1d0 - x)

  case ('Kbargg')
    val = -0.5d0 * PI**2 - 0.5d0 * log(1d0 - x)**2 - dilog(1d0 - x)
    val = 2d0 * CA * val

  case ('Ktgg')
    val = -log(1d0 - x)**2
    val = CA * val

  case default
    print *, 'Invalid xintPK type: ', trim(PKtype)
    stop

  end select
end function xint_h

       double precision function dilog(x)
       implicit double precision (a-z)
       parameter (pi6=1.644934066848226d+00)
       parameter (een=1.d+00)
       parameter (vier=0.25d+00)
       parameter (b2=+2.7777777777777778D-02)
       parameter (b3=-2.7777777777777778D-04)
       parameter (b4=+4.7241118669690098D-06)
       parameter (b5=-9.1857730746619641D-08)
       parameter (b6=+1.8978869988971001D-09)
       parameter (b7=-4.0647616451442256D-11)
       parameter (b8=+8.9216910204564523D-13)
1      if(x.lt.0.d0)go to 3
       if(x.gt.0.5d0)go to 4
       z=-dlog(1.d0-x)
7      z2=z*z
       dilog=z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b8+b7)+b6)+b5)+b4)+b3)+b2)+een)-z2*vier
       if(x.gt.een)dilog=-dilog-.5d0*u*u+2.d0*pi6
       return
3      if(x.gt.-een)go to 5
       y=een/(een-x)
       z=-dlog(een-y)
       z2=z*z
       u=dlog(y)
       dilog=z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b8+b7)+b6)+b5)+b4)+b3)+b2)+een)-z2*vier-u*(z+.5d0*u)-pi6
       return
4      if(x.ge.een)go to 10
       y=een-x
       z=-dlog(x)
6      u=dlog(y)
       z2=z*z
       dilog=-z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b8+b7)+b6)+b5)+b4)+b3)+b2)+een-u)+z2*vier+pi6
       if(x.gt.een)dilog=-dilog-.5d0*z*z+pi6*2.d0
       return
5      y=een/(een-x)
       z=-dlog(y)
       z2=z*z
       dilog=-z*(z2*(z2*(z2*(z2*(z2*(z2*(z2*b8+b7)+b6)+b5)+b4)+b3)+b2)+een)-z2*vier
       return
10     if(x.eq.een)go to 20
       xx=1.d0/x
       if(x.gt.2.d0)go to 11
       z=dlog(x)
       y=1.d0-xx
       go to 6
11     u=dlog(x)
       z=-dlog(1.d0-xx)
       go to 7
20     dilog=pi6
       return
       end function dilog

end module PK_functions_mod

