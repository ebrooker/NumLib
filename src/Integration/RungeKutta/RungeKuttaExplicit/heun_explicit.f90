!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Integation/RungeKutta/RungeKuttaExplicit/heun_explicit.f90
!!
!! SubModule containing Explicit Heun method, a 2nd order accurate explicit
!! Runge Kutta integration method.
!!
!!*****************************************************************************
SUBMODULE (rungekutta_explicit) heun_explicit

    REAL(rkp), PARAMETER :: one_half = 0.5_rkp

CONTAINS

    MODULE SUBROUTINE heun_explicit_1eqn (f, t, y0, y, h)
        !! Explicit Heun ODE integration method
        !! for a single ODE
        PROCEDURE(ode_single)    :: f
        REAL(rkp), INTENT(IN)    :: t
        REAL(rkp), INTENT(IN)    :: y0
        REAL(rkp), INTENT(OUT)   :: y
        REAL(rkp), INTENT(INOUT) :: h

        REAL(rkp) :: k1, k2

        k1 = h * f(t,     y0     )
        k2 = h * f(t + h, y0 + k1)

        y  = y0 + one_half*(k1 + k2)

    END SUBROUTINE heun_explicit_1eqn

    MODULE SUBROUTINE heun_explicit_neqn (f, t, y0, y, h)
        !! Explicit Heun ODE integration method
        !! for a system of ODEs
        PROCEDURE(ode_system)    :: f
        REAL(rkp), INTENT(IN)    :: t
        REAL(rkp), INTENT(IN)    :: y0(:)
        REAL(rkp), INTENT(OUT)   :: y(SIZE(y0,DIM=1))
        REAL(rkp), INTENT(INOUT) :: h

        REAL(rkp), DIMENSION(SIZE(y0,DIM=1)) :: k1, k2

        k1 = h * f(t,     y0     )
        k2 = h * f(t + h, y0 + k1)

        y  = y0 + one_half*(k1 + k2)

    END SUBROUTINE heun_explicit_neqn

END SUBMODULE heun_explicit