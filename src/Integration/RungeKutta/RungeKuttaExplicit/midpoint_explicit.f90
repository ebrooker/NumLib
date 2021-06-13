!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Integation/RungeKutta/RungeKuttaExplicit/midpoint_explicit.f90
!!
!! SubModule containing Explicit Midpoint method, a 2nd order accurate explicit
!! Runge Kutta integration method.
!!
!!*****************************************************************************
SUBMODULE (rungekutta_explicit) midpoint_explicit    

    REAL(rkp), PARAMETER :: one_half = 0.5_rkp

CONTAINS

    MODULE SUBROUTINE midpoint_explicit_1eqn (f, t, y0, y, h)
        !! Explicit Midpoint ODE integration method
        !! for a single ODE
        PROCEDURE(ode_single)    :: f
        REAL(rkp), INTENT(IN)    :: t
        REAL(rkp), INTENT(IN)    :: y0
        REAL(rkp), INTENT(OUT)   :: y
        REAL(rkp), INTENT(INOUT) :: h

        y = y0 + one_half*h*f(t,y0)
        y = y0 + h*f(t+one_half*h,y)

    END SUBROUTINE midpoint_explicit_1eqn

    MODULE SUBROUTINE midpoint_explicit_neqn (f, t, y0, y, h)
        !! Explicit Midpoint ODE integration method
        !! for a system of ODEs
        PROCEDURE(ode_system)    :: f
        REAL(rkp), INTENT(IN)    :: t
        REAL(rkp), INTENT(IN)    :: y0(:)
        REAL(rkp), INTENT(OUT)   :: y(SIZE(y0,DIM=1))
        REAL(rkp), INTENT(INOUT) :: h

        y = y0 + one_half*h*f(t,y0)
        y = y0 + h*f(t+one_half*h,y)

    END SUBROUTINE midpoint_explicit_neqn

END SUBMODULE midpoint_explicit