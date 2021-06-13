!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Integation/RungeKutta/RungeKuttaExplicit/forward_euler.f90
!!
!! Module containing Forward Euler method, a 1st order accurate explicit
!! Runge Kutta integration method. There is an implementation for both the
!! single ODE and system of ODE cases housed under a single interface name of
!! "forward_euler"
!!
!!*****************************************************************************
SUBMODULE (rungekutta_explicit) forward_euler

CONTAINS

    MODULE SUBROUTINE forward_euler_1eqn (f, t, y0, y, h)
        !! Forward Euler ODE integration method
        !! for a single ODE
        PROCEDURE(ode_single)    :: f
        REAL(rkp), INTENT(IN)    :: t
        REAL(rkp), INTENT(IN)    :: y0
        REAL(rkp), INTENT(OUT)   :: y
        REAL(rkp), INTENT(INOUT) :: h

        y = y0 + h * f(t,y0)

    END SUBROUTINE forward_euler_1eqn

    MODULE SUBROUTINE forward_euler_neqn (f, t, y0, y, h)
        !! Forward Euler ODE integration method
        !! for a system of ODEs
        PROCEDURE(ode_system)    :: f
        REAL(rkp), INTENT(IN)    :: t
        REAL(rkp), INTENT(IN)    :: y0(:)
        REAL(rkp), INTENT(OUT)   :: y(SIZE(y0,DIM=1))
        REAL(rkp), INTENT(INOUT) :: h

        y = y0 + h * f(t,y0)

    END SUBROUTINE forward_euler_neqn

END SUBMODULE forward_euler