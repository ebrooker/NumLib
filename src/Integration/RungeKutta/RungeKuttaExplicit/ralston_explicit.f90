!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Integation/RungeKutta/RungeKuttaExplicit/ralston_explicit.f90
!!
!! Module containing Explicit Ralston method, a 2nd order accurate explicit
!! Runge Kutta integration method.
!!
!!*****************************************************************************
SUBMODULE (rungekutta_explicit) ralston_explicit

    REAL(rkp), PARAMETER :: three     = 3.0_rkp
    REAL(rkp), PARAMETER :: two_thrds = 2.0_rkp/three
    REAL(rkp), PARAMETER :: one_forth = 0.25_rkp

CONTAINS

    MODULE SUBROUTINE ralston_explicit_1eqn (f, t, y0, y, h)
        !! Explicit Ralston ODE integration method
        !! for a single ODE
        PROCEDURE(ode_single)    :: f
        REAL(rkp), INTENT(IN)    :: t
        REAL(rkp), INTENT(IN)    :: y0
        REAL(rkp), INTENT(OUT)   :: y
        REAL(rkp), INTENT(INOUT) :: h

        REAL(rkp) :: k1, k2

        k1 = h * f(t,               y0               )
        k2 = h * f(t + two_thrds*h, y0 + two_thrds*k1)
        
        y  = y0 + one_forth*(k1 + three*k2)

    END SUBROUTINE ralston_explicit_1eqn

    MODULE SUBROUTINE ralston_explicit_neqn (f, t, y0, y, h)
        !! Explicit Ralston ODE integration method
        !! for a system of ODEs
        PROCEDURE(ode_system)    :: f
        REAL(rkp), INTENT(IN)    :: t
        REAL(rkp), INTENT(IN)    :: y0(:)
        REAL(rkp), INTENT(OUT)   :: y(SIZE(y0,DIM=1))
        REAL(rkp), INTENT(INOUT) :: h

        REAL(rkp), DIMENSION(SIZE(y0,DIM=1)) :: k1, k2

        k1 = h * f(t,               y0               )
        k2 = h * f(t + two_thrds*h, y0 + two_thrds*k1)
        
        y  = y0 + one_forth*(k1 + three*k2)

    END SUBROUTINE ralston_explicit_neqn

END SUBMODULE ralston_explicit