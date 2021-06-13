!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Integation/RungeKutta/RungeKuttaExplicit/rk4_explicit.f90
!!
!! Module containing Explicit RK4 method, a 4th order accurate explicit
!! Runge Kutta integration method.
!!
!!*****************************************************************************
SUBMODULE (rungekutta_explicit) rk4_explicit

    REAL(rkp), PARAMETER :: one_sixth = 1.0_rkp/6.0_rkp
    REAL(rkp), PARAMETER :: one_half  = 0.5_rkp
    REAL(rkp), PARAMETER :: two       = 2.0_rkp

CONTAINS

    MODULE SUBROUTINE rk4_explicit_1eqn (f, t, y0, y, h)
        !! Explicit RK4 ODE integration method
        !! for a single ODE
        PROCEDURE(ode_single)    :: f
        REAL(rkp), INTENT(IN)    :: t
        REAL(rkp), INTENT(IN)    :: y0
        REAL(rkp), INTENT(OUT)   :: y
        REAL(rkp), INTENT(INOUT) :: h

        REAL(rkp) :: k1, k2, k3, k4

        k1 = h * f(t,              y0              )
        k2 = h * f(t + one_half*h, y0 + one_half*k1)
        k3 = h * f(t + one_half*h, y0 + one_half*k2)
        k4 = h * f(t + h,          y0 + k3         )

        y  = y0 + one_sixth*(k1 + two*(k2+k3) + k4)

    END SUBROUTINE rk4_explicit_1eqn

    MODULE SUBROUTINE rk4_explicit_neqn (f, t, y0, y, h)
        !! Explicit RK4 ODE integration method
        !! for a system of ODEs
        PROCEDURE(ode_system)    :: f
        REAL(rkp), INTENT(IN)    :: t
        REAL(rkp), INTENT(IN)    :: y0(:)
        REAL(rkp), INTENT(OUT)   :: y(SIZE(y0,DIM=1))
        REAL(rkp), INTENT(INOUT) :: h

        REAL(rkp), DIMENSION(SIZE(y0,DIM=1)) :: k1,k2,k3,k4


        k1 = h * f(t,              y0              )
        k2 = h * f(t + one_half*h, y0 + one_half*k1)
        k3 = h * f(t + one_half*h, y0 + one_half*k2)
        k4 = h * f(t + h,          y0 + k3         )
        
        y  = y0 + one_sixth*(k1 + two*(k2+k3) + k4)

    END SUBROUTINE rk4_explicit_neqn

END SUBMODULE rk4_explicit