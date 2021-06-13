!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Integation/RungeKutta/RungeKuttaExplicit/ssprk3_explicit.f90
!!
!! Submodule containing the 3rd order accurate Strong Stability PReserving
!! Runge-Kutta (SSPRK3) integrator
!!
!!*****************************************************************************
SUBMODULE (rungekutta_explicit) ssprk3_explicit

    REAL(rkp), PARAMETER :: one_sixth = 1.0_rkp/6.0_rkp
    REAL(rkp), PARAMETER :: one_forth = 0.25_rkp
    REAL(rkp), PARAMETER :: two_thrds = 4.0_rkp * one_sixth
    REAL(rkp), PARAMETER :: one_half  = 0.5_rkp

CONTAINS

    MODULE SUBROUTINE ssprk3_explicit_1eqn (f, t, y0, y, h)
        !! Explicit SSPRK3 ODE integration method
        !! for a single ODE
        PROCEDURE(ode_single)    :: f
        REAL(rkp), INTENT(IN)    :: t
        REAL(rkp), INTENT(IN)    :: y0
        REAL(rkp), INTENT(OUT)   :: y
        REAL(rkp), INTENT(INOUT) :: h

        REAL(rkp) :: k1, k2, k3

        k1 = h * f(t,              y0                    )
        k2 = h * f(t +          h, y0 + k1               )
        k3 = h * f(t + one_half*h, y0 + one_forth*(k1+k2))

        y  = y0 + one_sixth*(k1 + k2) + two_thrds*k3

    END SUBROUTINE ssprk3_explicit_1eqn

    MODULE SUBROUTINE ssprk3_explicit_neqn (f, t, y0, y, h)
        !! Explicit SSPRK3 ODE integration method
        !! for a system of ODEs
        PROCEDURE(ode_system)    :: f
        REAL(rkp), INTENT(IN)    :: t
        REAL(rkp), INTENT(IN)    :: y0(:)
        REAL(rkp), INTENT(OUT)   :: y(SIZE(y0,DIM=1))
        REAL(rkp), INTENT(INOUT) :: h

        REAL(rkp), DIMENSION(SIZE(y0,DIM=1)) :: k1,k2,k3


        k1 = h * f(t,              y0                    )
        k2 = h * f(t +          h, y0 + k1               )
        k3 = h * f(t + one_half*h, y0 + one_forth*(k1+k2))

        y  = y0 + one_sixth*(k1 + k2) + two_thrds*k3

    END SUBROUTINE ssprk3_explicit_neqn

END SUBMODULE ssprk3_explicit