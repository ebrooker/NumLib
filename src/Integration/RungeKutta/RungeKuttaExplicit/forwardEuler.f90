!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Integation/RungeKutta/RungeKuttaExplicit/forwardEuler.f90
!!
!! Module containing Forward Euler method, the 1st order accurate explicit
!! Rungke Kutta integration method.
!!
!!*****************************************************************************
SUBMODULE (rungekutta_explicit) forwardEuler

CONTAINS

    MODULE PROCEDURE fwd_euler
        ! Forward Euler ODE integration method
        
        DO i = 1,neqn
            y(i) = y0(i) + h * f(t,y0(i))
        ENDDO

    END PROCEDURE fwd_euler


END SUBMODULE forwardEuler