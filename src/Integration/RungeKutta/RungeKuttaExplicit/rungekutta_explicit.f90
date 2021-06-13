!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Integation/RungeKutta/RungeKuttaExplicit/rungekutta_explicit.f90
!!
!! Module used to collate all of the explicit Runge Kutta integration routines
!! into one namespace for easier porting as a library code.
!!
!! For each formula, there are two implementations, one for a single ODE eqn
!! and one for a system of ODE eqns. These two implementations for each formula
!! are housed under an interface name of the formula. The actual implementation
!! of each version of a given formula is housed in a seperate submodule. This
!! module will only define the necessary interfaces for all of the methods.
!!
!! For example:
!!
!!     - forward_euler_1eqn: Single ODE implementation (1eqn = 1 equation)
!!     - forward_euler_neqn: System of ODEs implementation (neqn = N equations)
!!     - forward_euler_step: Interface name for above implementations
!!
!! Other explicit RK methods defined here:
!!     - midpoint_explicit_step: 2nd order, Midpoint method
!!     - heun_explicit_step:     2nd order, Heun's method
!!     - ralston_explicit_step:  2nd order, Ralston's method
!!     - ssprk3_explicit_step:   3rd order, Strong Stability Preserving RK
!!     - rk4_explicit_step:      4th order, generic RK4 method
!!
!!*****************************************************************************
MODULE rungekutta_explicit
USE kindSettings, ONLY : rkp, ilkp, iskp
USE functionInterface, ONLY : f_1d_intrfc, ode_single, ode_system
IMPLICIT NONE
PRIVATE

    PUBLIC :: forward_euler_step
    PUBLIC :: midpoint_explicit_step
    PUBLIC :: heun_explicit_step
    PUBLIC :: ralston_explicit_step
    PUBLIC :: ssprk3_explicit_step
    PUBLIC :: rk4_explicit_step

    INTERFACE forward_euler_step
        MODULE SUBROUTINE forward_euler_1eqn (f, t, y0, y, h)
            ! Forward Euler ODE integration method
            PROCEDURE(ode_single)    :: f
            REAL(rkp), INTENT(IN)    :: t
            REAL(rkp), INTENT(IN)    :: y0
            REAL(rkp), INTENT(OUT)   :: y
            REAL(rkp), INTENT(INOUT) :: h
        END SUBROUTINE forward_euler_1eqn

        MODULE SUBROUTINE forward_euler_neqn (f, t, y0, y, h)
            ! Forward Euler ODE integration method
            PROCEDURE(ode_system)    :: f
            REAL(rkp), INTENT(IN)    :: t
            REAL(rkp), INTENT(IN)    :: y0(:)
            REAL(rkp), INTENT(OUT)   :: y(SIZE(y0,DIM=1))
            REAL(rkp), INTENT(INOUT) :: h
        END SUBROUTINE forward_euler_neqn
    END INTERFACE forward_euler_step

    !!**********************************************************!!

    INTERFACE midpoint_explicit_step
        MODULE SUBROUTINE midpoint_explicit_1eqn (f, t, y0, y, h)
            ! Explicit Midpoint ODE integration method
            PROCEDURE(ode_single)    :: f
            REAL(rkp), INTENT(IN)    :: t
            REAL(rkp), INTENT(IN)    :: y0
            REAL(rkp), INTENT(OUT)   :: y
            REAL(rkp), INTENT(INOUT) :: h
        END SUBROUTINE midpoint_explicit_1eqn
        
        MODULE SUBROUTINE midpoint_explicit_neqn (f, t, y0, y, h)
            ! Explicit Midpoint ODE integration method
            PROCEDURE(ode_system)    :: f
            REAL(rkp), INTENT(IN)    :: t
            REAL(rkp), INTENT(IN)    :: y0(:)
            REAL(rkp), INTENT(OUT)   :: y(SIZE(y0,DIM=1))
            REAL(rkp), INTENT(INOUT) :: h
        END SUBROUTINE midpoint_explicit_neqn
    END INTERFACE midpoint_explicit_step

    !!**********************************************************!!

    INTERFACE ralston_explicit_step
        MODULE SUBROUTINE ralston_explicit_1eqn (f, t, y0, y, h)
            ! Explicit Ralston ODE integration method
            PROCEDURE(ode_single)    :: f
            REAL(rkp), INTENT(IN)    :: t
            REAL(rkp), INTENT(IN)    :: y0
            REAL(rkp), INTENT(OUT)   :: y
            REAL(rkp), INTENT(INOUT) :: h
        END SUBROUTINE ralston_explicit_1eqn
        
        MODULE SUBROUTINE ralston_explicit_neqn (f, t, y0, y, h)
            ! Explicit Ralston ODE integration method
            PROCEDURE(ode_system)    :: f
            REAL(rkp), INTENT(IN)    :: t
            REAL(rkp), INTENT(IN)    :: y0(:)
            REAL(rkp), INTENT(OUT)   :: y(SIZE(y0,DIM=1))
            REAL(rkp), INTENT(INOUT) :: h
        END SUBROUTINE ralston_explicit_neqn
    END INTERFACE ralston_explicit_step

    !!**********************************************************!!

    INTERFACE heun_explicit_step
        MODULE SUBROUTINE heun_explicit_1eqn (f, t, y0, y, h)
            ! Explicit Heun's ODE integration method
            PROCEDURE(ode_single)    :: f
            REAL(rkp), INTENT(IN)    :: t
            REAL(rkp), INTENT(IN)    :: y0
            REAL(rkp), INTENT(OUT)   :: y
            REAL(rkp), INTENT(INOUT) :: h
        END SUBROUTINE heun_explicit_1eqn
        
        MODULE SUBROUTINE heun_explicit_neqn (f, t, y0, y, h)
            ! Explicit Heun's ODE integration method
            PROCEDURE(ode_system)    :: f
            REAL(rkp), INTENT(IN)    :: t
            REAL(rkp), INTENT(IN)    :: y0(:)
            REAL(rkp), INTENT(OUT)   :: y(SIZE(y0,DIM=1))
            REAL(rkp), INTENT(INOUT) :: h
        END SUBROUTINE heun_explicit_neqn
    END INTERFACE heun_explicit_step

    !!**********************************************************!!

    INTERFACE ssprk3_explicit_step
        MODULE SUBROUTINE ssprk3_explicit_1eqn (f, t, y0, y, h)
            ! Explicit SSPRK3 ODE integration method
            PROCEDURE(ode_single)    :: f
            REAL(rkp), INTENT(IN)    :: t
            REAL(rkp), INTENT(IN)    :: y0
            REAL(rkp), INTENT(OUT)   :: y
            REAL(rkp), INTENT(INOUT) :: h
        END SUBROUTINE ssprk3_explicit_1eqn
        
        MODULE SUBROUTINE ssprk3_explicit_neqn (f, t, y0, y, h)
            ! Explicit SSPRK3 ODE integration method
            PROCEDURE(ode_system)    :: f
            REAL(rkp), INTENT(IN)    :: t
            REAL(rkp), INTENT(IN)    :: y0(:)
            REAL(rkp), INTENT(OUT)   :: y(SIZE(y0,DIM=1))
            REAL(rkp), INTENT(INOUT) :: h
        END SUBROUTINE ssprk3_explicit_neqn
    END INTERFACE ssprk3_explicit_step

    INTERFACE rk4_explicit_step
        MODULE SUBROUTINE rk4_explicit_1eqn (f, t, y0, y, h)
            ! Explicit RK4 ODE integration method
            PROCEDURE(ode_single)    :: f
            REAL(rkp), INTENT(IN)    :: t
            REAL(rkp), INTENT(IN)    :: y0
            REAL(rkp), INTENT(OUT)   :: y
            REAL(rkp), INTENT(INOUT) :: h
        END SUBROUTINE rk4_explicit_1eqn
        
        MODULE SUBROUTINE rk4_explicit_neqn (f, t, y0, y, h)
            ! Explicit RK4 ODE integration method
            PROCEDURE(ode_system)    :: f
            REAL(rkp), INTENT(IN)    :: t
            REAL(rkp), INTENT(IN)    :: y0(:)
            REAL(rkp), INTENT(OUT)   :: y(SIZE(y0,DIM=1))
            REAL(rkp), INTENT(INOUT) :: h
        END SUBROUTINE rk4_explicit_neqn
    END INTERFACE rk4_explicit_step



END MODULE rungekutta_explicit
