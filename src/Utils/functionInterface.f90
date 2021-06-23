!!*********************************************
!!
!! Author: Ezra Brooker
!! 2021
!! 
!! src/Utils/functionInterface.f90
!!
!! Module provides an abstract interface for
!! 1d simple functions of the form f(x) = y
!!
!! This allows for other procedures to accept
!! specific function handles as input arguments
!! as long as the generic argument handle is 
!! declared in the procedure as follows:
!!
!!     PROCEDURE(f_1d_intrfc) :: func_handl
!!
!! This adds a level of OOP to the parts of
!! NumLib that make use of this interface
!! such as the finite differentiation
!! mehods in the Differentiation subpackage
!!
!! Additional interfaces are provided for the
!! ODE integration subpackage to allow for the
!! integration of a single ODE and systems of
!! ODEs and is used to make a generic interface
!! for integrating both types
!!
!!*********************************************
MODULE functionInterface
IMPLICIT NONE
PRIVATE

    PUBLIC :: f_1d_intrfc

    ABSTRACT INTERFACE
        FUNCTION f_1d_intrfc(Xx) RESULT(Yy)
            USE kindSettings, ONLY : rkp
            REAL(rkp)                       :: Yy
            REAL(rkp), INTENT(IN)           :: Xx
        END FUNCTION f_1d_intrfc
    END INTERFACE

    ABSTRACT INTERFACE
        FUNCTION ode_single(x,y0) RESULT(y)
            USE kindSettings, ONLY : rkp
            REAL(rkp), INTENT(IN) :: x
            REAL(rkp), INTENT(IN) :: y0
            REAL(rkp)             :: y
        END FUNCTION ode_single
    END INTERFACE

    ABSTRACT INTERFACE
        FUNCTION ode_system(x,y0) RESULT(y)
            USE kindSettings, ONLY : rkp
            REAL(rkp), INTENT(IN) :: x
            REAL(rkp), INTENT(IN) :: y0(:)
            REAL(rkp)             :: y(SIZE(y0,DIM=1))
        END FUNCTION ode_system
    END INTERFACE

END MODULE functionInterface