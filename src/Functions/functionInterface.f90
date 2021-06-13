!!*********************************************
!!
!! Author: Ezra Brooker
!! 2021
!! 
!! src/Functions/functionInterface.f90
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

END MODULE functionInterface