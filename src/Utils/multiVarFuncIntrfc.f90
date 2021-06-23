!!*********************************************
!!
!! Author: Ezra Brooker
!! 2021
!! 
!! src/Utils/multiVarFuncIntrfc.f90
!!
!! Module provides an abstract interface for
!! nD simple functions of the form f(x(n)) = y(n)
!!
!! This allows for other procedures to accept
!! specific function handles as input arguments
!! as long as the generic argument handle is 
!! declared in the procedure as follows:
!!
!!     PROCEDURE(f_intrfc) :: func_handl
!!
!! This adds a level of OOP to the parts of
!! NumLib that make use of this interface
!! such as the global minimization methods
!! in the Optimization subpackage
!!
!!
!! Other abstract interfaces provided are 
!!
!!     gradf_intrfc
!!     hessf_intrfc
!!
!! these function similarly to f_intrfc,
!! but now adapted to allow nD gradients
!! and hessians
!!
!!*********************************************
MODULE  multiVarFuncIntrfc
IMPLICIT NONE
PRIVATE

    PUBLIC :: f_intrfc, gradf_intrfc, hessf_intrfc

    ABSTRACT INTERFACE
        FUNCTION f_intrfc(x) RESULT(f)
            USE kindSettings, ONLY : rkp
            REAL(rkp), INTENT(IN) :: x(:)
            REAL(rkp)             :: f
        END FUNCTION f_intrfc
    END INTERFACE

    ABSTRACT INTERFACE
        FUNCTION gradf_intrfc(x) RESULT(df)
            USE kindSettings, ONLY : rkp
            REAL(rkp), INTENT(IN) :: x(:)
            REAL(rkp)             :: df(SIZE(x),1)
        END FUNCTION gradf_intrfc
    END INTERFACE

    ABSTRACT INTERFACE
        FUNCTION hessf_intrfc(x) RESULT(d2f)
            USE kindSettings, ONLY : rkp
            REAL(rkp), INTENT(IN) :: x(:)
            REAL(rkp)             :: d2f(SIZE(x),SIZE(x))
        END FUNCTION hessf_intrfc
    END INTERFACE

END MODULE  multiVarFuncIntrfc