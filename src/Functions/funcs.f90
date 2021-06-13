!!*********************************************
!!
!! Author: Ezra Brooker
!! 2021
!! 
!! src/Functions/funcs.f90
!!
!! Module provides a number of function
!! expressions for use in the UNIT TESTS of 
!! NumLib. The 2D Rosenbrock equations are 
!! stored here.
!!
!!*********************************************
MODULE funcs
USE kindSettings, ONLY : rkp, ilkp, iskp
IMPLICIT NONE
PRIVATE

    PUBLIC :: f1, f2, fx, grad_fx

CONTAINS

    FUNCTION f1(x) RESULT(y)
        REAL(rkp)                       :: y
        REAL(rkp), INTENT(IN)           :: x

        y = 4.0*x**3 + 3.0*x**2 + 6.0*x + 7.0

	END FUNCTION f1


    FUNCTION f2(x) RESULT(y)
        REAL(rkp)                       :: y
        REAL(rkp), INTENT(IN)           :: x

        y = 5.0*COS(3.0*x)

	END FUNCTION f2

    FUNCTION fx(x) RESULT(f)
        ! compute 2D Rosenbrock equation
        REAL(rkp), INTENT(IN) :: x(:)
        REAL(rkp)             :: f
        f = 100.0 * (x(2) - x(1)**2)**2 + (1.0 - x(1))**2
    END FUNCTION fx


    FUNCTION grad_fx(x) RESULT(df)
        ! compute gradient of 2D Rosenbrock equation
        REAL(rkp), INTENT(IN) :: x(:)
        REAL(rkp)             :: df(SIZE(x,DIM=1),1)
        df(1,1) = 400.0 * (x(1)**3 - x(2)*x(1)) + 2.0 * (-1.0 + x(1))
        df(2,1) = 200.0 * (x(2) - x(1)**2)
    END FUNCTION grad_fx



END MODULE funcs