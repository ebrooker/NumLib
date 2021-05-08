MODULE funcs
USE kindSettings, ONLY : rkp, ilkp, iskp
IMPLICIT NONE
PRIVATE

    PUBLIC :: f1, f2

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

END MODULE funcs