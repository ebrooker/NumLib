MODULE funcs
USE kindSettings, ONLY : rkp, ilkp, iskp
IMPLICIT NONE
PRIVATE

	PUBLIC :: f1, f2

CONTAINS

	FUNCTION f1(x, args) RESULT(y)
		REAL(rkp)                       :: y
		REAL(rkp), INTENT(IN)           :: x
		REAL(rkp), OPTIONAL, INTENT(IN) :: args(:)

		y = 4.0*x**3 + 3.0*x**2 + 6.0*x + 7.0

	END FUNCTION f1


	FUNCTION f2(x,args) RESULT(y)
		REAL(rkp)                       :: y
		REAL(rkp), INTENT(IN)           :: x
		REAL(rkp), OPTIONAL, INTENT(IN) :: args(:)

		y = 5.0*COS(3.0*x) + args(1)

	END FUNCTION f2

END MODULE funcs