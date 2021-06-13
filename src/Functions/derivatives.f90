!!*********************************************
!!
!! Author: Ezra Brooker
!! 2021
!! 
!! src/Functions/derivatives.f90
!!
!! Module provides a number of analytical 
!! derivative expressions for use in the UNIT
!! TESTS of NumLib. There up to 4th order
!! derivative functions contained here.
!!
!!*********************************************
MODULE derivatives
USE kindSettings, ONLY : rkp, ilkp, iskp
IMPLICIT NONE
PRIVATE

    PUBLIC :: df1, df2, d2f2, d2f1, d3f1, d3f2, d4f1, d4f2

CONTAINS

	FUNCTION df1(x) RESULT(y)
		REAL(rkp)                       :: y
		REAL(rkp), INTENT(IN)           :: x

		y = 12.0*x**2 + 6.0*x + 6.0

	END FUNCTION df1

	FUNCTION d2f1(x) RESULT(y)
		REAL(rkp)                       :: y
		REAL(rkp), INTENT(IN)           :: x

		y = 24.0*x + 6.0

	END FUNCTION d2f1

	FUNCTION d3f1(x) RESULT(y)
		REAL(rkp)                       :: y
		REAL(rkp), INTENT(IN)           :: x

		y = 24.0

	END FUNCTION d3f1

	FUNCTION d4f1(x) RESULT(y)
		REAL(rkp)                       :: y
		REAL(rkp), INTENT(IN)           :: x

		y = 0.0

	END FUNCTION d4f1


	FUNCTION df2(x) RESULT(y)
		REAL(rkp)                       :: y
		REAL(rkp), INTENT(IN)           :: x

		y = -15.0*SIN(3.0*x)

	END FUNCTION df2

	FUNCTION d2f2(x) RESULT(y)
		REAL(rkp)                       :: y
		REAL(rkp), INTENT(IN)           :: x

		y = -45.0*COS(3.0*x)

	END FUNCTION d2f2

	FUNCTION d3f2(x) RESULT(y)
		REAL(rkp)                       :: y
		REAL(rkp), INTENT(IN)           :: x

		y = 135.0*SIN(3.0*x)

	END FUNCTION d3f2

	FUNCTION d4f2(x) RESULT(y)
		REAL(rkp)                       :: y
		REAL(rkp), INTENT(IN)           :: x

		y = 405.0*COS(3.0*x)

	END FUNCTION d4f2


END MODULE derivatives