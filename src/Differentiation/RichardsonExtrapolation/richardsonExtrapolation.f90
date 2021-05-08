MODULE richardsonExtrapolation
USE kindSettings, ONLY : rkp, ilkp, iskp
USE functionInterface, ONLY : f_1d_intrfc
IMPLICIT NONE
PRIVATE

	REAL(rkp), PARAMETER :: onethrd = 1.0/3.0

	PUBLIC :: richExtrap

CONTAINS

	FUNCTION richExtrap(d1,d2) result(y)
		REAL(rkp), INTENT(IN) :: d1, d2
		REAL(rkp)             :: y
		
		y = onethrd * (4.0*d2 - d1)

	END FUNCTION richExtrap


END MODULE richardsonExtrapolation