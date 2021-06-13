!!*******************************************************************************
!! Author: Ezra Brooker
!! 2021
!! 
!! NumLib/src/Differentiation/RichardsonExtrapolation/richardsonExtrapolation.f90
!!
!! Module contains the simple Richardson Extrapolation method for taking
!! two numerical derivative approximations made from differentation spacings
!! dx and dx/2 and creating a more approximate result from the two less accurate
!! approximations.
!!
!!*******************************************************************************
MODULE richardsonExtrapolation
USE kindSettings, ONLY : rkp, ilkp, iskp
USE functionInterface, ONLY : f_1d_intrfc
IMPLICIT NONE
PRIVATE

    REAL(rkp), PARAMETER :: onethrd = 1.0/3.0

    PUBLIC :: richExtrap

CONTAINS

    FUNCTION richExtrap(d1,d2) result(y)
        !!
        !! d1 (REAL) : derivative approximation using dx
        !! d2 (REAL) : derivative approximation using dx/2
        !! y(real)   : higher order accurate Richardson estimate
        !!             created using d1 and d2
        !!
        !! Richardson Extrapolation formula: (1/3) * ( 4*(df/dx) - (df/(dx/2)) )
        !!
        !!**********************************************************************
        REAL(rkp), INTENT(IN) :: d1, d2
        REAL(rkp)             :: y

        y = onethrd * (4.0*d2 - d1)

	END FUNCTION richExtrap

END MODULE richardsonExtrapolation