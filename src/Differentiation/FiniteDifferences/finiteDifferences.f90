!!*******************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Differentiation/FiniteDifferences/finteDifferences.f90
!!
!! Collate the finite difference modules under one
!! namespace: finiteDifferences
!!
!!*******************************************************************

MODULE finiteDifferences
    USE centralDifferences,  ONLY : cntrDiff
    USE forwardDifferences,  ONLY : fwrdDiff
    USE backwardDifferences, ONLY : bwrdDiff
    IMPLICIT NONE
    !!************************************************
    !!
    !! Central difference functions
    !!
    !!  - cntrDiff    (Driver method, ONLY PUBLIC)
    !!	- cntrDiff_d1 (1st derivative, 2nd/4th orders)
    !!	- cntrDiff_d2 (2nd derivative, 2nd/4th orders)
    !!	- cntrDiff_d3 (3rd derivative, 2nd/4th orders)
    !!	- cntrDiff_d4 (4th derivative, 2nd/4th orders)
    !!
    !! Forward difference functions
    !!
    !!  - fwrdDiff    (Driver method, ONLY PUBLIC)
    !!	- fwrdDiff_d1 (1st derivative, 1st/2nd orders)
    !!	- fwrdDiff_d2 (2nd derivative, 1st/2nd orders)
    !!	- fwrdDiff_d3 (3rd derivative, 1st/2nd orders)
    !!	- fwrdDiff_d4 (4th derivative, 1st/2nd orders)
    !!
    !! Backward difference functions
    !!
    !!  - bwrdDiff    (Driver method, ONLY PUBLIC)
    !!	- bwrdDiff_d1 (1st derivative, 1st/2nd orders)
    !!	- bwrdDiff_d2 (2nd derivative, 1st/2nd orders)
    !!	- bwrdDiff_d3 (3rd derivative, 1st/2nd orders)
    !!	- bwrdDiff_d4 (4th derivative, 1st/2nd orders)
    !!
    !!************************************************
END MODULE finiteDifferences