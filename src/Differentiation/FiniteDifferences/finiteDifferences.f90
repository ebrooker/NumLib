MODULE finiteDifferences
USE centralDifferences
USE forwardDifferences
USE backwardDifferences
IMPLICIT NONE
!!************************************************
!!
!! Central difference functions
!!
!!	- cntrDiff_d1 (1st derivative, 2nd/4th orders)
!!	- cntrDiff_d2 (2nd derivative, 2nd/4th orders)
!!	- cntrDiff_d3 (3rd derivative, 2nd/4th orders)
!!	- cntrDiff_d4 (4th derivative, 2nd/4th orders)
!!
!! Forward difference functions
!!
!!	- fwrdDiff_d1 (1st derivative, 1st/2nd orders)
!!	- fwrdDiff_d2 (2nd derivative, 1st/2nd orders)
!!	- fwrdDiff_d3 (3rd derivative, 1st/2nd orders)
!!	- fwrdDiff_d4 (4th derivative, 1st/2nd orders)
!!
!! Backward difference functions
!!
!!	- bwrdDiff_d1 (1st derivative, 1st/2nd orders)
!!	- bwrdDiff_d2 (2nd derivative, 1st/2nd orders)
!!	- bwrdDiff_d3 (3rd derivative, 1st/2nd orders)
!!	- bwrdDiff_d4 (4th derivative, 1st/2nd orders)
!!
!!************************************************
END MODULE finiteDifferences