!!*********************************************
!!
!! Author: Ezra Brooker
!! 2021
!! 
!! src/Utils/kindSettings.f90
!!
!! Defining specific integer and real
!! precision settings for the package
!!
!!*********************************************
MODULE kindSettings
USE, INTRINSIC :: ISO_FORTRAN_ENV
IMPLICIT NONE
PUBLIC
    INTEGER, PARAMETER :: rkp  = real64
    INTEGER, PARAMETER :: iskp = int16
    INTEGER, PARAMETER :: ilkp = int64
END MODULE kindSettings
