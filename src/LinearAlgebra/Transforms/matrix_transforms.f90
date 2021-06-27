!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/LinearAlgebra/matrix_systems_solvers.f90
!!
!! A collating module for a number of matrix transformation algorithms, such
!! as, but not limited to, LU decomposition, Gaussian Eliminaition, etc
!!
!!*****************************************************************************
MODULE matrix_transforms
USE kindSettings

IMPLICIT NONE
PRIVATE

    PUBLIC :: gaussian_elimination, gauss_jordan_inverse

    INTERFACE
    
        MODULE SUBROUTINE gaussian_elimination (A,U,b,y,form)
            REAL(rkp),                  INTENT(INOUT) :: A(:,:)
            REAL(rkp),        OPTIONAL, INTENT(OUT  ) :: U(:,:)
            REAL(rkp),        OPTIONAL, INTENT(INOUT) :: b(:)
            REAL(rkp),        OPTIONAL, INTENT(OUT  ) :: y(:)
            CHARACTER(LEN=*), OPTIONAL, INTENT(IN   ) :: form
        END SUBROUTINE gaussian_elimination

        MODULE SUBROUTINE gauss_jordan_inverse (A,B,tol)
            !! Uses Gaussian Elimination procedure in form='rref' mode
            REAL(rkp),           INTENT(IN ) :: A(:,:)
            REAL(rkp),           INTENT(OUT) :: B(:,:)
            REAL(rkp), OPTIONAL, INTENT(IN ) :: tol
        END SUBROUTINE gauss_jordan_inverse


    END INTERFACE

END MODULE matrix_transforms