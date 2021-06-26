!!******************************************************!!
!!
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/Optimization/Univariate/root_finding.f90
!!
!! Collates the different univariate (1D) root finding
!! methods implemented in this part of the Optimization
!! subpackage. It currently contains bracketing methods,
!! interpolation/approximation methods, and householder
!! methods that use derivative information.
!!
!!
!!******************************************************!!
MODULE root_finding_1d
USE kindSettings, ONLY : rkp, ilkp, iskp
USE functionInterface, only : f_1d_intrfc
IMPLICIT NONE
PRIVATE

    !! Bracketing methods
    PUBLIC :: bisection, brent_dekker, false_position, itp, ridder

    !! Interpolation/Approximation methods
    PUBLIC ::inv_quad_interp, muller, secant, steffenson

    !! Iterative Householder methods
    PUBLIC :: halley, newton_1d 



    INTERFACE !! Bracketing methods

        MODULE FUNCTION bisection(f,a,b,tol,nmax) result(c)
            PROCEDURE(f_1d_intrfc)    :: f
            REAL(rkp),     INTENT(IN) :: a,b,tol
            INTEGER(ilkp), INTENT(IN) :: nmax
            REAL(rkp)                 :: c
        END FUNCTION bisection

        MODULE FUNCTION brent_dekker(f,a,b,tol,nmax) result(xs)
            PROCEDURE(f_1d_intrfc)    :: f
            REAL(rkp),     INTENT(IN) :: a,b,tol
            INTEGER(ilkp), INTENT(IN) :: nmax
            REAL(rkp)                 :: xs
        END FUNCTION brent_dekker

        MODULE FUNCTION false_position(f,a,b,tol,nmax) result(c)
            PROCEDURE(f_1d_intrfc)    :: f
            REAL(rkp),     INTENT(IN) :: a,b,tol
            INTEGER(ilkp), INTENT(IN) :: nmax
            REAL(rkp)                 :: c
        END FUNCTION false_position

        MODULE FUNCTION itp(f,a,b,eps,kap1,kap2,n0,tol,nmax) result(c)
            PROCEDURE(f_1d_intrfc)    :: f
            REAL(rkp),     INTENT(IN) :: a,b,eps,kap1,kap2,n0
            REAL(rkp),     INTENT(IN) :: tol
            INTEGER(ilkp), INTENT(IN) :: nmax
            REAL(rkp)                 :: c   
        END FUNCTION itp

        MODULE FUNCTION ridder(f,x0,x2,tol,nmax) result(xs)
            PROCEDURE(f_1d_intrfc)    :: f
            REAL(rkp),     INTENT(IN) :: x0,x2,tol
            INTEGER(ilkp), INTENT(IN) :: nmax
            REAL(rkp)                 :: xs
        END FUNCTION ridder

    END INTERFACE !! End bracketing methods


    INTERFACE !! Interpolation/Finite Diff methods

        MODULE FUNCTION muller(f,x0,x1,x2,tol,nmax) RESULT (xs)
            PROCEDURE(f_1d_intrfc)    :: f
            REAL(rkp),     INTENT(IN) :: x0,x1,x2,tol
            INTEGER(ilkp), INTENT(IN) :: nmax
            REAL(rkp)                 :: xs
        END FUNCTION muller

        MODULE FUNCTION inv_quad_interp(f,x0,x1,x2,tol,nmax) RESULT(xs)
            PROCEDURE(f_1d_intrfc)    :: f
            REAL(rkp),     INTENT(IN) :: x0,x1,x2,tol
            INTEGER(ilkp), INTENT(IN) :: nmax
            REAL(rkp)                 :: xs
        END FUNCTION inv_quad_interp

        MODULE FUNCTION secant(f,x0,x1,tol,nmax) RESULT(xs)
            PROCEDURE(f_1d_intrfc)    :: f
            REAL(rkp),     INTENT(IN) :: x0,x1,tol
            INTEGER(ilkp), INTENT(IN) :: nmax
            REAL(rkp)                 :: xs
        END FUNCTION secant

        MODULE FUNCTION steffenson(f,x0,tol,nmax) result(xs)
            PROCEDURE(f_1d_intrfc)    :: f
            REAL(rkp),     INTENT(IN) :: x0,tol
            INTEGER(ilkp), INTENT(IN) :: nmax
            REAL(rkp)                 :: xs
        END FUNCTION steffenson

    END INTERFACE !! End interpolation methods


    INTERFACE !! Derivative bearing iterative methods

        MODULE FUNCTION halley(f,df,df2,x0,tol,nmax) result(xs)
            PROCEDURE(f_1d_intrfc)    :: f, df,df2
            REAL(rkp),     INTENT(IN) :: x0,tol
            INTEGER(ilkp), INTENT(IN) :: nmax
            REAL(rkp)                 :: xs
        END FUNCTION halley

        MODULE FUNCTION newton_1d(f,df,x0,tol,nmax) result(xs)
            PROCEDURE(f_1d_intrfc)    :: f, df
            REAL(rkp),     INTENT(IN) :: x0,tol
            INTEGER(ilkp), INTENT(IN) :: nmax
            REAL(rkp)                 :: xs
        END FUNCTION newton_1d

    END INTERFACE !! End iterative derivative methods

END MODULE root_finding_1d
