!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Optimization/Univariate/secant_submod.f90
!!
!! Secant method for root-finding in 1D
!!
!! Input arguments:
!!     f - PROCEDURE handle for the function f(x) for rooting
!!     x0 - initial guess
!!     x1 - next initial guess (ideally towards solution)
!!     tol - a numerical tolerance for changes in guess solution to
!!           serve as a DO WHILE exit criterion for finding root
!!     nmax - maximum number of integer counting iterations for
!!            DO WHILE search loop
!!
!! Output argument:
!!     xs - the best approximation of the root obtained
!!
!! The secant method uses two previous guesses to construct a linear 
!! approximation to the iterative update (essentially doing a 1st order
!! finite difference on the derivative of the Newton's method equation).
!!
!! The two previous root guesses will be used in function evaluations for
!! f(x_(k-2)), f(x_(k-1)), where k is the current iteration.
!!
!! These function evaluations are done prior to the iterative procedure, in
!! that way, only one additional function evaluation is required for each
!! iteration for the new guess. We have each iteration:
!!     - fx0 = fx1   (reassign k-1 guess func eval variable to k-2 variable)
!!     - fx1 = f(xs) (perform  f(new guess) eval and assign to k-1 variable)
!!
!!*****************************************************************************
SUBMODULE (root_finding_1d) secant_submod

CONTAINS

    MODULE FUNCTION secant(f,x0,x1,tol,nmax) RESULT(xs)
        PROCEDURE(f_1d_intrfc)    :: f
        REAL(rkp),     INTENT(IN) :: x0,x1,tol
        INTEGER(ilkp), INTENT(IN) :: nmax
        REAL(rkp)                 :: xs
        REAL(rkp)                 :: xt0,xt1,ft0,ft1
        REAL(rkp)                 :: relerr = HUGE(1.0_rkp)
        REAL(rkp)                 :: abserr = HUGE(1.0_rkp)        
        INTEGER(ilkp)             :: n = 0

        xt0 = x0
        xt1 = x1
        ft0 = f(xt0)

        DO WHILE (relerr .gt. tol .and. abserr .gt. tol .and. n .le. nmax)

            ft1    = f(xt1)
            xs     = xt1 - ft1*(xt1-xt0) / (ft1 - ft0)
            abserr = ABS(xs-xt1)
            relerr = abserr/xt1
            xt0    = xt1
            xt1    = xs
            ft0    = ft1
            n      = n + 1

        ENDDO

    END FUNCTION secant

END SUBMODULE secant_submod