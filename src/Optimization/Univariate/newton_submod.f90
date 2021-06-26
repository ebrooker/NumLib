!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Optimization/Univariate/newton_submod.f90
!!
!! Newton's iterative method for root-finding in 1D, a 2nd order Householder
!! root-finding method for univariate functions f(x) = 0
!!
!! Also known as the Newton-Raphson iterative method 
!!
!! Input arguments:
!!     f - PROCEDURE handle for the function f(x) for rooting
!!     df - PROCEDURE handle for first derivative function df(x) = f'
!!     x0 - initial guess
!!     tol - a numerical tolerance for changes in guess solution to
!!           serve as a DO WHILE exit criterion for finding root
!!     nmax - maximum number of integer counting iterations for
!!            DO WHILE search loop
!!
!! Output argument:
!!     xs - the best approximation of the root obtained
!!
!! Newton's method is one of the Householder methods that can be formulated 
!! for an aritrary number of derivatives. The more derivative degrees used the
!! the higher the order of accuracy. Newton's method uses the function and first
!! derivative to iteratively solve for the root of the original function. 
!! Because the method uses two levels of functions and dervatives, the method 
!! is 2nd order accurate/convergent.
!!
!!*****************************************************************************
SUBMODULE (root_finding_1d) newton_submod

CONTAINS

    MODULE FUNCTION newton_1d(f,df,x0,tol,nmax) result(xs)
        PROCEDURE(f_1d_intrfc)    :: f, df
        REAL(rkp),     INTENT(IN) :: x0,tol
        INTEGER(ilkp), INTENT(IN) :: nmax
        REAL(rkp)                 :: xs
        REAL(rkp)                 :: xt0
        REAL(rkp)                 :: relerr = HUGE(1.0_rkp)
        REAL(rkp)                 :: abserr = HUGE(1.0_rkp)
        INTEGER(ilkp)             :: n = 0

        xs = x0

        DO WHILE (relerr .gt. tol .and. abserr .gt. tol .and. n .le. nmax)

            xt0    = xs
            xs     = xt0 - f(xt0)/df(xt0)
            abserr = ABS(xs-xt0)
            relerr = abserr/xt0
            n      = n + 1

        ENDDO

    END FUNCTION newton_1d

END SUBMODULE newton_submod
