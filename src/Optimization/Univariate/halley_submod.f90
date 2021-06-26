!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Optimization/Univariate/halley_submod.f90
!!
!! Halley's iterative method for root-finding in 1D, a 3rd order Householder
!! root-finding method for univariate functions f(x) = 0
!!
!! Input arguments:
!!     f - PROCEDURE handle for the function f(x) for rooting
!!     df - PROCEDURE handle for first derivative function df(x) = f'
!!     df2 - PROCEDURE handle for second derivative function df(x) = f''
!!     x0 - initial guess
!!     tol - a numerical tolerance for changes in guess solution to
!!           serve as a DO WHILE exit criterion for finding root
!!     nmax - maximum number of integer counting iterations for
!!            DO WHILE search loop
!!
!! Output argument:
!!     xs - the best approximation of the root obtained
!!
!! Halley's method is one of the Householder methods that can be formulated 
!! for an aritrary number of derivatives. The more derivative degrees used the
!! the higher the order of accuracy. Halley's method uses the function, first
!! derivative, and second derivative to iteratively solve for the root of the
!! original function. Because the method uses three levels of functions and
!! dervatives, the method is 3rd order accurate/convergent. It is similar to
!! the Newton(-Raphson) method that is a second order Householder method using
!! only the function and its first derivative to iteratively solve for the
!! root of the function at f(x) = 0.
!!
!!*****************************************************************************
SUBMODULE (root_finding_1d) halley_submod

CONTAINS

    MODULE FUNCTION halley(f,df,df2,x0,tol,nmax) result(xs)
        PROCEDURE(f_1d_intrfc)    :: f, df,df2
        REAL(rkp),     INTENT(IN) :: x0,tol
        INTEGER(ilkp), INTENT(IN) :: nmax
        REAL(rkp)                 :: xs
        REAL(rkp)                 :: fx,fpx,xt0
        REAL(rkp)                 :: relerr = HUGE(1.0_rkp)
        REAL(rkp)                 :: abserr = HUGE(1.0_rkp)
        INTEGER(ilkp)             :: n = 0

        xs = x0

        DO WHILE (relerr .gt. tol .and. abserr .gt. tol .and. n .le. nmax)

            xt0    = xs
            fx     = f(xt0)
            fpx    = df(xt0) 
            xs     = xt0 - 2.0*fx*fpx / (2.0*fpx*fpx - fx*df2(xt0))
            abserr = ABS(xs-xt0)
            relerr = abserr/xt0
            n      = n + 1

        ENDDO

    END FUNCTION halley


END SUBMODULE halley_submod