!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Optimization/Univariate/muler_submod.f90
!!
!! Muller's interpolation method for root-finding in 1D
!!
!! Input arguments:
!!     f - PROCEDURE handle for the function f(x) for rooting
!!     x0 - initial guess
!!     x1 - next initial guess (ideally towards solution)
!!     x2 - final initial guess (ideally towards solution)
!!     tol - a numerical tolerance for changes in guess solution to
!!           serve as a DO WHILE exit criterion for finding root
!!     nmax - maximum number of integer counting iterations for
!!            DO WHILE search loop
!!
!! Output argument:
!!     xs - the best approximation of the root obtained
!!
!! Muller's method takes three initial guess for the root of a polynomial
!! f(x) and uses these three points to construct a parabola to estimate a new
!! guess estimating the nearest root of f(x) = 0. Once we obtain this estimate
!! we reset the initial guess: x0 = x1, x1 = x2, x2 = new guess. We proceed
!! iteratively until the new guess begins to change by a sufficiently small
!! enough tolerance value.
!!
!! This method is similar to the traditional secant root-finding method, but
!! uses three points to estimate an iterative update to the solution root, 
!! whereas secant method uses two points (a differnce between parabolic and
!! linear interpolation). This method makes use of solving a quadratic formula
!! with a parameter constructed from a sum of divided differences for the
!! previous three guesses for any given iteration.
!!
!! As such, this method makes use of the divided difference method, names
!! divide_diff.
!!
!! The three previous root guesses will be used in function evaluations for
!! f(x_(k-3)), f(x_(k-2)), f(x_(k-1)), where k is the current iteration.
!!
!! These function evaluations are done prior to the iterative procedure, in
!! that way, only one additional function evaluation is required for each
!! iteration for the new guess. We have each iteration:
!!     - fx0 = fx1   (reassign k-2 guess func eval variable to k-3 variable)
!!     - fx1 = fx2   (reassign k-1 guess func eval variable to k-2 variable)
!!     - fx2 = f(xs) (perform  f(new guess) eval and assign to k-1 variable)
!!
!!*****************************************************************************
SUBMODULE (root_finding_1d) muller_submod


CONTAINS

    MODULE FUNCTION muller(f,x0,x1,x2,tol,nmax) RESULT (xs)
        USE divided_differences, ONLY : divide_diff
        PROCEDURE(f_1d_intrfc)    :: f
        REAL(rkp),     INTENT(IN) :: x0,x1,x2,tol
        INTEGER(ilkp), INTENT(IN) :: nmax
        REAL(rkp)                 :: xs
        REAL(rkp)                 :: xt0,xt1,xt2,fx0,fx1,fx2
        REAL(rkp)                 :: w,fdd,neg,pos,dnm
        REAL(rkp)                 :: relerr = HUGE(1.0_rkp)
        REAL(rkp)                 :: abserr = HUGE(1.0_rkp)        
        INTEGER(ilkp)             :: n = 0

        !... set local guess variables to input guesses
        xt0 = x0
        xt1 = x1
        xt2 = x2

        !... obtain the initial three function evaluations
        fx0 = f(xt0)
        fx1 = f(xt1)
        fx2 = f(xt2)

        DO WHILE (relerr .gt. tol .and. abserr .gt. tol .and. n .le. nmax)

            !... Compute divided difference sums
            w = divided_difference([fx2,fx1],[xt2,xt1]) + &
                divided_difference([fx2,fx0],[xt2,xt0]) - &
                divided_difference([fx1,fx0],[xt1,xt0])

            !... compute larger divided difference for all 3 guesses
            fdd = divided_difference([fx2,fx1,fx0],[xt2,xt1,xt0])

            !... Compute both quadratic solutions and find largest ABS(denom)
            pos = w + SQRT(w**2 - 4.0_rkp*fx2*fdd)
            neg = w - SQRT(w**2 - 4.0_rkp*fx2*fdd)
            dnm = MAX(ABS(pos),ABS(neg))

            !... Update the solution and guess variables
            xs     = xt2 - 2.0_rkp*fx2/dnm
            abserr = ABS(xs-xt2)
            relerr = abserr/xt2
            xt0    = xt1
            xt1    = xt2
            xt2    = xs
            fx0    = fx1
            fx1    = fx2
            fx2    = f(xs)
            n      = n + 1

        ENDDO

    END FUNCTION muller


END SUBMODULE muller_submod
