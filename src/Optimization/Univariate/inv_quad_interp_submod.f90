!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Optimization/Univariate/inv_quad_interp_submod.f90
!!
!! Inverse quadratic interpolation method for root-finding in 1D
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
!! Inverse quadratic interpolation uses three previous solution guesses to
!! construct a parabolic/quadratic estimate to the solution guess in an
!! iterative procedure. Similar to Muller's method in concept, but does not
!! make use of divided difference tables to weight the solution.
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
SUBMODULE (root_finding_1d) inv_quad_interp_submod

CONTAINS

    MODULE FUNCTION inv_quad_interp(f,x0,x1,x2,tol,nmax) RESULT(xs)
        PROCEDURE(f_1d_intrfc)    :: f
        REAL(rkp),     INTENT(IN) :: x0,x1,x2,tol
        INTEGER(ilkp), INTENT(IN) :: nmax
        REAL(rkp)                 :: xs
        REAL(rkp)                 :: xt0,xt1,xt2,ft0,ft1,ft2
        REAL(rkp)                 :: relerr = HUGE(1.0_rkp)
        REAL(rkp)                 :: abserr = HUGE(1.0_rkp)        
        INTEGER(ilkp)             :: n = 0

        xt0 = x0
        xt1 = x1
        xt2 = x2
        ft0 = f(xt0)
        ft1 = f(xt1)
        ft2 = f(xt2)

        DO WHILE (relerr .gt. tol .and. abserr .gt. tol .and. n .le. nmax)

            xs     = (xt0*ft1*ft2) / ((ft0-ft1)*(ft0-ft2)) + &
                     (xt1*ft0*ft2) / ((ft1-ft0)*(ft1-ft2)) + &
                     (xt2*ft0*ft1) / ((ft2-ft0)*(ft2-ft1))

            abserr = ABS(xs-xt2)
            relerr = abserr/xt2
            xt0    = xt1
            xt1    = xt2
            xt2    = xs
            ft0    = ft1
            ft1    = ft2
            ft2    = f(xt2)
            n      = n + 1

        ENDDO

    END FUNCTION inv_quad_interp

END SUBMODULE inv_quad_interp_submod