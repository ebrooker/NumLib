!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Optimization/Univariate/steffenson_submod.f90
!!
!! Steffenson's method for root-finding in 1D
!!
!! Input arguments:
!!     f - PROCEDURE handle for the function f(x) for rooting
!!     x0 - initial guess
!!     tol - a numerical tolerance for changes in guess solution to
!!           serve as a DO WHILE exit criterion for finding root
!!     nmax - maximum number of integer counting iterations for
!!            DO WHILE search loop
!!
!! Output argument:
!!     xs - the best approximation of the root obtained
!!
!! Steffenson's method uses a single previous guess to construct two
!! function evals, fx = f(x0) and gx = (f(x0 + fx)/fx) - 1; after this, 
!! the algorithm proceeds much like Newton's method except gx is used to
!! approximate the derivative evaluation seen in Newton's method.
!!
!!*****************************************************************************
SUBMODULE (root_finding_1d) steffenson_submod

CONTAINS

    MODULE FUNCTION steffenson(f,x0,tol,nmax) result(xs)
        PROCEDURE(f_1d_intrfc)    :: f
        REAL(rkp),     INTENT(IN) :: x0,tol
        INTEGER(ilkp), INTENT(IN) :: nmax
        REAL(rkp)                 :: xs
        REAL(rkp)                 :: fx,gx,xt0
        REAL(rkp)                 :: relerr = HUGE(1.0_rkp)
        REAL(rkp)                 :: abserr = HUGE(1.0_rkp)
        INTEGER(ilkp)             :: n = 0

        xs = x0

        DO WHILE (relerr .gt. tol .and. abserr .gt. tol .and. n .le. nmax)

            xt0    = xs
            fx     = f(xt0)
            gx     = f(xt0+fx)/fx - 1.0_rkp
            xs     = xt0 - (fx/gx)
            abserr = ABS(xs-xt0)
            relerr = abserr/xt0
            n      = n + 1

        ENDDO

    END FUNCTION steffenson

END SUBMODULE steffenson_submod
