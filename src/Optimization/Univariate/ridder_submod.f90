!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Optimization/Univariate/ridder_submod.f90
!!
!! Ridder's bracketing method for root-finding in 1D
!!
!! Input arguments:
!!     f - PROCEDURE handle for the function f(x) for rooting
!!     x0 - left hand of initial bracket search area
!!     x2 - right hand of initial bracket search area
!!     tol - a numerical tolerance for changes in bracker size to
!!           serve as a DO WHILE exit criterion for finding root
!!     nmax - maximum number of integer counting iterations for
!!            DO WHILE search loop
!!
!! Output argument:
!!     xs - the best approximation of the root obtained
!!
!! Ridder's method takes initial domain bracket to search for the 
!! root of f(x) provided as input. See the wikipedia page on
!! Ridder's method for more details. In gist, we
!!
!! Once the space between a and b is sufficiently small, we exit
!! and expect to receive the correct root... or, we hit n==nmax and
!! exit due to running out of iterations, likely have a hard
!! problem to solve or a bad initial bracket interval.
!!
!!*****************************************************************************
SUBMODULE (root_finding_1d) ridder_submod

CONTAINS

    MODULE FUNCTION ridder(f,x0,x2,tol,nmax) result(xs)
        PROCEDURE(f_1d_intrfc)    :: f
        REAL(rkp),     INTENT(IN) :: x0,x2,tol
        INTEGER(ilkp), INTENT(IN) :: nmax
        REAL(rkp)                 :: xs
        REAL(rkp)                 :: xt0,xt1,xt2
        REAL(rkp)                 :: fx0,fx1,fx2,fxs
        REAL(rkp)                 :: w,fdd,neg,pos,dnm
        REAL(rkp)                 :: err = HUGE(1.0_rkp)
        INTEGER(ilkp)             :: n = 0

        fx0 = f(x0)
        fx2 = f(x2)

        IF ( fx0*fx2 .ge. 0.0_rkp) STOP "Root not contained within initial bracket"

        IF ( x0 > x2 ) THEN
            !... Swap bracket bounds if left point is greater than right point
            xt0 = x2
            xt2 = x0
            fx1 = fx0
            fx0 = fx2
            fx2 = fx1
        ELSE
            xt0 = x0
            xt2 = x2
        END IF

        DO WHILE ( err .gt. tol .and. n .lt. nmax )
 
            xt1 = 0.5_rkp*(xt0+xt2)
            fx1 = f(xt1)
            xs  = xt1 + (xt1-xt0) * &
            (SIGN(fx0,1.0_rkp)*fx1)/SQRT(fx1*fx1 - fx0*fx2)

            fxs = f(xs)

            IF ( fxs .gt. 0.0_rkp ) THEN
                xt2 = xs
                fx2 = fxs
            ELSEIF ( fxs .lt. 0.0_rkp ) THEN
                xt0 = xs
                fx0 = fxs
            ELSE
                RETURN
            ENDIF
            
            IF (fx1*fxs .lt. 0.0_rkp) THEN
                xt0 = MIN(xt1,xs)
                fx0 = f(xt0)
                xt2 = MAX(xt1,xs)
                fx2 = f(xt2)
            ENDIF
            err = 0.5_rkp*(xt2-xt0)

        END DO

    END FUNCTION ridder

END SUBMODULE ridder_submod
