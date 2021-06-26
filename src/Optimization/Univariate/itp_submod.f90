!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Optimization/Univariate/itp_submod.f90
!!
!! ITP bracketing method for root-finding in 1D
!!
!! Input arguments:
!!     f - PROCEDURE handle for the function f(x) for rooting
!!     a - left hand of initial bracket search area
!!     b - right hand of initial bracket search area
!!     eps - special parameter
!!     kap1 - special parameter
!!     kap2 - special parameter
!!     n0 - special parameter
!!     tol - a numerical tolerance for changes in bracker size to
!!           serve as a DO WHILE exit criterion for finding root
!!     nmax - maximum number of integer counting iterations for
!!            DO WHILE search loop
!!
!! Output argument:
!!     c - the best approximation of the root obtained
!!
!! ITP takes initial domain bracket to search for the root of f(x)
!! provided as input and special parameters. See the wikipedia page
!! on the ITP method for more details. In gist, we combine the 
!! new bracket point evaluations of BISECTION and FALSE POSITION to
!! improve a final estimator for the root.
!!
!! Once the space between a and b is sufficiently small, we exit
!! and expect to receive the correct root... or, we hit n==nmax and
!! exit due to running out of iterations, likely have a hard
!! problem to solve or a bad initial bracket interval.
!!
!!*****************************************************************************
SUBMODULE (root_finding_1d) itp_submod

CONTAINS

    MODULE FUNCTION itp(f,a,b,eps,kap1,kap2,n0,tol,nmax) result(c)
        PROCEDURE(f_1d_intrfc)    :: f
        REAL(rkp),     INTENT(IN) :: a,b,eps,kap1,kap2,n0
        REAL(rkp),     INTENT(IN) :: tol
        INTEGER(ilkp), INTENT(IN) :: nmax
        REAL(rkp)                 :: c
        REAL(rkp)                 :: cl,cr,cf,ct,citp
        REAL(rkp)                 :: sigma,delta,rhok,n12
        REAL(rkp)                 :: yl,yr,yc,yitp
        REAL(rkp)                 :: err = HUGE(1.0_rkp)
        INTEGER(ilkp)             :: n = 0

        yl = f(a)
        yr = f(b)

        IF ( yl*yr .ge. 0.0_rkp) STOP "Root not contained within initial bracket"

        IF ( a > b ) THEN
            cl = b
            cr = a
            yc = yl
            yl = yr
            yr = yc
        ELSE
            cl = a
            cr = b
        END IF

        n12 = log10((b-a)/(2.0_rkp*eps))/log10(2.0_rkp)

        DO WHILE ( err .gt. 2.0_rkp*eps .and. n .lt. nmax )

            !... Calculate bisection and regula falsi points
            c  = 0.50_rkp*(cl+cr)
            cf = (cl*yr - b*yl) / (yr - yl)

            !... Truncate by perturbing estimator towards center
            sigma = SIGN(c-cf, c)
            delta = MIN( kap1*ABS(b-a)**kap2, ABS(c-cf) )
            ct    = cf + sigma*delta

            !... Project the estimator to minmax interval for citp = c - sigma*rhok
            rhok = MIN( (eps*(2.0_rkp)**(n12+n0-n) - (cr-cl)*0.5_rkp), ABS(ct-c) )
            citp = c - sigma*rhok
            yitp = f(citp)

            !... Update interval
            IF( yitp .gt. 0.0_rkp ) THEN
                cr = citp
                yr = yitp
            ELSEIF  (yitp .lt. 0.0_rkp ) THEN
                cl = citp
                yl = yitp
            ELSE
                cl = citp
                cr = citp
            ENDIF
            err = cr-cl
            n   = n + 1

        END DO

        c = (cl+cr)*0.5_rkp

    END FUNCTION itp


END SUBMODULE itp_submod