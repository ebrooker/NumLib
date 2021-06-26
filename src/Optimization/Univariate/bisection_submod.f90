!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Optimization/Univariate
!!
!! BISECTION, traditional bracketing method for root-finding in 1D
!!
!! Input arguments:
!!     f - PROCEDURE handle for the function f(x) for rooting
!!     a - left hand of initial bracket search area
!!     b - right hand of initial bracket search area
!!     tol - a numerical tolerance for changes in bracker size to
!!           serve as a DO WHILE exit criterion for finding root
!!     nmax - maximum number of integer counting iterations for
!!            DO WHILE search loop
!!
!! Output argument:
!!     c - the best approximation of the root obtained
!!
!! BISECTION takes an initial domain bracket to search for the root
!! of f(x) provided as input. The method evalutes a midpoint "c"
!! and determines which side of the root it is on. The bracket is
!! adjust to so that a=c if c is left of root and b=c otherwise.
!!
!! Once the space between a and b is sufficiently small, we exit
!! and expect to receive the correct root... or, we hit n==nmax and
!! exit due to running out of iterations, likely have a hard
!! problem to solve or a bad initial bracket interval.
!!
!!*****************************************************************************
SUBMODULE (root_finding_1d) bisection_submod

CONTAINS

    MODULE FUNCTION bisection(f,a,b,tol,nmax) result(c)
        PROCEDURE(f_1d_intrfc)    :: f
        REAL(rkp),     INTENT(IN) :: a,b,tol
        INTEGER(ilkp), INTENT(IN) :: nmax
        REAL(rkp)                 :: c
        REAL(rkp)                 :: cl,cr,yl,yc
        REAL(rkp)                 :: err = HUGE(1.0_rkp)
        INTEGER(ilkp)             :: n = 0

        IF ( f(a)*f(b) .ge. 0.0_rkp) STOP "Root not contained within initial bracket"

        IF ( a > b ) THEN
            cl = b
            cr = a
        ELSE
            cl = a
            cr = b
        END IF

        yl = f(cl)

        DO WHILE (err .gt. tol .and. n .le. nmax)
            
            c  = 0.5*(cr+cl)
            yc = f(c)
            IF (yc .eq. 0.0_rkp) RETURN

            IF ( yc*yl .ge. 0.0_rkp ) THEN
                cl = c
                yl = yc
            ELSE
                cr = c
            ENDIF
            err = 0.5e0*(cr-cl)
            n   = n + 1

        ENDDO

    END FUNCTION bisection


END SUBMODULE bisection_submod
