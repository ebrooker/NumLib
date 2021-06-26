!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Optimization/Univariate/false_position_submod
!!
!! FALSE POSITION or REGULA FALSI bracketing method for 1D
!! root-finding
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
!! REGULA FALSI takes initial domain bracket to search for the root
!! of f(x) provided as input. See the wikipedia page on the False
!! Position method for more details. In gist, it is exactly the 
!! BISECTION method, but using a linear interpolant, instead of a
!! simple midpoint to determine the new bracket point.
!!
!! Once the space between a and b is sufficiently small, we exit
!! and expect to receive the correct root... or, we hit n==nmax and
!! exit due to running out of iterations, likely have a hard
!! problem to solve or a bad initial bracket interval.
!!
!!*****************************************************************************
SUBMODULE (root_finding_1d) false_position_submod

CONTAINS

    MODULE FUNCTION false_position(f,a,b,tol,nmax) result(c)
        PROCEDURE(f_1d_intrfc)    :: f
        REAL(rkp),     INTENT(IN) :: a,b,tol
        INTEGER(ilkp), INTENT(IN) :: nmax
        REAL(rkp)                 :: c
        REAL(rkp)                 :: cl,cr,yl,yr,yc
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

        DO WHILE (err .gt. tol .and. n .le. nmax)
            
            c  = (cl*yr - cr*yl) / (yr - yl)
            yc = f(c)

            IF (yc .eq. 0.0_rkp) RETURN

            IF ( yc*yl .ge. 0.0_rkp) THEN
                cl = c
                yl = yc
            ELSE
                cr = c
                yr = yc
            ENDIF
            err = (cr-cl)
            n   = n + 1

        ENDDO

    END FUNCTION false_position

END SUBMODULE false_position_submod
