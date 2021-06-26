!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Optimization/Univariate/brent_dekker_submod.f90
!!
!! BRENT or BRENT-DEKKER bracketing method for root-finding in 1D
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
!!     xs - the best approximation of the root obtained
!!
!! BRENT-DEKKER takes initial domain bracket to search for the root
!! of f(x) provided as input. See the wikipedia page on the Brent
!! root finding method for more details. In gist, we combine the 
!! different parts of secant, inverse quadratic interpolation,
!! bisection, and false position methods to create a superior
!! method in terms of successful root finding.
!!
!! Once the space between a and b is sufficiently small, we exit
!! and expect to receive the correct root... or, we hit n==nmax and
!! exit due to running out of iterations, likely have a hard
!! problem to solve or a bad initial bracket interval.
!!
!!*****************************************************************************
SUBMODULE (root_finding_1d) brent_dekker_submod

CONTAINS

    MODULE FUNCTION brent_dekker(f,a,b,tol,nmax) result(xs)
        PROCEDURE(f_1d_intrfc)    :: f
        REAL(rkp),     INTENT(IN) :: a,b,tol
        INTEGER(ilkp), INTENT(IN) :: nmax
        REAL(rkp)                 :: xs
        REAL(rkp)                 :: xa,xb,xc,xd,xm
        REAL(rkp)                 :: fa,fb,fc,fs
        REAL(rkp)                 :: err = HUGE(1.0_rkp)
        INTEGER(ilkp)             :: mf = 1, n = 0
        LOGICAL                   :: mk(5)

        fa = f(a)
        fb = f(b)

        IF ( fa*fb .ge. 0.0_rkp ) STOP "Root not contained within bracket"

        !... Switch the brackets if needed to make them left to right in
        !... function space
        IF ( ABS(fb) .le. ABS(fa) ) THEN
            xa = b
            xb = a
            fc = fa
            fa = fb
            fb = fc
        ELSE
            xa = a
            xb = b
        END IF

        xc = xa

        DO WHILE ( err .gt. tol .and. n .lt. nmax )

            fc = f(xc)
            
            IF ( fa .ne. fc .and. fb .ne. fc ) THEN
                !... use inverse quadratic interpolation
                xs = (xa*fb*fc) / ((fa-fb)*(fa-fc)) + &
                     (xb*fa*fc) / ((fb-fa)*(fb-fc)) + &
                     (xc*fa*fb) / ((fc-fa)*(fc-fb))
            ELSE
                !... use secant method
                xs = xb - fb*(xb-xa) / (fb-fa)
            ENDIF

            !... optimize the bracketing
            xm = 0.25_rkp*(3.0_rkp*xa+xb)

            !... set the 5 masking conditions
            mk(1) = ( xs .lt. xm  .or.           xs .gt. xb                 )
            mk(2) = ( mf .eq. 1   .and.  ABS(xs-xb) .ge. ABS(xb-xc)*0.5_rkp )
            mk(3) = ( mf .eq. 0   .and.  ABS(xs-xb) .ge. ABS(xc-xd)*0.5_rkp )
            mk(4) = ( mf .eq. 1   .and.  ABS(xb-xc) .lt. ABS(tol)           )
            mk(5) = ( mf .eq. 0   .and.  ABS(xc-xd) .lt. ABS(tol)           )

            !... check if any mask conditions are true
            IF ( mk(1) .or. mk(2) .or. mk(3) .or. mk(4) .or. mk(5) ) THEN
                !... use bisection method
                xs = 0.5_rkp * (xb+xa)
                mf = 1
            ELSE
                mf = 0
            ENDIF

            fs = f(xs)
            xd = xc
            xc = xb

            IF (fa*fs .lt. 0.0_rkp) THEN
                xb = xs
                fb = fs
            ELSE
                xa = xs
                fa = fs
            END IF

            IF ( ABS(fb) .le. ABS(fa) ) THEN
                !... swap bracket bounds
                xm = xa
                xa = xb
                xb = xm

                !... swap function evals of bracket bounds
                xm = fa
                fa = fb
                fb = xm
            END IF

            IF ( fb .eq. 0.0_rkp ) THEN
                xs = xb
                RETURN
            ENDIF

            IF ( fs .eq. 0.0_rkp) RETURN
            
            err = xb - xa
            n   = n + 1

        END DO 

    END FUNCTION brent_dekker

END SUBMODULE brent_dekker_submod
