!!********************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Differentiation/FiniteDifferences/forwardDifferences.f90
!!
!! Implementations of Forward Finite Difference methods for numerical
!! differentiation.
!!
!! Example of 1st order accurate formula for a 1st order derivative:
!!        f'(x) = (f(x+h) - f(x)) / h
!!
!! There are Forward Diff methods for up to 4th derivatives with a 1st
!! and 2nd order accurate formula for each order of derivative.
!!
!! The specific method is appended with a "_d<INT>", where INT is the
!! integer number representing the order of the derivative the method
!! is designed to approximate.
!!
!! The procedure "fwrdDiff" is the driver differentiator for all of
!! the Forward Diff methods. It functions by passing in the usual args
!! of the, function (f), derivative point (x), estimation spacing (h),
!! with two additional args, the order of derivative (deriv) and the
!! order accuracy method desired (order)
!!
!!********************************************************************
MODULE forwardDifferences
USE kindSettings, ONLY : rkp, ilkp, iskp
USE functionInterface, ONLY : f_1d_intrfc
IMPLICIT NONE
PRIVATE

    PUBLIC  :: fwrdDiff
    PRIVATE :: fwrdDiff_d1, fwrdDiff_d2, fwrdDiff_d3, fwrdDiff_d4

CONTAINS

    FUNCTION fwrdDiff(f,x,h,deriv,order) RESULT(y)
        !!
        !! f (REAL): Function handle (designed with a specific) abstract
        !!           interface to allow 1D functions of f(x) form pass
        !!           through
        !!
        !! x (REAL): numerical point where derivative is take at
        !!
        !! h (REAL): finite difference spacing for the numerical estimation
        !!
        !! derv (INT): Order of the derivative passed as an integer, i.e.
        !!             deriv=1 for 1st order derivative, deriv=2 for 2nd
        !!             order, etc
        !!
        !! order (INT): Order of accuracy method desired for numerical 
        !!              derivative, i.e., order=1 selects the 1st order
        !!              accurate method, order=2 selects 2nd order accurate
        !!              method (slightly more computationally expensive)
        !!
        !!*****************************************************************
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: deriv,order

        INTEGER   :: ideriv, iord

        IF (PRESENT(deriv)) THEN
          ideriv = deriv
        ELSE
          ideriv = 1 ! default 1st order derivative
          WRITE(*,*) "[Forward Differentiation Warning] &
          Approximating 1st order derivative by default"
          WRITE(*,*) '    Set "deriv" to a value between &
          1 and 4 in "fwrdDiff" argument list to disable &
          this warning'
        ENDIF

        IF (PRESENT(order)) THEN
          iord = order
        ELSE
          iord = 1 ! default 1st order accuracy
          WRITE(*,*) "[Forward Differentiation Warning] &
          Using 1st order accurate method by default"
          WRITE(*,*) '    Set "order" to a value of either &
          1 or 2 in "fwrdDiff" argument list to disable &
          this warning'
        ENDIF

        ! Use "deriv" to select what order of derivative
        ! is to be numerically approximated (up to f^(4)(x))
        SELECT CASE (ideriv)
            CASE (1)
                ! 1st order derivative
                y = fwrdDiff_d1(f,x,h,iord)
            CASE (2)
                ! 2nd order derivative
                y = fwrdDiff_d2(f,x,h,iord)
            CASE (3)
                ! 3rd order derivative
                y = fwrdDiff_d3(f,x,h,iord)
            CASE (4)
                ! 4th order derivative
                y = fwrdDiff_d4(f,x,h,iord)
            CASE default
                ! Abort due to error
                write(*,*) "[Forward Differentiation Error] &
                Cannot numerically differentiate derivative of order:", deriv
                stop "!!*** ABORTING ***!!"
        END SELECT
    
    END FUNCTION fwrdDiff


    FUNCTION fwrdDiff_d1(f,x,h,order) RESULT(y)
        !!
        !! 1st order derivative method
        !! Input arguments match driver method,
        !! see "fwrdDiff" procedure for details
        !!
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order

        SELECT CASE (order)

            CASE (1)
                ! 1st order accurate formula
                y = (f(x+h) - f(x)) / h

            CASE (2)
                ! 2nd order accurate formula
                y = 0.5 * ( 4.0*f(x+h    ) &
                          -     f(x+2.0*h) &
                          - 3.0*f(x      ) ) / h

            CASE default
                ! Abort due to error
                write(*,*) "[Forward Differentiation Error] &
                Cannot use differentiation formula of accuracy order:", order
                stop "!!*** ABORTING ***!!"

        END SELECT

    END FUNCTION fwrdDiff_d1


    FUNCTION fwrdDiff_d2(f,x,h,order) RESULT(y)
        !!
        !! 2nd order derivative method
        !! Input arguments match driver method,
        !! see "fwrdDiff" procedure for details
        !!
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order

        SELECT CASE (order)

            CASE (1)
                ! 1st order accurate formula
                y = (     f(x+2.0*h) &
                    - 2.0*f(x+h    ) &
                    +     f(x      ) ) / h**2

            CASE (2)
                ! 2nd order accurate formula
                y = ( 2.0*f(x      ) &
                    - 5.0*f(x+h    ) &
                    + 4.0*f(x+2.0*h) &
                    -     f(x+3.0*h) ) / h**2


            CASE default
                ! Abort due to error
                write(*,*) "[Forward Differentiation Error] &
                Cannot use differentiation formula of accuracy order:", order
                stop "!!*** ABORTING ***!!"

        END SELECT

    END FUNCTION fwrdDiff_d2

    FUNCTION fwrdDiff_d3(f,x,h,order) RESULT(y)
        !!
        !! 3rd order derivative method
        !! Input arguments match driver method,
        !! see "fwrdDiff" procedure for details
        !!
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order

        SELECT CASE (order)

            CASE (1)
                ! 1st order accurate formula
                y = (      f(x+3.0*h) &
                     - 3.0*f(x+2.0*h) &
                     + 3.0*f(x+h    ) &
                     -     f(x      ) ) / h**3

            CASE (2)
                ! 2nd order accurate formula
                y = (-  2.5*f(x      ) &
                     +  9.0*f(x+h    ) &
                     - 12.0*f(x+2.0*h) &
                     +  7.0*f(x+3.0*h) &
                     -  1.5*f(x+4.0*h) ) / h**3


            CASE default
                ! Abort due to error
                write(*,*) "[Forward Differentiation Error] &
                Cannot use differentiation formula of accuracy order:", order
                stop "!!*** ABORTING ***!!"

        END SELECT

    END FUNCTION fwrdDiff_d3

    FUNCTION fwrdDiff_d4(f,x,h,order) RESULT(y)
        !!
        !! 4th order derivative method
        !! Input arguments match driver method,
        !! see "fwrdDiff" procedure for details
        !!
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order

        SELECT CASE (order)

            CASE (1)
                ! 1st order accurate formula
                y = (      f(x+4.0*h) &
                     - 4.0*f(x+3.0*h) &
                     + 6.0*f(x+2.0*h) &
                     - 4.0*f(x+h    ) &
                     +     f(x      ) ) / h**4

            CASE (2)
                ! 2nd order accurate formula
                y = (+  3.0*f(x      ) &
                     - 14.0*f(x+h    ) &
                     + 26.0*f(x+2.0*h) &
                     - 24.0*f(x+3.0*h) &
                     + 11.0*f(x+4.0*h) &
                     -  2.0*f(x+5.0*h) ) / h**4


            CASE default
                ! Abort due to error
                write(*,*) "[Forward Differentiation Error] &
                Cannot use differentiation formula of accuracy order:", order
                stop "!!*** ABORTING ***!!"

        END SELECT

    END FUNCTION fwrdDiff_d4

END MODULE forwardDifferences