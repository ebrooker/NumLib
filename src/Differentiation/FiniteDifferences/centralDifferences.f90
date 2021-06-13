!!********************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Differentiation/FiniteDifferences/centralDifferences.f90
!!
!! Implementations of Central Finite Difference methods for numerical
!! differentiation.
!!
!! Example of 1st order accurate formula for a 1st order derivative:
!!        f'(x) = (f(x+h) - f(x-h)) / 2h
!!
!! There are Central Diff methods for up to 4th derivatives with a 2nd
!! and 4th order accurate formula for each order of derivative.
!!
!! The specific method is appended with a "_d<INT>", where INT is the
!! integer number representing the order of the derivative the method
!! is designed to approximate.
!!
!! The procedure "bwrdDiff" is the driver differentiator for all of
!! the Central Diff methods. It functions by passing in the usual args
!! of the, function (f), derivative point (x), estimation spacing (h),
!! with two additional args, the order of derivative (deriv) and the
!! order accuracy method desired (order)
!!
!!********************************************************************
MODULE centralDifferences
USE kindSettings, ONLY : rkp, ilkp, iskp
USE functionInterface, ONLY : f_1d_intrfc
IMPLICIT NONE
PRIVATE

    ! Some constants used for simplicity sake
    REAL(rkp), PARAMETER :: twothrd = 2.0/3.0
    REAL(rkp), PARAMETER :: forthrd = 4.0/3.0
    REAL(rkp), PARAMETER :: fivehlf = 5.0/2.0
    REAL(rkp), PARAMETER :: oneeght = 1.0/8.0
    
    REAL(rkp), PARAMETER :: onetwlf = 1.0/12.0
    
    REAL(rkp), PARAMETER :: thrteht = 13.0/8.0

    PUBLIC  :: cntrDiff
    PRIVATE :: cntrDiff_d1, cntrDiff_d2, cntrDiff_d3, cntrDiff_d4

CONTAINS


    FUNCTION cntrDiff(f,x,h,deriv,order) RESULT(y)
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
          WRITE(*,*) "[Central Differentiation Warning] &
          Approximating 1st order derivative by default"
          WRITE(*,*) '    Set "deriv" to a value between &
          1 and 4 in "cntrDiff" argument list to disable &
          this warning'
        ENDIF

        IF (PRESENT(order)) THEN
          iord = order
        ELSE
          iord = 1 ! default 1st order accuracy
          WRITE(*,*) "[Central Differentiation Warning] &
          Using 2ndt order accurate method by default"
          WRITE(*,*) '    Set "order" to a value of either &
          2 or 4 in "cntrDiff" argument list to disable &
          this warning'
        ENDIF

        ! Use "deriv" to select what order of derivative
        ! is to be numerically approximated (up to f^(4)(x))
        SELECT CASE (ideriv)
            CASE (1)
                ! 1st order derivative
                y = cntrDiff_d1(f,x,h,iord)
            CASE (2)
                ! 2nd order derivative
                y = cntrDiff_d2(f,x,h,iord)
            CASE (3)
                ! 3rd order derivative
                y = cntrDiff_d3(f,x,h,iord)
            CASE (4)
                ! 4th order derivative
                y = cntrDiff_d4(f,x,h,iord)
            CASE default
                ! Abort due to error
                write(*,*) "[Central Differentiation Error] &
                Cannot numerically differentiate derivative of order:", deriv
                stop "!!*** ABORTING ***!!"
        END SELECT
    
    END FUNCTION cntrDiff

    FUNCTION cntrDiff_d1(f,x,h,order) RESULT(y)
        !!
        !! 1st order derivative method
        !! Input arguments match driver method,
        !! see "cntrDiff" procedure for details
        !!    
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order

        SELECT CASE (order)

            CASE (2)
                ! 2nd order accurate formula
                y = 0.5 * (f(x+h) - f(x-h)) / h

            CASE (4)
                ! 4th order accurate formula
                y = ( twothrd*(f(x+h    ) - f(x-h    )) &
                    + onetwlf*(f(x-2.0*h) - f(x+2.0*h)) ) / h

            CASE default
                ! Abort due to error
                write(*,*) "[Central Differentiation Error] &
                Cannot use differentiation formula of accuracy order:", order
                stop "!!*** ABORTING ***!!"

        END SELECT

    END FUNCTION cntrDiff_d1


    FUNCTION cntrDiff_d2(f,x,h,order) RESULT(y)
        !!
        !! 2nd order derivative method
        !! Input arguments match driver method,
        !! see "cntrDiff" procedure for details
        !!    
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order

        SELECT CASE (order)

            CASE (2)
                ! 2nd order accurate formula
                y = (f(x+h) - 2.0*f(x) + f(x-h)) / h**2

            CASE (4)
                ! 4th order accurate formula
                y = (- onetwlf * (f(x+2.0*h) + f(x-2.0*h)) &
                     + forthrd * (f(x+h    ) + f(x-h    )) &
                     - fivehlf *  f(x      )) / h**2

            CASE default
                ! Abort due to error
                write(*,*) "[Central Differentiation Error] &
                Cannot use differentiation formula of accuracy order:", order
                stop "!!*** ABORTING ***!!"

        END SELECT

    END FUNCTION cntrDiff_d2

    FUNCTION cntrDiff_d3(f,x,h,order) RESULT(y)
        !!
        !! 3rd order derivative method
        !! Input arguments match driver method,
        !! see "cntrDiff" procedure for details
        !!    
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order

        SELECT CASE (order)

            CASE (2)
                ! 2nd order accurate formula
                y = (- 0.5*f(x-2.0*h) &
                     +     f(x-h    ) &
                     -     f(x+h    ) &
                     + 0.5*f(x+2.0*h) ) / h**3

            CASE (4)
                ! 4th order accurate formula
                y = (  oneeght*(f(x-3.0*h) - f(x+3.0*h)) &
                     +          f(x+2.0*h) - f(x-2.0*h)  &
                     + thrteht*(f(x-h    ) - f(x+h    )) ) / h**3


            CASE default
                ! Abort due to error
                write(*,*) "[Central Differentiation Error] &
                Cannot use differentiation formula of accuracy order:", order
                stop "!!*** ABORTING ***!!"

        END SELECT

    END FUNCTION cntrDiff_d3

    FUNCTION cntrDiff_d4(f,x,h,order) RESULT(y)
        !!
        !! 4th order derivative method
        !! Input arguments match driver method,
        !! see "cntrDiff" procedure for details
        !!    
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order

        SELECT CASE (order)

            CASE (2)
                ! 2nd order accurate formula
                y = (      f(x-2.0*h) &
                     - 4.0*f(x-h    ) &
                     + 6.0*f(x      ) &
                     - 4.0*f(x+h    ) &
                     +     f(x+2.0*h) ) / h**4

            CASE (4)
                ! 4th order accurate formula
                y = (  56.0* f(x      )               &
                     -      (f(x+3.0*h) + f(x-3.0*h)) &
                     + 12.0*(f(x+2.0*h) + f(x-2.0*h)) &
                     - 39.0*(f(x+h    ) + f(x-h    )) ) / (6.0*h**4)

            CASE default
                ! Abort due to error
                write(*,*) "[Central Differentiation Error] &
                Cannot use differentiation formula of accuracy order:", order
                stop "!!*** ABORTING ***!!"

        END SELECT

    END FUNCTION cntrDiff_d4

END MODULE centralDifferences
