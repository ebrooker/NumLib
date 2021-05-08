MODULE centralDifferences
USE kindSettings, ONLY : rkp, ilkp, iskp
USE functionInterface, ONLY : f_1d_intrfc
IMPLICIT NONE
PRIVATE

    
    REAL(rkp), PARAMETER :: twothrd = 2.0/3.0
    REAL(rkp), PARAMETER :: forthrd = 4.0/3.0
    REAL(rkp), PARAMETER :: fivehlf = 5.0/2.0
    REAL(rkp), PARAMETER :: oneeght = 1.0/8.0
    
    REAL(rkp), PARAMETER :: onetwlf = 1.0/12.0
    
    REAL(rkp), PARAMETER :: thrteht = 13.0/8.0

    PUBLIC    :: cntrDiff, cntrDiff_d1, cntrDiff_d2, cntrDiff_d3, cntrDiff_d4

CONTAINS


    FUNCTION cntrDiff(f,x,h,deriv,order) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: deriv,order

        INTEGER   :: ideriv, iord
        ideriv = -1
        IF (PRESENT(deriv)) ideriv = deriv

        iord = -1
        IF (PRESENT(order)) iord   = order

        SELECT CASE (ideriv)
            CASE (1)
                y = cntrDiff_d1(f,x,h,iord)
            CASE (2)
                y = cntrDiff_d2(f,x,h,iord)
            CASE (3)
                y = cntrDiff_d3(f,x,h,iord)
            CASE (4)
                y = cntrDiff_d4(f,x,h,iord)
            CASE default
                y = cntrDiff_d1(f,x,h,iord)
        END SELECT
    
    END FUNCTION cntrDiff

    FUNCTION cntrDiff_d1(f,x,h,order) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order

        INTEGER ::          iord = -1
        IF (PRESENT(order)) iord = order

        SELECT CASE (iord)

        CASE (2)
            y = 0.5 * (f(x+h) - f(x-h)) / h

        CASE (4)
            y = ( twothrd*(f(x+h    ) - f(x-h    )) &
                + onetwlf*(f(x-2.0*h) - f(x+2.0*h)) ) / h

        CASE default
            y = 0.5 * (f(x+h) - f(x-h)) / h

        END SELECT

    END FUNCTION cntrDiff_d1


    FUNCTION cntrDiff_d2(f,x,h,order) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order

        INTEGER ::          iord = -1
        IF (PRESENT(order)) iord = order

        SELECT CASE (iord)

        CASE (2)
            y = (f(x+h) - 2.0*f(x) + f(x-h)) / h**2

        CASE (4)
            y = (- onetwlf * (f(x+2.0*h) + f(x-2.0*h)) &
                 + forthrd * (f(x+h    ) + f(x-h    )) &
                 - fivehlf *  f(x      )) / h**2

        CASE default
            y = (f(x+h) - 2.0*f(x) + f(x-h)) / h**2

        END SELECT

    END FUNCTION cntrDiff_d2

    FUNCTION cntrDiff_d3(f,x,h,order) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order

        INTEGER ::          iord = -1
        IF (PRESENT(order)) iord = order

        SELECT CASE (iord)

        CASE (2)
            y = (- 0.5*f(x-2.0*h) &
                 +     f(x-h    ) &
                 -     f(x+h    ) &
                 + 0.5*f(x+2.0*h) ) / h**3

        CASE (4)
            y = (  oneeght*(f(x-3.0*h) - f(x+3.0*h)) &
                 +          f(x+2.0*h) - f(x-2.0*h)  &
                 + thrteht*(f(x-h    ) - f(x+h    )) ) / h**3


        CASE default
            y = (- 0.5*f(x-2.0*h) &
                 +     f(x-h    ) &
                 -     f(x+h    ) &
                 + 0.5*f(x+2.0*h) ) / h**3

        END SELECT

    END FUNCTION cntrDiff_d3

    FUNCTION cntrDiff_d4(f,x,h,order) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order

        INTEGER ::          iord = -1
        IF (PRESENT(order)) iord = order

        SELECT CASE (iord)

        CASE (2)
            y = (      f(x-2.0*h) &
                 - 4.0*f(x-h    ) &
                 + 6.0*f(x      ) &
                 - 4.0*f(x+h    ) &
                 +     f(x+2.0*h) ) / h**4

        CASE (4)
            y = (  56.0* f(x      )               &
                 -      (f(x+3.0*h) + f(x-3.0*h)) &
                 + 12.0*(f(x+2.0*h) + f(x-2.0*h)) &
                 - 39.0*(f(x+h    ) + f(x-h    )) ) / (6.0*h**4)

        CASE default
            y = (      f(x-2.0*h) &
                 - 4.0*f(x-h    ) &
                 + 6.0*f(x      ) &
                 - 4.0*f(x+h    ) &
                 +     f(x+2.0*h) ) / h**4

        END SELECT

    END FUNCTION cntrDiff_d4

END MODULE centralDifferences
