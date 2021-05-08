MODULE backwardDifferences
USE kindSettings, ONLY : rkp, ilkp, iskp
USE functionInterface, ONLY : f_1d_intrfc
IMPLICIT NONE
PRIVATE

    PUBLIC :: bwrdDiff, bwrdDiff_d1, bwrdDiff_d2, bwrdDiff_d3, bwrdDiff_d4

CONTAINS


    FUNCTION bwrdDiff(f,x,h,deriv,order) RESULT(y)
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
                y = bwrdDiff_d1(f,x,h,iord)
            CASE (2)
                y = bwrdDiff_d2(f,x,h,iord)
            CASE (3)
                y = bwrdDiff_d3(f,x,h,iord)
            CASE (4)
                y = bwrdDiff_d4(f,x,h,iord)
            CASE default
                y = bwrdDiff_d1(f,x,h,iord)
        END SELECT
    
    END FUNCTION bwrdDiff


    FUNCTION bwrdDiff_d1(f,x,h,order) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order

        INTEGER          :: iord = -1
        IF (PRESENT(order)) iord = order

        SELECT CASE (iord)

        CASE (1)
            y = (f(x) - f(x-h)) / h

        CASE (2)
            y = 0.5 * ( 3.0*f(x      ) &
                      - 4.0*f(x-h    ) &
                      +     f(x-2.0*h) ) / h

        CASE default
            y = (f(x) - f(x-h)) / h

        END SELECT

    END FUNCTION bwrdDiff_d1


    FUNCTION bwrdDiff_d2(f,x,h,order) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order

        INTEGER          :: iord = -1
        IF (PRESENT(order)) iord = order

        SELECT CASE (iord)

        CASE (1)
            y = (     f(x-2.0*h) &
                - 2.0*f(x-h    ) &
                +     f(x      ) ) / h**2

        CASE (2)
            y = ( 2.0*f(x      ) &
                - 5.0*f(x-h    ) &
                + 4.0*f(x-2.0*h) &
                -     f(x-3.0*h) ) / h**2

        CASE default
            y = (     f(x-2.0*h) &
                - 2.0*f(x-h    ) &
                +     f(x      ) ) / h**2

        END SELECT

    END FUNCTION bwrdDiff_d2

    FUNCTION bwrdDiff_d3(f,x,h,order) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order

        INTEGER          :: iord = -1
        IF (PRESENT(order)) iord = order

        SELECT CASE (iord)

        CASE (1)
            y = (-     f(x-3.0*h) &
                 + 3.0*f(x-2.0*h) &
                 - 3.0*f(x-h    ) &
                 +     f(x      ) ) / h**3

        CASE (2)
            y = (+  2.5*f(x      ) &
                 -  9.0*f(x-h    ) &
                 + 12.0*f(x-2.0*h) &
                 -  7.0*f(x-3.0*h) &
                 +  1.5*f(x-4.0*h) ) / h**3


        CASE default
            y = (-     f(x-3.0*h) &
                 + 3.0*f(x-2.0*h) &
                 - 3.0*f(x-h    ) &
                 +     f(x      ) ) / h**3

        END SELECT

    END FUNCTION bwrdDiff_d3

    FUNCTION bwrdDiff_d4(f,x,h,order) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order

        INTEGER          :: iord = -1
        IF (PRESENT(order)) iord = order

        SELECT CASE (iord)

        CASE (1)
            y = (      f(x-4.0*h) &
                 - 4.0*f(x-3.0*h) &
                 + 6.0*f(x-2.0*h) &
                 - 4.0*f(x-h    ) &
                 +     f(x      ) ) / h**4

        CASE (2)
            y = (+  3.0*f(x      ) &
                 - 14.0*f(x-h    ) &
                 + 26.0*f(x-2.0*h) &
                 - 24.0*f(x-3.0*h) &
                 + 11.0*f(x-4.0*h) &
                 -  2.0*f(x-5.0*h) ) / h**4


        CASE default
            y = (      f(x-4.0*h) &
                 - 4.0*f(x-3.0*h) &
                 + 6.0*f(x-2.0*h) &
                 - 4.0*f(x-h    ) &
                 +     f(x      ) ) / h**4

        END SELECT

    END FUNCTION bwrdDiff_d4

END MODULE backwardDifferences