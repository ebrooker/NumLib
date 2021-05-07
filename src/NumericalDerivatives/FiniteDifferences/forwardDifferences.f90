MODULE forwardDifferences
USE kindSettings, ONLY : rkp, ilkp, iskp
USE functionInterface, ONLY : f_1d_intrfc
IMPLICIT NONE
PRIVATE

    PUBLIC :: fwrdDiff_d1, fwrdDiff_d2, fwrdDiff_d3, fwrdDiff_d4

CONTAINS

    FUNCTION fwrdDiff_d1(f,x,h,order,args) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order
        REAL(rkp), OPTIONAL, INTENT(IN) :: args(:)

        INTEGER   :: iord

        iord = -1
        IF (PRESENT(order)) iord = order

        SELECT CASE (iord)

        CASE (1)
            y = (f(x+h,args) - f(x,args)) / h

        CASE (2)
            y = 0.5 * ( 4.0*f(x+h,args) &
                      - f(x+2.0*h,args) &
                      - 3.0*f(x,  args) ) / h

        CASE default
            y = (f(x+h,args) - f(x,args)) / h

        END SELECT

    END FUNCTION fwrdDiff_d1


    FUNCTION fwrdDiff_d2(f,x,h,order,args) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order
        REAL(rkp), OPTIONAL, INTENT(IN) :: args(:)

        INTEGER   :: iord

        iord = -1
        IF (PRESENT(order)) iord = order

        SELECT CASE (iord)

        CASE (1)
            y = (     f(x+2.0*h,args) &
                - 2.0*f(x+h,    args) &
                +     f(x,      args) ) / h**2

        CASE (2)
            y = ( 2.0*f(x,      args) &
                - 5.0*f(x+h,    args) &
                + 4.0*f(x+2.0*h,args) &
                -     f(x+3.0*h,args) ) / h**2


        CASE default
            y = (     f(x+2.0*h,args) &
                - 2.0*f(x+h,    args) &
                +     f(x,      args) ) / h**2

        END SELECT

    END FUNCTION fwrdDiff_d2

    FUNCTION fwrdDiff_d3(f,x,h,order,args) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order
        REAL(rkp), OPTIONAL, INTENT(IN) :: args(:)

        INTEGER   :: iord

        iord = -1
        IF (PRESENT(order)) iord = order

        SELECT CASE (iord)

        CASE (1)
            y = (      f(x+3.0*h,args) &
                 - 3.0*f(x+2.0*h,args) &
                 + 3.0*f(x+h,    args) &
                 -     f(x,      args) ) / h**3

        CASE (2)
            y = (-  2.5*f(x,      args) &
                 +  9.0*f(x+h,    args) &
                 - 12.0*f(x+2.0*h,args) &
                 +  7.0*f(x+3.0*h,args) &
                 -  1.5*f(x+4.0*h,args) ) / h**3


        CASE default
            y = (      f(x+3.0*h,args) &
                 - 3.0*f(x+2.0*h,args) &
                 + 3.0*f(x+h,    args) &
                 -     f(x,      args) ) / h**3

        END SELECT

    END FUNCTION fwrdDiff_d3

    FUNCTION fwrdDiff_d4(f,x,h,order,args) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order
        REAL(rkp), OPTIONAL, INTENT(IN) :: args(:)

        INTEGER   :: iord

        iord = -1
        IF (PRESENT(order)) iord = order

        SELECT CASE (iord)

        CASE (1)
            y = (      f(x+4.0*h,args) &
                 - 4.0*f(x+3.0*h,args) &
                 + 6.0*f(x+2.0*h,args) &
                 - 4.0*f(x+h,    args) &
                 +     f(x,      args) ) / h**4

        CASE (2)
            y = (+  3.0*f(x,      args) &
                 - 14.0*f(x+h,    args) &
                 + 26.0*f(x+2.0*h,args) &
                 - 24.0*f(x+3.0*h,args) &
                 + 11.0*f(x+4.0*h,args) &
                 -  2.0*f(x+5.0*h,args) ) / h**4


        CASE default
            y = (      f(x+4.0*h,args) &
                 - 4.0*f(x+3.0*h,args) &
                 + 6.0*f(x+2.0*h,args) &
                 - 4.0*f(x+h,    args) &
                 +     f(x,      args) ) / h**4

        END SELECT

    END FUNCTION fwrdDiff_d4

END MODULE forwardDifferences