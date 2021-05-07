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

    PUBLIC    :: cntrDiff_d1, cntrDiff_d2, cntrDiff_d3, cntrDiff_d4

CONTAINS

    FUNCTION cntrDiff_d1(f,x,h,order,args) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order
        REAL(rkp), OPTIONAL, INTENT(IN) :: args(:)

        REAL(rkp) :: numer
        INTEGER   :: iord

        iord = -1
        IF (PRESENT(order)) iord = order

        SELECT CASE (iord)

        CASE (2)
            y = 0.5 * (f(x+h,args) - f(x-h,args)) / h

        CASE (4)
            y = ( twothrd*(f(x+h,args    ) - f(x-h,args    )) &
                + onetwlf*(f(x-2.0*h,args) - f(x+2.0*h,args)) ) / h

        CASE default
            y = 0.5 * (f(x+h,args) - f(x-h,args)) / h

        END SELECT

    END FUNCTION cntrDiff_d1


    FUNCTION cntrDiff_d2(f,x,h,order,args) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order
        REAL(rkp), OPTIONAL, INTENT(IN) :: args(:)

        REAL(rkp) :: numer
        INTEGER   :: iord

        iord = -1
        IF (PRESENT(order)) iord = order

        SELECT CASE (iord)

        CASE (2)
            y = (f(x+h,args) - 2.0*f(x,args) + f(x-h,args)) / h**2

        CASE (4)
            y = (- onetwlf * (f(x+2.0*h,args) + f(x-2.0*h,args)) &
                 + forthrd * (f(x+h,    args) + f(x-h,    args)) &
                 - fivehlf *  f(x,      args)) / h**2

        CASE default
            y = (f(x+h,args) - 2.0*f(x,args) + f(x-h,args)) / h**2

        END SELECT

    END FUNCTION cntrDiff_d2

    FUNCTION cntrDiff_d3(f,x,h,order,args) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order
        REAL(rkp), OPTIONAL, INTENT(IN) :: args(:)

        INTEGER   :: iord

        iord = -1
        IF (PRESENT(order)) iord = order

        SELECT CASE (iord)

        CASE (2)
            y = (- 0.5*f(x-2.0*h,args) &
                 +     f(x-h,    args) &
                 -     f(x+h,    args) &
                 + 0.5*f(x+2.0*h,args) ) / h**3

        CASE (4)
            y = (  oneeght*(f(x-3.0*h,args) - f(x+3.0*h,args)) &
                 +          f(x+2.0*h,args) - f(x-2.0*h,args)  &
                 + thrteht*(f(x-h,    args) - f(x+h,    args)) ) / h**3


        CASE default
            y = (- 0.5*f(x-2.0*h,args) &
                 +     f(x-h,    args) &
                 -     f(x+h,    args) &
                 + 0.5*f(x+2.0*h,args) ) / h**3

        END SELECT

    END FUNCTION cntrDiff_d3

    FUNCTION cntrDiff_d4(f,x,h,order,args) RESULT(y)
        REAL(rkp)                       :: y
        PROCEDURE(f_1d_intrfc)          :: f
        REAL(rkp),           INTENT(IN) :: x,h
        INTEGER,   OPTIONAL, INTENT(IN) :: order
        REAL(rkp), OPTIONAL, INTENT(IN) :: args(:)

        INTEGER   :: iord

        iord = -1
        IF (PRESENT(order)) iord = order

        SELECT CASE (iord)

        CASE (2)
            y = (      f(x-2.0*h,args) &
                 - 4.0*f(x-h,    args) &
                 + 6.0*f(x,      args) &
                 - 4.0*f(x+h,    args) &
                 +     f(x+2.0*h,args) ) / h**4

        CASE (4)
            y = (  56.0* f(x,      args)                    &
                 -      (f(x+3.0*h,args) + f(x-3.0*h,args)) &
                 + 12.0*(f(x+2.0*h,args) + f(x-2.0*h,args)) &
                 - 39.0*(f(x+h,    args) + f(x-h,    args)) ) / (6.0*h**4)

        CASE default
            y = (      f(x-2.0*h,args) &
                 - 4.0*f(x-h,    args) &
                 + 6.0*f(x,      args) &
                 - 4.0*f(x+h,    args) &
                 +     f(x+2.0*h,args) ) / h**4

        END SELECT

    END FUNCTION cntrDiff_d4

END MODULE centralDifferences
