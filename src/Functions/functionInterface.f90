MODULE functionInterface
USE kindSettings, ONLY : rkp, ilkp, iskp
IMPLICIT NONE
! PRIVATE

    PUBLIC :: f_1d_noargs, f_1d_intrfc

    ABSTRACT INTERFACE
        FUNCTION f_1d_intrfc(Xx,Aargs) RESULT(Yy)
            USE kindSettings, ONLY : rkp, ilkp, iskp
            REAL(rkp)                       :: Yy
            REAL(rkp), INTENT(IN)           :: Xx
            REAL(rkp), OPTIONAL, INTENT(IN) :: Aargs(:)
        END FUNCTION f_1d_intrfc

        FUNCTION f_1d_noargs(Xx) RESULT(Yy)
            USE kindSettings, ONLY : rkp, ilkp, iskp
            REAL(rkp)                       :: Yy
            REAL(rkp), INTENT(IN)           :: Xx
        END FUNCTION f_1d_noargs

    END INTERFACE

END MODULE functionInterface