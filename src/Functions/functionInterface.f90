MODULE functionInterface
USE kindSettings, ONLY : rkp
IMPLICIT NONE
PRIVATE

    PUBLIC :: f_1d_intrfc

    ABSTRACT INTERFACE
        FUNCTION f_1d_intrfc(Xx) RESULT(Yy)
            USE kindSettings, ONLY : rkp
            REAL(rkp)                       :: Yy
            REAL(rkp), INTENT(IN)           :: Xx
        END FUNCTION f_1d_intrfc
    END INTERFACE

END MODULE functionInterface