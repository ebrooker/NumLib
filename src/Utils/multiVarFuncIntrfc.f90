MODULE  multiVarFuncIntrfc
IMPLICIT NONE
PRIVATE

    PUBLIC :: f_intrfc, gradf_intrfc, hessf_intrfc

    ABSTRACT INTERFACE
        FUNCTION f_intrfc(x) RESULT(f)
            USE kindSettings, ONLY : rkp
            REAL(rkp), INTENT(IN) :: x(:)
            REAL(rkp)             :: f
        END FUNCTION f_intrfc
    END INTERFACE

    ABSTRACT INTERFACE
        FUNCTION gradf_intrfc(x) RESULT(df)
            USE kindSettings, ONLY : rkp
            REAL(rkp), INTENT(IN) :: x(:)
            REAL(rkp)             :: df(SIZE(x),1)
        END FUNCTION gradf_intrfc
    END INTERFACE

    ABSTRACT INTERFACE
        FUNCTION hessf_intrfc(x) RESULT(d2f)
            USE kindSettings, ONLY : rkp
            REAL(rkp), INTENT(IN) :: x(:)
            REAL(rkp)             :: d2f(SIZE(x),SIZE(x))
        END FUNCTION hessf_intrfc
    END INTERFACE

END MODULE  multiVarFuncIntrfc