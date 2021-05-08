PROGRAM test
USE kindSettings, ONLY : rkp
USE funcs
USE derivatives
USE proctorClass, ONLY : proctor_1d
IMPLICIT NONE
    
    TYPE(proctor_1d) :: proctor
    REAL(rkp) :: x=2.0, h=0.01

    CALL proctor%init(fin=f2,dfin=df2,ideriv=1,xin=x,hin=h)
    CALL proctor%dtest()

    CALL proctor%set_f(fin=f2,dfin=d2f2,ideriv=2)
    CALL proctor%dtest()

    CALL proctor%set_f(fin=f2,dfin=d3f2,ideriv=3)
    CALL proctor%dtest()

    CALL proctor%set_f(fin=f2,dfin=d4f2,ideriv=4)
    CALL proctor%dtest()

END PROGRAM test