PROGRAM test
USE kindSettings, ONLY : rkp
USE funcs
USE derivatives
USE proctorDerivative, ONLY : proctor_d1f
USE proctorOptimize,   ONLY : proctorOpt
IMPLICIT NONE
    
    TYPE(proctor_d1f) :: proctorD
    TYPE(proctorOpt ) :: proctorOp
    REAL(rkp) :: x=2.0, h=0.01
    REAL(rkp) :: B(2,2), iB(2,2), x0(2,1)
    INTEGER   :: m=2

    x0 = RESHAPE([-3.,-4.], SHAPE(x0))
    B  = TRANSPOSE(RESHAPE([1.,0.,0.,1.], SHAPE(B)))
    iB = B


    CALL proctorD%init(fin=f2,dfin=df2,ideriv=1,xin=x,hin=h)
    CALL proctorD%dtest()

    ! CALL proctorD%set_f(fin=f2,dfin=d2f2,ideriv=2)
    ! CALL proctorD%dtest()

    ! CALL proctorD%set_f(fin=f2,dfin=d3f2,ideriv=3)
    ! CALL proctorD%dtest()

    ! CALL proctorD%set_f(fin=f2,dfin=d4f2,ideriv=4)
    ! CALL proctorD%dtest()

    CALL proctorOp%init(fx,grad_fx,m,x0,B,iB)
    CALL proctorOp%test_steepDescent()
    CALL proctorOp%set_x0(x0)
    CALL proctorOp%test_BGFS()

END PROGRAM test