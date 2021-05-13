MODULE lineSearch
USE kindSettings, ONLY : rkp
USE linearAlgebra, ONLY : outerProd, innerProd, shermanMorrison
USE multiVarFuncIntrfc, ONLY : f_intrfc, gradf_intrfc
IMPLICIT NONE
PRIVATE

    PUBLIC :: BGFS
    PUBLIC :: steepDescent
    PUBLIC :: backtrack

CONTAINS

    SUBROUTINE steepDescent(gradf,x,p)
        PROCEDURE(gradf_intrfc) :: gradf
        REAL(rkp),  INTENT(IN ) :: x(:,:)
        REAL(rkp),  INTENT(OUT) :: p(:,:)
        
        p = -gradf(x(:,1))

    END SUBROUTINE steepDescent


    SUBROUTINE BGFS(gradf,first,x,xi,p,B,iB)
        PROCEDURE(gradf_intrfc)  :: gradf
        LOGICAL,   INTENT(IN   ) :: first
        REAL(rkp), INTENT(IN   ) :: x(:,:), xi(:,:)
        REAL(rkp), INTENT(OUT  ) :: p(:,:)
        REAL(rkp), INTENT(INOUT) :: B(:,:), iB(:,:)
        REAL(rkp), ALLOCATABLE   :: s(:,:), u(:,:), v(:,:)
        REAL(rkp)                :: alpha, beta
        REAL(rkp)                :: sTu, isTu, t1
        INTEGER                  :: m

        ! If first iteration, use initial iB and exit
        IF (first) THEN
            p = - MATMUL(iB,gradf(x(:,1)))
            RETURN
        ENDIF

        ! ALLOCATE local helper (column) vectors
        m = SIZE(x,DIM=1)
        IF(ALLOCATED(s)) DEALLOCATE(s)
        IF(ALLOCATED(u)) DEALLOCATE(u)
        IF(ALLOCATED(v)) DEALLOCATE(v)
        ALLOCATE(s(m,1),source=0.0_rkp)
        ALLOCATE(u(m,1),source=0.0_rkp)
        ALLOCATE(v(m,1),source=0.0_rkp)

        ! Compute helper vectors and scalars
        s = x - xi
        u = gradf(x(:,1)) - gradf(xi(:,1))
        v = MATMUL(B,s)
        alpha =  1.0 / DOT_PRODUCT(u(:,1),s(:,1))
        beta  = -1.0 / DOT_PRODUCT(s(:,1),v(:,1))

        ! Update approximation for Hessian
        B = B + alpha * MATMUL(u,TRANSPOSE(u)) &
              +  beta * MATMUL(v,TRANSPOSE(v))

        ! Invert approx Hess with Sherman Morrison formula
        ! formula is applied twice to get solution we seel
        sTu  = DOT_PRODUCT(s(:,1),u(:,1))
        isTu = 1.0/sTu
        t1   = sTu + INNERPROD(MATMUL(TRANSPOSE(u),iB),u)
        iB = iB - ( MATMUL( MATMUL(iB,u), TRANSPOSE(s) ) &
                  + MATMUL( MATMUL(s,TRANSPOSE(u)), iB ) ) * isTu

        iB = iB + MATMUL(s,TRANSPOSE(s)) * t1 * (isTu**2)

        ! Calculate solution descent direction
        p = - MATMUL(iB,gradf(x(:,1)))

    END SUBROUTINE BGFS


    SUBROUTINE backtrack(f,gradf,x,p,alpha0,rho,gamma,alpha)
        PROCEDURE(gradf_intrfc) :: gradf
        PROCEDURE(f_intrfc    ) :: f
        REAL(rkp),  INTENT(IN ) :: x(:,:), p(:,:)
        REAL(rkp),  INTENT(IN ) :: alpha0, rho, gamma
        REAL(rkp),  INTENT(OUT) :: alpha
        REAL(rkp)               :: f1,f2,f2a,f2bp

        alpha = alpha0
        f1    = f(x(:,1)+alpha*p(:,1))
        f2a   = f(x(:,1))
        f2bp  = INNERPROD(gradf(x(:,1)),p)
        f2    = f2a + gamma*alpha*f2bp

        DO WHILE (f1 > f2)
            alpha = alpha*rho
            f1    = f(x(:,1) + alpha*p(:,1))
            f2a   = f(x(:,1))
            f2bp  = INNERPROD(p,gradf(x(:,1)))
            f2    = f2a + gamma*alpha*f2bp
        ENDDO

    END SUBROUTINE backtrack


END MODULE lineSearch
