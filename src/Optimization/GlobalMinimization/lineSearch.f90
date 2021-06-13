!!********************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Optimization/GlobalMinimization/lineSearch.f90
!!
!! Module used for line search algorithm in global minimization.
!!     
!!       - backtrack: Backtracking procedure for line search step size
!!                   calculation
!!
!!       - steepDescent: Steepest Gradient Descent method
!!
!!       - BFGS: Broyden-Fletcher-Goldfarb-Shanno algorithm
!!
!!********************************************************************
MODULE lineSearch
USE kindSettings, ONLY : rkp
USE linearAlgebra, ONLY : outerProd, innerProd, shermanMorrison
USE multiVarFuncIntrfc, ONLY : f_intrfc, gradf_intrfc
IMPLICIT NONE
PRIVATE

    PUBLIC :: BFGS
    PUBLIC :: steepDescent
    PUBLIC :: backtrack

CONTAINS

    SUBROUTINE steepDescent(gradf,x,p)
        !!
        !! Steepest Gradient Descent
        !! 
        !! gradf : function handle for gradient using
        !!         prodecure interface "gradf_intrfc"
        !! 
        !! x (REAL) : shape (:,:) row/column vector input
        !!            for gradient function
        !!
        !! p (REAL) : same shape as x, this is resulting
        !!            values from grad(x) prodedure and is
        !!            the direction we iterat towards solution
        !!
        PROCEDURE(gradf_intrfc) :: gradf
        REAL(rkp),  INTENT(IN ) :: x(:,:)
        REAL(rkp),  INTENT(OUT) :: p(:,:)
        
        p = -gradf(x(:,1))

    END SUBROUTINE steepDescent


    SUBROUTINE BFGS(gradf,first,x,xi,p,B,iB)
        !!
        !! Broyden-Fletcher-Goldfarb-Shanno method
        !! 
        !! gradf : function handle for gradient using
        !!         prodecure interface "gradf_intrfc"
        !! 
        !! first (LOGICAL) : True on first procedure call
        !!                   for initializing inverse of
        !!                   Hessian as identity matrix
        !!
        !! x (REAL) : shape (n,1) row/column vector input
        !!            for gradient function and the current
        !!            solution guess
        !!
        !! xi (REAL) : old guess of solution
        !!
        !! p (REAL) : same shape as x, this is resulting
        !!            values from grad(x) prodedure
        !!
        !! B (REAL) : shape (n,n) matrix that serves as 
        !!            an approximation to the Hessian of f
        !!
        !! iB (REAL) : same shape as B, inverse of B, it is
        !!             computed using doubly applied Sherman-
        !!             Morrison formula
        !!
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
        ! formula is applied twice to get solution we see
        sTu  = DOT_PRODUCT(s(:,1),u(:,1))
        isTu = 1.0/sTu
        t1   = sTu + INNERPROD(MATMUL(TRANSPOSE(u),iB),u)
        iB = iB - ( MATMUL( MATMUL(iB,u), TRANSPOSE(s) ) &
                  + MATMUL( MATMUL(s,TRANSPOSE(u)), iB ) ) * isTu

        iB = iB + MATMUL(s,TRANSPOSE(s)) * t1 * (isTu**2)

        ! Calculate solution descent direction
        p = - MATMUL(iB,gradf(x(:,1)))

    END SUBROUTINE BFGS


    SUBROUTINE backtrack(f,gradf,x,p,alpha0,rho,gamma,alpha)
        !!
        !! Backtracking method used as part of a line search method
        !! to solve multivariate unconstrained optimization problems.
        !!
        !! f : function handle using a multidimensional procedure interface
        !!
        !! gradf: Same as above, but for the gradient of the function
        !!
        !! x (REAL) : shape (n,1) technically, it is the initial guess of 
        !! the solution
        !!
        !! p (REAL) : shape(x), it is the descent direction towards solution
        !!
        !! alpha0 (REAL) : constant tunable parameter
        !!
        !! rho  (REAL) : constant tunable parameter
        !!
        !! gamma (REAL) : constant tunable parameter
        !!
        !! alpha (REAL) : the variable stepsize that we need to calculate
        !!                for advancing the line search solution at each
        !!                iteration
        !!
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
