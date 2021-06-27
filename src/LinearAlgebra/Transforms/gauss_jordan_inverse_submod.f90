!!******************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/LinearAlgebra/gauss_jordan_inverse_submod.f90
!!
!! The Gauss-Jordan Inverse is an extension of the Gaussian Elimination
!! transformation used to determine the inverse of a matrix, if it exists.
!!
!! Given an NxN matrix A, one extends the matrix to an Nx2N matrix of the form
!! [A|I] where I is the identity matrix. Then one proceeds with the Gaussian
!! elimination procedure. The goal is to reduce the A side of the new block
!! matrix to the identity matrix, which will lead to the I side of the block
!! matrix transforming into a new matrix B which is the inverse of A.
!!
!! [A|I] --> Gaussian Elimination --> [I|B], B = A^(-1)
!!
!! If the A side of [A|I] cannot be reduced to I, then A is not invertible
!!
!! INPUT:
!!     A - Shape (n,n) matrix to be transformed using Gauss-Jordan Elimination
!!     B - Shape (n,n) matrix that we store the inverted form of A in from the 
!!         transformed block matrix [A|I] --> [I|B] where B = A^(-1)
!!
!! OPTIONAL INPUT:
!!    tol - A numerical tolerance for determining if the mean of the trace of
!!          the left side block of the transformed [A|I] block matrix has been
!!          reduced to the Identity Matrix. A numerical tolerance is used since
!!          this is a floating point comparison and very small roundoff will
!!          likely occur. This parameter is optional if the user would like to
!!          increase or decrease the tolerance variable "eps" that is used.
!!
!!******************************************************************************
SUBMODULE (matrix_transforms) gauss_jordan_inverse_submod

CONTAINS

    MODULE SUBROUTINE gauss_jordan_inverse (A,B,tol)
        REAL(rkp),           INTENT(IN ) :: A(:,:)
        REAL(rkp),           INTENT(OUT) :: B(:,:)
        REAL(rkp), OPTIONAL, INTENT(IN ) :: tol

        INTEGER(ilkp)          :: i,j,k,n,n2
        REAL(rkp)              :: x, eps=1.0e-6
        REAL(rkp), ALLOCATABLE :: AI(:,:)

        n = SIZE(A,DIM=1)
        IF ( n .ne. SIZE(A,DIM=2) ) STOP &
        "[Gauss-Jordan Elimination ERROR] A is not a square matrix"

        n2 = 2*n

        !... adjust tolerance if provided
        IF ( PRESENT(tol) ) eps = tol

        !... allocate temporary arrays
        IF ( ALLOCATED(AI ) ) DEALLOCATE(AI )

        ALLOCATE( AI(n,n2) )

        !... extend A to [A|I] augmented block matrix
        AI = 0.0_rkp
        DO i = 1,n
            AI(i,1:n) = A(i,:)
            AI(i,n+i) = 1.0_rkp
        ENDDO

        !... Use Gaussian Elimination in RREF mode
        CALL gaussian_elimination(AI,form='rref')

        !... Get trace of left side block to check if left block is ~I
        x = 0.0_rkp
        DO i = 1,n
            x = x + AI(i,i)
        ENDDO

        !... Mean of trace should be close to 1.0 within some specified tolerance
        x = x/REAL(n,kind=rkp)
        IF (1.0_rkp - eps .gt. x .or. x .gt. 1.0_rkp + eps) stop &
        "[Gauss-Jordon Elimination ERROR] Matrix A is not invertible"

        !... Take right side block from transformed [A|I] --> [I|B] to get B
        B = 0.0_rkp
        DO i = 1,n
            B(i,:) = AI(i,n+1:n2)
        ENDDO

        !... FREE YOUR MEMORY
        IF ( ALLOCATED(AI ) ) DEALLOCATE(AI )

    END SUBROUTINE gauss_jordan_inverse

END SUBMODULE gauss_jordan_inverse_submod
