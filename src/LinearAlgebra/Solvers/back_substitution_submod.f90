!!*********************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/LinearAlgebra/back_substitution_submod.f90
!!
!! Back Substitution algorithm for taking an upper right triangular matrix 
!! and the RHS column vector of Ux=y system of equations and solving for x
!!
!! This is often used in conjunction with Gaussian Elimination to get a
!! matrix A from Ax=y into an upper right triangular form before performing
!! back substitution
!!
!! INPUT:
!!     U - Shape (m,n) upper right triangular matrix, possibly from Gaussian
!!         Elimination procedure. U is also an augmented matrix, holding the 
!!         RHS of the Ux=y equation in the rightmost column of U.
!!     b - Shape(m) column vector holding the RHS of the linear system
!!
!! OPTIONAL OUTPUT:
!!     x - Shape (m) column vector solution to Ux=y that we compute with back
!!         substitution algorithm.
!!
!!*********************************************************************************
SUBMODULE (matrix_systems_solvers) back_substitution_submod

CONTAINS

    MODULE SUBROUTINE back_substitution (U,b,x)
        REAL(rkp), INTENT(IN ) :: U(:,:)
        REAL(rkp), INTENT(IN ) :: b(:)
        REAL(rkp), INTENT(OUT) :: x(:)

        INTEGER(ilkp) :: i,j,m,n
        REAL(rkp)     :: y

        m = SIZE(U,DIM=1)
        n = SIZE(U,DIM=2)
        x = 0.0_rkp

        DO i = m, 1, -1

            x(i) = b(i) 

            DO j = i+1,m
                x(i) = x(i) - U(i,j)*x(j)
            ENDDO

            x(i) = x(i) / U(i,i)

        ENDDO

    END SUBROUTINE back_substitution

END SUBMODULE back_substitution_submod