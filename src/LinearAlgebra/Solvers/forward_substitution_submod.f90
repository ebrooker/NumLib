!!*********************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/LinearAlgebra/forward_substitution_submod.f90
!!
!! Back Substitution algorithm for taking an upper right triangular matrix 
!! and the RHS column vector of Ly=b system of equations and solving for y
!! from the primary LUx=b equation derived from a factored form of Ax=b
!!
!! This is often used in conjunction with LU type decompositions to get a
!! matrix A from Ax=b into a lower right triangular form before performing
!! forward substitution
!!
!! INPUT:
!!     L - Shape (m,n) upper right triangular matrix, possibly from Gaussian
!!         Elimination procedure. L is also an augmented matrix, holding the 
!!         RHS of the Ly=b equation in the rightmost column of L.
!!     b - Shape(m) column vector holding the RHS of the linear system
!!
!! OPTIONAL OUTPUT:
!!     y - Shape (m) column vector solution to Ly=b that we compute with forward
!!         substitution algorithm.
!!
!!*********************************************************************************
SUBMODULE (matrix_systems_solvers) forward_substitution_submod

CONTAINS

    MODULE SUBROUTINE forward_substitution (L,b,y)
        REAL(rkp), INTENT(IN ) :: L(:,:)
        REAL(rkp), INTENT(IN ) :: b(:)
        REAL(rkp), INTENT(OUT) :: y(:)

        INTEGER(ilkp) :: i,j,m,n
        REAL(rkp)     :: y

        m = SIZE(L,DIM=1)
        n = SIZE(L,DIM=2)
        y = 0.0_rkp

        y(1) = b(1) / L(1,1)

        DO i = 2,m

            y(i) = b(i) 

            !... compute dot product
            DO j = 1,i-1
                y(i) = y(i) - L(i,j)*y(j)
            ENDDO

            y(i) = y(i) / L(i,i)

        ENDDO

    END SUBROUTINE forward_substitution

END SUBMODULE forward_substitution_submod