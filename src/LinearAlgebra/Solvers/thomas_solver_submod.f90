!!*********************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/LinearAlgebra/thomas_solver_submod.f90
!!
!! Thomas' Algorithm for solving Ax=b linear systems where A is a tridiagonal
!! matrix.
!!
!! Implemented here is a version that can either save memory management by
!! overwriting the input arrays or allocating memory for temporary arrays if
!! user does not want to lose original array data. A simple logical flag will
!! toggle this on or off (overwrite_arrays = .true. or .false.)
!!
!! Structure of linear system being solved,
!!
!! ----                      ----  --    --     --    --
!! | b(1) c(1)             0    |  | x(1) |     | d(1) |
!! | a(2) b(2) c(2)             |  | x(2) |     | d(2) |
!! |      a(3) b(3) ....        |  | x(3) |  =  | d(3) |
!! |           .... .... c(n-1) |  | .... |     | .... |
!! | 0              a(n) b(n)   |  | x(n) |     | d(n) |
!! ----                      ----  --    --     --    --
!!
!! Where we have Ax=b form and "A" is a tridiagonal matrix. The vectors a,b,c are
!! derived from the subdiagonal, main diagonal, and superdiagonal of the matrix
!! "A", and the vector d is derived from the column vector "b". We are looking to
!! find a solution for the vector "x" in the linear system.
!!
!! For the SUBROUTINE implemented here, thomas_solver, we have the following
!! in/out arguments:
!!
!! INPUT:
!!     a - Subdiagonal vector from matrix "A"
!!     c - Superdiagonal vector from matrix "A"
!!
!! INPUT/OUTPUT:
!!     b - Main diagonal vector from matrix "A", can be optionally overwritten
!!     d - the right-hand side column vector "b", can be optionally overwritten
!!
!! OUTPUT:
!!     x - the column vector whose solution is determined by Thomas' Algorithm
!!
!! OPTIONAL INPUT:
!!     overwrite_arrays - logical flag that turns on/off overwriting the arrays of
!!                        b and d given as inputs to the subroutine. If .true. the 
!!                        arrays are overwritten to reduce memory usage. If .false.
!!                        temporary arrays are allocated to preserve the values
!!                        that are stored in arrays b and d
!!
!!*********************************************************************************
SUBMODULE (matrix_systems_solvers) thomas_solver_submod

CONTAINS

    MODULE SUBROUTINE thomas_solver (a,b,c,d,x,overwrite_arrays)
        REAL(rkp),         INTENT(IN   ) :: a(:)
        REAL(rkp),         INTENT(INOUT) :: b(:)
        REAL(rkp),         INTENT(IN   ) :: c(:)
        REAL(rkp),         INTENT(INOUT) :: d(:)
        REAL(rkp),         INTENT(OUT  ) :: x(:)
        LOGICAL, OPTIONAL, INTENT(IN   ) :: overwrite_arrays

        !! local variables
        INTEGER(ilkp)          :: i,n
        REAL(rkp)              :: w
        REAL(rkp), ALLOCATABLE :: ck(:),dk(:)
        LOGICAL                :: flag = .false.

        n = SIZE(x, DIM=1)

        IF ( PRESENT(overwrite_arrays) ) flag = overwrite_arrays

        IF ( flag ) THEN

            !... Forward Sweep
            DO i = 2,n
                
                w    = a(i) / b(i-1)
                b(i) = b(i) - w*c(i-1)
                d(i) = d(i) - w*d(i-1)

            ENDDO

            !... Back Substitution
            x(n) = d(n) / b(n)
            DO i = n-1, 1

                x(i) = (d(i) - c(i)*x(i+1)) / b(i)

            ENDDO

        ELSE !... not overwriting the input vectors b and d

            !... allocate temporary arrays for Thomas solver
            IF ( ALLOCATED(ck) ) DEALLOCATE(ck)
            IF ( ALLOCATED(dk) ) DEALLOCATE(dk)
            ALLOCATE(ck(n))
            ALLOCATE(dk(n))

            !... Forward Sweep
            ck(1) = c(i) / b(i)
            dk(1) = d(i) / b(i)
            DO i = 2,n

                ck(i) = c(i) / (b(i) - a(i)*ck(i-1))
                dk(i) = (d(i) - a(i)*dk(i-1)) / (b(i) - a(i)*ck(i-1))

            ENDDO

            !... Back Substitution
            x(n) = dk(n)
            DO i = n-1,1

                x(i) = dk(i) - ck(i)*x(i+1)

            ENDDO

            DEALLOCATE(ck)
            DEALLOCATE(dk)

        ENDIF

    END SUBROUTINE thomas_solver

END SUBMODULE thomas_solver_submod
