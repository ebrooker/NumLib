!!*****************************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/LinearAlgebra/matrix_systems_solvers.f90
!!
!! A collating module for a number of solvers for Ax=b and similar types of
!! systems of equations that are solved using different matrix algorithms.
!!
!! The interfaces are defined here and the actual implementations are realized
!! in submodule files for each method. Some examples of methods are the Thomas
!! algorithm for solving Ax=b where A is a tridiagonal matrix
!!
!!*****************************************************************************
MODULE matrix_systems_solvers
USE kindSettings
IMPLICIT NONE
PRIVATE

    PUBLIC :: thomas_solver
    PUBLIC :: back_substitution

    INTERFACE
    
        MODULE SUBROUTINE thomas_solver (a,b,c,d,x,overwrite_arrays)
            REAL(rkp),         INTENT(IN   ) :: a(:)
            REAL(rkp),         INTENT(INOUT) :: b(:)
            REAL(rkp),         INTENT(IN   ) :: c(:)
            REAL(rkp),         INTENT(INOUT) :: d(:)
            REAL(rkp),         INTENT(OUT  ) :: x(:)
            LOGICAL, OPTIONAL, INTENT(IN   ) :: overwrite_arrays
        END SUBROUTINE thomas_solver

        MODULE SUBROUTINE back_substitution (U,b,x)
            REAL(rkp), INTENT(IN ) :: U(:,:)
            REAL(rkp), INTENT(IN ) :: b(:)
            REAL(rkp), INTENT(OUT) :: x(:)
        END SUBROUTINE back_substitution

    END INTERFACE

END MODULE matrix_systems_solvers