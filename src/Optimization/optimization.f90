!!********************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/Optimization/optimization.f90
!!
!! Module used to collate all of the optimzation routines and subunits
!! into one namespace for easier porting as a library code.
!!
!!********************************************************************
MODULE optimization
    USE lineSearch
    IMPLICIT NONE
    !!
    !! Unconstrained GlobalMinmization Methods:
    !!      - lineSearch: Global Minimization routine w/ backtracking
    !!          -- steepestDescent (multidimensional optimizer)
    !!          -- BFGS: Broyden-Fletcher-Goldfarb-Shanno algorithm
    !!                   use Hessian approximations
    !!
    !!  
    !!
    !!
END MODULE optimization