!!********************************************************!!
!!
!! Ezra Brooker (2021)
!! 
!! NumLib/src/Optimization/Univariate/divided_differences.f90
!!
!! A form of Newton's Divided Difference 
!! method                                
!!                                       
!! Input Arguments:                      
!!     x - spatial points where function is evaluated                  
!!     y - function evaluation results   
!! (BOTH are length N arrays)            
!!                                       
!! Outpute Argument:                     
!!     ys - scalar divided difference    
!!********************************************************!!
MODULE divided_differences
use kindSettings
IMPLICIT NONE

CONTAINS

    FUNCTION divide_diff(x,y) RESULT (ys)
        REAL(rkp),  INTENT(IN) :: x(:),y(:)
        REAL(rkp)              :: ys
        REAL(rkp)              :: denom
        INTEGER(ilkp)          :: j,k

        ys = 0.0_rkp

        DO j = 1, SIZE(x,DIM=1)
            denom = x(j) - x(j+1)
            DO k = j+2, j+SIZE(x,DIM=1)-1
                denom = denom * (x(j)-x(k))
            END DO
            ys = ys + y(j)/denom
        END DO

    END FUNCTION divide_diff

END MODULE divided_differences