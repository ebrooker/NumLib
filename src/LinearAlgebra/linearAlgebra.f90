!!********************************************************************
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/LinearAlgebra/linearAlgebra.f90
!!
!! Module currently holds a few wrappers for vector and matrix
!! multiplication methods. Included are the following:
!!
!!     - outerprod: An interface wrapper for two methods that can 
!!                  calculate outer products:
!!                    -- outerprod_1d: two vectors of shape (n)
!!                    -- outerprod_2d: two matrices of shape (n,m)
!!                                      including row/column vectors 
!!
!!     - innerprod: An interface wrapper for two methods that can 
!!                  calculate inner products:
!!                    -- innerprod_1d: two vectors of shape (n)
!!                    -- innerprod_2d: two matrices of shape (n,m)
!!                                     including row/column vectors 
!!
!!
!!********************************************************************
MODULE linearAlgebra
USE kindSettings, ONLY : rkp
IMPLICIT NONE
PRIVATE

    PUBLIC :: outerprod, innerprod


    !! outer product interface wrapper
    INTERFACE outerprod
        PROCEDURE :: outerprod_1d, outerprod_2d
    END INTERFACE outerprod

    !! inner product interface wrapper
    INTERFACE innerprod
        PROCEDURE :: innerprod_1d, innerprod_2d
    END INTERFACE innerprod


CONTAINS

    FUNCTION outerprod_1d(vec1, vec2) RESULT(mat)
        ! Vector 1d array version of outer product multiply
        ! Results in a matrix with shape (n1,n2) from the
        ! vectors with shape vec1(n1) and vec2(n2)
        REAL(rkp), INTENT(IN) :: vec1(:), vec2(:)
        REAL(rkp)             :: mat(SIZE(vec1),SIZE(vec2))
        INTEGER               :: i

        DO i = 1, SIZE(vec1)
            mat(i,:) = vec1(i) * vec2
        ENDDO 

    END FUNCTION outerprod_1d

    FUNCTION outerprod_2d(vec1, vec2) RESULT(mat)
        ! "Vector" 2d array version of outer product multiply
        ! Results in a matrix with shape (n1,n2) from the
        ! vectors with shape vec1(n1,m1>=1) and vec2(n2,m2>=1)
        REAL(rkp), INTENT(IN) :: vec1(:,:), vec2(:,:)
        REAL(rkp)             :: mat(SIZE(vec1,DIM=1),SIZE(vec2,DIM=1))
        INTEGER               :: i

        mat = MATMUL(vec1,vec2)

    END FUNCTION outerprod_2d

    FUNCTION innerprod_1d(vec1, vec2) RESULT(scalar)
        ! "Vector" 1d array version of inner product multiply
        ! Results in a scalar and uses DOT_PRODUCT intrinic
        ! Procedure exists for interchanging vectors of shape
        ! (n) and (n,1) without worrying about proper indexing
        REAL(rkp), INTENT(IN) :: vec1(:), vec2(:)
        REAL(rkp)             :: scalar
        
        scalar = DOT_PRODUCT(vec1,vec2)
    
    END FUNCTION innerprod_1d

    FUNCTION innerprod_2d(vec1, vec2) RESULT(scalar)
        ! "Vector" 2d array version of inner product multiply
        ! Results in a scalar and uses DOT_PRODUCT intrinic
        ! Procedure exists for interchanging vectors of shape
        ! (n) and (n,1) without worrying about proper indexing
        REAL(rkp), INTENT(IN) :: vec1(:,:), vec2(:,:)
        REAL(rkp)             :: scalar
    
        scalar = DOT_PRODUCT(vec1(:,1),vec2(:,1))
    
    END FUNCTION innerprod_2d

END MODULE linearAlgebra