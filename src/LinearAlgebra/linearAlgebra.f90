MODULE linearAlgebra
USE kindSettings, ONLY : rkp
IMPLICIT NONE
PRIVATE

    PUBLIC :: outerprod, innerprod, shermanMorrison

    INTERFACE outerprod
        PROCEDURE :: outerprod_1d, outerprod_2d
    END INTERFACE outerprod

    INTERFACE innerprod
        PROCEDURE :: innerprod_1d, innerprod_2d
    END INTERFACE innerprod


CONTAINS

    FUNCTION outerprod_1d(vec1, vec2) RESULT(mat)
        REAL(rkp), INTENT(IN) :: vec1(:), vec2(:)
        REAL(rkp)             :: mat(SIZE(vec1),SIZE(vec2))
        INTEGER               :: i

        DO i = 1, SIZE(vec1)
            mat(i,:) = vec1(i) * vec2
        ENDDO 

    END FUNCTION outerprod_1d

    FUNCTION outerprod_2d(vec1, vec2) RESULT(mat)
        REAL(rkp), INTENT(IN) :: vec1(:,:), vec2(:,:)
        REAL(rkp)             :: mat(SIZE(vec1,DIM=1),SIZE(vec2,DIM=1))
        INTEGER               :: i

        mat = MATMUL(vec1,vec2)

    END FUNCTION outerprod_2d

    FUNCTION innerprod_1d(vec1, vec2) RESULT(scalar)
        REAL(rkp), INTENT(IN) :: vec1(:), vec2(:)
        REAL(rkp)             :: scalar
        
        scalar = DOT_PRODUCT(vec1,vec2)
    
    END FUNCTION innerprod_1d

    FUNCTION innerprod_2d(vec1, vec2) RESULT(scalar)
        REAL(rkp), INTENT(IN) :: vec1(:,:), vec2(:,:)
        REAL(rkp)             :: scalar
    
        scalar = DOT_PRODUCT(vec1(:,1),vec2(:,1))
    
    END FUNCTION innerprod_2d

    SUBROUTINE shermanMorrison(s,u,v,iB)
        REAL(rkp), INTENT(IN   ) :: u(:,:), v(:,:), s(:,:)
        REAL(rkp), INTENT(INOUT) :: iB(:,:)
        REAL(rkp)                :: sTu, isTu, t1

        sTu  = DOT_PRODUCT(s(:,1),u(:,1))
        isTu = 1.0/sTu
        t1   = sTu + INNERPROD(MATMUL(TRANSPOSE(u),iB),u)

        iB = iB - ( MATMUL( MATMUL(iB,u), TRANSPOSE(s) ) &
                  + MATMUL( MATMUL(s,TRANSPOSE(u)), iB ) ) * isTu

        iB = iB + MATMUL(s,TRANSPOSE(s)) * t1 * (isTu**2)

    END SUBROUTINE shermanMorrison

END MODULE linearAlgebra