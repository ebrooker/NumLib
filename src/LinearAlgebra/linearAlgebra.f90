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

END MODULE linearAlgebra