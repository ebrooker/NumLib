!!******************************************************************************
!!
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/LinearAlgebra/lu_decomposition_submod.f90
!!
!!
!!******************************************************************************
SUBMODULE (matrix_transforms) lu_decomposition_submod


CONTAINS

    MODULE SUBROUTINE lu_decomposition (A,L,U,P,Q,D,S)
        REAL(rkp),                  INTENT(IN   ) :: A(:,:)
        REAL(rkp),                  INTENT(OUT  ) :: L(:,:)
        REAL(rkp),                  INTENT(OUT  ) :: U(:,:)
        REAL(rkp),        OPTIONAL, INTENT(INOUT) :: P(:,:)
        REAL(rkp),        OPTIONAL, INTENT(INOUT) :: Q(:,:)
        REAL(rkp),        OPTIONAL, INTENT(INOUT) :: D(:,:)
        REAL(rkp),        OPTIONAL, INTENT(INOUT) :: S(:,:)

        !... local variables
        INTEGER(ilkp)                 :: i,j,k,n,imax,qmax(2),irow,icol
        REAL(rkp)                     :: total
        CHARACTER(LEN=:), ALLOCATABLE :: method

        !... Get size and check for matrix squareness first
        n = SIZE(A,DIM=1)
        IF ( n .ne. SIZE(A,DIM=2) ) STOP &
        "[LU Decomposition ERROR] Matrix A must be a square matrix."

        IF ( PRESENT(P) ) P(:,:) = 0.0_rkp
        IF ( PRESENT(Q) ) Q(:,:) = 0.0_rkp
        IF ( PRESENT(D) ) D(:,:) = 0.0_rkp
        IF ( PRESENT(S) ) S(:,:) = 0.0_rkp

        method = 'LU'
        IF ( PRESENT(P) ) method = 'LUP'
        IF ( PRESENT(Q) ) method = 'LUPQ'

        U(:,:) = A(:,:)
        L(:,:) = 0.0e0_rkp
        DO i = 1,n
            IF ( PRESENT(P) ) P(i,i) = 1.0_rkp
            IF ( PRESENT(Q) ) Q(i,i) = 1.0_rkp
        ENDDO

        IF ( method .eq. 'LUP' ) THEN ! LUP decomp

            DO k = 1,n-1

                !... Find index of row with maxval for column i
                imax = MAXLOC(ABS(U(k:n,k)), DIM=1) + k - 1 

                !... Check principal diag element, abort if zero
                IF ( U(k,imax) .eq. 0.0e0_rkp ) THEN
                    WRITE(*,'(a,2i3, es12.4)') "A(i,j) = ", i, j, U(k,imax)
                    STOP "[Gaussian Elimination ERROR] A is singular"
                ENDIF

                !... perform the row swaps in L, U, and P
                IF ( imax .ne. k ) THEN

                    CALL swap_slices(U,imax,k,mode='row')
                    CALL swap_slices(P,imax,k,mode='row')

                ENDIF

                !... update U and L
                CALL lu(U,k)

            ENDDO !... end k-loop

        ELSEIF ( method .eq. 'LUPQ' ) THEN !... LUPQ decomp

            DO k = 1,n-1

                !... Find index row and column index of max(abs(A))
                qmax = MAXLOC(ABS(U(k:n,k:n)))
                irow = qmax(1) + k - 1
                icol = qmax(2) + k - 1

                !... Check principal diag element, abort if zero
                IF ( U(irow,icol) .eq. 0.0e0_rkp ) THEN
                    WRITE(*,'(a,2i3, es12.4)') "A(i,j) = ", i, j, U(irow,icol)
                    STOP "[Gaussian Elimination ERROR] A is singular"
                ENDIF

                !... perform the row and column swaps in L, U, and P
                IF ( irow .ne. k ) THEN
                    
                    !... Swap rows first with irow and k
                    CALL swap_slices(U,irow,k,mode='row')
                    CALL swap_slices(P,irow,k,mode='row')

                ENDIF

                IF ( icol .ne. k ) THEN

                    !... Swap columns next with icol and k
                    CALL swap_slices(U,icol,k,mode='column')
                    CALL swap_slices(Q,icol,k,mode='column')

                ENDIF

                !... Update U and L
                CALL lu(U,k)

            ENDDO !... end k-loop

        ELSE !... basic LU decomp

            DO k = 1, n

                !... Update U and L naively
                CALL lu(U,k)

            ENDDO !... end k-loop

        ENDIF

        !... Extract L from U
        DO i = 1,n
            L(i,i)     = U(i,i) / U(i,i)
            L(i,1:i-1) = U(i,1:i-1)
            U(i,1:i-1) = 0.0_rkp
        ENDDO

        !... Get diagonal matrix for A=LDU and adjust U
        IF ( PRESENT(D) ) THEN
            DO i = 1,n
                D(i,i) = U(i,i)
                U(i,:) = U(i,:) / D(i,i)
            ENDDO
        ENDIF

    END SUBROUTINE lu_decomposition


    SUBROUTINE swap_slices(mat,i,k,mode)
        !! Helper subroutine that swaps two 1D slices given given slice index i and k
        !! Can swap either rows or columns based on mode setting ("row" or "column")
        REAL(rkp),                  INTENT(INOUT) :: mat(:,:)
        INTEGER(ilkp),              INTENT(IN   ) :: i,k
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN   ) :: mode

        REAL(rkp), ALLOCATABLE :: temp(:)

        IF ( ALLOCATED(temp) ) DEALLOCATE(temp)
        ALLOCATE(temp(SIZE(mat,DIM=1)))

        IF ( mode .eq. 'row' ) THEN
            temp(:)  = mat(i,:)
            mat(i,:) = mat(k,:)
            mat(k,:) = temp(:)

        ELSEIF ( mode .eq. 'column' ) THEN
            temp(:)  = mat(:,i)
            mat(:,i) = mat(:,k)
            mat(:,k) = temp(:)

        ELSE
            STOP "[Swap Array Slices ERROR] Unrecognized swap mode"
        
        ENDIF

        IF ( ALLOCATED(temp) ) DEALLOCATE(temp)

    END SUBROUTINE swap_slices

    SUBROUTINE lu (U,k)
        !! Helper subroutine that does the actual computations
        REAL(rkp),     INTENT(INOUT) :: U(:,:)
        INTEGER(ilkp), INTENT(IN   ) :: k
        INTEGER(ilkp)                :: i,j

        DO i = k+1,SIZE(U,DIM=1)

            !... modify subdominant column elements in L
            U(i,k) = U(i,k) / U(k,k)

            !... modify superdominat column elements in U
            DO j = k+1,SIZE(U,DIM=2)
                U(i,j) = U(i,j) - U(i,k)*U(k,j)
            ENDDO

        ENDDO

    END SUBROUTINE lu


END SUBMODULE lu_decomposition_submod

