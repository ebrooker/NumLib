!!******************************************************************************
!!
!! Author: Ezra Brooker
!! 2021
!!
!! NumLib/src/LinearAlgebra/gaussian_elimination_submod.f90
!!
!! Performs Gaussian Elimination on a given matrix of shape (m,n) using partial
!! pivoting to obtain an upper triangular matrix U. This new matrix can also be
!! put into either Row Echelon Form or Reduced Row Echelon Form with an optional
!! argument character string flag.
!!
!! The procedure will overwrite the input matrix unless user gives an optional
!! second matrix, then the solution will be stored in that second argument, and
!! thus preserving the original input matrix.
!!
!! The use of partial pivoting, i.e. selecting the row with the maximum
!! absolute value in the kth-column, is implemented for maximal stability.
!! Normal pivoting with the maximum non-absolute value is less stable, and
!! no pivoting, i.e. no row swapping, is the least stable version.
!!
!! INPUT/OUTPUT:
!!     A - Shape (m,n) matrix to be transformed using Gaussian Elimination; is
!!         overwritten unless a separate solution matrix variable is provided.
!!         This can either be a typical square matrix or the augmented matrix
!!         [A|b] from Ax=b
!!
!! OPTIONAL INPUT:
!!     form - Character string flag that can either be form="ref" or form="rref"
!!            and is used to set flags for what row operations are performed 
!!            during the elimination procedure.
!!
!! OPTIONAL INPUT/OUTPUT:
!!     b - Shape(m) column vector to be transformed into REF, similarly to A.
!!         It is an optional parameter and will be overwritten unless another
!!         1D array is passed as the saving array (array "y")
!!
!! OPTIONAL OUTPUT:
!!     U - Shape (m,n) matrix that will house the solution to the Gaussian
!!         Elimination performed on A. If this OPTIONAL argument is provided,
!!         then the data of matrix A will be preserved after procedure. This
!!         must be same shape as A.
!!     y - Shape(m) column vector that will hold the transformed state of b in
!!         order to preserve data in b.
!!
!!******************************************************************************
SUBMODULE (matrix_transforms) gaussian_elimination_submod

CONTAINS

    MODULE SUBROUTINE gaussian_elimination (A,U,b,y,form)
        REAL(rkp),                  INTENT(INOUT) :: A(:,:)
        REAL(rkp),        OPTIONAL, INTENT(OUT  ) :: U(:,:)
        REAL(rkp),        OPTIONAL, INTENT(INOUT) :: b(:)
        REAL(rkp),        OPTIONAL, INTENT(OUT  ) :: y(:)
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN   ) :: form

        !! Local variables
        INTEGER(ilkp)          :: i,j,k,m,n,imax,irow,flag
        REAL(rkp)              :: x
        REAL(rkp), ALLOCATABLE :: row(:)
        LOGICAL                :: use_U, use_b, use_y

        !... Get logical flags for optional arguments
        use_U = PRESENT(U)
        use_b = PRESENT(b)
        use_y = PRESENT(y)

        !... Determine if base Gaussian form, REF, or RREF
        flag = 0
        IF ( PRESENT(form) ) THEN
            IF ( form .eq. 'ref'  ) flag = 1
            IF ( form .eq. 'rref' ) flag = 2
        ENDIF

        m = SIZE(A,DIM=1)
        n = SIZE(A,DIM=2)

        IF ( ALLOCATED(row) ) DEALLOCATE(row)
        ALLOCATE(row(n))

        !... Save data into y if being used based on A or b
        IF ( use_y .and. use_b ) THEN
            y(:) = b(:)
        ELSEIF ( use_y .and. .not. use_b ) THEN
            y(:) = A(:,n)
        ENDIF


        IF ( use_U ) THEN !... Preserve input matrix A

            !... Copy matrix A to matrix U
            U(:,:) = A(:,:)

            !... Begin Gaussian elimination with partial pivoting
            DO k = 1,m

                !... Find index of row with maxval for column i
                imax = MAXLOC(ABS(U(:,k)), DIM=1)

                !... Check principal diag element, abort if zero
                IF ( U(k,imax) .eq. 0.0_rkp ) THEN
                    WRITE(*,'(a,2i3, es12.4)') "A(i,j) = ", i, j, U(k,imax)
                    STOP "[Gaussian Elimination ERROR] A is singular"
                ENDIF

                !... Perform row swap between imax and k
                IF ( imax .ne. k ) THEN
                    row(:)    = U(k,:)
                    U(k,:)    = U(imax,:)
                    U(imax,:) = row(:)
                ENDIF

                !... Normalize if REF or RREF
                IF ( flag .gt. 0 ) THEN
                    U(k,:) = U(k,:) / U(k,k)
                    IF ( use_b .and. .not. use_y ) b(k) = b(k) / U(k,k)
                    IF ( use_b .and.       use_y ) y(k) = y(k) / U(k,k)
                ENDIF

                !... Set row loop range for RREF or not
                irow = k + 1
                IF ( flag .eq. 2 ) irow = 1

                DO i = irow,m

                    IF ( i .eq. k ) CYCLE

                    x = U(i,k) / U(k,k)

                    !... If using some combo of 1D arrays for RHS                    
                    IF ( use_b .and. .not. use_y ) b(i) = b(i) - x*b(k)
                    IF ( use_b .and.       use_y ) y(i) = y(i) - x*y(k)

                    !... subtract multiple from corresponding row elements
                    DO j = k+1,n
                            U(i,j) = U(i,j) - x*U(k,j)
                    ENDDO
                    U(i,k) = 0.0_rkp

                ENDDO !... end i-loop

            ENDDO !... end k-loop

        ELSE !... Overwrite input matrix A            

            !... Begin Gaussian elimination with  partial pivoting
            DO k = 1,m

                !... Find index of row with maxval for column i
                imax = MAXLOC(ABS(A(:,k)), DIM=1)

                !... Check principal diag element, abort if zero
                IF ( A(k,imax) .eq. 0.0_rkp ) THEN
                    WRITE(*,'(a,2i3, es12.4)') "A(i,j) = ", i, j, A(k,imax)
                    STOP "[Gaussian Elimination ERROR] A is singular"
                ENDIF

                !... Perform row swap between imax and k
                IF ( imax .ne. k ) THEN
                    row(:)    = A(k,:)
                    A(k,:)    = A(imax,:)
                    A(imax,:) = row(:)
                ENDIF

                !... Normalize if REF or RREF
                IF ( flag .gt. 0 ) THEN
                    A(k,:) = A(k,:) / A(k,k)
                    IF ( use_b .and. .not. use_y ) b(k) = b(k) / A(k,k)
                    IF ( use_b .and.       use_y ) y(k) = y(k) / A(k,k)
                ENDIF

                !... Set row loop range for RREF or not
                irow = k + 1
                IF ( flag .eq. 2 ) irow = 1

                DO i = irow,m

                    IF ( i .eq. k ) CYCLE

                    x = A(i,k) / A(k,k)

                    !... If using some combo of 1D arrays for RHS
                    IF ( use_b .and. .not. use_y ) b(i) = b(i) - x*b(k)
                    IF ( use_b .and.       use_y ) y(i) = y(i) - x*y(k)

                    !... subtract multiple from corresponding row elements
                    DO j = k+1,n
                            A(i,j) = A(i,j) - x*A(k,j)
                    ENDDO
                    A(i,k) = 0.0_rkp

                ENDDO !... end i-loop

            ENDDO !... end k-loop

        ENDIF

        IF ( ALLOCATED(row) ) DEALLOCATE(row)

        !... If using y array and not b array then decide which matrix to obtain RHS
        IF ( use_y .and. .not. use_b .and. use_U ) THEN
            y(:) = U(:,n)
        ELSEIF( use_y .and. .not. use_b .and. .not. use_U ) THEN
            y(:) = A(:,n)
        ENDIF

    END SUBROUTINE gaussian_elimination

END SUBMODULE gaussian_elimination_submod
