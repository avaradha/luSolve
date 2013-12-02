!-----------------------------------------------------------------------
! ASSIGNMENT - 6
!
! Aravindh Varadharaju
! Email ID: avaradha@indiana.edu
! Department of Computer Science
! Indiana University
!
! Perform LU factorization with partial pivoting on the matrix B 
! 
! Parameters:
!   B       - Matrix in column-major order
!   m       - Number of rows in B
!   n       - Number of columns in B
!   ldb     - Leading dimension of matrix B
!   ipiv    - Vector of length of m containing the Pivot indices
!   errflag - -k means argument number k is wrong
!             0 means success and everything is OK
!           - k means a near-zero divisor was found on row k
!   Note:
!       B can be a rectangular matrix where m =/= n
!
!   Had to refer to Professor Bramley's Matlab code LU4.m for the 
!   approach taken to swap pivot rows. Earlier, I was swapping the rows
!   with pivot value in one step. Once, the complete row was swapped 
!   DTRSV was then changing the values above the diagonal elements which
!   we did not want. 
!   Professor was actually swapping values within the column, doing the
!   transformations and at the end swapping values within the next column.
!
!-----------------------------------------------------------------------
SUBROUTINE LU4(B, ldb, m, n, ipiv, errflag)
        
    INTEGER :: ldb, m, n, errflag, incx, incy
    INTEGER :: ipiv(n)
    DOUBLE PRECISION :: B(ldb,n)
    DOUBLE PRECISION, PARAMETER :: ONE = 1.0d0, ZERO = 0.0d0
    DOUBLE PRECISION :: alpha = -1.0, SFAC
    DOUBLE PRECISION :: small = 1000*EPSILON(ONE)
    CHARACTER*1 UPLO, TRANS, DIAG
    INTEGER, DIMENSION(1) :: pivot

    INTRINSIC ABS       ! Computes the absolute value
    INTRINSIC MAXLOC    ! Returns the location of the maximum value of
                        ! elements in a specified dimension of an array

    ! BLAS 1 Routine 
    EXTERNAL DSWAP      ! Swap a vector with another, in our case, used
                        ! to swap the rows of the matrix
    EXTERNAL DSCAL      ! Used to scale the elements of the given 
                        ! vector
    
    ! BLAS 2 Routines
    EXTERNAL DGEMV      ! Matrix-Vector Product
    EXTERNAL DTRSV      ! Solves a system of linear equations whose
                        ! co-efficients are in a triangular matrix
    
    UPLO    = 'L'       ! Upper or Lower Triangular Matrix
    TRANS   = 'N'       ! Transposed?
    DIAG    = 'N'       ! Unit Diagonal         
    incx    = 1
    incy    = 1
    
    errflag = 0         ! No error to begin with

!   WRITE(*,*) ("LU4: A on entry..")
!   CALL PrintArrayD(0,0,m,n,B,ldb,'B')

    DO k = 1,n
        
        !Apply any previous interchanges to k-th column
        DO j=1,k-1
            IF (ipiv(j) /= j) THEN
                ! CALL DSWAP(n, x, incx, y, incy)
                CALL DSWAP(j,B(j,k),ldb,B(ipiv(j),k),ldb)
            END IF
        END DO  
                                
        ! Compute elements 1:k-1 of the k-th column
        
        ! CALL DTRSV(uplo, trans, diag, n, a, lda, x, incx)
        CALL DTRSV(UPLO,TRANS,DIAG,k-1,B,ldb,B(1,k),incx)

        ! Update the elements k:n of k-th column
        
        ! CALL DGEMV(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
        CALL DGEMV(TRANS,m-k+1,k-1,alpha,B(k,1),ldb,B(1,k),incx,ONE,B(k,k),incy)
                        
        ! If the diagonal element is small then set the error flag
        ! and return to the caller.
        IF (ABS(B(k,k)) < small) THEN
            errflag = k
            RETURN
        END IF

        ! Find the pivot element in the current column
        pivot = MAXLOC(ABS(B(k:m,k)))
        ipiv(k) = pivot(1) + k - 1
        
        ! If the Pivot element is small, then set the error flag
        ! and return to the caller.

        IF (ABS(pivot(1)) < small) THEN
            errflag = k
            RETURN
        END IF

        ! Swap the current row with the row which has the 
        ! pivot element. If the pivot is in the current row
        ! do not swap
        IF (ipiv(k) /= k) THEN
            ! CALL DSWAP(n, x, incx, y, incy)
            CALL DSWAP(k,B(k,1),ldb,B(ipiv(k),1),ldb)
        END IF

        ! Scale the elements 
        IF (k < n .OR. m > n) THEN
            SFAC=(1/B(k,k))
            ! CALL DSCAL(n, a, x, incx)
            CALL DSCAL(m-k,SFAC,B(k+1,k),1)
            !B(k+1:m,k) = B(k+1:m,k)/B(k,k)
        END IF
    END DO
!   WRITE(*,*) ("LU4: A on exit")
!   CALL PrintArrayD(0,0,m,n,B,ldb,'B')
END SUBROUTINE LU4
