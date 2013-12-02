SUBROUTINE LU8(A, lda, n, ipiv, errflag, nb)
!------------------------------------------------------------------------------------
!
! Perform LU factorization with partial (row) pivoting on the matrix A contained in
! the array A. A is n x n with a declared leading dimension of lda. This version is
! based on matrix-matrix products, and calls BLAS-3 (BLAS level 3) routines. The
! variable errflag is returned by LU4, and may need to be adjusted to global array
! values before it is returned to the caller.
!
! Variables are:
!
! n     The number of rows = columns of A
! A     On entry, the m by n *matrix* to be factored is stored in the *array* A.
!       On exit, the LU factors are stored in the array. L is unit lower 
!       triangular, so its diagonal entries are not stored explicitly
! lda   Declared leading dimension of the array A
! ipiv  The pivot indices; for 1 <= i <= min(m,n), row i of the matrix was 
!       interchanged with row ipiv(i)
! errflag  = 0: success
!          < 0: if errflag = -k, the k-th argument had an illegal value
!          > 0: if errflag =  k, U(k,k) is too small to rely upon.
! errflag can be passed to LU4, but its value may need to be adjusted to global
! indexing before returning it to the program calling LU8.
!
! Credit:
!   The idea of LU8 is kind of based on Professor Bramley's LU8.m
!   Had to run Professor's code in Matlab to understand the logic
!   behind the code... 
!
!------------------------------------------------------------------------------------
    
    INTEGER, INTENT(IN):: lda, n, nb
    INTEGER, INTENT(OUT):: errflag
    INTEGER, INTENT(OUT):: ipiv(n)
    INTEGER :: i = 1, j = 1, k = 1, cols, nob, l, tmp
    DOUBLE PRECISION, INTENT(INOUT):: A(lda,n)
    CHARACTER*1 UPLO, TRANS, DIAG

    DOUBLE PRECISION, PARAMETER:: ONE   = 1.0d0, ZERO = 0.0d0
   
    DOUBLE PRECISION, PARAMETER:: small = 1000*EPSILON(ONE)
    
    DOUBLE PRECISION :: alpha = -1.0, beta = 1.0

    INTRINSIC MIN       ! Returns the minimum value of the elements

    ! BLAS 1 Routine 
    EXTERNAL DSWAP      ! Swap a vector with another, in our case, used
                        ! to swap the rows of the matrix
    
    ! BLAS 3 Routines
    EXTERNAL DGEMM      ! Computes a scalar-matrix-matrix product and
                        ! adds the result to a scalar-matrix product
    EXTERNAL DTRSM      ! Solves a Matrix equation (one matrix operand 
                        ! is triangular
    
    UPLO    = 'L'       ! Upper or Lower Triangular Matrix
    TRANS   = 'N'       ! Transposed?
    DIAG    = 'N'       ! Unit Diagonal         
    incx    = 1
    incy    = 1

    errflag = 0         ! No error to begin with
    
    ! Without this loop ipiv was containing garbage
    DO i = 1,n
        ipiv(i) = i
    END DO
 
    ! My 4 X 4 sample array to start with block size of 2
    !
    ! 2 -1  5  1
    ! 3  2  2 -6
    ! 1  3  3 -1 
    ! 5  2 -3  3
    
!   WRITE (*,*) ("LU8: Initialized array ....")
!   CALL PrintArrayD(0,0,n,n,A,lda,'A')
 
    IF (nb>n) THEN
        ! If the block size is greater than n call LU4 directly
        CALL LU4(A, lda, n, n, ipiv, errflag)
    ELSE
        DO k = 1,n,nb
            nob = MIN(n-k+1,nb) ! blocks
            noc = MIN(n,k+nob-1) ! columns
            
            ! CALL DGEMM(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            CALL DGEMM(TRANS,TRANS,n-k+1,nob,k-1,alpha,A(k,1),lda,A(1,k),lda,beta,A(k,k),lda)
            
            !LU factorize the current block 
            ! CALL LU4(B, ldb, m, n, ipiv, errflag)
            CALL LU4(A(k,k),lda,n-k+1,nob,ipiv(k),errflag)
        
            !Adjust the pivot values
            DO i = k,noc
                ipiv(i) = ipiv(i)+k-1
            END DO

            !  5  1                          -3  3
            !  2 -6  needs to be changed to   3 -1 
            !  3 -1                           2  6
            ! -3  3                           5  1
            
            DO l = k,k+nob-1
                IF (tmp /= l) THEN
                    ! CALL DSWAP(n, x, incx, y, incy)
                    CALL DSWAP(k-1, A(l,1),lda,A(tmp,1),lda)
                END IF
            END DO

            DO l = k,k+nob-1
                tmp = ipiv(l)
                IF (tmp /= l) THEN
                    ! CALL DSWAP(n, x, incx, y, incy)
                    CALL DSWAP(n-k-nob+1, A(l,l+nob),lda,A(tmp,l+nob),lda)
                END IF
            END DO
                        
            IF (k+nob <= n) THEN
                ! CALL DGEMM(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
                CALL DGEMM(TRANS,TRANS,nob,n-k-nob+1,k-1,alpha,A(k,1),lda,a(1,j+nob),lda,beta,A(k,k+nob),lda)
                DIAG='U'
                ! CALL DTRSM (side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
                CALL DTRSM('L',UPLO,TRANS,DIAG,nob,n-k-nob,alpha,A(k,k),lda,A(k,k+nob),lda)
            END IF
        END DO
    END IF
!   WRITE (*,*) ("LU8: Output array A ....")
!   CALL PrintArrayD(0,0,n,n,A,lda,'A')
    RETURN
END SUBROUTINE LU8
