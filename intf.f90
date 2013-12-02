
subroutine LU4(B, ldb, m, n, ipiv, errflag)
!------------------------------------------------------------------------------------
!
! Perform LU factorization with partial (row) pivoting on the matrix B contained in
! the array B. B is m x n with a declared leading dimension of ldb. This version is
! based on matrix-vector products, and call BLAS-2 (BLAS level 2) routines. The
! parameter "small" determines when partial pivoting has failed from the absolute
! value of the diagonal element B(k,k) being too small to safely divide by. When
! that happens errflag is set to k and LU4 bails out. The matrix B can be
! rectangular (viz., m /= n)
!
!------------------------------------------------------------------------------------
    integer:: ldb, m, n, errflag
    integer:: ipiv(n)
    double precision:: B(ldb,n)
    double precision, parameter:: ONE = 1.0d0, ZERO = 0.0d0
    double precision:: small = 1000*epsilon(ONE)
end subroutine LU4

!=============================================================================

subroutine LU8(A, lda, n, ipiv, errflag, nb)
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
!------------------------------------------------------------------------------------
    integer, intent(in):: lda, n, nb
    integer, intent(out):: errflag
    integer, intent(out):: ipiv(n)
    double precision, intent(inout):: A(lda,n)

    double precision, parameter:: ONE   = 1.0d0, ZERO = 0.0d0
    double precision, parameter:: small = 1000*epsilon(ONE)

end subroutine LU8

