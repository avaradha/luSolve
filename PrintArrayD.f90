!-----------------------------------------------------------------------

! Credit: The PrintArrayD subroutine taken from examples in Intel's MKL
!		  libraries. A cool set of lines to print a matrix in Upper or 
!		  lower triangular format or just the complete matrix

!-----------------------------------------------------------------------
		SUBROUTINE PrintArrayD(flag1, flag2, m, n, a, lda, name)
				INTEGER flag1, flag2, m, n, lda
				CHARACTER*1 name
				DOUBLE PRECISION a(lda,*)
				INTEGER i, j
			 
				IF (flag1.eq.0) THEN
					PRINT 100, name, name, lda
				ELSE
					PRINT 101, name
				END IF
			  
				IF (flag2.eq.0) THEN
					DO i=1,m
						PRINT 110, (a(i,j),j=1,n)
					END DO
				ELSE IF (flag2.eq.1) THEN
					DO i=1, m
						GOTO (10,20,30,40,50) i
		  10       		PRINT 110, (a(i,j),j=i,m)
						GOTO  60
		  20       		PRINT 120, (a(i,j),j=i,m)
						GOTO 60
		  30       		PRINT 130, (a(i,j),j=i,m)
						GOTO 60
		  40       		PRINT 140, (a(i,j),j=i,m)
						GOTO 60
		  50       		PRINT 150, (a(i,j),j=i,m)
		  60       		CONTINUE
					END DO
				ELSE IF (flag2.eq.-1) THEN
					DO i=1, m
					 PRINT 110, (a(i,j),j=1,i)
					END DO
				END IF

		 100  	FORMAT(7x,'ARRAY ',a1,'   LD',a1,'=',i1)
		 101	FORMAT(7x,'ARRAY ',a1)
		 110  	FORMAT(9x,10(f8.3,2x))
		 120  	FORMAT(19x,10(f8.3,2x))
		 130  	FORMAT(29x,10(f8.3,2x))
		 140  	FORMAT(39x,10(f8.3,2x))
		 150	FORMAT(49x,10(f8.3,2x))
				RETURN
		END SUBROUTINE
