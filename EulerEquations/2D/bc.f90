 SUBROUTINE BC(U)
 use dat, only : N, M
 IMPLICIT NONE
 REAL(8) :: U(-2:M+3,-2:N+3)
 INTEGER :: I,J
!  Non-Reflective bondary condition

 DO I = 1, N
 DO J = 0, 3
	U(I,1-J) = U(I,M-J)
	U(I,M+J) = U(I,1+J)
 END DO
 END DO
 
 DO J = -2, M+3
 DO I = 0, 3
	U(1-I,J) = U(N-I,J)
	U(N+I,J) = U(1+I,J)
 END DO
 END DO
 
END SUBROUTINE

