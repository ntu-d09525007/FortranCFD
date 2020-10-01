 SUBROUTINE BC(U)
 use dat, only : N
 IMPLICIT NONE
 REAL(8) :: U(-2:N+3)

!  Non-Reflective bondary condition
 U(1)= U(2)
 U(0)= U(1)
 U(-1)= U(1)
 U(-2)= U(1)
 U(N)= U(N-1)
 U(N+1)= U(N)
 U(N+2)= U(N)
 U(N+3)= U(N)

END SUBROUTINE

SUBROUTINE BC_R(U)
USE DAT, ONLY : N
REAL(8) :: U(-2:N+3)

U(1)=0.0D0
U(N)=0.0D0

U(0)=-U(2)
U(-1)=-U(3)
U(-2)=-U(4)

U(N+1)=-U(N-1)
U(N+2)=-U(N-2)
U(N+3)=-U(N-3)

CALL BC(U)

END SUBROUtine

