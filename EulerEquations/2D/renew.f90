SUBROUTINE RENEW()
USE DAT
IMPLICIT NONE
INTEGER :: I

 CALL BC(E)
 CALL BC(RHO)
 CALL BC(IX)
 CALL BC(IY)
 
!$omp parallel do
DO J = -2, M+3
DO I = -2, N+3
	U(I,J) = IX(I,J)/RHO(I,J)
	V(I,J) = IY(I,J)/RHO(I,J)
	P(I,J) = (gamma-1.0d0)*(E(I,J)-0.5d0*RHO(I,J)*(U(I,J)**2.0d0+V(I,J)**2.0D0))
	A(I,J) = dsqrt( gamma*P(I,J)/rho(I,J) )
	H(i,j) = (E(I,j)+P(I,j))/rho(I,j)
END DO
END DO
!$omp end parallel do

 CALL BC(P)
 CALL BC(U)
 CALL BC(V)
 CALL BC(A)
 CALL BC(H)

END SUBROUTINE
