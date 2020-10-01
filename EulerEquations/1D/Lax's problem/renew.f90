SUBROUTINE RENEW()
USE DAT
IMPLICIT NONE
INTEGER :: I

 CALL BC(E)
 CALL BC(RHO)
 CALL BC_R(IX)

!$omp parallel do
DO I = -2, N+3
	U(I) =  IX(I)/RHO(I)
	P(I) = (gamma-1.0d0)*(E(I)-0.5d0*RHO(I)*U(I)**2.0d0 )
	H(I) = (E(i)+P(i))/rho(i)
	A(I) = dsqrt( gamma*P(I)/rho(I) )
END DO
!$omp end parallel do

 CALL BC(P)
 CALL BC_R(U)
 call BC(H)
 CALL BC(A)

END SUBROUTINE
