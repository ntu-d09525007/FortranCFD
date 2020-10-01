SUBROUTINE MPLS()
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA
USE DUMMY_DATA
IMPLICIT NONE
INTEGER :: I,J
  
 CALL MPLS_SOURCE(S0)
 
 !$omp parallel do
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
   PHI(i,j) = PHI(i,j) + s0(i,j)
 ENDDO
 enddo
 !$omp end parallel do
 
 CALL BC_LS()
 CALL MPLS_SOURCE(S0)
 
 !$omp parallel do
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
   PHI(i,j) = PHI(i,j) + s0(i,j)
 ENDDO
 enddo
 !$omp end parallel do
 
  CALL BC_LS()
 
END SUBROUTINE

SUBROUTINE MPLS_RK3()
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA
USE DUMMY_DATA
IMPLICIT NONE
INTEGER :: I,J

 CALL MPLS_SOURCE(S0)
 
 !$omp parallel do
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
   PHI(i,j) = PHI(i,j) + s0(i,j)
 ENDDO
 enddo
 !$omp end parallel do  
 
 CALL BC_LS()
 
 CALL MPLS_SOURCE(S1)

 !$omp parallel do
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
   PHI(i,j) = PHI(i,j) + (-3.0_DP*s0(i,j)+s1(i,j))/4.0_DP
 ENDDO
 enddo
 !$omp end parallel do  

 CALL BC_LS()
 
 CALL MPLS_SOURCE(S2)
 
 !$omp parallel do
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
   PHI(i,j)=PHI(i,j) + (-s0(i,j)-s1(i,j)+8.0_DP*s2(i,j))/12.0_DP
 enddo
 ENDDO
 !$omp end parallel do
 
 CALL BC_LS()
 
END SUBROUTINE

SUBROUTINE MPLS_SOURCE(S)
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA
IMPLICIT NONE
INTEGER :: I,J
REAL(DP) :: NX,NY,FS
REAL(DP),DIMENSION(-2:NODE_X+3,-2:NODE_Y+3) :: S

 MASS_LS = 0.0_DP
 FS = 0.0_DP

 CALL HEAVY_F(PHI) 
 
 !$OMP PARALLEL DO PRIVATE(NX,NY), REDUCTION(+:MASS_LS,FS)
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
 
	NX = 0.5_DP*(PHI(I+1,J)-PHI(I-1,J))/DX
    NY = 0.5_DP*(PHI(I,J+1)-PHI(I,J-1))/DY
	
	GRAD(I,J) = DSQRT(NX**2+NY**2)
	
	MASS_LS = MASS_LS + HEAVY(I,J)*(HEAVY(I,J) + RATIO_RHO*(1.0-HEAVY(I,J)))
	
	FS = FS + DELTA(I,J)**2*GRAD(I,J)*( 2.0_DP*(1.0_DP-RATIO_RHO)*HEAVY(I,J)+RATIO_RHO )

 ENDDO
 ENDDO
 !$OMP END PARALLEL DO
 
 FS=(IMASS_LS-MASS_LS)/FS
 
 !$OMP PARALLEL DO 
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
	S(I,J) = FS*DELTA(I,J)*GRAD(I,J)
 ENDDO
 ENDDO
 !$OMP END PARALLEL DO
 
END SUBROUTINE
