SUBROUTINE SOLVE_IBM()
USE PROBLEM_DEF
IMPLICIT NONE

 IF( IBM_SOLVER .EQ. 0) RETURN
 CALL SOLID_POS()
 CALL VELOCITY_CORRECTION()

END SUBROUTINE

SUBROUTINE SOLID_POS()
USE PRECISION
USE PROBLEM_DEF
USE IBM_DATA
IMPLICIT NONE
INTEGER :: I,J,K
 
 SOLID_Z = SOLID_Z + SOLID_W*DT

CALL FIND_SOLID

END SUBROUTINE

SUBROUTINE FIND_SOLID()
USE PRECISION
USE PROBLEM_DEF
USE IBM_DATA
IMPLICIT NONE
INTEGER :: I,J,K
INTEGER :: II,JJ,KK
REAL(DP) :: XS,YS,ZS

 IF( IBM_SOLVER .EQ. 0) RETURN

 !$OMP PARALLEL DO PRIVATE(II,JJ,KK,XS,YS,ZS)
 DO K = 1, NODE_Z
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
 
 	SOLID_VOF(I,J,K) = 0.0
	
	DO II = 1, UNIT_GRID2
    DO JJ = 1, UNIT_GRID2
    DO KK = 1, UNIT_GRID2
	
		XS = 0.5_DP*(X(I-1)+X(I)) + REAL(II,KIND=DP)*DX/REAL(UNIT_GRID2,KIND=DP)
    	YS = 0.5_DP*(Y(J-1)+Y(J)) + REAL(JJ,KIND=DP)*DY/REAL(UNIT_GRID2,KIND=DP) 
    	ZS = 0.5_DP*(Z(K-1)+Z(K)) + REAL(KK,KIND=DP)*DZ/REAL(UNIT_GRID2,KIND=DP)
		
		IF( ABS(XS-SOLID_X) <= 0.5_DP*DX .AND. ZS>=SOLID_Z )THEN
			SOLID_VOF(I,J,K) = SOLID_VOF(I,J,K) + 1.0
		ENDIF
		
	ENDDO
	ENDDO
	ENDDO
	
	SOLID_VOF(I,J,K) = SOLID_VOF(I,J,K) / REAL(UNIT_GRID2**3,KIND=DP)
	
 ENDDO
 ENDDO
 ENDDO
 !$OMP END PARALLEL DO

END SUBROUTINE

SUBROUTINE VELOCITY_CORRECTION()
USE PRECISION
USE PROBLEM_DEF
USE IBM_DATA
USE FLUID_PROPERTIES
IMPLICIT NONE
INTEGER :: I,J,K
REAL(DP) :: S

 !$omp parallel do private(s)
 do k = 1, node_z-1
 do j = 1, node_y
 do i = 1, node_x
	 s = 0.5_DP * (SOLID_VOF(i,j,k+1)+SOLID_VOF(i,j,k)) 
	 w(i,j,k) = w(i,j,k)*(1.0_DP-s) + s*SOLID_W
 enddo
 enddo
 enddo
 !$omp end parallel do
 
  !$omp parallel do private(s)
 do k = 1, node_z
 do j = 1, node_y-1
 do i = 1, node_x
	 s = 0.5_DP * (SOLID_VOF(i,j+1,k)+SOLID_VOF(i,j,k)) 
	 V(i,j,k) = V(i,j,k)*(1.0_DP-s) + s*SOLID_V
 enddo
 enddo
 enddo
 !$omp end parallel do
 
  !$omp parallel do private(s)
 do k = 1, node_z
 do j = 1, node_y
 do i = 1, node_x-1
	 s = 0.5_DP * (SOLID_VOF(i+1,j,k)+SOLID_VOF(i,j,k)) 
	 U(i,j,k) = U(i,j,k)*(1.0_DP-s) + s*SOLID_U
 enddo
 enddo
 enddo
 !$omp end parallel do

END SUBROUTINE