SUBROUTINE POISSON_CG_SOLVER
USE PRECISION
USE PROBLEM_DEF
USE FLUID_PROPERTIES, ONLY : P, P_TMP, RHS
USE DUMMY_DATA, ONLY : A,B,C,F,G
IMPLICIT NONE
INTEGER :: I,J,K, INITER
REAL(DP) :: BB, AP2, LAM, MU

 INITER = 0

 ! A - R
 ! B - P
 ! C - S
 ! F - A*P
 ! G - AT*A*P

 CALL CG_DOT(RHS,RHS,BB)
 CALL CG_RES(P,RHS,A)
 CALL CG_AX(A,B)            ! Calculate P
 CALL CG_COPY(C,B)          ! Copy P to S

 DO
 	 
 	 INITER = INITER + 1
 	 	 
 	 CALL CG_AX(B,F)          ! Calculate A*P
 	 CALL CG_AX(F,G)          ! Calculate AT*A*P
 	 
 	 CALL CG_DOT(F,F,AP2)     
 	 CALL CG_DOT(B,C,LAM)     
 	 LAM = LAM / AP2

 	 CALL CG_UPDATE_P(P,P_TMP,LAM,B)
 	 CALL CG_SOL_STAB(P)
 	 CALL CG_UPDATE(A,-LAM,F)
 	 CALL CG_UPDATE(C,-LAM,G)
 	 
 	 CALL CG_DOT(G,C,MU)  
 	 MU = - MU / AP2
 	 
 	 CALL CG_UPDATE_2(B,MU,C)
 	 
 	 CALL CG_DOT(A,A,CG_ERROR)  
 	 LINF_P = DSQRT( CG_ERROR / BB )
 	 
 	 IF( LINF_P < TOL_P )EXIT
   IF( CG_ERROR < 1.0E-5 )EXIT
   
 	 IF( MOD(INITER,5000).EQ.0 )THEN
 	 	 WRITE(*,*)"Residual:",LINF_P
 	 endif
 	 
 	 if( initer > PRESS_MAX_IT )then
 	 	 write(2,*)"Pressure can not converge, stop at ",Linf_P
 	 	 stop
 	 endif
 	
 ENDDO
 
 CALL BC3D(P)
 
END SUBROUTINE

SUBROUTINE CG_RES(X,B,R)
USE PRECISION
USE PROBLEM_DEF, only : node_x,node_y,node_z
IMPLICIT NONE
REAL(DP), DIMENSION(-2:NODE_X+3,-2:NODE_Y+3,-2:NODE_Z+3) :: X,B,R
INTEGER :: I,J,K
REAL(DP):: CC,CR,CL,CF,CB,CU,CD,DX2

 DX2 = DX*DX

 !$OMP PARALLEL DO PRIVATE(CR,CL,CF,CB,CU,CD,CC)
 DO K = 1, NODE_Z
 DO J = 1, NODE_Y 
 DO I = 1, NODE_X
 	
 	 CR = 1.0_DP / DX2
 	 CL = CR
 	 CF = CR
 	 CB = CR
 	 CU = CR
 	 CD = CR
 	 CC = -(CR+CL+CU+CD+CF+CB)
 	 
 	 IF(I.EQ.1)THEN
 	 	 CC=CC+CL
 	 	 CL=0.0
 	 ENDIF
 	 
 	 IF(I.EQ.NODE_X)THEN
 	 	 CC=CC+CR
 	 	 CR=0.0
 	 ENDIF
 	 
 	 IF(J.EQ.1)THEN
 	 	 CC=CC+CB
 	 	 CB=0.0
 	 ENDIF
 	 
 	 IF(J.EQ.NODE_Y)THEN
 	 	 CC=CC+CF
 	 	 CF=0.0
 	 ENDIF
 	 
 	 IF(K.EQ.1)THEN
 	 	 CC=CC+CD
 	 	 CD=0.0
 	 ENDIF
 	 
 	 IF(K.EQ.NODE_Z)THEN
 	 	 CC=CC+CU
 	 	 CU=0.0
 	 ENDIF
 	 
 	 R(I,J,K) = B(I,J,K) - CC*X(I,J,K) &
 	                    &- CR*X(I+1,J,K) - CL*X(I-1,J,K) &
 	                    &- CF*X(I,J+1,K) - CB*X(I,J-1,K) &
 	                    &- CU*X(I,J,K+1) - CD*X(I,J,K-1)
 	   
 ENDDO
 ENDDO
 ENDDO
 !$OMP END PARALLEL DO


END SUBROUTINE

SUBROUTINE CG_AX(X,SOL)
USE PRECISION
USE PROBLEM_DEF, only : node_x,node_y,node_z
IMPLICIT NONE
REAL(DP), DIMENSION(-2:NODE_X+3,-2:NODE_Y+3,-2:NODE_Z+3) :: X, SOL
INTEGER :: I,J,K
REAL(DP):: CC,CR,CL,CF,CB,CU,CD,DX2

 DX2 = DX*DX

 !$OMP PARALLEL DO PRIVATE(CR,CL,CF,CB,CU,CD,CC)
 DO K = 1, NODE_Z
 DO J = 1, NODE_Y 
 DO I = 1, NODE_X
 	
 	 CR = 1.0_DP / DX2
 	 CL = CR
 	 CF = CR
 	 CB = CR
 	 CU = CR
 	 CD = CR
 	 CC = -(CR+CL+CU+CD+CF+CB)
 	 
 	 IF(I.EQ.1)THEN
 	 	 CC=CC+CL
 	 	 CL=0.0
 	 ENDIF
 	 
 	 IF(I.EQ.NODE_X)THEN
 	 	 CC=CC+CR
 	 	 CR=0.0
 	 ENDIF
 	 
 	 IF(J.EQ.1)THEN
 	 	 CC=CC+CB
 	 	 CB=0.0
 	 ENDIF
 	 
 	 IF(J.EQ.NODE_Y)THEN
 	 	 CC=CC+CF
 	 	 CF=0.0
 	 ENDIF
 	 
 	 IF(K.EQ.1)THEN
 	 	 CC=CC+CD
 	 	 CD=0.0
 	 ENDIF
 	 
 	 IF(K.EQ.NODE_Z)THEN
 	 	 CC=CC+CU
 	 	 CU=0.0
 	 ENDIF
 	 
 	 SOL(I,J,K) = CC*X(I,J,K) &
 	           &+ CR*X(I+1,J,K) + CL*X(I-1,J,K) &
 	           &+ CF*X(I,J+1,K) + CB*X(I,J-1,K) &
 	           &+ CU*X(I,J,K+1) + CD*X(I,J,K-1)
 	   
 ENDDO
 ENDDO
 ENDDO
 !$OMP END PARALLEL DO
 
END SUBROUTINE

SUBROUTINE CG_DOT(A,B,SOL)
USE PRECISION
USE PROBLEM_DEF, only : node_x,node_y,node_z
IMPLICIT NONE
REAL(DP), DIMENSION(-2:NODE_X+3,-2:NODE_Y+3,-2:NODE_Z+3) :: A,B
REAL(DP) :: SOL
INTEGER :: I,J,K

 SOL=0.0
 
 !$OMP PARALLEL DO REDUCTION(+:SOL)
 DO K = 1, NODE_Z
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
   SOL = SOL + A(I,J,K)*B(I,J,K)	
 ENDDO
 ENDDO
 ENDDO
 !$OMP END PARALLEL DO

END SUBROUTINE

SUBROUTINE CG_COPY(A,B)
USE PRECISION
USE PROBLEM_DEF, only : node_x,node_y,node_z
IMPLICIT NONE
REAL(DP), DIMENSION(-2:NODE_X+3,-2:NODE_Y+3,-2:NODE_Z+3) :: A, B
INTEGER :: I,J,K

 !$OMP PARALLEL DO
 DO K = -2, NODE_Z+3
 DO J = -2, NODE_Y+3
 DO I = -2, NODE_X+3
 	 A(I,J,K) = B(I,J,K)
 ENDDO
 ENDDO
 ENDDO
 !$OMP END PARALLEL DO
 
END SUBROUTINE

SUBROUTINE CG_UPDATE(X,LAM,S)
USE PRECISION
USE PROBLEM_DEF, only : node_x,node_y,node_z
IMPLICIT NONE
REAL(DP), DIMENSION(-2:NODE_X+3,-2:NODE_Y+3,-2:NODE_Z+3) :: X,S
INTEGER :: I,J,K
REAL(DP),INTENT(IN) :: LAM

 !$OMP PARALLEL DO
 DO K = 1, NODE_Z
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
   X(I,J,K) = X(I,J,K) + LAM*S(I,J,K)
 ENDDO
 ENDDO
 ENDDO
 !$OMP END PARALLEL DO
 
END SUBROUTINE

SUBROUTINE CG_UPDATE_2(X,MU,S)
USE PRECISION
USE PROBLEM_DEF, only : node_x,node_y,node_z
IMPLICIT NONE
REAL(DP), DIMENSION(-2:NODE_X+3,-2:NODE_Y+3,-2:NODE_Z+3) :: X,S
INTEGER :: I,J,K
REAL(DP),INTENT(IN) :: MU

 !$OMP PARALLEL DO
 DO K = 1, NODE_Z
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
   X(I,J,K) = S(I,J,K) + MU*X(i,J,K)
 ENDDO
 ENDDO
 ENDDO
 !$OMP END PARALLEL DO
 
END SUBROUTINE

SUBROUTINE CG_UPDATE_P(X,Y,LAM,S)
USE PRECISION
USE PROBLEM_DEF, only : node_x,node_y,node_z
IMPLICIT NONE
REAL(DP), DIMENSION(-2:NODE_X+3,-2:NODE_Y+3,-2:NODE_Z+3) :: X,Y,S
INTEGER :: I,J,K
REAL(DP),INTENT(IN) :: LAM

 LINF_P = 0.0

 !$OMP PARALLEL DO REDUCTION(MAX:LINF_P)
 DO K = 1, NODE_Z
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
 	 Y(I,J,K) = LAM*S(I,J,K)
   X(I,J,K) = X(I,J,K) + Y(I,J,K)
   LINF_P = MAX(LINF_P,ABS(Y(I,J,K)))
 ENDDO
 ENDDO
 ENDDO
 !$OMP END PARALLEL DO
 
END SUBROUTINE

SUBROUTINE CG_SOL_STAB(X)
USE PRECISION
USE PROBLEM_DEF, only : node_x,node_y,node_z
IMPLICIT NONE
REAL(DP), DIMENSION(-2:NODE_X+3,-2:NODE_Y+3,-2:NODE_Z+3) :: X
INTEGER :: I,J,K
REAL(DP) :: S

 S = 0.0
 
 !$OMP PARALLEL DO REDUCTION(+:S)
 DO K = 1, NODE_Z
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
   	S = S + X(I,J,K)
 ENDDO
 ENDDO
 ENDDO
 !$OMP END PARALLEL DO
 
 S = S / REAL(NODE_X*NODE_Y*NODE_Z,KIND=DP)
  
 !$OMP PARALLEL DO 
 DO K = 1, NODE_Z
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
   	X(I,J,K) = X(I,J,K) - S
 ENDDO
 ENDDO
 ENDDO
 !$OMP END PARALLEL DO
   
END SUBROUTINE