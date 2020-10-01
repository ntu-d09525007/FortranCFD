SUBROUTINE POISSON_SOLVER()
USE PRECISION
USE PROBLEM_DEF
USE FLUID_PROPERTIES, ONLY : P, P_TMP, RHS
IMPLICIT NONE
INTEGER :: I,J,K, INITER
REAL(DP) :: SUMP,CC,CR,CL,CF,CB,CU,CD,DX2


 INITER = 0
 DX2 = DX*DX
 
 DO
 	
 	 INITER = INITER + 1
 	 
 	 SUMP = 0.0
 	 
 	 !$OMP PARALLEL DO
 	 DO K = 1, NODE_Z
 	 DO J = 1, NODE_Y
 	 DO I = 1, NODE_X
 	   	P_TMP(I,J,K) = P(I,J,K)
 	 ENDDO
 	 ENDDO
 	 ENDDO
 	 !$OMP END PARALLEL DO
 	 
 	 !$OMP PARALLEL DO PRIVATE(CR,CL,CF,CB,CU,CD,CC), REDUCTION(+:SUMP)
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
 	 	 
 	   P(I,J,K) = RHS(I,J,K) - CR*P(I+1,J,K) - CL*P(I-1,J,K) &
 	                         - CF*P(I,J+1,K) - CB*P(I,J-1,K) &
 	                         - CU*P(I,J,K+1) - CD*P(I,J,K-1)
 	                         
 	   P(I,J,K) = P(I,J,K) / CC	
 	   
 	   P(i,j,k) = P(i,j,k) * 1.5 - 0.5 * P_tmp(i,j,k)
 	   
 	   SUMP = SUMP + P(I,J,K)
 	    	 	 
 	 ENDDO
 	 ENDDO
 	 ENDDO
 	 !$OMP END PARALLEL DO
 	 
 	 SUMP = SUMP / REAL(NODE_X*NODE_Y*NODE_Z,KIND=DP)
 	
 	 LINF_P = 0.0
 	 
 	 !$OMP PARALLEL DO REDUCTION(MAX:LINF_P)
 	 DO K = 1, NODE_Z
 	 DO J = 1, NODE_Y
 	 DO I = 1, NODE_X
 	    P(I,J,K) = P(I,J,K) - SUMP 	
 	    LINF_P = MAX(LINF_P,ABS(P(I,J,K)-P_TMP(I,J,K)))
 	 ENDDO
 	 ENDDO
 	 ENDDO
 	 !$OMP END PARALLEL DO
 	 
 	 IF( LINF_P < TOL_P )EXIT
 	 
 	 IF( INITER > PRESS_MAX_IT )THEN
 	 	 WRITE(2,*)"Pressure can not converge, ",Linf_P
 	   STOP
 	 ENDIF
 	 
 	 if (mod(initer,5000)==0)then 
      WRITE(*,*) "AMAXP=" , LINF_P
   end if
 	 
 ENDDO

END SUBROUTINE

