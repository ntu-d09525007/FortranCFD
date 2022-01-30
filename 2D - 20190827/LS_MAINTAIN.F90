SUBROUTINE LS_MAINTAIN()
USE PRECISION
USE LS_DATA
USE VOF_DATA
USE DUMMY_DATA
USE PROBLEM_DEF
IMPLICIT NONE
INTEGER :: I,J


  IF(INTERFACE_METHOD==0)THEN
  	
    IF(MOD(ITER,5)==0)CALL LS_REDISTANCE(PHI,1) 
  	
  ELSE IF( INTERFACE_METHOD==2 )THEN
  	
    IF( MOD(ITER,10)==0 )THEN
    	
  	    !$OMP PARALLEL DO
  	    DO J = -2, NODE_Y+3
  	    DO I = -2, NODE_X+3
  	    	PHI_V(I,J) = 2.0_DP*VOF(I,J)-1.0_DP
  	    END DO
  	    END DO
  	    !$OMP END PARALLEL DO
  	    
  	    CALL LS_REDISTANCE(PHI_V,2)

  	    !$OMP PARALLEL DO
  	    DO J = -2, NODE_Y+3
  	    DO I = -2, NODE_X+3
  	    	IF( ABS(PHI(I,J))<1.5*INTERFACE_WIDTH )PHI(I,J)=PHI_V(I,J)
  	    END DO
  	    END DO
  	    !$OMP END PARALLEL DO
  	    
  	    !IF(METHOD_CNT==1)THEN
 			CALL LS_REDISTANCE(PHI,1)
		!ELSE
		!    CALL LS_REDISTANCE(PHI,3)
		!ENDIF
		
  	    
    ELSE IF( MOD(ITER,5)==0 )THEN
  	
        !CALL LS_REDISTANCE(PHI,1)
        
    END IF
  	      	
  ELSE IF( INTERFACE_METHOD==3 )THEN

  	    !$OMP PARALLEL DO
  	    DO J = -2, NODE_Y+3
  	    DO I = -2, NODE_X+3
  	    	PHI(I,J) = 2.0_DP*VOF(I,J)-1.0_DP
  	    END DO
  	    END DO
  	    !$OMP END PARALLEL DO
  	    
  	    CALL LS_REDISTANCE(PHI,2)  	
  	
  END IF
  
END SUBROUTINE
