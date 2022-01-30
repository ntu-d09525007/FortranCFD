SUBROUTINE QUICK()
IMPLICIT NONE
 
  CALL VELA_Q
  CALL VELB_Q

END SUBROUTINE

SUBROUTINE QUICK_CELL_FACE(FI)
USE PRECISION
USE PROBLEM_DEF
USE DUMMY_DATA, ONLY : FP,FM,GP,GM
IMPLICIT NONE
INTEGER :: I,J
REAL(DP), DIMENSION(-2:NODE_X+3,-2:NODE_Y+3) :: FI

 !$OMP PARALLEL DO 
 DO J = 0, NODE_Y
 DO I = 0, NODE_X
   	FP(I,J) = (-FI(I-1,J)+6.0*FI(I,J)+3.0*FI(I+1,J))/8.0
   	FM(I,J) = (-FI(I+2,J)+6.0*FI(I+1,J)+3.0*FI(I,J))/8.0
   	
   	GP(I,J) = (-FI(I,J-1)+6.0*FI(I,J)+3.0*FI(I,J+1))/8.0
   	GM(I,J) = (-FI(I,J+2)+6.0*FI(I,J+1)+3.0*FI(I,J))/8.0
 ENDDO
 ENDDO
 !$OMP END PARALLEL DO
 
END SUBROUTINE

SUBROUTINE QUICK_U()
USE PROBLEM_DEF
USE PRECISION
USE DUMMY_DATA, ONLY : F,G,FP,FM,GP,GM
USE FLUID_PROPERTIES, ONLY : U,V
IMPLICIT NONE
INTEGER :: I,J
REAL(DP) :: UU,VV

 CALL QUICK_CELL_FACE(U)
  
 !$OMP PARALLEL DO PRIVATE(UU,VV)
 DO J = 0, NODE_Y
 DO I = 0, NODE_X
 	 
 	 UU = 0.5*(U(I,J)+U(I+1,J))
 	 VV = 0.5*(V(I,J)+V(I+1,J))
 	 
 	 IF(UU>=0.0)THEN
 	 	  F(I,J) = UU*FP(I,J)
 	 ELSE
 	 	  F(I,J) = UU*FM(I,J)
 	 ENDIF

 	 IF(VV>=0.0)THEN
 	 	  G(I,J) = VV*GP(I,J)
 	 ELSE
 	 	  G(I,J) = VV*GM(I,J)
 	 ENDIF
 	  	  	 
 ENDDO
 ENDDO
 !$OMP END PARALLEL DO
  	
END SUBROUTINE

SUBROUTINE QUICK_V()
USE PROBLEM_DEF
USE PRECISION
USE DUMMY_DATA, ONLY : F,G,FP,FM,GP,GM
USE FLUID_PROPERTIES, ONLY : U,V
IMPLICIT NONE
INTEGER :: I,J
REAL(DP) :: UU,VV

 CALL QUICK_CELL_FACE(V)
  
 !$OMP PARALLEL DO PRIVATE(UU,VV)
 DO J = 0, NODE_Y
 DO I = 0, NODE_X
 	 
 	 UU = 0.5*(U(I,J)+U(I,J+1))
 	 VV = 0.5*(V(I,J)+V(I,J+1))
 	 
 	 IF(UU>=0.0)THEN
 	 	  F(I,J) = UU*FP(I,J)
 	 ELSE
 	 	  F(I,J) = UU*FM(I,J)
 	 ENDIF

 	 IF(VV>=0.0)THEN
 	 	  G(I,J) = VV*GP(I,J)
 	 ELSE
 	 	  G(I,J) = VV*GM(I,J)
 	 ENDIF

 ENDDO
 ENDDO
 !$OMP END PARALLEL DO
  	
END SUBROUTINE

SUBROUTINE VELA_Q()
USE PRECISION
USE PROBLEM_DEF
USE FLUID_PROPERTIES
USE LS_DATA
use vof_data
USE DUMMY_DATA, ONLY : F,G
IMPLICIT NONE
INTEGER :: I,J
REAL(DP) :: DX2,DY2,RE2X,RE2Y
REAL(DP) :: CM,CR,CC,DETF,CURVK,FINX,A1,A2,A3,A4,A5,A6,A7,TMP
    
  DX2=DX*DX
  DY2=DY*DY

  RE2X=RE*DX2
  RE2Y=RE*DY2
  
  CALL QUICK_U
  
 !$omp parallel do private(cm,cr,cc,detf,curvk,finx,a1,a2,a3,a4,a5,a6,a7,TMP)
 DO J = 1, NODE_Y
 DO I = 1, NODE_X-1

    CM=(AMU_OLD(I,J)+AMU_OLD(I+1,J))
    CR=(RHO_OLD(I,J)+RHO_OLD(I+1,J))
    CC=CM/CR
    
    !detf = 0.5_DP*(DELTA(i,j,k)+DELTA(i+1,j,k))
    TMP = HEAV(0.5_DP*(PHI(I,J)+PHI(I+1,J)),DETF)
    curvk= 0.5_DP*(curv(i,j)+curv(i+1,j))
    finx = 0.5_DP*(NORMAL_X(i,j)+NORMAL_X(i+1,j))
    
    IF(INTERFACE_METHOD.eq.1)THEN
    	DETF=1.0
    	FINX=(VOF(I+1,J)-VOF(I,J))/DX
    ENDIF

    A1=1.0_DP/DX*(F(I,J)-F(I-1,J))
    
    A2=1.0_DP/DX*(G(I,J)-G(I,J-1))
       
    A4=1.0_DP/RE2X*CC*(U(I+1,J)-2.0_DP*U(I,J)+U(I-1,J))
   
    A5=1.0_DP/RE2Y*CC*(U(I,J+1)-2.0_DP*U(I,J)+U(I,J-1))
   	
    if(ST_FORCE==0)THEN
    	A7=0.0_dp
    ELSE IF(ST_FORCE==1)THEN
    	A7=2.0_DP/(CR*WE)*curvk*detf*finx
    ELSE IF(ST_FORCE==2)THEN
    	A7=2.0_DP/(CR*WE)*curvk*detf*finx*2.0_DP*TMP
    ENDIF
    
    US(I,J)=A1+A2-A4-A5+A7

  ENDDO
  ENDDO
  !$omp end parallel do

  END SUBROUTINE VELA_Q   
  
 SUBROUTINE VELB_Q()
 USE PRECISION
 USE PROBLEM_DEF
 USE FLUID_PROPERTIES
 USE LS_DATA
 use vof_data
 USE DUMMY_DATA, ONLY : F,G
 IMPLICIT NONE
 INTEGER :: I,J
 REAL(DP) :: DX2,DY2,RE2X,RE2Y
 REAL(DP) :: CM,CC,CR,DETF,FINY,CURVK,B1,B2,B3,B4,B5,B6,B7,TMP
 
  DX2=DX*DX
  DY2=DY*DY

  RE2X=RE*DX2
  RE2Y=RE*DY2
 
  CALL QUICK_V
 
 !$omp parallel do private(cm,cc,cr,detf,finy,curvk,b1,b2,b3,b4,b5,b6,b7,TMP)
 DO J = 1, NODE_Y-1
 DO I = 1, NODE_X
    
    CM=(AMU_OLD(I,J)+AMU_OLD(I,J+1))
    CR=(RHO_OLD(I,J)+RHO_OLD(I,J+1))
    CC=CM/CR
    
    !detf = 0.5_DP*(DELTA(i,j,k)+DELTA(i,j+1,k))
    TMP = HEAV(0.5_DP*(PHI(I,J)+PHI(I,J+1)),DETF)
    curvk= 0.5_DP*(curv(i,j)+curv(i,j+1))
    finy = 0.5_DP*(NORMAL_Y(i,j)+NORMAL_Y(i,j+1))   

    IF(INTERFACE_METHOD.eq.1)THEN
      DETF=1.0
      FINY=(VOF(I,J+1)-VOF(I,J))/DX
    ENDIF
    
    B1=1.0_DP/DX*(F(I,J)-F(I-1,J))
    
    B2=1.0_DP/DY*(G(I,J)-G(I,J-1))
   
    B4=1.0_DP/RE2X*CC*(V(I+1,J)-2.0_DP*V(I,J)+V(I-1,J))
   
    B5=1.0_DP/RE2Y*CC*(V(I,J+1)-2.0_DP*V(I,J)+V(I,J-1))
    
    if(ST_FORCE==0)THEN
    	B7=0.0_dp
    ELSE IF(ST_FORCE==1)THEN
    	B7=2.0/(CR*WE)*curvk*detf*finy
    ELSE IF(ST_FORCE==2)THEN
    	B7=2.0/(CR*WE)*curvk*detf*finy*2.0_DP*TMP
    ENDIF

    VS(I,J)=B1+B2-B4-B5+B7
	
	IF(G_FORCE==1)VS(i,j)=VS(i,j)+1.0/FR**2

  ENDDO
  ENDDO
  !$omp end parallel do

  END SUBROUTINE VELB_Q 