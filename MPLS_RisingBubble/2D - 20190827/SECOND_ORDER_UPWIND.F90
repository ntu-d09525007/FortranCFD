SUBROUTINE SECOND_ORDER_UPWIND()
IMPLICIT NONE

 CALL VELA()
 CALL VELB()

END SUBROUTINE

SUBROUTINE VELA()
USE PRECISION
USE PROBLEM_DEF
USE FLUID_PROPERTIES
USE LS_DATA
use vof_data
IMPLICIT NONE
INTEGER :: I,J
REAL(DP) :: DX2,DY2,DZ2,RE2X,RE2Y,RE2Z
REAL(DP) :: UP,UM,VV,VP,VM,CM,CR,CC,DETF,CURVK,FINX,A1,A2,A3,A4,A5,A6,A7,TMP
    
  DX2=DX*DX
  DY2=DY*DY

  RE2X=RE*DX2
  RE2Y=RE*DY2

 
 !$omp parallel do private(up,um,vv,vp,vm,cm,cr,cc,detf,curvk,finx,a1,a2,a3,a4,a5,a6,a7,TMP)
 DO J = 1, NODE_Y
 DO I = 1, NODE_X-1

    UP=0.5_DP*(U(I,J)+ABS(U(I,J)))
    UM=0.5_DP*(U(I,J)-ABS(U(I,J)))
    
    VV=0.25_DP*(V(I+1,J)+V(I,J)+V(I+1,J-1)+V(I,J-1))
    VP=0.5_DP*(VV+ABS(VV))
    VM=0.5_DP*(VV-ABS(VV))
    
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
 
    A1=0.5_DP/DX*(UP*(3.0_DP*U(I,J)-4.0_DP*U(I-1,J)+U(I-2,J)) &
                 +UM*(-U(I+2,J)+4.0_DP*U(I+1,J)-3.0_DP*U(I,J)))
   
    A2=0.5_DP/DY*(VP*(3.0_DP*U(I,J)-4.0_DP*U(I,J-1)+U(I,J-2)) &
                 +VM*(-U(I,J+2)+4.0_DP*U(I,J+1)-3.0_DP*U(I,J)))
      
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

  END SUBROUTINE VELA
  
  
  
 SUBROUTINE VELB()
 USE PRECISION
 USE PROBLEM_DEF
 USE FLUID_PROPERTIES
 USE LS_DATA
 use vof_data
 IMPLICIT NONE
 INTEGER :: I,J,K
 REAL(DP) :: DX2,DY2,DZ2,RE2X,RE2Y,RE2Z
 REAL(DP) :: UU,UP,UM,VP,VM,CM,CC,CR,DETF,FINY,CURVK,B1,B2,B3,B4,B5,B6,B7,TMP
 
  DX2=DX*DX
  DY2=DY*DY
  RE2X=RE*DX2
  RE2Y=RE*DY2
 
 !$omp parallel do private(UU,up,um,vp,vm,cm,cc,cr,detf,finy,curvk,b1,b2,b3,b4,b5,b6,b7,TMP)
 DO J = 1, NODE_Y-1
 DO I = 1, NODE_X

    UU=0.25_DP*(U(I,J+1)+U(I,J)+U(I-1,J+1)+U(I-1,J))
    UP=0.5_DP*(UU+ABS(UU))
    UM=0.5_DP*(UU-ABS(UU))
    
    VP=0.5_DP*(V(I,J)+ABS(V(I,J)))
    VM=0.5_DP*(V(I,J)-ABS(V(I,J)))
    
    CM=(AMU_OLD(I,J)+AMU_OLD(I,J+1))
    CR=(RHO_OLD(I,J)+RHO_OLD(I,J+1))
    CC=CM/CR
    
    !detf = 0.5_DP*(DELTA(i,j)+DELTA(i,j+1))
	TMP = HEAV(0.5_DP*(PHI(I,J)+PHI(I,J+1)),DETF)
    curvk= 0.5_DP*(curv(i,j)+curv(i,j+1))
    finy = 0.5_DP*(NORMAL_Y(i,j)+NORMAL_Y(i,j+1))   

    B1=0.5_DP/DX*(UP*(3.0_DP*V(I,J)-4.0_DP*V(I-1,J)+V(I-2,J)) &
                 +UM*(-V(I+2,J)+4.0_DP*V(I+1,J)-3.0_DP*V(I,J)))
   
    B2=0.5_DP/DY*(VP*(3.0_DP*V(I,J)-4.0_DP*V(I,J-1)+V(I,J-2)) &
                 +VM*(-V(I,J+2)+4.0_DP*V(I,J+1)-3.0_DP*V(I,J)))
   
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

  END SUBROUTINE VELB 
  