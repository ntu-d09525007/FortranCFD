SUBROUTINE LS_REDISTANCE(TARGET,BTN)
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA
USE DUMMY_DATA
IMPLICIT NONE
INTEGER :: BTN,I,J
REAL(DP),DIMENSION(-2:NODE_X+3,-2:NODE_Y+3) :: TARGET
REAL(DP) :: AMX, RED_TIME, ERR, TIMESTOP

IF( BTN==0 .AND. BTN==2 )THEN
	
  !$OMP PARALLEL DO
  DO J = 1, NODE_Y
  DO I = 1, NODE_X
	  SGN(I,J) = TARGET(I,J)
  END DO
  END DO
  !$OMP END PARALLEL DO
  	
  CALL GRADFI(TARGET)
  
  AMX = MAXVAL(GRAD(1:NODE_X,1:NODE_Y))
  
  !$OMP PARALLEL DO
  DO J = 1, NODE_Y
  DO I = 1, NODE_X
	  TARGET(I,J) = TARGET(I,J) / AMX
  END DO
  END DO
  !$OMP END PARALLEL DO
  
END IF

  CALL sgnf(TARGET) 
  
  RED_TIME = 0.0_DP
 
  SELECT CASE(BTN)
	CASE(0)
	  TIMESTOP=0.5*MAX(XL,YL)
	CASE(1)
	  TIMESTOP=2.0*INTERFACE_WIDTH
	CASE(2)
	  TIMESTOP=2.0*INTERFACE_WIDTH
	CASE(3)
	  TIMESTOP=2.0*INTERFACE_WIDTH
  END SELECT
  
  DO
  	
  	RED_TIME = RED_TIME + RDT
  	
    !$OMP PARALLEL DO
    DO J = 1, NODE_Y
    DO I = 1, NODE_X
	    PHI_TMP(I,J) = TARGET(I,J)
    END DO
    END DO
    !$OMP END PARALLEL DO
  	
    CALL LS_REDIS(TARGET,BTN)
  	
    CALL FIND_ERROR(ERR,PHI_TMP,TARGET,2)
    
  	IF(RED_TIME>TIMESTOP)EXIT
	
  END DO

END SUBROUTINE

SUBROUTINE LS_REDIS(TARGET,BTN)
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA
USE DUMMY_DATA
IMPLICIT NONE
INTEGER :: BTN,I,J
REAL(DP),DIMENSION(-2:NODE_X+3,-2:NODE_Y+3) :: TARGET

 CALL GRADFI(TARGET)
 CALL RIS_MASS(TARGET,BTN)
 CALL LS_REDIS_SORCE(S0)
 
 !$omp parallel do
 do j=1,NODE_Y
 do i=1,NODE_X
   TARGET(i,j)=TARGET(i,j)-RDT*S0(i,j)
 enddo
 enddo
 !$omp end parallel do
 
 call bc3d(TARGET)

 CALL GRADFI(TARGET)
 CALL RIS_MASS(TARGET,BTN)
 CALL LS_REDIS_SORCE(S1)
 
 !$omp parallel do
 do j=1,NODE_Y
 do i=1,NODE_X
   TARGET(i,j)=TARGET(i,j)-RDT/4.0_DP*(-3.0_DP*S0(i,j)+S1(i,j))
 enddo
 enddo
 !$omp end parallel do
 
 call bc3d(TARGET)

 CALL GRADFI(TARGET)
 CALL RIS_MASS(TARGET,BTN)
 CALL LS_REDIS_SORCE(S2)
 
 !$omp parallel do
 do j=1,NODE_Y
 do i=1,NODE_X
   TARGET(i,j)=TARGET(i,j)-RDT/12.0_DP*(-S0(i,j)-S1(i,j)+8.0_DP*S2(i,j))
 enddo
 enddo
 !$omp end parallel do
 
 call bc3d(TARGET)
   
END SUBROUTINE

SUBROUTINE LS_REDIS_SORCE(S)
USE PRECISION
USE PROBLEM_DEF
USE DUMMY_DATA
USE LS_DATA
IMPLICIT NONE
INTEGER :: I,J
REAL(DP),DIMENSION(-2:NODE_X+3,-2:NODE_Y+3) :: S

 !$omp parallel do
 do j=1,NODE_Y
 do i=1,NODE_X
     S(i,j)=(SGN(i,j)-C(i,j))*GRAD(i,j)-SGN(i,j)
 enddo
 enddo
 !$omp end parallel do

END SUBROUTINE

SUBROUTINE RIS_MASS(TARGET,BTN)
USE PRECISION
USE PROBLEM_DEF
USE DUMMY_DATA
USE LS_DATA
IMPLICIT NONE
INTEGER :: I,J,BTN
REAL(DP),DIMENSION(-2:NODE_X+3,-2:NODE_Y+3) :: TARGET
REAL(DP) :: SUMA, SUMB

 IF( BTN==0 .and. btn==2 )THEN
 	
    !$omp parallel do
 	do j=1,NODE_Y
	do i=1,NODE_X
       C(I,J) = 0.0_DP
    enddo
    enddo
    !$omp end parallel do
 	
 ELSE
 	
    CALL HEAVY_F(TARGET)
    
    !$omp parallel do
 	do j=1,NODE_Y
	do i=1,NODE_X
       A(i,j)= (SGN(i,j)*DELTA(i,j))*(GRAD(i,j)-1.0_dp)
       B(i,j)= GRAD(i,j)*DELTA(i,j)**2
    enddo
    enddo
    !$omp end parallel do
      
    call bc3d(a)
    call bc3d(b) 	
    
	 !$omp parallel do private(suma,sumb)
     do j=1,NODE_Y
	 do i=1,NODE_X

	    suma=16.0_dp*a(i,j)+sum(a(i-1:i+1,j-1:j+1))
	  	sumb=16.0_dp*b(i,j)+sum(b(i-1:i+1,j-1:j+1))

      if(sumb /=0.0_dp ) then
		    c(i,j)=suma/sumb
      else
		    c(i,j)=0.0_dp
	  end if
      
      c(i,j) = c(i,j)*DELTA(i,j)
      
	  enddo
	  enddo
	  !$omp end parallel do  
  
 END IF

END SUBROUTINE

SUBROUTINE GRADFI(TARGET)
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA
USE DUMMY_DATA
IMPLICIT NONE
INTEGER :: I,J,K
REAL(DP),DIMENSION(-2:NODE_X+3,-2:NODE_Y+3) :: TARGET
     
  CALL FV_DERI(TARGET)
  
  !$OMP PARALLEL DO
  DO J = 1, NODE_Y
  DO I = 1, NODE_X
  	GRAD(I,J) = gradf(up(I,J),um(I,J),vp(I,J),vm(I,J),SGN(I,J))  	 
  END DO
  END DO
  !$OMP END PARALLEL DO

END SUBROUTINE

subroutine FV_DERI(TARGET)
USE PRECISION
USE PROBLEM_DEF
USE DUMMY_DATA
IMPLICIT NONE
INTEGER :: I,J
REAL(DP),DIMENSION(-2:NODE_X+3,-2:NODE_Y+3) :: TARGET

!$OMP PARALLEL DO
DO J = 1, NODE_Y
DO I = 1, NODE_X
  F(I,J) = (TARGET(I,J)-TARGET(I-1,J))/DX
  G(I,J) = (TARGET(I,J)-TARGET(I,J-1))/DY
 END DO
 END DO
!$OMP END PARALLEL DO

 CALL BC3D(F)
 CALL BC3D(G)

!$OMP PARALLEL DO
 DO J = 1, NODE_Y
  CALL weno_js(F(-2:NODE_X+3,j),UP(-2:NODE_X+3,j),UM(-2:NODE_X+3,j),NODE_X)
 END DO
 !$OMP END PARALLEL DO
 
 !$OMP PARALLEL DO
 DO I = 1, NODE_X
  CALL weno_js(G(i,-2:NODE_Y+3),VP(i,-2:NODE_Y+3),VM(i,-2:NODE_Y+3),NODE_Y)
 END DO
!$OMP END PARALLEL DO


end subroutine

function gradf(up,um,vp,vm,fi0)
USE PRECISION
IMPLICIT NONE
REAL(KIND=DP) :: GRADF    
real(KIND=DP)  :: up,um,vp,vm,upm,upp,umm,ump,vpm,vpp,vmm,vmp,fi0

  upm=-MIN(up,0.0_DP)
  upp= MAX(up,0.0_DP)
  umm=-MIN(um,0.0_DP)
  ump= MAX(um,0.0_DP)
  vpm=-MIN(vp,0.0_DP)
  vpp= MAX(vp,0.0_DP)
  vmm=-MIN(vm,0.0_DP)
  vmp= MAX(vm,0.0_DP)
  
if (fi0 >=0.0_DP ) then
  gradf=sqrt(MAX(upm,ump)**2+MAX(vpm,vmp)**2)
else
  gradf=sqrt(MAX(upp,umm)**2+MAX(vpp,vmm)**2)
endif

end FUNCTION 

SUBROUTINE HEAVY_F(TARGET)
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA
IMPLICIT NONE
INTEGER :: I,J,K
REAL(DP),DIMENSION(-2:NODE_X+3,-2:NODE_Y+3) :: TARGET

 !$OMP PARALLEL DO
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
   HEAVY(I,J) = HEAV(TARGET(I,J),DELTA(I,J))    
 END DO
 END DO
 !$OMP END PARALLEL DO

END SUBROUTINE

FUNCTION HEAV(X,HP)  
USE PRECISION
USE LS_DATA
IMPLICIT NONE     
REAL(KIND=DP) :: HEAV, X, HP

IF(X > INTERFACE_WIDTH) THEN
  HEAV=1.0_DP
  HP=0.0_DP
ELSE IF(X < -INTERFACE_WIDTH) THEN
  HEAV=0.0_DP
  HP=0.0_DP
ELSE
  HEAV=0.5_DP*(1.0_DP+X/INTERFACE_WIDTH+1.0_DP/PI*SIN(PI*X/INTERFACE_WIDTH))
  HP=0.5_DP*(1.0_DP/INTERFACE_WIDTH+1.0_DP/INTERFACE_WIDTH*COS(PI*X/INTERFACE_WIDTH))
ENDIF

END FUNCTION HEAV 

FUNCTION SGNN(X)  
USE PRECISION  
IMPLICIT NONE     
REAL(KIND=DP) :: SGNN,X,HP,H

H=HEAV(X,HP)
SGNN=2.0_DP*(H-0.5_DP)

END FUNCTION SGNN

subroutine sgnf(TARGET) 
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA
IMPLICIT NONE
INTEGER :: I,J
REAL(DP),DIMENSION(-2:NODE_X+3,-2:NODE_Y+3) :: TARGET

 !$OMP PARALLEL DO
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
     SGN(i,j)=SGNN(TARGET(I,J))
 Enddo
 ENDDO
 !$omp end parallel do

end subroutine sgnf