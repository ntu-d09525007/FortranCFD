SUBROUTINE GS_PROJ_SOLVER()
use problem_def
IMPLICIT NONE

 CALL Adams_Bashforth()
 CALL GS_RHSP()
 CALL GS_PRESS_POISSON()
 CALL SOLVE_MOMENTUM()

END SUBROUTINE

SUBROUTINE GS_RHSP()
USE PRECISION
USE PROBLEM_DEF
USE FLUID_PROPERTIES
IMPLICIT NONE
INTEGER :: I,J,K

!$omp parallel do
DO  K=0,NODE_Z+1
DO  J=0,NODE_Y+1
  U1(0,J,K)=0.0_DP
  U1(NODE_X,J,K)=0.0_DP
ENDDO
ENDDO
!$omp end parallel do

!$omp parallel do
  DO K=0,NODE_Z+1
  DO I=0,NODE_X+1
   V1(I,0,K)=0.0_DP
   V1(I,NODE_Y,K)=0.0_DP
  ENDDO
  ENDDO
!$omp end parallel do
 
!$omp parallel do
  DO J=0,NODE_Y+1
  DO I=0,NODE_X+1
   W1(I,J,0)=0.0_DP
   W1(I,J,NODE_Z)=0.0_DP
  ENDDO
  ENDDO
!$omp end parallel do

!$omp parallel do
 DO K = 1, NODE_Z
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
   RHS(I,J,K)=(U_OLD(I,J,K)/DT-U_OLD(I-1,J,K)/DT-U1(I,J,K)+U1(I-1,J,K))/DX &
             +(V_OLD(I,J,K)/DT-V_OLD(I,J-1,K)/DT-V1(I,J,K)+V1(I,J-1,K))/DY &
             +(W_OLD(I,J,K)/DT-W_OLD(I,J,K-1)/DT-W1(I,J,K)+W1(I,J,K-1))/DZ
  ENDDO
  ENDDO
  ENDDO
!$omp end parallel do
  

END SUBROUTINE

SUBROUTINE GS_PRESS_POISSON()
USE PRECISION
USE PROBLEM_DEF
USE FLUID_PROPERTIES
USE DUMMY_DATA, ONLY : FP, FM, GP, GM, HP, HM, C
IMPLICIT NONE
INTEGER :: I,J,K,ICOUNTP
REAL(DP) :: RR,RL,RF,RB,RU,RD,CL,CR,CF,CB,CD,CC,CU,DX2,DY2,DZ2
REAL(DP) :: POLD,PCAL,SUMP

  DX2=DX*DX
  DY2=DY*DY
  DZ2=DZ*DZ 

!$omp parallel do private(rr,rl,rf,rb,ru,rd,cl,cr,cf,cb,cd,cc,cu,pold,pcal)
DO K=1,NODE_Z
DO J=1,NODE_Y
DO I=1,NODE_X

  RR=0.5*(RHO(I,J,K)+RHO(I+1,J,K))
  RL=0.5*(RHO(I,J,K)+RHO(I-1,J,K))
  RF=0.5*(RHO(I,J,K)+RHO(I,J+1,K))
  RB=0.5*(RHO(I,J,K)+RHO(I,J-1,K))
  RU=0.5*(RHO(I,J,K)+RHO(I,J,K+1))
  RD=0.5*(RHO(I,J,K)+RHO(I,J,K-1))
  
  CL=1.0/(DX2*RL)
  CR=1.0/(DX2*RR)
  CF=1.0/(DY2*RF)
  CB=1.0/(DY2*RB)
  CU=1.0/(DZ2*RU)
  CD=1.0/(DZ2*RD)
  
  CC=-(CL+CR+CF+CB+CU+CD)

  IF(I.EQ.1) THEN
    CC=CC+CL
    CL=0.0_DP
  ENDIF

  IF(I.EQ.NODE_X) THEN
    CC=CC+CR
    CR=0.0_DP
  ENDIF

  IF(J.EQ.1) THEN
    CC=CC+CB
    CB=0.0_DP
  ENDIF

  IF(J.EQ.NODE_Y) THEN
    CC=CC+CF
    CF=0.0_DP
  ENDIF

  IF(K.EQ.1) THEN
    CC=CC+CD
    CD=0.0_DP
  ENDIF

  IF(K.EQ.NODE_Z) THEN
    CC=CC+CU
    CU=0.0_DP
  ENDIF
  
  FP(I,J,K) = CR
  FM(I,J,K) = CL
  
  GP(I,J,K) = CF
  GM(I,J,K) = CB
  
  HP(I,J,K) = CU
  HM(I,J,K) = CD
  
  C(I,J,K) = CC
  
ENDDO
ENDDO
ENDDO
!$omp end parallel do
  
  ICOUNTP=0
  10 LINF_P = 0.0
  ICOUNTP=ICOUNTP+1
  
!$omp parallel do PRIVATE(PCAL)
DO K=1,NODE_Z
DO J=1,NODE_Y
DO I=1,NODE_X

  P_TMP(i,j,k) = P(i,j,k)

  PCAL=RHS(I,J,K)-FM(I,J,K)*P(I-1,J,K)-FP(I,J,K)*P(I+1,J,K) &
                 -GM(I,J,K)*P(I,J-1,K)-GP(I,J,K)*P(I,J+1,K) &
                 -HM(I,J,K)*P(I,J,K-1)-HP(I,J,K)*P(I,J,K+1)
				 
  PCAL=PCAL/C(I,J,K)
  
  P(I,J,K) = PCAL
  
ENDDO
ENDDO
ENDDO
!$omp end parallel do

sump = 0.0
  
!$OMP PARALLEL DO reduction(+:SUMP)
 DO K = 1, NODE_Z
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
  P(I,J,K) = PRESS_OMEGA*P(I,J,K)+(1.0_DP-PRESS_OMEGA)*P_TMP(i,j,k)
  SUMP  = SUMP + P(I,J,K)
 END DO
 END DO
 END DO
!$OMP END PARALLEL DO

SUMP = SUMP / REAL(NODE_X*NODE_Y*NODE_Z,KIND=DP)
!SUMP = SUMP / REAL(NODE_X*NODE_Y,KIND=DP)
  
!$OMP PARALLEL DO reduction(max:linf_p)
 DO K = 1, NODE_Z
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
  P(I,J,K) = P(I,J,K) - SUMP
  LINF_P = max(LINF_P,abs(P(i,j,k)-P_TMP(i,j,k)))
 END DO
 END DO
 END DO
!$OMP END PARALLEL DO


if (mod(ICOUNTP,5000)==0)then 
  WRITE(*,*) "AMAXP=" , LINF_P
end if

IF(LINF_P.LT.TOL_P) GOTO 150
IF(ICOUNTP.GT.PRESS_MAX_IT) GOTO 100
GOTO 10
100 WRITE(2,110) LINF_P
110 FORMAT(3X,'PRESSURE CAN NOT CONVERGENCE, MAX. ERROR IS :',F10.6 &
      /3X,'STOP THE CALCULATION')
STOP
150 CONTINUE

END SUBROUTINE