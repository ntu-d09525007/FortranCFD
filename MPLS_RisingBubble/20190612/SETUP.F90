SUBROUTINE SETUP()
USE PRECISION
USE PROBLEM_DEF
USE FLUID_PROPERTIES
USE LS_DATA
USE VOF_DATA
USE IBM_DATA
IMPLICIT NONE
INTEGER :: I,J,K,II,JJ,KK
REAL(DP) :: XS,YS,ZS
REAL(DP) :: IMF, B, DEL, DELX
REAL(DP) :: OH, height, tmp

 Xstart =  0.0_DP
 Xend   =  4.0_DP
 
 Ystart =  0.0_DP
 Yend   =  4.0_DP  
 
 Zstart =  0.0_DP
 Zend   =  8.0_DP

XL = XEND-XSTART
YL = YEND-YSTART
ZL = ZEND-ZSTART

UNIT_GRID = 35
UNIT_GRID2= 20

NODE_X = UNIT_GRID * (XEND-XSTART)
NODE_Y = UNIT_GRID * (YEND-YSTART)
NODE_Z = UNIT_GRID * (ZEND-ZSTART)

DX = 1.0_DP / REAL(UNIT_GRID,KIND=DP)
DY = DX
DZ = DX

DT = 0.005_DP*DX
RDT = 0.5_DP*DX
INTERFACE_WIDTH = 1.5*DX

PRESS_OMEGA = 1.5_DP
PRESS_MAX_IT = 5000000

REC_MASS = 0

! Tolerance of convervence

TOL_P = 1.0E-7_DP
TOL_M = 1.0E-5_DP

! Enviornment Parameters

RHO_L = 995.6_DP
RHO_G = 1.177_DP

MU_L = 0.7972E-3
MU_G = 1.846E-5

GRAVITY = 9.81
SURFACE_TENSION = 0.072_DP

! Characteristic Value

 CHAR_LENGTH = 0.05715
 CHAR_VELOCITY = DSQRT(GRAVITY*CHAR_LENGTH)
 CHAR_TIME = CHAR_LENGTH / CHAR_VELOCITY

 TIME_TO_STOP = 1.0
 TIME_TO_PLOT = 0.1

! Dimensionless Parameters

 RE = RHO_L * CHAR_VELOCITY * CHAR_LENGTH / MU_L
 WE = RHO_L * CHAR_VELOCITY**2 * CHAR_LENGTH / SURFACE_TENSION
 FR = CHAR_VELOCITY / DSQRT(GRAVITY*CHAR_LENGTH)
 RATIO_RHO = RHO_G / RHO_L
 RATIO_AMU = MU_G / MU_L

 RATIO_RHO = 0.001
 RATIO_AMU = 0.01
 Re   = 67.27_DP
 Fr   = 1.0_dp  
 WE   = 16.0_DP
 
! Numerical Setup

!===================
! VOF SOLVER
!
! 1 - THINC/WLIC(LS)
! 2 - THINC/WLIC(VOF)
! 3 - THINC
! 4 - THINC-SW
!
!====================

INTERFACE_METHOD = 0
IBM_SOLVER = 0
VEL_BC = 0
ST_FORCE = 1
G_FORCE = 1

LS_SOLVER = 10
VOF_SOLVER = 1

SOLID_MOTION = 0
VELOCITY_INTERPOLATION = 1

 CALL TO_FILE
 
! Data Initialization

CALL ALLOCATE_ARRAY

X(0) = XSTART
X(NODE_X+1) = XEND

Y(0) = YSTART
Y(NODE_Y+1) = YEND

Z(0) = ZSTART
Z(NODE_Z+1) = ZEND

!$OMP PARALLEL DO
DO I = 1, NODE_X
 X(I) = XSTART + (REAL(I,KIND=DP)-0.5_DP)*DX
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
DO J = 1, NODE_Y
 Y(J) = YSTART + (REAL(J,KIND=DP)-0.5_DP)*DY
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
DO K = 1, NODE_Z
 Z(K) = ZSTART + (REAL(K,KIND=DP)-0.5_DP)*DZ
END DO
!$OMP END PARALLEL DO

 CALL READ_DATA
 IF(READ_OPER==1)GOTO 10

 HEIGHT = 11.0
 
!WRITE(*,*)'H=?'
!READ(*,*)HEIGHT

!$OMP PARALLEL DO PRIVATE(XS,YS,ZS,II,JJ,KK)
DO K = 1, NODE_Z
DO J = 1, NODE_Y
DO I = 1, NODE_X

  U(I,J,K) = 0.0_DP
  V(I,J,K) = 0.0_DP
  W(I,J,K) = 0.0_DP
  P(I,J,K) = 0.0_DP
  VOF(I,J,K) = 0.0_DP
  VVOF(I,J,K) = 0.0_DP
  SOLID_VOF(I,J,K) = 0.0

  DO II = 1, UNIT_GRID2
  DO JJ = 1, UNIT_GRID2
  DO KK = 1, UNIT_GRID2
  
    XS = 0.5_DP*(X(I-1)+X(I)) + REAL(II,KIND=DP)*DX/REAL(UNIT_GRID2,KIND=DP)
    YS = 0.5_DP*(Y(J-1)+Y(J)) + REAL(JJ,KIND=DP)*DY/REAL(UNIT_GRID2,KIND=DP) 
    ZS = 0.5_DP*(Z(K-1)+Z(K)) + REAL(KK,KIND=DP)*DZ/REAL(UNIT_GRID2,KIND=DP)

     if( sqrt((xs-2.0)**2+(ys-2.0)**2+(zs-1.0)**2)-0.5 < 0.0 )then
          vof(i,j,k)=vof(i,j,k)
     else if( sqrt((xs-2.0)**2+(ys-2.0)**2+(zs-2.5)**2)-0.5 < 0.0 )then
          vof(i,j,k)=vof(i,j,k)
     else
        vof(i,j,k)=vof(i,j,k)+1.0                       
     end if 


  END DO
  END DO
  END DO

  VOF(I,J,K)=VOF(I,J,K)/REAL(UNIT_GRID2**3,KIND=DP)
  
  XS = X(I)
  YS = Y(J)
  ZS = Z(K)

   PHI1(i,j,k) = -sqrt((xs-2.0)**2+(ys-2.0)**2+(zs-1.0)**2)+0.5
   phi2(i,j,k) = -sqrt((xs-2.0)**2+(ys-2.0)**2+(zs-2.5)**2)+0.5

   phi(i,j,k) = max(phi1(i,j,k), phi2(i,j,k))
  
END DO
END DO
END DO
!$OMP END PARALLEL DO

 CALL BC3D(PHI)
 CALL BC3D(VOF)
 CALL BC3D(U)
 CALL BC3D(V)
 CALL BC3D(W)
 
 CALL LS_REDISTANCE(PHI,1)
 
 PLOT_INT=0
 TIME=ITER*DT

 CALL PLOT

10 CALL AMURHO
 CALL STORE_DATA 

 CALL SOLVE_NS
 CALL STORE_DATA
 CALL FIND_MASS

 CALL CPU_START_REC

END SUBROUTINE

