module dat
implicit none
integer :: m, n
real(8) :: dx, dy, dt, gamma
REAL(8),ALLOCATABLE,DIMENSION(:) :: X, Y
REAL(8),ALLOCATABLE,DIMENSION(:,:) :: RHO, IX, IY, P, E, U, V, H, A
REAL(8),ALLOCATABLE,DIMENSION(:,:,:) :: S0, S1, S2, Q
REAL(8),ALLOCATABLE,DIMENSION(:,:,:,:) :: FS, FSR, FSL, FF, FFR, FFL, QR, QL
REAL(8),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: QA, iQA, DF, DS
end module

PROGRAM MAIN
use dat
IMPLICIT NONE
INTEGER :: I, j, ITER, plt, met
REAL(8) :: xstart, xend, ystart, yend, TIME, pltime, time2stop
REAL(8) :: B, R
Character(20) :: name

N = 64
M = 64
N = N + 1
M = M + 1

name = "OCRWENO-LD"
met = 3

XSTART = -5.0d0
XEND =  5.0d0

YSTART = -5.0D0
YEND = 5.0D0

DX = (XEND-XSTART)/REAL(N-1,8)
DY = (YEND-YSTART)/REAL(M-1,8)
DT = 0.01d0*DX

gamma = 1.4d0

allocate( X(N), Y(M) )
ALLOCATE( RHO(-2:N+3,-2:M+3), IX(-2:N+3,-2:M+3), IY(-2:N+3,-2:M+3), P(-2:N+3,-2:M+3), E(-2:N+3,-2:M+3) )
ALLOCATE( U(-2:N+3,-2:M+3), V(-2:N+3,-2:M+3), A(-2:N+3,-2:M+3), H(-2:N+3,-2:M+3) )
ALLOCATE( S0(-2:N+3,-2:M+3,4), S1(-2:N+3,-2:M+3,4), S2(-2:N+3,-2:M+3,4) )
ALLOCATE(  Q(-2:N+3,-2:M+3,4), QR(2,-2:N+3,-2:M+3,4), QL(2,-2:N+3,-2:M+3,4) )
ALLOCATE( FS(2,-2:N+3,-2:M+3,4), FSR(2,-2:N+3,-2:M+3,4), FSL(2,-2:N+3,-2:M+3,4) )
ALLOCATE( FF(2,-2:N+3,-2:M+3,4), FFR(2,-2:N+3,-2:M+3,4), FFL(2,-2:N+3,-2:M+3,4) )
ALLOCATE( QA(2,-2:N+3,-2:M+3,4,4), IQA(2,-2:N+3,-2:M+3,4,4), DF(2,-2:N+3,-2:M+3,4,4), DS(2,-2:N+3,-2:M+3,4,4) )

DO I = 1, N
	X(I) = XSTART + REAL(I-1,8)*DX 
END DO

DO J = 1, M
	Y(J) = XSTART + REAL(J-1,8)*DY
END DO

! Isentropic vortex convection

B = 0.5D0

!$OMP PARALLEL DO PRIVATE(R)
DO J = 1, M
DO I = 1, N
	
	R = DSQRT( X(I)**2.0D0 + Y(J)**2.0D0 )
	
	RHO(I,J) = ( 1.0D0 - (GAMMA-1.0D0)*B**2.0D0*DEXP(1.0-R**2.0D0) / ( 8.0D0*GAMMA*DACOS(-1.0D0)**2.0D0 ) )**(1.0D0/(GAMMA-1.0D0))
	P(I,J) = RHO(I,J)**GAMMA
	U(I,J) = 0.1D0 - B/(2.0*DACOS(-1.0D0))*DEXP(0.5D0*(1.0D0-R**2.0D0))*Y(J)
	V(I,J) = 0.0D0 + B/(2.0*DACOS(-1.0D0))*DEXP(0.5D0*(1.0D0-R**2.0D0))*X(I)
	
	IX(I,J) = RHO(I,J)*U(I,J)
	IY(I,J) = RHO(I,J)*V(I,J)
	E(I,J)  = P(I,J) / (GAMMA-1.0D0) + 0.5D0*RHO(I,J)*(U(I,J)**2.0D0+V(I,J)**2.0D0)
	
END DO
END DO
!$OMP END PARALLEL DO
 
call RENEW
 

TIME = 0.0d0
ITER = 0
time2stop=200.0D0

plt=0
pltime=0.5

call plot(plt,Name)

DO
   
	ITER = ITER + 1
	TIME = TIME + DT
	
	write(*,*)time
   
	!IF( mod(iter,2)==0 )THEN
		CALL SOLVE_1D(0)
	!ELSE
	!	CALL SOLVE_1D(0)
    !ENDIF
   
	if( abs(time-plt*pltime)<dt )then
		call plot(plt,Name)
	endif
   
	if( time>time2stop )exit
   
END DO

 CONTAINS
 
 INCLUDE 'bc.f90'
 include 'plot.f90'
 include 'renew.f90'
 include 'solve_2d.f90'
 include 'tdma.f90'
 include 'scheme.f90'
 include 'scheme_weno.f90'
 
END PROGRAM
