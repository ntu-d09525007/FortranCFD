module dat
implicit none
integer :: n
real(8) :: dx, dt, gamma
REAL(8),ALLOCATABLE,DIMENSION(:) :: X
REAL(8),ALLOCATABLE,DIMENSION(:) :: RHO, IX, P, E, U, H, A
REAL(8),ALLOCATABLE,DIMENSION(:) :: RHO_OLD, IX_OLD, E_OLD, U_OLD, P_OLD, H_OLD, A_OLD
REAL(8),ALLOCATABLE,DIMENSION(:,:) :: S0,S1,S2
REAL(8),ALLOCATABLE,DIMENSION(:,:) :: FS, FSR, FSL, FF, FFR, FFL, QR, QL
real(8),allocatable,dimension(:,:) :: FDD, FH
REAL(8),ALLOCATABLE,DIMENSION(:,:,:) :: Q, IQ, D
end module

PROGRAM MAIN
use dat
IMPLICIT NONE
INTEGER :: I, ITER, plt, met
REAL(8) :: xstart, xend, TIME, pltime, time2stop
Character(20) :: name

N = 200
N = N + 1

name = "OCRWENO-LD"
met = 3


XSTART = 0.0d0
XEND = 1.0d0

DX = (XEND-XSTART)/REAL(N-1,8)
DT = 0.005d0*DX

gamma = 1.4d0

allocate( X(N) )
ALLOCATE( RHO(-2:N+3), IX(-2:N+3), P(-2:N+3), E(-2:N+3), U(-2:N+3), H(-2:N+3), A(-2:N+3) )
ALLOCATE( RHO_OLD(-2:N+3), IX_old(-2:N+3), P_OLD(-2:N+3), E_OLD(-2:N+3), U_OLD(-2:N+3), H_OLD(-2:N+3), A_OLD(-2:N+3) )
ALLOCATE( S0(-2:N+3,3), S1(-2:N+3,3), S2(-2:N+3,3) )
ALLOCATE( FS(-2:N+3,3), FSR(-2:N+3,3), FSL(-2:N+3,3) )
ALLOCATE( FF(-2:N+3,3), FFR(-2:N+3,3), FFL(-2:N+3,3) ) 
allocate( FDD(-2:N+3,3), FH(-2:N+3,3) )
ALLOCATE( QR(-2:N+3,3), QL(-2:N+3,3), D(-2:N+3,3,3) )
ALLOCATE( Q(-2:N+3,3,3), IQ(-2:N+3,3,3) )

DO I = 1, N
	X(I) = XSTART + REAL(I-1,8)*DX 
END DO

! Blast-waves interaction

DO I = 1, N
	IF( X(I) > 0.5D0 )THEN
        RHO(I) = 0.125D0
        U(I) = 0.0D0
        P(I) = 0.1D0
    ELSE
        RHO(I) = 1.0D0
        U(I) = 0.0D0
        P(I) = 1.0D0
    ENDIF
    E(I) = P(I)/(gamma-1.0d0) + 0.5d0*RHO(I)*U(I)**2.0d0
    IX(I) = RHO(I)*U(I)
END DO
 
call RENEW
 
TIME = 0.0d0
ITER = 0
time2stop=0.14d0

plt=0
pltime=time2stop

call plot(plt,X,RHO,N,Name)

DO
   
   ITER = ITER + 1
   TIME = TIME + DT
   
   IF( mod(iter,3)==0 )THEN
	 CALL SOLVE_1D(0)
   ELSE
     CALL SOLVE_1D(met)
   ENDIF
   
   if( abs(time-plt*pltime)<dt )then
       call plot(plt,X,RHO,N,Name)
   endif
   
   if( time>time2stop )exit
   
END DO

 CONTAINS
 
 INCLUDE 'bc.f90'
 include 'plot.f90'
 include 'renew.f90'
 include 'solve_1d.f90'
 include 'tdma.f90'
 include 'scheme.f90'
 include 'scheme_weno.f90'
 
END PROGRAM
