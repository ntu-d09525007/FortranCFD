module dat
implicit none
integer :: n
real(8) :: dx, dt, gamma
REAL(8),ALLOCATABLE,DIMENSION(:) :: X
REAL(8),ALLOCATABLE,DIMENSION(:) :: RHO, IX, P, E, U, H, A
REAL(8),ALLOCATABLE,DIMENSION(:) :: RHO_OLD, IX_OLD, E_OLD, U_OLD, P_OLD, H_OLD, A_OLD
REAL(8),ALLOCATABLE,DIMENSION(:) :: S0,S1,S2
REAL(8),ALLOCATABLE,DIMENSION(:) :: SR, SL
REAL(8),ALLOCATABLE,DIMENSION(:) :: F, FH, FR, FL, GR, GL, HR, HL
end module

PROGRAM MAIN
use dat
IMPLICIT NONE
INTEGER :: I, ITER, plt, met
REAL(8) :: xstart, xend, TIME, pltime
Character(20) :: name

N = 2000
N = N + 1

name = "Exact_"
met = 0

XSTART = 0.0d0
XEND = 1.0d0

DX = (XEND-XSTART)/REAL(N-1,8)
DT = 0.001d0*DX

gamma = 1.4d0

allocate( X(N) )
ALLOCATE( RHO(-2:N+3), IX(-2:N+3), P(-2:N+3), E(-2:N+3), U(-2:N+3), H(-2:N+3), A(-2:N+3) )
ALLOCATE( RHO_OLD(-2:N+3), IX_old(-2:N+3), P_OLD(-2:N+3), E_OLD(-2:N+3), U_OLD(-2:N+3), H_OLD(-2:N+3), A_OLD(-2:N+3) )
ALLOCATE( S0(-2:N+3), S1(-2:N+3), S2(-2:N+3) )
ALLOCATE( SR(-2:N+3), SL(-2:N+3) )
ALLOCATE( F(-2:N+3), FH(-2:N+3), FR(-2:N+3), FL(-2:N+3), GR(-2:N+3), GL(-2:N+3), HR(-2:N+3), HL(-2:N+3)  )


DO I = 1, N
	X(I) = XSTART + REAL(I-1,8)*DX 
END DO

! Blast-waves interaction

DO I = 1, N

    rho(i) = 1.0d0
    u(i)  = 0.0d0

	IF( X(I) < 0.1D0 )THEN
		P(I) = 1000.0D0
    else if (x(i)>0.9d0 )then
		P(I) = 100.0D0
    else
        P(i) = 0.001d0
	ENDIF


    E(I) = P(I)/(gamma-1.0d0) + 0.5d0*RHO(I)*U(I)**2.0d0
    IX(I) = RHO(I)*U(I)
END DO
 
call RENEW
 
TIME = 0.0;ITER = 0

plt=0
pltime=0.038

 call plot(plt,X,RHO,N,Name)


DO
   
   ITER = ITER + 1
   TIME = TIME + DT
   
   ! 1 -- WENO_JS
   ! 2 -- WENO_Z
   ! 3 -- WENO_M
   ! 4 -- DRPCRWENO_LD
   
   CALL SOLVE_1D(met)

   if( abs(time-plt*pltime)<dt )then
       call plot(plt,X,RHO,N,Name)
   endif
   
   if( time>0.038 )exit
   
END DO

 CONTAINS
 
 INCLUDE 'bc.f90'
 include 'e.f90'
 include 'ix.f90'
 include 'plot.f90'
 include 'renew.f90'
 include 'rho.f90'
 include 'rk3.f90'
 include 'solve_1d.f90'
 include 'sorce.f90'
 include 'tdma.f90'
 include 'scheme.f90'
 include 'scheme_weno.f90'
 
END PROGRAM
