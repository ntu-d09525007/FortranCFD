SUBROUTINE SOLVE_1D(BTN)
use dat
IMPLICIT NONE
INTEGER :: BTN, i
	
	CALL FAST_SLOW_CHARA(S0,BTN)
	
	!$omp parallel do 
	do i = 1, n
		rho(i) = rho(i) + s0(i,1)*dt
		ix(i) = ix(i) + s0(i,2)*dt
		e(i) = e(i) + s0(i,3)*dt
	end do
	!$omp end parallel do 
	
	CALL RENEW()

	CALL FAST_SLOW_CHARA(S1,BTN)
	
	!$omp parallel do 
	do i = 1, n
		rho(i) = rho(i) + DT/4.0*(-3.0*S0(i,1)+S1(i,1))
		ix(i) = ix(i) + DT/4.0*(-3.0*S0(i,2)+S1(i,2))
		e(i) = e(i) + DT/4.0*(-3.0*S0(i,3)+S1(i,3))
	end do
	!$omp end parallel do 
	
	CALL RENEW()
	
	CALL FAST_SLOW_CHARA(S2,BTN)
	
	!$omp parallel do 
	do i = 1, n
		rho(i) = rho(i) + DT/12.0*(-S0(i,1)-S1(i,1)+8.0*S2(i,1))
		ix(i) = ix(i) + DT/12.0*(-S0(i,2)-S1(i,2)+8.0*S2(i,2))
		e(i) = e(i) + DT/12.0*(-S0(i,3)-S1(i,3)+8.0*S2(i,3))
	end do
	!$omp end parallel do 
	
	CALL RENEW()
   
END SUBROUTINE

SUBROUTINE JACOB_DECOM
USE DAT
IMPLICIT NONE
INTEGER :: I
REAL(8) :: CON

!$OMP PARALLEL DO
DO I = -2, N+3
 
	Q(I,1,1) = 1.0D0
	Q(I,1,2) = RHO(I) / (2.0D0*A(I))
	Q(I,1,3) = -RHO(I) / (2.0D0*A(I))
	
	Q(I,2,1) = U(I)
	Q(I,2,2) = RHO(I)/(2.0D0*A(I)) * (U(I)+A(I))
	Q(I,2,3) = -RHO(I)/(2.0D0*A(I)) * (U(I)-A(I))
	
	Q(I,3,1) = U(I)**2.0D0/2.0D0 
	Q(I,3,2) = RHO(I)/(2.0D0*A(i))*( U(I)**2.0D0/2.0D0 + A(I)**2.0D0/(GAMMA-1.0D0) + A(I)*U(I) )
	Q(I,3,3) = -RHO(I)/(2.0D0*A(I))*( U(I)**2.0D0/2.0D0 + A(I)**2.0D0/(GAMMA-1.0D0) - A(I)*U(I) )
	
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(CON)
DO I = -2, N+3

	CON = (GAMMA-1.0D0)/(RHO(I)*A(I))
	
	IQ(I,1,1) = CON*RHO(I)/A(I)*(-U(I)**2.0D0/2.0D0+A(I)**2.0D0/(GAMMA-1.0D0))
	IQ(I,1,2) = CON*RHO(I)*U(I)/A(I)
	IQ(I,1,3) = -CON*RHO(I)/A(I)
	
	IQ(I,2,1) = CON*( U(I)**2.0D0 / 2.0D0 - A(I)*U(I)/(GAMMA-1.0D0) )
	IQ(I,2,2) = CON*( -U(I) + A(I) / (GAMMA-1.0D0) )
	IQ(I,2,3) = CON
	
	IQ(I,3,1) = CON*( -U(I)**2.0D0 / 2.0D0 - A(I)*U(I)/(GAMMA-1.0D0) )
	IQ(I,3,2) = CON*( U(I) + A(I) / (GAMMA-1.0D0) )
	IQ(I,3,3) = -CON
	
END DO
!$OMP END PARALLEL DO

END SUBROUTINE

SUBROUTINE FAST_SLOW_CHARA(src,btn)
USE DAT
IMPLICIT NONE
INTEGER :: I,btn
real(8), dimension(-2:N+3,3) :: src

CALL JACOB_DECOM

!$OMP PARALLEL DO
DO I = -2, N+3

	FS(I,1) = (GAMMA-1.0D0)/GAMMA * RHO(I)*U(I)
	FS(I,2) = (GAMMA-1.0D0)/GAMMA * RHO(I)*U(I)**2.0D0
	FS(I,3) = 0.5D0*(GAMMA-1.0D0)/GAMMA* RHO(I)*U(I)**3.0D0
	
	FF(I,1) = RHO(I)*U(I) / GAMMA
	FF(I,2) = RHO(I)*U(I)**2.0D0/GAMMA + P(I)
	FF(I,3) = (E(I)+P(I))*U(I)-0.5D0*(GAMMA-1.0D0)/GAMMA*RHO(I)*U(I)**3.0D0

	
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
DO I = 0, N
	D(I,1,1) = MAX( ABS(U(I)), ABS(U(I+1)) )
	D(I,1,2) = 0.0D0
	D(I,1,3) = 0.0D0
	
	D(I,2,1) = 0.0D0
	D(I,2,2) = MAX( ABS(U(I))+A(I), ABS(U(I+1))+A(I) )
	D(I,2,3) = 0.0D0
	
	D(I,3,1) = 0.0D0
	D(I,3,2) = 0.0D0
	D(I,3,3) = MAX( ABS(U(I))+A(I), ABS(U(I+1))+A(I) )
END DO
!$OMP END PARALLEL DO

call find_numer_flux(btn)
	
!$omp parallel do
do i = 0, N
	FDD(i,:) = matmul(matmul(matmul(Q(i,:,:),D(i,:,:)),IQ(i,:,:)),QR(i,:)-QL(i,:))
	 FH(i,:) = 0.5d0*(FFR(i,:)+FFL(i,:)) + 0.5d0*(FSR(i,:)+FSL(i,:)) - 0.5d0*FDD(i,:)
	FDD(i,:) = FS(i,:) + FF(i,:)
end do
!$omp end parallel do

!CALL limiter(FDD(:,1),FH(:,1),N)
!CALL limiter(FDD(:,2),FH(:,2),N)
!CALL limiter(FDD(:,3),FH(:,3),N)

!$omp parallel do
do i = 1, n
	src(i,:) = -(FH(i,:)-FH(i-1,:))/dx
end do
!$omp end parallel do 

END SUBROUTINE

subroutine limiter(f,fh,n)
implicit none
real(8), dimension(-2:N+3) :: f,fh,r
integer :: i,n
real(8) :: rd, rr

rd=0.75d0

!$omp parallel do
do i = 0, n
    R(i) = ( abs(2.0d0*(f(i+1)-f(i))*f(i)-f(i-1)) + 1.0d-13 )/ ( (f(i+1)-f(i))**2.0d0 + (f(i)-f(i-1))**2.0d0 + 1.0d-13 )
end do
!$omp end parallel do

!$omp parallel do private(rr)
do i = 0, n
    rr = min(r(i),r(i+1))
    rr = 0.5d0 * ( 1.0d0 + dtanh(3.0d0*(rr-rd)/max(rd,abs(rr-rd)))/dtanh(3.0d0) )
    rr = 1.0d0 - rr
    !fh(i) = ( 1.0d0-rr ) * fh(i) + rr * (-f(i-1)+7.0d0*f(i)+7.0d0*f(i+1)-f(i+2))/12.0d0
	fh(i) = ( 1.0d0-rr ) * fh(i) + rr * 0.5d0*( f(i) + f(I+1) )
end do
!$omp end parallel do

end subroutine

subroutine find_numer_flux(btn)
use dat
implicit none
integer :: i, btn

if( btn==0 )then

    CALL weno_js(FS(:,1),FSR(:,1),FSL(:,1),N)
	CALL weno_js(FS(:,2),FSR(:,2),FSL(:,2),N)
	CALL weno_js(FS(:,3),FSR(:,3),FSL(:,3),N)
	
    CALL weno_js(FF(:,1),FFR(:,1),FFL(:,1),N)
	CALL weno_js(FF(:,2),FFR(:,2),FFL(:,2),N)
	CALL weno_js(FF(:,3),FFR(:,3),FFL(:,3),N)

	call weno_js(rho,QR(:,1),QL(:,1),N)
	call weno_js(ix,QR(:,2),QL(:,2),N)
	call weno_js(e,QR(:,3),QL(:,3),N)
	
else if (btn==1) then

	CALL crweno(FS(:,1),FSR(:,1),FSL(:,1),N)
	CALL crweno(FS(:,2),FSR(:,2),FSL(:,2),N)
	CALL crweno(FS(:,3),FSR(:,3),FSL(:,3),N)
	
    CALL crweno(FF(:,1),FFR(:,1),FFL(:,1),N)
	CALL crweno(FF(:,2),FFR(:,2),FFL(:,2),N)
	CALL crweno(FF(:,3),FFR(:,3),FFL(:,3),N)

	call crweno(rho,QR(:,1),QL(:,1),N)
	call crweno(ix,QR(:,2),QL(:,2),N)
	call crweno(e,QR(:,3),QL(:,3),N)

else if (btn==2) then

	CALL crweno_LD(FS(:,1),FSR(:,1),FSL(:,1),N)
	CALL crweno_LD(FS(:,2),FSR(:,2),FSL(:,2),N)
	CALL crweno_LD(FS(:,3),FSR(:,3),FSL(:,3),N)
	
    CALL crweno_LD(FF(:,1),FFR(:,1),FFL(:,1),N)
	CALL crweno_LD(FF(:,2),FFR(:,2),FFL(:,2),N)
	CALL crweno_LD(FF(:,3),FFR(:,3),FFL(:,3),N)

	call crweno_LD(rho,QR(:,1),QL(:,1),N)
	call crweno_LD(ix,QR(:,2),QL(:,2),N)
	call crweno_LD(e,QR(:,3),QL(:,3),N)
	
else if (btn==3) then

    CALL ocrweno_LD(FS(:,1),FSR(:,1),FSL(:,1),N)
	CALL ocrweno_LD(FS(:,2),FSR(:,2),FSL(:,2),N)
	CALL ocrweno_LD(FS(:,3),FSR(:,3),FSL(:,3),N)
	
    CALL ocrweno_LD(FF(:,1),FFR(:,1),FFL(:,1),N)
	CALL ocrweno_LD(FF(:,2),FFR(:,2),FFL(:,2),N)
	CALL ocrweno_LD(FF(:,3),FFR(:,3),FFL(:,3),N)

	call ocrweno_LD(rho,QR(:,1),QL(:,1),N)
	call ocrweno_LD(ix,QR(:,2),QL(:,2),N)
	call ocrweno_LD(e,QR(:,3),QL(:,3),N)

end if


end subroutine

