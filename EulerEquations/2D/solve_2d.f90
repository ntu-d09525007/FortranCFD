SUBROUTINE SOLVE_1D(BTN)
use dat
IMPLICIT NONE
INTEGER :: BTN, I, J

	CALL CHARA_PART(S0,0)
   
	!$OMP PARALLEL DO 
	DO J = 1, M
	DO I = 1, N
		RHO(I,J) = RHO(I,J) + S0(I,J,1)*DT
		IX(I,J) = IX(I,J) + S0(I,J,2)*DT
		IY(I,J) = IY(I,J) + S0(I,J,3)*DT
		E(I,J) = E(I,J) + S0(I,J,4)*DT
	END DO
	END DO
	!$OMP END PARALLEL DO 
	
	CALL RENEW()

	CALL CHARA_PART(S1,0)
   
	!$OMP PARALLEL DO 
	DO J = 1, M
	DO I = 1, N
		RHO(I,J) = RHO(I,J) + DT * (-3.0D0*S0(I,J,1)+S1(I,J,1))/4.0D0
		IX(I,J) = IX(I,J) + DT * (-3.0D0*S0(I,J,2)+S1(I,J,2))/4.0D0
		IY(I,J) = IY(I,J) + DT * (-3.0D0*S0(I,J,3)+S1(I,J,3))/4.0D0
		E(I,J) = E(I,J) + DT * (-3.0D0*S0(I,J,4)+S1(I,J,4))/4.0D0
	END DO
	END DO
	!$OMP END PARALLEL DO 
	
	CALL RENEW()
	
	CALL CHARA_PART(S2,0)
   
	!$OMP PARALLEL DO 
	DO J = 1, M
	DO I = 1, N
		RHO(I,J) = RHO(I,J) + DT * ( -S0(I,J,1)-S1(I,J,1)+8.0D0*S2(I,J,1) )/12.0D0
		IX(I,J) = IX(I,J) + DT * ( -S0(I,J,2)-S1(I,J,2)+8.0D0*S2(I,J,2) )/12.0D0
		IY(I,J) = IY(I,J) + DT * ( -S0(I,J,3)-S1(I,J,3)+8.0D0*S2(I,J,3) )/12.0D0
		E(I,J) = E(I,J) + DT * ( -S0(I,J,4)-S1(I,J,4)+8.0D0*S2(I,J,4) )/12.0D0
	END DO
	END DO
	!$OMP END PARALLEL DO 
	
	CALL RENEW()
	
END SUBROUTINE

SUBROUTINE JACOB_DECOM
USE DAT
IMPLICIT NONE
INTEGER :: I,J,ii,jj
REAL(8),DIMENSION(4) :: R1,R2,R3,R4,L1,L2,L3,L4

!$OMP PARALLEL DO PRIVATE(R1,R2,R3,R4,L1,L2,L3,L4)
DO J = -2, M+3
DO I = -2, N+3

	CALL EIGEN_SYSTEM(R1,R2,R3,R4,L1,L2,L3,L4,1.0D0,0.0D0,A(I,J),U(I,J),V(I,J),H(I,J),DF(1,I,J,:,:),DS(1,I,J,:,:))	
	
	 QA(1,I,J,:,1) = R1;  QA(1,I,J,:,2) = R2;  QA(1,I,J,:,3) = R3;  QA(1,I,J,:,4) = R4
	IQA(1,I,J,1,:) = L1; IQA(1,I,J,2,:) = L2; IQA(1,I,J,3,:) = L3; IQA(1,I,J,4,:) = L4
	
	CALL EIGEN_SYSTEM(R1,R2,R3,R4,L1,L2,L3,L4,0.0D0,1.0D0,A(I,J),U(I,J),V(I,J),H(I,J),DF(2,I,J,:,:),DS(2,I,J,:,:))	
	
	 QA(2,I,J,:,1) = R1;  QA(2,I,J,:,2) = R2;  QA(2,I,J,:,3) = R3;  QA(2,I,J,:,4) = R4
	IQA(2,I,J,1,:) = L1; IQA(2,I,J,2,:) = L2; IQA(2,I,J,3,:) = L3; IQA(2,I,J,4,:) = L4
	
END DO
END DO
!$OMP END PARALLEL DO

END SUBROUTINE

SUBROUTINE EIGEN_SYSTEM(R1,R2,R3,R4,L1,L2,L3,L4,NX,NY,A,U,V,H,DF,DS)
USE DAT, ONLY : GAMMA
IMPLICIT NONE
REAL(8) :: NX,NY
REAL(8),DIMENSION(4),INTENT(OUT) :: R1,R2,R3,R4,L1,L2,L3,L4
REAL(8),DIMENSION(4,4),INTENT(OUT) :: DF, DS
REAL(8),INTENT(IN) :: A,U,V,H
REAL(8) :: VN,EK,AA,A2,G

VN = U*NX+V*NY
EK = 0.5D0*(U*U+V*V)
AA=A*A
A2=A*A*2.0D0
G=GAMMA-1.0D0

R1 = (/ 1.0D0, U-A*NX, V-A*NY,    H-A*VN /)
R2 = (/ 1.0D0,      U,      V,        EK /)
R3 = (/ 1.0D0, U+A*NX, V+A*NY,    H+A*VN /)
R4 = (/ 0.0D0,     NY,    -NX, U*NY-V*NX /)

L1 = (/ (G*EK+A*VN)/A2, (-G*U-A*NX)/A2,  (-G*V-A*NY)/A2,  G/A2 /)
L2 = (/   (AA-G*EK)/AA,         G*U/AA,          G*V/AA, -G/AA /)
L3 = (/ (G*EK-A*VN)/A2, (-G*U+A*NX)/A2,  (-G*V+A*NY)/A2,  G/A2 /)
L4 = (/   (V-VN*NY)/NX,             NY,             -NX, 0.0D0 /)

DF(:,1) = (/  VN-A, 0.0D0, 0.0D0, 0.0D0 /)
DF(:,2) = (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /)
DF(:,3) = (/ 0.0D0, 0.0D0,  VN+A, 0.0D0 /)
DF(:,4) = (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /)

DS(:,1) = (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /)
DS(:,2) = (/ 0.0D0,    VN, 0.0D0, 0.0D0 /)
DS(:,3) = (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /)
DS(:,4) = (/ 0.0D0, 0.0D0, 0.0D0,   VN  /)
	
END SUBROUTINE

SUBROUTINE CHARA_PART(S,BTN)
USE DAT
IMPLICIT NONE
INTEGER :: I,J,K,BTN
REAL(8) :: VV1,VV2,S(-2:N+3,-2:M+3,4)

CALL JACOB_DECOM

!$OMP PARALLEL DO 
DO J = 1, M
DO I = 1, N
	Q(I,J,:) = (/ RHO(I,J), IX(I,J), IY(I,J), E(I,J) /)
	DO K = 1, 2
		FS(K,I,J,:) = MATMUL(MATMUL(MATMUL(QA(K,I,J,:,:),DS(K,I,J,:,:)),IQA(K,I,J,:,:)),Q(I,J,:))
		FF(K,I,J,:) = MATMUL(MATMUL(MATMUL(QA(K,I,J,:,:),DF(K,I,J,:,:)),IQA(K,I,J,:,:)),Q(I,J,:))
	END DO
END DO
END DO
!$OMP END PARALLEL DO 

CALL FIND_RIGHT_LEFT

!$omp parallel do PRIVATE(VV1,vv2)
do j = 0, M
do i = 0, N
	
	VV1 = max( abs(u(i+1,J))+a(i+1,J) , abs(u(i,J))+a(i,J) ) 
	VV2 = MAX( ABS(U(I+1,J)), ABS(U(I,J)) )
	
	DF(1,i,j,:,1) = (/   VV1, 0.0D0, 0.0D0, 0.0D0 /)
	DF(1,i,j,:,2) = (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /)
	DF(1,i,j,:,3) = (/ 0.0D0, 0.0D0,   VV1, 0.0D0 /)
	DF(1,i,j,:,4) = (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /)

	DS(1,i,j,:,1) = (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /)
	DS(1,i,j,:,2) = (/ 0.0D0,   VV2, 0.0D0, 0.0D0 /)
	DS(1,i,j,:,3) = (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /)
	DS(1,i,j,:,4) = (/ 0.0D0, 0.0D0, 0.0D0,   VV2 /)
	
	FF(1,I,J,:) = 0.5D0*( FFR(1,I,J,:)+FFL(1,I,J,:) - MATMUL(MATMUL(MATMUL(QA(1,I,J,:,:),DF(1,I,J,:,:)),IQA(1,I,J,:,:)),QR(1,I,J,:)-QL(1,I,J,:)) )
	FS(1,I,J,:) = 0.5D0*( FSR(1,I,J,:)+FSL(1,I,J,:) - MATMUL(MATMUL(MATMUL(QA(1,I,J,:,:),DS(1,I,J,:,:)),IQA(1,I,J,:,:)),QR(1,I,J,:)-QL(1,I,J,:)) )

	VV1 = max( abs(V(i,J+1))+A(i,J+1) , abs(V(i,J))+A(i,J) )
	VV2 = MAX( ABS(V(I,J+1)), ABS(V(I,J)) )
	
	DF(2,i,j,:,1) = (/   VV1, 0.0D0, 0.0D0, 0.0D0 /)
	DF(2,i,j,:,2) = (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /)
	DF(2,i,j,:,3) = (/ 0.0D0, 0.0D0,   VV1, 0.0D0 /)
	DF(2,i,j,:,4) = (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /)

	DS(2,i,j,:,1) = (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /)
	DS(2,i,j,:,2) = (/ 0.0D0,   VV2, 0.0D0, 0.0D0 /)
	DS(2,i,j,:,3) = (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0 /)
	DS(2,i,j,:,4) = (/ 0.0D0, 0.0D0, 0.0D0,   VV2 /)
	
	FF(2,I,J,:) = 0.5D0*( FFR(2,I,J,:)+FFL(2,I,J,:) - MATMUL(MATMUL(MATMUL(QA(2,I,J,:,:),DF(2,I,J,:,:)),IQA(2,I,J,:,:)),QR(2,I,J,:)-QL(2,I,J,:)) )
	FS(2,I,J,:) = 0.5D0*( FSR(2,I,J,:)+FSL(2,I,J,:) - MATMUL(MATMUL(MATMUL(QA(2,I,J,:,:),DS(2,I,J,:,:)),IQA(2,I,J,:,:)),QR(2,I,J,:)-QL(2,I,J,:)) )
	
end do
end do
!$omp end parallel do 

!$omp parallel do 
do j = 1, M
do i = 1, N
	s(i,j,:) =  - ( ff(1,i,j,:)-ff(1,i-1,j,:) ) /dx  - ( fs(1,i,j,:)-fs(1,i-1,j,:) ) /dx
	s(i,j,:) = s(i,j,:) - ( ff(2,i,j,:)-ff(2,i,j-1,:) ) /dy - ( fs(2,i,j,:)-fs(2,i,j-1,:) ) /dy
end do
end do 
!$omp end parallel do 


END SUBROUTINE

SUBROUTINE FIND_RIGHT_LEFT
USE DAT
IMPLICIT NONE
INTEGER :: I,J,K

!$OMP PARALLEL DO
DO I = 1, 4
DO J = 1, 2
	CALL BC(FS(J,:,:,I))
	CALL BC(FF(J,:,:,I))
END DO
CALL BC(Q(:,:,I))
END DO
!$OMP END PARALLEL DO 

DO K = 1, 4
DO J = 0, M
	CALL WENO_JS(FS(1,:,J,K),FSR(1,:,J,K),FSL(1,:,J,K),N)
	CALL WENO_JS(FF(1,:,J,K),FFR(1,:,J,K),FFL(1,:,J,K),N)
	CALL WENO_JS(   Q(:,J,K), QR(1,:,J,K), QL(1,:,J,K),N)
END DO 
END DO

DO K = 1, 4
DO I = 0, N
	CALL WENO_JS(FS(2,I,:,K),FSR(2,I,:,K),FSL(2,I,:,K),M)
	CALL WENO_JS(FF(2,I,:,K),FFR(2,I,:,K),FFL(2,I,:,K),M)
	CALL WENO_JS(   Q(I,:,K), QR(2,I,:,K), QL(2,I,:,K),M)
END DO 
END DO

END SUBROUTINE





