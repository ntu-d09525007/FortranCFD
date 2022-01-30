subroutine SOLVE_LS()
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA

 IF(INTERFACE_METHOD==1)RETURN
 
 CALL LS_RK3_SOLVER

 CALL LS_MAINTAIN
   	  
end subroutine

SUBROUTINE LS_RK3_SOLVER()
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA
USE DUMMY_DATA
IMPLICIT NONE
INTEGER :: I,J

 IF(INTERFACE_METHOD==3)RETURN
 
 !CALL CELL_FACE_VELOCITY()
 CALL VORTEX_VELOCITY()
  
 !$omp parallel do
 DO J = -2, NODE_Y+3
 DO I = -2, NODE_X+3
	PHI_OLD(I,J) = PHI(I,J)
 ENDDO
 enddo
 !$omp end parallel do  
 
 CALL LS_sorce(S0)
 
 !$omp parallel do
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
   PHI(i,j) = PHI(i,j) + dt*s0(i,j)
 ENDDO
 enddo
 !$omp end parallel do  
 
 CALL BC_LS()
 
 CALL LS_sorce(S1)

 !$omp parallel do
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
   PHI(i,j) = PHI(i,j) + dt/4.0_DP*(-3.0_DP*s0(i,j)+s1(i,j))
 ENDDO
 enddo
 !$omp end parallel do  

 CALL BC_LS()
 
 CALL LS_sorce(S2)
 
 !$omp parallel do
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
   PHI(i,j)=PHI(i,j) + dt/12.0_DP*(-s0(i,j)-s1(i,j)+8.0_DP*s2(i,j))
 enddo
 ENDDO
 !$omp end parallel do
 
 CALL BC_LS()
  
END SUBROUTINE

SUBROUTINE CELL_FACE_VELOCITY()
USE PRECISION
USE FLUID_PROPERTIES
USE PROBLEM_DEF
IMPLICIT NONE
INTEGER :: I,J

!$OMP PARALLEL DO
DO J = 1, NODE_Y
DO I = 1, NODE_X
	UH(I,J) = 0.5_DP*(U(I,J)+U(I-1,J))
	VH(I,J) = 0.5_DP*(V(I,J)+V(I,J-1))
END DO
END DO
!$OMP END PARALLEL DO

 CALL BC3D(UH)
 CALL BC3D(VH)

END SUBROUTINE

subroutine LS_sorce(S)
USE PRECISION
USE PROBLEM_DEF
USE DUMMY_DATA
USE FLUID_PROPERTIES
USE LS_DATA
IMPLICIT NONE
INTEGER :: I,J
REAL(DP),DIMENSION(-2:NODE_X+3,-2:NODE_Y+3) :: S
REAL(DP) :: TMP, FX, FY

if( LS_SOLVER < 10 )then

!$OMP PARALLEL DO
DO J = -2, NODE_Y+3
DO I = -2, NODE_X+3
 up(I,J) = 0.5*(uH(I,J)+abs(uH(I,J)))*PHI(I,J)
 um(I,J) = 0.5*(uH(I,J)-abs(uH(I,J)))*PHI(I,J)
 vp(I,J) = 0.5*(vH(I,J)+abs(vH(I,J)))*PHI(I,J)
 vm(I,J) = 0.5*(vH(I,J)-abs(vH(I,J)))*PHI(I,J)              
end do
END DO               
!$OMP END PARALLEL DO
                  
 !$OMP PARALLEL DO
 DO J = 1, NODE_Y 
 	call split_sorce(um(:,j),up(:,j),fp(:,j),fm(:,j),NODE_X,LS_SOLVER) 	
 end do
 !$OMP END PARALLEL DO
 
 !$OMP PARALLEL DO
 DO I = 1, NODE_X
 	call split_sorce(vm(i,:),vp(i,:),gp(i,:),gm(i,:),NODE_Y,LS_SOLVER) 
 END DO
 !$OMP END PARALLEL DO

!$OMP PARALLEL DO
DO J = -2, NODE_Y+3
DO I = -2, NODE_X+3
 F(I,J) = FP(I,J)  +  FM(I,J)
 G(I,J) = GP(I,J)  +  GM(I,J)             
end do
END DO               
!$OMP END PARALLEL DO

 !$omp parallel do 
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
 	 s(i,j) = -(F(i,j)-F(i-1,j))/dx -(G(i,j)-G(i,j-1))/dy
 end do
 end do
 !$omp end parallel do
 
else if( LS_SOLVER==10 )then
	
 !$OMP PARALLEL DO
 DO J = 1, NODE_Y
 	call UCCD(NODE_X,dx,UH(1:NODE_X,j),PHI(1:NODE_X,j),F(1:NODE_X,j))
 end do
 !$OMP END PARALLEL DO
 	
 !$OMP PARALLEL DO
 DO I = 1, NODE_X
 	call UCCD(NODE_Y,dy,VH(i,1:NODE_Y),PHI(i,1:NODE_Y),G(i,1:NODE_Y))
 end do
 !$OMP END PARALLEL DO

 !$omp parallel do
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
 	 s(i,j) = -( F(i,j)*UH(i,j) + G(i,j)*VH(i,j) )	 
 end do
 end do
 !$omp end parallel do
 	
end if
  
end subroutine

subroutine split_sorce(f,g,fp,gm,n,btn)
USE PRECISION
implicit none
integer :: n,btn
real(kind=dp),dimension(-2:n+3) :: f,g,fp,fm,gp,gm

if( btn==0 )then
 call weno_js(f,fp,fm,N)
 call weno_js(g,gp,gm,N)
else if( btn==1 )then
 call weno_z(f,fp,fm,N)
 call weno_z(g,gp,gm,N)	
else if( btn==2 )then
 call weno_m(f,fp,fm,N)
 call weno_m(g,gp,gm,N)
else if( btn==3 )then
 call ocrweno(f,fp,fm,N)
 call ocrweno(g,gp,gm,N) 
else if( btn==4 )then
 call ocrweno_LD(f,fp,fm,N)
 call ocrweno_LD(g,gp,gm,N)
end if

end subroutine

SUBROUTINE MPLS_CORRECTION
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA
IMPLICIT NONE
INTEGER :: I,J,K
REAL(DP) :: NX,NY,FS
REAL(DP) :: MASOLD, VOLOLD, TMP


DO K = 1, 1

 VOL_LS = 0.0
 VOLOLD = 0.0
 FS=0.0
 
 CALL HEAVY_F(PHI) 
 
 !$OMP PARALLEL DO REDUCTION(+:VOL_LS,FS,VOLOLD), PRIVATE(NX,NY)
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
 
   NX = 0.5_DP*(PHI(I+1,J)-PHI(I-1,J))/DX
   NY = 0.5_DP*(PHI(I,J+1)-PHI(I,J-1))/DY
   
   GRAD(I,J) = DSQRT(NX**2+NY**2)
   
   FS = FS + DELTA(I,J)**2*GRAD(I,J)
   
   VOL_LS = VOL_LS + HEAVY(I,J)
   VOLOLD = VOLOLD + HEAV(PHI_OLD(I,J),TMP)
   
 ENDDO
 ENDDO
 !$OMP END PARALLEL DO
 
 !FS=(VOLOLD-VOL_LS)/FS
 FS=(IVOL_LS-VOL_LS)/FS
 
 !$OMP PARALLEL DO 
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
   PHI(I,J) = PHI(I,J) + FS*DELTA(I,J)*GRAD(I,J)
 ENDDO
 ENDDO
 !$OMP END PARALLEL DO 
 
 CALL BC_LS()

ENDDO

END SUBROUTINE

