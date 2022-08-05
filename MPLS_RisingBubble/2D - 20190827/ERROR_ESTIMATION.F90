SUBROUTINE ERROR_ESTIMATION()
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA
USE VOF_DATA
IMPLICIT NONE
INTEGER :: I,J
REAL(DP) :: TMPA, TMPB, TMPC


 ! VOL_LS = 0.0_DP
 ! VOL_VOF = 0.0_DP
 
 ! CALL HEAVY_F(PHI)
 
 ! !$OMP PARALLEL DO REDUCTION(+:VOL_LS,VOL_VOF)
 ! DO J = 1, NODE_Y
 ! DO I = 1, NODE_X
 !   VOL_LS = VOL_LS + HEAVY(I,J)
 !   VOL_VOF = VOL_VOF + VOF(I,J)
 ! END DO
 ! END DO
 ! !$OMP END PARALLEL DO
 
 ! EM_LS = (IVOL_LS - VOL_LS)/IVOL_LS
 ! EM_VOF = (IVOL_VOF - VOL_VOF)/IVOL_vof
 
 EI_LS = 0.0_DP
 ! EI_VOF = 0.0_DP
 
 !$OMP PARALLEL DO REDUCTION(+:EI_LS,EI_VOF), PRIVATE(TMPA,TMPB,TMPC)
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
 
   TMPA = max(HEAV(phi(I,J),TMPC),0.5)
   TMPB = max(HEAV(PHI_I(I,J),TMPC),0.5)
   
   !IF(VOF(I,J)>=0.5)TMPA=1.0_DP
   !IF(VOF_I(I,J)>=0.5)TMPB=1.0_DP
   
   ! EI_VOF = EI_VOF + ABS(TMPA-TMPB)
 
   TMPA = MAX(0.0, HEAV(PHI(I,J),TMPC)-0.5)+0.5
   TMPB = MAX(0.0, HEAV(PHI_I(I,J),TMPC)-0.5)+0.5
   
   TMPA = HEAV(PHI(I,J),TMPC)
   TMPB = HEAV(PHI_I(I,J),TMPC)

   EI_LS = EI_LS + ABS(TMPA-TMPB)!/(0.15*2*ACOS(-1.0D0))
   
 END DO
 END DO
 !$OMP END PARALLEL DO
 
 ! EI_VOF = EI_VOF / REAL(NODE_X*NODE_Y,KIND=DP)
 EI_LS = EI_LS / REAL(NODE_X*NODE_Y,KIND=DP)
 
 ! IF( ITER_CNT.EQ.1)THEN
 !   IF( METHOD_CNT.EQ.1 )THEN
 !  	 WRITE(*,*)"Mass Preserving LS"
 !   ELSE IF( METHOD_CNT.EQ.2)THEN
 !    WRITE(*,*)"Classical LS"
 !   ENDIF
 ! ENDIF
 
 WRITE(*,*)"============================================="
 WRITE(*,'(A25,I4,A2,I4)')"Grid number:",NODE_X,"x",NODE_Y
 WRITE(*,'(A,F3.1)')"T=",TIME_TO_STOP
   
 ! goto 100 
 
 ! IF( INTERFACE_METHOD.EQ.0 )THEN
   ! WRITE(*,'(A25,2ES15.4)')"Volume loss:",em_ls
   ! WRITE(*,'(A25,2ES15.4)')"Averaged loss:",em_ls_a
   write(*,'(A25,2ES15.4)')"Interface error:",ei_ls
 ! ELSE IF( INTERFACE_METHOD.EQ.2)THEN
   ! WRITE(*,'(A25,2ES15.4)')"Volume loss:",EM_VOF
   ! WRITE(*,'(A25,2ES15.4)')"Averaged loss:",EM_VOF_A
   ! write(*,'(A25,2ES15.4)')"Interface error:",EI_VOF
 ! END IF
 
 ! 100 WRITE(*,'(A25,2ES15.4)')"CPU time:",cpu_cost
 
END SUBROUTINE