SUBROUTINE ERROR_ESTIMATION()
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA
USE VOF_DATA
IMPLICIT NONE
INTEGER :: I,J
REAL(DP) :: TMPA, TMPB, TMPC

 MASS_LS = 0.0_DP

 CALL HEAVY_F(PHI)
 
 EI_LS = 0.0_DP
 
 !$OMP PARALLEL DO REDUCTION(+:EI_LS,EI_VOF), PRIVATE(TMPA,TMPB,TMPC)
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
 
   TMPA = HEAV(PHI(I,J),TMPC)    
   TMPB = HEAV(PHI_I(I,J),TMPC)  
   
   EI_LS = EI_LS + ABS(TMPA-TMPB)
   
 END DO
 END DO
 !$OMP END PARALLEL DO

 EI_LS = EI_LS / REAL(NODE_X*NODE_Y,KIND=DP)
 

 IF( GRID_CNT .EQ. 1)THEN
  WRITE(*,*)"============================================="
  WRITE(*,'(A,F8.5)')"Density ratio:", RATIO_rho
  write(*,*)""
 ENDIF
 
  WRITE(*,'(A25,I5,A2,I5)')"Grid number:",NODE_X,"x",NODE_Y
  WRITE(*,'(A25,2ES15.4)')"Averaged loss:",em_ls_a
  write(*,'(A25,2ES15.4)')"Maximum loss:",EM_MAX
  write(*,'(A25,2ES15.4)')"Interface error:",ei_ls

 
 WRITE(*,'(A25,2ES15.4)')"CPU time:",cpu_cost
 write(*,*)""
 
END SUBROUTINE