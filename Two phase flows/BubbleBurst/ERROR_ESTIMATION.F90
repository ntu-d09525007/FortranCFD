SUBROUTINE ERROR_ESTIMATION()
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA
USE VOF_DATA
IMPLICIT NONE
INTEGER :: I,J,K
REAL(DP) :: TMPA, TMPB, TMPC
 
 !101 IF( GRID_CNT .EQ. 1 )THEN
   WRITE(77,*)"============================================="
   !WRITE(*,'(A,F8.5)')"Density ratio:",RATIO_RHO
   !WRITE(*,*)"-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
   !IF(METHOD_CNT.EQ.2)THEN
   !	WRITE(*,*)"Classical LS"
   !ELSE
   !	WRITE(*,*)"Mass-preserving LS"
   !ENDIF
   !WRITE(*,*)"-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
 !ENDIF
  
  WRITE(77,'(A25,I5,A2,I5,A2,I5)')"Grid number:",NODE_X,"x",NODE_Y,"x",NODE_Z
  !WRITE(*,*)""
 
 !IF( INTERFACE_METHOD.EQ.0  )THEN
   !WRITE(*,'(A25,2ES15.4)')"Volume loss:",em_ls
   WRITE(77,'(A25,2ES15.4)')"Averaged loss:",em_ls_a/time_to_stop
   !write(*,'(A25,2ES15.4)')"Interface error:",ei_ls
 !ELSE 
   !WRITE(*,'(A25,2ES15.4)')"Volume loss:",EM_VOF
   !WRITE(*,'(A25,2ES15.4)')"Averaged loss:",EM_VOF_A
   !write(*,'(A25,2ES15.4)')"Interface error:",EI_VOF
 !END IF
   !WRITE(*,*)""
   WRITE(77,'(A,F10.2)')"CPU TIME:",CPU_COST
   WRITE(77,*)""
   
END SUBROUTINE
