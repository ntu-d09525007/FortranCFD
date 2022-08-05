SUBROUTINE PRINT_SCREEN()
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA
USE VOF_DATA
USE IBM_DATA
IMPLICIT NONE

    IF( MOD(iter,5) == 0 ) THEN
       WRITE(*,*) "============================="
       WRITE(*,'(A15,I7)'  ) "Step=",iter 
       WRITE(*,'(A15,F8.5)') "Time=",time 
       WRITE(*,'(A15,ES13.5)') "DIV=",div 
       WRITE(*,'(A15,ES13.5)') "L2_U=",L2_U 
       WRITE(*,'(A15,ES13.5)') "L2_V=",L2_V 
       WRITE(*,'(A15,ES13.5)') "Linf_P=",LINF_P 
       !WRITE(*,'(A15,ES13.5)') "< mass,VOF >",(1.0-MASS_VOF/IMASS_VOF)*100
	   !WRITE(*,'(A15,ES13.5)') "< vol ,VOF >",(1.0-VOL_VOF/IVOL_VOF)*100
       WRITE(*,'(A15,ES13.5)') "< mass,LS >",(1.0-MASS_LS/IMASS_LS)*100
	   WRITE(*,'(A15,ES13.5)') "< vol ,LS >",(1.0-VOL_LS/IVOL_LS)*100
       if( IBM_SOLVER==1 )then
         write(*,*)"----------------------------"
         write(*,'(A10,3A18)')" ","X","Y"
         write(*,'(A10,3ES18.5)') "Position",SOLID_X,SOLID_Y
         write(*,'(A10,3ES18.5)') "Velocity",SOLID_U,SOLID_V
         write(*,*)"----------------------------"
         write(*,'(A10,3ES18.5)') "ANGLE",SOLID_TX,SOLID_TY
         write(*,'(A10,3ES18.5)') "ROTATION",SOLID_WX,SOLID_WY
      end if
      WRITE(*,*) "============================="
    ENDIF 


END SUBROUTINE
