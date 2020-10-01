SUBROUTINE PRINT_SCREEN()
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA
USE VOF_DATA
USE IBM_DATA
USE FLUID_PROPERTIES
IMPLICIT NONE
    
    IF( MOD(iter,5) == 0 ) THEN
       WRITE(*,*) "============================="
       WRITE(*,'(A15,I7)'  ) "Step=",iter 
       WRITE(*,'(A15,F8.5)') "Time=",time 
       WRITE(*,'(A15,ES13.5)') "DIV=",div 
	   
       WRITE(*,'(A15,ES13.5)')"L2_velocity=",max(l2_u,L2_v,l2_w)
	   
       !WRITE(*,'(A15,ES13.5)') "Linf_U=",Linf_U 
       !WRITE(*,'(A15,ES13.5)') "Linf_V=",Linf_V 
       !WRITE(*,'(A15,ES13.5)') "Linf_W=",Linf_W 
	   
	   WRITE(*,'(A15,ES13.5)')"Linf_velocity=",max(linf_u,Linf_v,linf_w)
	   
       WRITE(*,'(A15,ES13.5)')"Linf_P=",LINF_P 
	   
       WRITE(*,'(A15,ES13.5)') "Loss(%) VOF=",(1.0-MASS_VOF/IMASS_VOF)*100
       WRITE(*,'(A15,ES13.5)') " ",(1.0-VOL_VOF/IVOL_VOF)*100
       WRITE(*,'(A15,ES13.5)') "Mass loss(%):",(1.0-MASS_LS/IMASS_LS)*100
       WRITE(*,'(A15,ES13.5)') "Vol. loss(%):",(1.0-VOL_LS/IVOL_LS)*100
       !WRITE(*,'(A15,ES13.5)') "momentum.x :", MASS_u
       !write(*,'(a15,es13.5)') "momentum.y :", mass_v
       !write(*,'(a15,es13.5)') "momentum.z :", mass_w
       if( IBM_SOLVER==1 )then
         write(*,*)"----------------------------"
         write(*,'(A10,3A18)')" ","X","Y","Z"
         write(*,'(A10,3ES18.5)') "Position",SOLID_X,SOLID_Y,SOLID_Z
         write(*,'(A10,3ES18.5)') "Velocity",SOLID_U,SOLID_V,SOLID_W 
         write(*,*)"----------------------------"
         write(*,'(A10,3ES18.5)') "ANGLE",SOLID_TX,SOLID_TY,SOLID_TZ
         write(*,'(A10,3ES18.5)') "ROTATION",SOLID_WX,SOLID_WY,SOLID_WZ 
      end if
      WRITE(*,*) "============================="
    ENDIF 

END SUBROUTINE

