SUBROUTINE BC_LS()
USE LS_DATA
IMPLICIT NONE

 CALL BC3D(PHI)

END SUBROUTINE

subroutine bc3d(fi)   
USE PRECISION
USE PROBLEM_DEF
IMPLICIT NONE
INTEGER :: I,J
REAL(DP),DIMENSION(-2:NODE_X+3,-2:NODE_Y+3) :: FI

!$omp parallel do
 DO J = 1, NODE_Y

   fi(0,J)=fi(1,J)
   fi(-1,J)=fi(1,J)
   fi(-2,J)=fi(1,J)
   
   fi(NODE_X+1,J)=fi(NODE_X,J)
   fi(NODE_X+2,J)=fi(NODE_X,J)
   fi(NODE_X+3,J)=fi(NODE_X,J)
ENDDO
!$omp end parallel do

!$omp parallel do
DO I =-2, NODE_X+3
	    
   fi(I,0)=fi(I,1)
   fi(I,-1)=fi(I,1)
   fi(I,-2)=fi(I,1)
   

   fi(I,NODE_Y+1)=fi(I,NODE_Y)
   fi(I,NODE_Y+2)=fi(I,NODE_Y)
   fi(I,NODE_Y+3)=fi(I,NODE_Y)   
     
ENDDO
!$omp end parallel do


end subroutine bc3d

subroutine bcu(fi)
USE PRECISION
USE PROBLEM_DEF
IMPLICIT NONE
INTEGER :: I,J
REAL(DP),DIMENSION(-2:NODE_X+3,-2:NODE_Y+3) :: FI

if(VEL_BC==0) then 

    !$omp parallel do
    DO J=1,NODE_Y
       fi(0,J)=0.0_DP
       fi(-1,J)= -fi(1,J)
       fi(-2,J)= -fi(2,J)
       fi(NODE_X,J)=0.0_DP
       fi(NODE_X+1,j)=-fi(NODE_X-1,J)
       fi(NODE_X+2,J)=-fi(NODE_X-2,J)
       fi(NODE_X+3,J)=-fi(NODE_X-3,J)     
    ENDDO
    !$omp end parallel do
    	  
    !$omp parallel do
    DO I=-2,NODE_X+3

       fi(I,0)=-fi(I,1)
       fi(I,-1)=-fi(I,2)
       fi(I,-2)=-fi(I,3)
       fi(I,NODE_Y+1)=-fi(I,NODE_Y)
       fi(I,NODE_Y+2)=-fi(I,NODE_Y-1)
       fi(I,NODE_Y+3)=-fi(I,NODE_Y-2)     
    ENDDO
    !$omp end parallel do

else
	  	
   !$omp parallel do
   DO J=1,NODE_Y
       fi(0,J)=0.0_DP
       fi(-1,J)=-fi(1,J)
       fi(-2,J)=-fi(2,J)
       fi(NODE_X,J)=0.0_DP
       fi(NODE_X+1,J)=-fi(NODE_X-1,J)
       fi(NODE_X+2,J)=-fi(NODE_X-2,J)
       fi(NODE_X+3,J)=-fi(NODE_X-3,J)     
    ENDDO
   !$omp end parallel do

    !$omp parallel do
    DO I=-2,NODE_X+3
       fi(I,0)=fi(I,1)
       fi(I,-1)=fi(I,2)
       fi(I,-2)=fi(I,3)
       fi(I,NODE_Y+1)=fi(I,NODE_Y)
       fi(I,NODE_Y+2)=fi(I,NODE_Y-1)
       fi(I,NODE_Y+3)=fi(I,NODE_Y-2)     
    ENDDO
	!$omp end parallel do
	  
endif 
	  
end subroutine bcu
	  
subroutine bcv(fi)
USE PRECISION
USE PROBLEM_DEF
IMPLICIT NONE
INTEGER :: I,J
REAL(DP),DIMENSION(-2:NODE_X+3,-2:NODE_Y+3) :: FI


if(VEL_BC==0) then 

 !$omp parallel do
  DO I=1,NODE_X
     fi(I,0)=0.0_DP
     fi(I,-1)=-fi(I,1)
     fi(I,-2)=-fi(I,2)
     fi(i,NODE_Y)=0.0_DP
     fi(I,NODE_Y+1)=-fi(I,NODE_Y-1)
     fi(I,NODE_Y+2)=-fi(I,NODE_Y-2)
     fi(I,NODE_Y+3)=-fi(I,NODE_Y-3)     
  ENDDO
	!$omp end parallel do

  !$omp parallel do
  DO J=-2,NODE_Y+3
     fi(0,J)=-fi(1,J)
     fi(-1,J)=-fi(2,J)
     fi(-2,J)=-fi(3,J)
     fi(NODE_X+1,J)=-fi(NODE_X,J)
     fi(NODE_X+2,J)=-fi(NODE_X-1,J)
     fi(NODE_X+3,J)=-fi(NODE_X-2,J)     
  ENDDO
	!$omp end parallel do

	  else

  !$omp parallel do
   DO I=1,NODE_X
      fi(I,0)=0.0_DP
      fi(I,-1)=-fi(I,1)
      fi(I,-2)=-fi(I,2)
	  fi(i,NODE_Y)=0.0_DP
      fi(I,NODE_Y+1)=-fi(I,NODE_Y-1)
      fi(I,NODE_Y+2)=-fi(I,NODE_Y-2)
      fi(I,NODE_Y+3)=-fi(I,NODE_Y-3)     
   ENDDO
	!$omp end parallel do

  !$omp parallel do 
   DO J=-2,NODE_Y+3
      fi(0,J)=fi(1,J)
      fi(-1,J)=fi(2,J)
      fi(-2,J)=fi(3,J)
      fi(NODE_X+1,J)=fi(NODE_X,J)
      fi(NODE_X+2,J)=fi(NODE_X-1,J)
      fi(NODE_X+3,J)=fi(NODE_X-2,J)     
   ENDDO
	!$omp end parallel do

endif

end subroutine bcv
