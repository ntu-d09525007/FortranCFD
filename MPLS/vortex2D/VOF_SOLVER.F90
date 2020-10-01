subroutine SOLVE_VOF()
USE PRECISION
USE PROBLEM_DEF
USE VOF_DATA
IMPLICIT NONE
 
 if( INTERFACE_METHOD==0 )RETURN 
 
 CALL WLIC_COF
 
 CALL SPLIT_X
 CALL SPLIT_Y
  
end subroutine

SUBROUTINE WLIC_COF
USE PRECISION
USE PROBLEM_DEF
USE LS_DATA
USE VOF_DATA
IMPLICIT NONE
INTEGER :: I,J
REAL(DP) :: A,NX,NY

!$OMP PARALLEL DO
DO J = -2, NODE_Y+3
DO I = -2, NODE_X+3
  VOF_OLD(I,J) = VOF(I,J)
ENDDO
ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO 
DO J = 0, NODE_Y
DO I = 0, NODE_X

  NORMAL_X(I,J) = 0.5*(PHI(I+1,J)-PHI(I-1,J))/DX
  NORMAL_Y(I,J) = 0.5*(PHI(I,J+1)-PHI(I,J-1))/DY
 		 
  A = DSQRT( NORMAL_X(I,J)**2 + NORMAL_Y(I,J)**2 )
	
  VOF_SX(I,J) = ABS(NORMAL_X(I,J)) / A
  VOF_SY(I,J) = ABS(NORMAL_Y(I,J)) / A
  
  NORMAL_X(I,J) = NORMAL_X(I,J) / A
  NORMAL_Y(I,J) = NORMAL_Y(I,J) / A
  
END DO
END DO
!$OMP END PARALLEL DO

END SUBROUTINE

subroutine SPLIT_X
USE PRECISION
USE PROBLEM_DEF
USE FLUID_PROPERTIES
USE LS_DATA
USE VOF_DATA
USE DUMMY_DATA
implicit none
integer :: i,j,k,ii,ib,BTN
real(kind=dp) :: a1,a3,xc,beta,alpha,isgn,a4,a5
 
 !$omp parallel do private( a1,a3,a4,a5,xc,alpha,isgn,ib,ii,beta )
 do j = 0, NODE_Y
 do i = 0, NODE_X
 	
 	if( u(i,j)>0.0 )then
 		ii = i
 		isgn = 1.0
 	else
 	  ii = i+1
 	  isgn = 0.0
 	end if
 	
 	if( VOF_SOLVER==0 )then
 		beta = 2.3_dp
 	else if( VOF_SOLVER==1 .or. VOF_SOLVER==2  )then
 	  beta = 2.3*ABS(NORMAL_X(ii,j)) + 0.01 
  end if


 	if( VOF(ii,j) > 1.0-ePS .or. VOF(ii,j) < ePS ) then
 		
 		F(i,j) = VOF(ii,j) * u(i,j) * dt
 	
 	else
 	  
 	  ib = max(1,ii-1)
 	  
 	  if( VOF(ib,j) < VOF(ii+1,j) )then
 	  	alpha = 1.0
 	  else
 	  	alpha = -1.0
 	  end if
 	  
 	a1 = Exp( beta * ( 2.0 * VOF(ii,j) - 1.0 ) / alpha )
    a3 = Exp(beta)
    xc = 0.5 / beta * Log( ( a3 * a3 - a1 * a3 ) / ( a1 * a3 - 1.0 ) )
    a4 = Cosh( beta * ( isgn - U(i,j) * dt / dx - xc ) )
    a5 = Cosh( beta * ( isgn - xc ) )
 	  
    F(i,j) = 0.5 * ( U(i,j) * dt - alpha * Dx / beta * Log( a4 / a5 ) )

    if(VOF_SOLVER==0 .or. VOF_SOLVER==2) F(i,j) = F(i,j) * VOF_SX(ii,j) + VOF(ii,j) * U(i,j) * dt * (1.0 - VOF_Sx(ii,j))
    
  end if 
 	 	
 end do
 end do
 !$omp end parallel do
 
 !$omp parallel do 
 do j = 1, NODE_Y
 do i = 1, NODE_X
    VOF(i,j) = VOF(i,j) - ( F(i,j) - F(i-1,j) ) / Dx + VOF_OLD(I,J)*DT*(U(I,J)-U(I-1,J))/DX 
 End Do
 end do
 !$omp end parallel do
  
 
end subroutine

subroutine SPLIT_Y
USE PRECISION
USE PROBLEM_DEF
USE FLUID_PROPERTIES
USE LS_DATA
USE VOF_DATA
USE DUMMY_DATA
implicit none
integer :: i,j,ii,ib,BTN
real(kind=dp) :: a1,a3,xc,beta,alpha,isgn,a4,a5

 
 !$omp parallel do private( a1,a3,a4,a5,xc,alpha,isgn,ib,ii,beta )
 do j = 0, NODE_Y
 do i = 0, NODE_X
 	
 	if( V(i,j)>0.0 )then
 		ii = j
 		isgn = 1.0
 	else
 	  ii = j+1
 	  isgn = 0.0
 	end if
 	
 	
 	if( VOF_SOLVER==0 )then
 		beta = 2.3_dp
 	else if( VOF_SOLVER==1 .or. VOF_SOLVER==2  )then
 	  beta = 2.3*ABS(NORMAL_Y(i,II)) + 0.01 
  end if

 	if( VOF(i,ii) > 1.0-EPS .or. VOF(i,ii) < EPS ) then
 		
 		G(i,j) = VOF(i,ii) * V(i,j) * dt
 	
 	else
 	  
 	  ib = max(1,ii-1)
 	  
 	  if( VOF(i,ib) < VOF(i,ii+1) )then
 	  	alpha = 1.0
 	  else
 	  	alpha = -1.0
 	  end if
 	  
 	  a1 = Exp( beta * ( 2.0 * VOF(i,ii) - 1.0 ) / alpha )
    a3 = Exp(beta)
    xc = 0.5 / beta * Log( ( a3 * a3 - a1 * a3 ) / ( a1 * a3 - 1.0 ) )
    a4 = Cosh( beta * ( isgn - V(i,j) * dt / dx - xc ) )
    a5 = Cosh( beta * ( isgn - xc ) )
 	  
    G(i,j) = 0.5 * ( V(i,j) * dt - alpha * Dx / beta * Log( a4 / a5 ) )

    if(VOF_SOLVER==0 .or. VOF_SOLVER==2) G(i,j) = G(i,j) * VOF_SY(i,ii) + VOF(i,ii) * V(i,j) * dt * (1.0 - VOF_SY(i,ii))
    
  end if 

 end do
 end do
 !$omp end parallel do
 
 !$omp parallel do 
 do j = 1, NODE_Y
 do i = 1, NODE_X
    VOF(i,j) = VOF(i,j) - ( G(i,j) - G(i,j-1) ) / Dx + VOF_OLD(I,J)*DT*(V(I,J)-V(I,J-1))/DY
	  VOF(I,J) = MIN(MAX(0.0,VOF(I,J)),1.0)
 End Do
 end do
 !$omp end parallel do
 
 CALL BC3D(VOF)
 

end subroutine
