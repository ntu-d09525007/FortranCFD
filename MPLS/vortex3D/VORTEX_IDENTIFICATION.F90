SUBROUTINE VORTEX_IDENTIFICATION()
USE PRECISION
USE PROBLEM_DEF
USE FLUID_PROPERTIES, only : U,V,W,uh,vh,wh
USE VORTEX_DATA
IMPLICIT NONE
INTEGER :: I,J,K
REAL(DP) :: UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,R
REAL(DP) :: L1, L2, L3

!$OMP PARALLEL DO PRIVATE(UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,R,L1,L2,L3)
 DO K = 1, NODE_Z
 DO J = 1, NODE_Y
 DO I = 1, NODE_X
	
	Ux = 0.5*(uh(i+1,j,k)-uh(i-1,j,k))/dx
	Uy = 0.5*(uh(i,j+1,k)-uh(i,j-1,k))/dx
	Uz = 0.5*(uh(i,j,k+1)-uh(i,j,k-1))/dx

	Vx = 0.5*(vh(i+1,j,k)-vh(i-1,j,k))/dx
	Vy = 0.5*(vh(i,j+1,k)-vh(i,j-1,k))/dx
	Vz = 0.5*(vh(i,j,k+1)-vh(i,j,k-1))/dx
	
	Wx = 0.5*(wh(i+1,j,k)-wh(i-1,j,k))/dx
	Wy = 0.5*(wh(i,j+1,k)-wh(i,j-1,k))/dx
	Wz = 0.5*(wh(i,j,k+1)-wh(i,j,k-1))/dx
	
    Call eigenvalue_of_tensor(Ux,Uy,Uz,Vx,Vy,Vz,Wx,Wy,Wz,L1,L2,L3)
  
    R = Ux*Vy*Wz - Ux*Vz*Wy - Uy*Vx*Wz + Uy*Vz*Wx + Uz*Vx*Wy - Uz*Vy*Wx
			
	Q_CRI(I,J,K) = -0.5*(L1+L2+L3)
	
	DELTA_CRI(i,j,k) = (Q_CRI(I,J,K)/3.0)**3.0 + (R/2.0)**2.0
	
	LAM2_CRI(I,J,K) = L1+L2+L3 - MAX(MAX(L1,L2),L3) - MIN(MIN(L1,L2),L3)
	
	Vor_x(i,j,k) = Wy-Vz
	Vor_y(i,j,k) = Uz-Wx
	Vor_z(i,j,k) = Vx-Uy
	
	!Lamb_x(i,j,k) = - Vh(i,j,k)*(Uy - Vx) - Wh(i,j,k)*(Uz - Wx)
	!Lamb_y(i,j,k) =   Uh(i,j,k)*(Uy - Vx) - Wh(i,j,k)*(Vz - Wy)
	!Lamb_z(i,j,k) =   Uh(i,j,k)*(Uz - Wx) + Vh(i,j,k)*(Vz - Wy)
	
	Lamb_x(i,j,k) = VH(I,J,K)*VOR_Z(I,J,K) - WH(I,J,K)*VOR_Y(I,J,K)
	LAMB_Y(I,J,K) = WH(I,J,K)*VOR_X(I,J,K) - UH(I,J,K)*VOR_Z(I,J,K)
	LAMB_Z(I,J,K) = UH(I,J,K)*VOR_Y(I,J,K) - VH(I,J,K)*VOR_X(I,J,K)
	
	FLEX_PROD(I,J,K) =  - UH(I,J,K)*(UH(I+1,J,K)-2.0*UH(I,J,K)+UH(I-1,J,K))/DX**2  &
					    - VH(I,J,K)*(VH(I,J+1,K)-2.0*VH(I,J,K)+VH(I,J-1,K))/DY**2  &
					    - WH(I,J,K)*(WH(I,J,K+1)-2.0*WH(I,J,K)+WH(I,J,K-1))/DZ**2
						       
	ENTROPY(I,J,K) = Vor_x(i,j,k)**2 + Vor_y(i,j,k)**2 + Vor_z(i,j,k)**2
	
	LAMB_DIV(I,J,K) = FLEX_PROD(I,J,K) - ENTROPY(I,J,K) 
	
	HELICITY(I,J,K) = UH(I,J,K)*VOR_X(I,J,K) + VH(I,J,K)*VOR_Y(I,J,K) + WH(I,J,K)*VOR_Z(I,J,K)
		
 ENDDO
 ENDDO
 ENDDO
 !$OMP END PARALLEL DO
 
END SUBROUTINE

SUBROUTINE eigenvalue_of_tensor(Ux,Uy,Uz,Vx,Vy,Vz,Wx,Wy,Wz,L1,L2,L3)
USE PRECISION
implicit none
real(dp) :: Ux,Uy,Uz,Vx,Vy,Vz,Wx,Wy,Wz
real(dp) :: a, b, c, d, e, f, p, q
real(dp) :: cfB, cfC, cfD
real(dp) :: theta, L1, L2, L3

 a = Ux**2 + Uy*Vx + Uz*Wx
 b = (Ux*Uy)/2 + (Ux*Vx)/2 + (Uy*Vy)/2 + (Vx*Vy)/2 + (Uz*Wy)/2 + (Vz*Wx)/2
 c = (Ux*Uz)/2 + (Uy*Vz)/2 + (Ux*Wx)/2 + (Uz*Wz)/2 + (Vx*Wy)/2 + (Wx*Wz)/2
 d = Vy**2 + Uy*Vx + Vz*Wy
 e = (Uz*Vx)/2 + (Uy*Wx)/2 + (Vy*Vz)/2 + (Vy*Wy)/2 + (Vz*Wz)/2 + (Wy*Wz)/2
 f = Wz**2 + Uz*Wx + Vz*Wy
 
 cfB = - (a + d + f)
 cfC = - (b**2 + c**2 + e**2 - a*d - a*f - d*f)
 cfD = f*b**2 - 2*b*c*e + d*c**2 + a*e**2 - a*d*f
 
 p = cfC - cfB**2/3.0
 q = cfD - cfB*cfC/3.0 + 2.0*cfB**3/27.0
 
 theta = dacos(-0.5*q / ( dsqrt((-p/3.0)**3) + eps ))
 
 L1 = -cfB/3.0 + dsqrt(-p/3.0)* 2.0*dcos(theta/3.0)
 L2 = -cfB/3.0 - dsqrt(-p/3.0)*( dcos(theta/3.0) + dsqrt(3.0_8)*dsin(theta/3.0) )
 L3 = -cfB/3.0 - dsqrt(-p/3.0)*( dcos(theta/3.0) - dsqrt(3.0_8)*dsin(theta/3.0) )
 
end subroutine
