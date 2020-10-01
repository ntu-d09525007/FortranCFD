SUBROUTINE VORTEX_VELOCITY()
USE PRECISION
USE PROBLEM_DEF
USE FLUID_PROPERTIES
IMPLICIT NONE
INTEGER :: I,J,K
REAL(DP) :: XS,YS

 !$OMP PARALLEL DO PRIVATE(XS,YS)
 DO J = 0, NODE_Y
 DO I = 0, NODE_X
 
  xs=0.5*(x(i)+x(i+1))
  ys=0.5*(y(j)+y(j+1))

  u(i,j) =  sin(pi*xs)**2 * sin(2.0_dp*pi*y(j)) * cos(pi*time/TIME_TO_STOP)
  v(i,j) = -sin(pi*ys)**2 * sin(2.0_dp*pi*x(i)) * cos(pi*time/TIME_TO_STOP)

  xs=X(I)
  ys=Y(J)

  uh(i,j) =  sin(pi*xs)**2 * sin(2.0_dp*pi*ys) * cos(pi*time/TIME_TO_STOP)
  vh(i,j) = -sin(pi*ys)**2 * sin(2.0_dp*pi*xs) * cos(pi*time/TIME_TO_STOP)
   
 END DO
 END DO
 !$OMP END PARALLEL DO
 
 CALL BC3D(UH)
 CALL BC3D(VH)
 CALL BC3D(U)
 CALL BC3D(V)
 

END SUBROUTINE

SUBROUTINE double_vortex_velocity()
USE PRECISION
USE PROBLEM_DEF
USE FLUID_PROPERTIES
IMPLICIT NONE
INTEGER :: I,J,K
REAL(DP) :: XS,YS

 !$OMP PARALLEL DO PRIVATE(XS,YS)
 DO J = 0, NODE_Y
 DO I = 0, NODE_X
 
  xs=0.5*(x(i)+x(i+1))
  ys=0.5*(y(j)+y(j+1))

  u(i,j) = dsin(4.0*pi*xs/xl) * dsin(4.0*pi*y(j)/yl) * cos(pi*time/TIME_TO_STOP)
  v(i,j) = dcos(4.0*pi*x(i)/xl) * dcos(4.0*pi*ys/yl) * cos(pi*time/TIME_TO_STOP)

  xs=X(I)
  ys=Y(J)

  uh(i,j) = dsin(4.0*pi*xs/xl) * dsin(4.0*pi*ys/yl) * cos(pi*time/TIME_TO_STOP)
  vh(i,j) = dcos(4.0*pi*xs/xl) * dcos(4.0*pi*ys/yl) * cos(pi*time/TIME_TO_STOP)
   
 END DO
 END DO
 !$OMP END PARALLEL DO
 
 CALL BC3D(UH)
 CALL BC3D(VH)
 CALL BC3D(U)
 CALL BC3D(V)
 

END SUBROUTINE
