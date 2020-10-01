SUBROUTINE PLOT(ITER,X,RHO,N,nam)
IMPLICIT NONE
INTEGER :: N, I, ITER
REAL(8) :: X(N),RHO(-2:N+3)
CHARACTER(3) :: NAME
Character(*) :: nam 

 WRITE(NAME,'(I3.3)')ITER
 OPEN( UNIT=11, FILE=trim(nam)//Name//'_.plt' )
 write(11,*)'variables = "x" "rho" '
 
DO I = 1, N
	write(11,*)x(i),rho(i)
END DO
 close(11)
 
 iter=iter+1
 
end subroutine


