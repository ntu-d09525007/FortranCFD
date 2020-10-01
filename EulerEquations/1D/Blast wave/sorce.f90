SUBROUTINE SORCE(Src,Ud,BTN)
use dat
IMPLICIT NONE
INTEGER :: I, BTN
REAL(8) :: Src(-2:N+3),Ud(-2:N+3)

CALL BC(F)
 
if( btn==0 ) then
    CALL weno_JS(F,FR,FL,N)
else if (btn==1)then
    CALL CRWENO(F,FR,FL,N)
else IF(BTN==2)then
    CALL CRWENO_LD(F,FR,FL,N)
ELSE
   CALL oCRWENO_LD(F,FR,FL,N)
end if

call limiter(f,fr,fl,n)

!$omp parallel do
DO I = 0, N

	IF( SL(I) >= 0.0D0 )THEN
		FH(I) = FL(I)
	ELSE IF ( SR(I) <= 0.0D0 )THEN
		FH(I) = FR(I)
	ELSE
		FH(I) = ( SR(I)*FL(I)-SL(I)*FR(I)+SL(I)*SR(I)*( Ud(i+1)-Ud(I) ) )/ ( SR(I)-SL(I) )
	END IF 
	
 END DO
 !$omp end parallel do
 
 DO I = 1, N
 	Src(I) = -(FH(I)-FH(I-1))/DX
 END DO
 
END SUBROUTINE
