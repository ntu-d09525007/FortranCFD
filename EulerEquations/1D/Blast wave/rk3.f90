 SUBROUTINE RK3_1(Fq)
 use dat
 IMPLICIT NONE
 INTEGER :: I
 REAL(8) :: Fq(-2:N+3)
 
 !$omp parallel do
 DO I = 1, N
 	Fq(I) = Fq(I) + DT*S0(I)
 END DO
 !$omp end parallel do
 
 END SUBROUTINE
 
 SUBROUTINE RK3_2(Fq)
 use dat
 IMPLICIT NONE
 INTEGER :: I
 REAL(8) :: Fq(-2:N+3)
 
 !$omp parallel do
 DO I = 1, N
 	Fq(i) = Fq(i) + DT/4.0*(-3.0*S0(i)+S1(i))
 END DO
 !$omp end parallel do
 
 END SUBROUTINE
 
 SUBROUTINE RK3_3(Fq)
 use dat
 IMPLICIT NONE
 INTEGER :: I
 REAL(8) :: Fq(-2:N+3)
 
 !$omp parallel do
 DO I = 1, N
 	Fq(i) = Fq(i) + DT/12.0*(-S0(i)-S1(i)+8.0*S2(i))
 END DO
 !$omp end parallel do
 
 END SUBROUTINE
