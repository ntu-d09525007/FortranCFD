MODULE SRKCCD
use my_precision
USE MATRIX_SOLVER
IMPLICIT NONE

TYPE SRKCCD_1D
INTEGER :: IS, IE
real(kind=ap) :: DX, DT
real(kind=ap),ALLOCATABLE,DIMENSION(:) :: F, C3
real(kind=ap),ALLOCATABLE,DIMENSION(:) :: S, SS 
real(kind=ap),ALLOCATABLE,DIMENSION(:,:) :: A, AA, B, BB
CONTAINS
PROCEDURE INIT => SRKCCD_INITIALIZED
PROCEDURE ASSIGN_MATRIX_up => SRKCCD_ASSIGN_MATRIX_up
PROCEDURE ASSIGN_MATRIX_dn => SRKCCD_ASSIGN_MATRIX_dn
PROCEDURE ASSIGN_SRC_up => SRKCCD_ASSIGN_SOURCE_up
PROCEDURE ASSIGN_SRC_dn => SRKCCD_ASSIGN_SOURCE_dn
PROCEDURE SOLVE => SRKCCD_SOLVE_1D
END TYPE SRKCCD_1D

CONTAINS

FUNCTION srkccd_C3(X) RESULT(Y)
IMPLICIT NONE
real(kind=ap) :: a,b,c,d,e
real(kind=ap) :: X, Y

 A= 0.0001671_AP
 B= 7.43943_AP
 C= 0.840798_AP
 D= 0.00084915_AP
 E= -0.159194_AP
 
 Y= A*X**B+C*X**D+E
 
END FUNCTION

SUBROUTINE SRKCCD_SOLVE_1D(THIS,U,F,FX,FXX)
IMPLICIT NONE
CLASS(SRKCCD_1D) :: THIS
real(kind=ap),DIMENSION(THIS%IS-3:THIS%IE+3) :: U,F,FX
real(kind=ap),DIMENSION(THIS%IS-3:THIS%IE+3),optional :: FXX
INTEGER :: I

 !$OMP PARALLEL DO
 DO I = THIS%IS-3, THIS%IE+3
   THIS%F(I) = F(I)
   THIS%C3(I) = srkccd_C3(ABS(U(I))*THIS%DT/THIS%DX)
 ENDDO
 !$OMP END PARALLEL DO
 
 CALL THIS%ASSIGN_MATRIX_up
 CALL THIS%ASSIGN_SRC_up
 CALL TWIN_DEC( THIS%A, THIS%B, THIS%AA, THIS%BB, THIS%IS, THIS%IE )
 CALL TWIN_BKS( THIS%A, THIS%B, THIS%AA, THIS%BB, THIS%S, THIS%SS, THIS%IS, THIS%IE )
 
 !$OMP PARALLEL DO
 DO I = THIS%IS-3, THIS%IE+3 
 	FX(I) = THIS%S(I)
 ENDDO
 !$OMP END PARALLEL DO
 
 IF( PRESENT(FXX) )THEN
    !$OMP PARALLEL DO
 	DO I = THIS%IS-3, THIS%IE+3 
		FXX(I) = THIS%SS(I)
 	ENDDO 
	!$OMP END PARALLEL DO
 ENDIF
 
 CALL THIS%ASSIGN_MATRIX_dn
 CALL THIS%ASSIGN_SRC_dn
 CALL TWIN_DEC( THIS%A, THIS%B, THIS%AA, THIS%BB, THIS%IS, THIS%IE )
 CALL TWIN_BKS( THIS%A, THIS%B, THIS%AA, THIS%BB, THIS%S, THIS%SS, THIS%IS, THIS%IE )

 !$OMP PARALLEL DO
 DO I = THIS%IS-3, THIS%IE+3 
    if( u(i)<0.0_ap )then
 	  FX(I) = THIS%S(I)
	endif
 ENDDO
 !$OMP END PARALLEL DO
 
 IF( PRESENT(FXX) )THEN
    !$OMP PARALLEL DO
 	DO I = THIS%IS-3, THIS%IE+3 
		FXX(I) = 0.5_ap*(fxx(i)+THIS%SS(I))
 	ENDDO 
	!$OMP END PARALLEL DO
ENDIF
	
END SUBROUTINE


SUBROUTINE SRKCCD_INITIALIZED(THIS,IS,IE,DX,DT)
IMPLICIT NONE
CLASS(SRKCCD_1D) :: THIS
INTEGER :: IS, IE
real(kind=ap) :: DX, dt

THIS%IS = IS; THIS%IE = IE; THIS%DX = DX; THIS%DT = DT

ALLOCATE( THIS%F(IS-3:IE+3), THIS%C3(IS-3:IE+3) )
ALLOCATE( THIS%S(IS-3:IE+3), THIS%SS(IS-3:IE+3) )
ALLOCATE( THIS%A(3,IS-3:IE+3), THIS%AA(3,IS-3:IE+3), THIS%B(3,IS-3:IE+3), THIS%BB(3,IS-3:IE+3) )

END SUBROUTINE

SUBROUTINE SRKCCD_ASSIGN_MATRIX_up(THIS)
IMPLICIT NONE
CLASS(SRKCCD_1D) :: THIS
REAL(KIND=AP) :: A1, A3, B1, B3
REAL(KIND=AP) :: AA1, AA3, BB1, BB3
REAL(KIND=AP) :: C3
INTEGER :: I

  AA1 = -9.0_ap / 8.0_ap
  AA3 =  9.0_ap / 8.0_ap
  BB1 = -1.0_ap / 8.0_ap
  BB3 = -1.0_ap / 8.0_ap

  !For I = 1
  !I
   !For U_X
   this%A(2,THIS%IS-3) = 1.0_ap
   this%B(2,THIS%IS-3) = 0.0_ap
   !For U_XX
   this%AA(2,THIS%IS-3) = 0.0_ap
   this%BB(2,THIS%IS-3) = 1.0_ap
  !I+1
   !For U_X
   this%A(3,THIS%IS-3) = 2.0_ap
   this%B(3,THIS%IS-3) = -this%Dx
   !For U_XX
   this%AA(3,THIS%IS-3) = -2.5_ap / this%Dx
   this%BB(3,THIS%IS-3) =  8.5_ap

! ==============================================

  !$omp parallel do private(A1,A3,B1,B3)
  DO I = this%IS-2 , this%Ie+2
  
    A1 = 131.0_AP / 128.0_AP - 5.0_AP * this%c3(i) / 8.0_AP
	A3 = -19.0_AP / 128.0_AP + 5.0_AP * this%c3(i) / 8.0_AP
	
	B1 = 23.0_AP / 128.0_AP - this%c3(i) / 8.0_AP
	B3 =  7.0_AP / 128.0_AP - this%c3(i) / 8.0_AP
	
   !For U_X
    this%A(1,I) = A1
	this%A(2,I) = 1.0_ap
	this%A(3,I) = A3
	
    this%B(1,I) = B1 * this%Dx
	this%B(2,I) = 0.0_AP
	this%B(3,I) = B3 * this%Dx
	
   !For U_XX
    this%AA(1,I) = AA1 / this%Dx
	this%AA(2,I) = 0.0_ap
	this%AA(3,I) = AA3 / this%Dx
	
    this%BB(1,I) = BB1
	this%BB(2,I) = 1.0_ap
	this%BB(3,I) = BB3
	
  END DO
  !$omp end parallel do
  
  ! ==============================================

   !For U_X
   this%A(1,this%ie+3) = 2.0_ap
   this%B(1,this%ie+3) = this%Dx
   !For U_XX
   this%AA(1,this%ie+3) = 2.5_ap / THIS%Dx
   this%BB(1,this%ie+3) = 8.5_ap

   !For U_X
   this%A(2,this%ie+3) = 1.0_ap
   this%B(2,this%ie+3) = 0.0_ap
   !For U_XX
   this%AA(2,this%ie+3) = 0.0_ap
   this%BB(2,this%ie+3) = 1.0_ap  
   
END SUBROUTINE

SUBROUTINE SRKCCD_ASSIGN_MATRIX_dn(THIS)
IMPLICIT NONE
CLASS(SRKCCD_1D) :: THIS
REAL(KIND=AP) :: A1, A3, B1, B3
REAL(KIND=AP) :: AA1, AA3, BB1, BB3
INTEGER :: I

  AA1 = -9.0_ap / 8.0_ap
  AA3 =  9.0_ap / 8.0_ap
  BB1 = -1.0_ap / 8.0_ap
  BB3 = -1.0_ap / 8.0_ap

  !For I = 1
  !I
   !For U_X
   this%A(2,THIS%IS-3) = 1.0_ap
   this%B(2,THIS%IS-3) = 0.0_ap
   !For U_XX
   this%AA(2,THIS%IS-3) = 0.0_ap
   this%BB(2,THIS%IS-3) = 1.0_ap
  !I+1
   !For U_X
   this%A(3,THIS%IS-3) = 2.0_ap
   this%B(3,THIS%IS-3) = -this%Dx
   !For U_XX
   this%AA(3,THIS%IS-3) = -2.5_ap / this%Dx
   this%BB(3,THIS%IS-3) =  8.5_ap

! ============================================== 

  !For I = 2 , N-1
  
  !$omp parallel do private(A1,A3,B1,B3)
  DO I = this%IS-2 , this%IE+2
  
    A1 = 131.0_AP / 128.0_AP - 5.0_AP * this%c3(i) / 8.0_AP
	A3 = -19.0_AP / 128.0_AP + 5.0_AP * this%c3(i) / 8.0_AP
	
	B1 = 23.0_AP / 128.0_AP - this%c3(i) / 8.0_AP
	B3 =  7.0_AP / 128.0_AP - this%c3(i) / 8.0_AP
	
   !For U_X
    this%A(1,I) = A3
	this%A(2,I) = 1.0_ap
	this%A(3,I) = A1
	
    this%B(1,I) = -B3 * this%Dx
	this%B(2,I) = 0.0_AP
	this%B(3,I) = -B1 * this%Dx
	
   !For U_XX
    this%AA(1,I) = AA1 / this%Dx
	this%AA(2,I) = 0.0_ap
	this%AA(3,I) = AA3 / this%Dx
	
    this%BB(1,I) = BB1
	this%BB(2,I) = 1.0_ap
	this%BB(3,I) = BB3
	
  END DO
  !$omp end parallel do
  
  ! ==============================================

   !For U_X
   this%A(1,this%ie+3) = 2.0_ap
   this%B(1,this%ie+3) = this%Dx
   !For U_XX
   this%AA(1,this%ie+3) = 2.5_ap / THIS%Dx
   this%BB(1,this%ie+3) = 8.5_ap
   !For U_X
   this%A(2,this%ie+3) = 1.0_ap
   this%B(2,this%ie+3) = 0.0_ap
   !For U_XX
   this%AA(2,this%ie+3) = 0.0_ap
   this%BB(2,this%ie+3) = 1.0_ap  
   
END SUBROUTINE


SUBROUTINE SRKCCD_ASSIGN_SOURCE_up(THIS)
IMPLICIT NONE
CLASS(SRKCCD_1D) :: THIS
real(kind=ap) :: C1 , C2 , C3, C4
real(kind=ap) :: CC1 , CC2 , CC3
INTEGER:: I

 !For uxx eq.
  CC1 =  3.0_ap
  CC2 = -6.0_ap
  CC3 =  3.0_ap
 
 
  THIS%S(THIS%IS-3)  = (-3.5_ap * THIS%F(THIS%IS-3) + 4.0_ap * THIS%F(THIS%IS-2) - 0.5_ap * THIS%F(THIS%IS-1)) / THIS%Dx

  THIS%SS(THIS%IS-3) = ( 34.0_ap / 3.0_ap * THIS%F(THIS%IS-3) - 83.0_ap / 4.0_ap * THIS%F(THIS%IS-2) &
          &  +10.0_ap * THIS%F(THIS%IS-1) -7.0_ap/12.0_ap * THIS%F(THIS%IS) ) / (THIS%Dx**2)

  !For I = 2 , N-1
  !$omp parallel do private(c3,c2,c1,c4)
  DO I = THIS%IS-2 , THIS%IE+2
     
    C3 = THIS%C3(I)
    C2 = 15.0_AP / 8.0_AP - 2.0_AP*C3
  	C1 = -15.0_AP / 8.0_AP + C3
    
  	!U_X
    THIS%S(I)  = (C1*THIS%F(I-1) + C2*THIS%F(I) + C3*THIS%F(I+1)  ) / THIS%Dx
    !U_XX
    THIS%SS(I) = (CC1*THIS%F(I-1) + CC2*THIS%F(I) + CC3*THIS%F(I+1)) / (THIS%Dx**2)
  END DO
  !$omp end parallel do

  !For I = N
   !U_X
   THIS%S(THIS%IE+3)  = -(-3.5_ap * THIS%F(THIS%IE+3) + 4.0_ap * THIS%F(THIS%IE+2) - 0.5_ap * THIS%F(THIS%IE+1)) / THIS%Dx
   !U_XX
   THIS%SS(THIS%IE+3) = ( 34.0_ap / 3.0_ap * THIS%F(THIS%IE+3) - 83.0_ap / 4.0_ap * THIS%F(THIS%IE+2) &
         &  +10.0_ap * THIS%F(THIS%IE+1) -7.0_ap/12.0_ap * THIS%F(THIS%IE) ) / (THIS%Dx**2)


END SUBROUTINE
   
SUBROUTINE SRKCCD_ASSIGN_SOURCE_dn(THIS)
IMPLICIT NONE
CLASS(SRKCCD_1D) :: THIS
real(kind=ap) :: C1 , C2 , C3, C4
real(kind=ap) :: CC1 , CC2 , CC3
INTEGER:: I

 !For uxx eq.
  CC1 =  3.0_ap
  CC2 = -6.0_ap
  CC3 =  3.0_ap
 
 
  THIS%S(THIS%IS-3)  = (-3.5_ap * THIS%F(THIS%IS-3) + 4.0_ap * THIS%F(THIS%IS-2) - 0.5_ap * THIS%F(THIS%IS-1)) / THIS%Dx

  THIS%SS(THIS%IS-3) = ( 34.0_ap / 3.0_ap * THIS%F(THIS%IS-3) - 83.0_ap / 4.0_ap * THIS%F(THIS%IS-2) &
          &  +10.0_ap * THIS%F(THIS%IS-1) -7.0_ap/12.0_ap * THIS%F(THIS%IS) ) / (THIS%Dx**2)

  !For I = 2 , N-1
  !$omp parallel do private(c3,c2,c1)
  DO I = THIS%IS-2 , THIS%IE+2
  
    C3 = THIS%C3(I)
    C2 = 15.0_AP / 8.0_AP - 2.0_AP*C3
  	C1 = -15.0_AP / 8.0_AP + C3
    
  	!U_X
    THIS%S(I)  = -(C3*THIS%F(I-1) + C2*THIS%F(I) + C1*THIS%F(I+1) ) / THIS%Dx
    !U_XX
    THIS%SS(I) = (CC1*THIS%F(I-1) + CC2*THIS%F(I) + CC3*THIS%F(I+1)) / (THIS%Dx**2)
  END DO
  !$omp end parallel do
  
   !U_X
   THIS%S(THIS%IE+3)  = -(-3.5_ap * THIS%F(THIS%IE+3) + 4.0_ap * THIS%F(THIS%IE+2) - 0.5_ap * THIS%F(THIS%IE+1)) / THIS%Dx
   !U_XX
   THIS%SS(THIS%IE+3) = ( 34.0_ap / 3.0_ap * THIS%F(THIS%IE+3) - 83.0_ap / 4.0_ap * THIS%F(THIS%IE+2) &
         &  +10.0_ap * THIS%F(THIS%IE+1) -7.0_ap/12.0_ap * THIS%F(THIS%IE) ) / (THIS%Dx**2)


END SUBROUTINE
	 
END MODULE