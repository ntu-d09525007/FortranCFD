MODULE MATRIX_SOLVER
IMPLICIT NONE

CONTAINS

SUBROUTINE TWIN_DEC( A, B, AA, BB, IS, IE )
IMPLICIT NONE
real(8), DIMENSION(3,IS:IE) :: A,  B
real(8), DIMENSION(3,IS:IE) :: AA, BB
!--------------
INTEGER ::  I, IS, IE
real(8) ::  den , SC1 , SC2
 
 A (1,IS) = 0.0_8
 AA(1,IS) = 0.0_8
 A (3,IE) = 0.0_8
 AA(3,IE) = 0.0_8

 B (1,IS) = 0.0_8
 BB(1,IS) = 0.0_8
 B (3,IE) = 0.0_8
 BB(3,IE) = 0.0_8

 DO i = IS+1, IE
   Den = A(2,i-1)*BB(2,i-1) - AA(2,i-1)*B(2,i-1)
   SC1 = -( AA(2,i-1)*B(1,i) - A(1,i)*BB(2,i-1) )/Den
   SC2 =  (  B(1,i)*A(2,i-1) - B(2,i-1)*A(1,i)  )/Den

   A(1,i) = SC1
   A(2,i) = A(2,i) - ( SC1*A(3,i-1) + SC2*AA(3,i-1) )

   B(1,i) = SC2
   B(2,i) = B(2,i) - ( SC1*B(3,i-1) + SC2*BB(3,i-1) )
  !=========================================
   SC1 = -( AA(2,i-1)*BB(1,i) - AA(1,i)*BB(2,i-1) )/Den
   SC2 =  ( BB(1,i)*A(2,i-1)  -  B(2,i-1)*AA(1,i) )/Den

   AA(1,i) = SC1
   AA(2,i) = AA(2,i) - (SC1*A(3,i-1) + SC2*AA(3,i-1))

   BB(1,i) = SC2
   BB(2,i) = BB(2,i) - (SC1*B(3,i-1) + SC2*BB(3,i-1))
 END DO

END SUBROUTINE TWIN_DEC

SUBROUTINE TWIN_BKS( A, B, AA, BB, S, SS, IS, IE )
IMPLICIT NONE
real(8), DIMENSION(3,IS:IE) :: A,  B
real(8), DIMENSION(3,IS:IE) :: AA, BB
real(8), DIMENSION(IS:IE)   :: S , SS
!--------------
INTEGER :: I, IS, IE
real(8) :: Den, SolS, SolSS

 DO i = IS+1, IE
   S(i)  = S(i)  -(  A(1,i)*S(i-1) +  B(1,i)*SS(i-1) )
   SS(i) = SS(i) -( AA(1,i)*S(i-1) + BB(1,i)*SS(i-1) )
 END DO

 Den   = A(2,IE)*BB(2,IE) - AA(2,IE)*B(2,IE)
 SolS  = -( B(2,IE)*SS(IE) - BB(2,IE)*S(IE) )/Den
 SolSS =  ( A(2,IE)*SS(IE) - AA(2,IE)*S(IE) )/Den
 S(IE)  = SolS
 SS(IE) = SolSS

 DO i = IE-1, IS, -1
   S(i)  = S(i)  - (  A(3,i)*SolS +  B(3,i)*SolSS )
   SS(i) = SS(i) - ( AA(3,i)*SolS + BB(3,i)*SolSS )

   Den   = A(2,i)*BB(2,i) - AA(2,i)*B(2,i)
   SolS  = -( B(2,i)*SS(i) - BB(2,i)*S(i))/Den
   SolSS =  ( A(2,i)*SS(i) - AA(2,i)*S(i))/Den
   S(i)  = SolS
   SS(i) = SolSS
 END DO

END SUBROUTINE TWIN_BKS

END MODULE MATRIX_SOLVER

MODULE SRKCCD_solvers
USE MATRIX_SOLVER
IMPLICIT NONE

TYPE SRKCCD_roots
INTEGER :: IS, IE
real(8) :: DX, DT
real(8),ALLOCATABLE,DIMENSION(:) :: F, C3
real(8),ALLOCATABLE,DIMENSION(:) :: S, SS 
real(8),ALLOCATABLE,DIMENSION(:,:) :: A, AA, B, BB
CONTAINS
PROCEDURE alloc => SRKCCD_roots_alloc
PROCEDURE ASSIGN_MATRIX_up => SRKCCD_roots_ASSIGN_MATRIX_up
PROCEDURE ASSIGN_MATRIX_dn => SRKCCD_roots_ASSIGN_MATRIX_dn
PROCEDURE ASSIGN_SRC_up => SRKCCD_roots_ASSIGN_SOURCE_up
PROCEDURE ASSIGN_SRC_dn => SRKCCD_roots_ASSIGN_SOURCE_dn
PROCEDURE SOLVE => SRKCCD_roots_SOLVE
END TYPE SRKCCD_roots

type srkccd_solver
type(srkccd_roots) :: x, y
contains
procedure alloc => srkccd_solver_alloc
end type srkccd_solver

CONTAINS

FUNCTION SRKCCD_roots_c3(X) RESULT(Y)
IMPLICIT NONE
real(8) :: a,b,c,d,e, y
real(8),intent(in) :: X

 A= 0.0001671_8
 B= 7.43943_8
 C= 0.840798_8
 D= 0.00084915_8
 E= -0.159194_8
 
 Y= A*x**B+C*x**D+E
 
END FUNCTION

SUBROUTINE SRKCCD_roots_SOLVE(THIS,U,F,FX,FXX)
IMPLICIT NONE
CLASS(SRKCCD_roots) :: THIS
real(8),DIMENSION(THIS%IS:THIS%IE) :: U,F,FX
real(8),DIMENSION(THIS%IS:THIS%IE),optional :: FXX
INTEGER :: I

 DO I = THIS%IS, THIS%IE
   THIS%F(I) = F(I)
   THIS%C3(I) = SRKCCD_roots_c3(ABS(U(I))*THIS%DT/THIS%DX)
 ENDDO
 
 CALL THIS%ASSIGN_MATRIX_up
 CALL THIS%ASSIGN_SRC_up
 CALL TWIN_DEC( THIS%A, THIS%B, THIS%AA, THIS%BB, THIS%IS, THIS%IE )
 CALL TWIN_BKS( THIS%A, THIS%B, THIS%AA, THIS%BB, THIS%S, THIS%SS, THIS%IS, THIS%IE )
 
 DO I = THIS%IS, THIS%IE
 	FX(I) = THIS%S(I)
 ENDDO
 
 IF( PRESENT(FXX) )THEN
 	DO I = THIS%IS, THIS%IE
		FXX(I) = THIS%SS(I)
 	ENDDO 
 ENDIF
 
 CALL THIS%ASSIGN_MATRIX_dn
 CALL THIS%ASSIGN_SRC_dn
 CALL TWIN_DEC( THIS%A, THIS%B, THIS%AA, THIS%BB, THIS%IS, THIS%IE )
 CALL TWIN_BKS( THIS%A, THIS%B, THIS%AA, THIS%BB, THIS%S, THIS%SS, THIS%IS, THIS%IE )

 DO I = THIS%IS, THIS%IE
    if( u(i)<0.0_8 )then
 	  FX(I) = THIS%S(I)
	endif
 ENDDO
 
 IF( PRESENT(FXX) )THEN
 	DO I = THIS%IS, THIS%IE
		FXX(I) = 0.5_8*(fxx(i)+THIS%SS(I))
 	ENDDO 
ENDIF
	
END SUBROUTINE

SUBROUTINE SRKCCD_roots_alloc(THIS,IS,IE,DX,DT)
IMPLICIT NONE
CLASS(SRKCCD_roots) :: THIS
INTEGER, intent(in) :: IS, IE
real(8), intent(in) :: DX, dt

THIS%IS = IS; THIS%IE = IE; THIS%DX = DX; THIS%DT = DT

ALLOCATE( THIS%F(IS:IE), THIS%C3(IS:IE) )
ALLOCATE( THIS%S(IS:IE), THIS%SS(IS:IE) )
ALLOCATE( THIS%A(3,IS:IE), THIS%AA(3,IS:IE), THIS%B(3,IS:IE), THIS%BB(3,IS:IE) )

END SUBROUTINE

SUBROUTINE SRKCCD_roots_ASSIGN_MATRIX_up(THIS)
IMPLICIT NONE
CLASS(SRKCCD_roots) :: THIS
REAL(8) :: A1, A3, B1, B3
REAL(8) :: AA1, AA3, BB1, BB3
INTEGER :: I

  AA1 = -9.0_8 / 8.0_8
  AA3 =  9.0_8 / 8.0_8
  BB1 = -1.0_8 / 8.0_8
  BB3 = -1.0_8 / 8.0_8

  !For I = 1
  !I
   !For U_X
   this%A(2,THIS%IS) = 1.0_8
   this%B(2,THIS%IS) = 0.0_8
   !For U_XX
   this%AA(2,THIS%IS) = 0.0_8
   this%BB(2,THIS%IS) = 1.0_8
  !I+1
   !For U_X
   this%A(3,THIS%IS) = 2.0_8
   this%B(3,THIS%IS) = -this%Dx
   !For U_XX
   this%AA(3,THIS%IS) = -2.5_8 / this%Dx
   this%BB(3,THIS%IS) =  8.5_8

  !For I = 2 , N-1
  
  !$omp parallel do private(A1,A3,B1,B3)
  DO I = this%IS+1 , this%IE-1
  
    A1 =  131.0_8 / 128.0_8 - 5.0_8 * THIS%C3(I) / 8.0_8
	A3 = -19.0_8 / 128.0_8 + 5.0_8 * THIS%C3(I) / 8.0_8
	
	B1 = 23.0_8 / 128.0_8 - THIS%C3(I) / 8.0_8
	B3 = 7.0_8 / 128.0_8 - THIS%C3(I) / 8.0_8
	
   !For U_X
    this%A(1,I) = A1
	this%A(2,I) = 1.0_8
	this%A(3,I) = A3
	
    this%B(1,I) = B1 * this%Dx
	this%B(2,I) = 0.0_8
	this%B(3,I) = B3 * this%Dx
	
   !For U_XX
    this%AA(1,I) = AA1 / this%Dx
	this%AA(2,I) = 0.0_8
	this%AA(3,I) = AA3 / this%Dx
	
    this%BB(1,I) = BB1
	this%BB(2,I) = 1.0_8
	this%BB(3,I) = BB3
	
  END DO
  !$omp end parallel do

  !For I = N
  !I-1
   !For U_X
   this%A(1,this%ie) = 2.0_8
   this%B(1,this%ie) = this%Dx
   !For U_XX
   this%AA(1,this%ie) = 2.5_8 / THIS%Dx
   this%BB(1,this%ie) = 8.5_8
  !I
   !For U_X
   this%A(2,this%ie) = 1.0_8
   this%B(2,this%ie) = 0.0_8
   !For U_XX
   this%AA(2,this%ie) = 0.0_8
   this%BB(2,this%ie) = 1.0_8  
   
END SUBROUTINE

SUBROUTINE SRKCCD_roots_ASSIGN_MATRIX_dn(THIS)
IMPLICIT NONE
CLASS(SRKCCD_roots) :: THIS
REAL(8) :: A1, A3, B1, B3
REAL(8) :: AA1, AA3, BB1, BB3
INTEGER :: I

  AA1 = -9.0_8 / 8.0_8
  AA3 =  9.0_8 / 8.0_8
  BB1 = -1.0_8 / 8.0_8
  BB3 = -1.0_8 / 8.0_8

  !For I = 1
  !I
   !For U_X
   this%A(2,THIS%IS) = 1.0_8
   this%B(2,THIS%IS) = 0.0_8
   !For U_XX
   this%AA(2,THIS%IS) = 0.0_8
   this%BB(2,THIS%IS) = 1.0_8
  !I+1
   !For U_X
   this%A(3,THIS%IS) = 2.0_8
   this%B(3,THIS%IS) = -this%Dx
   !For U_XX
   this%AA(3,THIS%IS) = -2.5_8 / this%Dx
   this%BB(3,THIS%IS) =  8.5_8

  !For I = 2 , N-1
  
  !$omp parallel do private(A1,A3,B1,B3)
  DO I = this%IS+1 , this%IE-1
  
    A1 =  131.0_8 / 128.0_8 - 5.0_8 * THIS%C3(I) / 8.0_8
	A3 = -19.0_8 / 128.0_8 + 5.0_8 * THIS%C3(I) / 8.0_8
	
	B1 = 23.0_8 / 128.0_8 - THIS%C3(I) / 8.0_8
	B3 = 7.0_8 / 128.0_8 - THIS%C3(I) / 8.0_8
	
   !For U_X
    this%A(1,I) = A3
	this%A(2,I) = 1.0_8
	this%A(3,I) = A1
	
    this%B(1,I) = -B3 * this%Dx
	this%B(2,I) = 0.0_8
	this%B(3,I) = -B1 * this%Dx
	
   !For U_XX
    this%AA(1,I) = AA1 / this%Dx
	this%AA(2,I) = 0.0_8
	this%AA(3,I) = AA3 / this%Dx
	
    this%BB(1,I) = BB1
	this%BB(2,I) = 1.0_8
	this%BB(3,I) = BB3
	
  END DO
  !$omp end parallel do

  !For I = N
  !I-1
   !For U_X
   this%A(1,this%ie) = 2.0_8
   this%B(1,this%ie) = this%Dx
   !For U_XX
   this%AA(1,this%ie) = 2.5_8 / THIS%Dx
   this%BB(1,this%ie) = 8.5_8
  !I
   !For U_X
   this%A(2,this%ie) = 1.0_8
   this%B(2,this%ie) = 0.0_8
   !For U_XX
   this%AA(2,this%ie) = 0.0_8
   this%BB(2,this%ie) = 1.0_8  
   
END SUBROUTINE


SUBROUTINE SRKCCD_roots_ASSIGN_SOURCE_up(THIS)
IMPLICIT NONE
CLASS(SRKCCD_roots) :: THIS
real(8) :: C1 , C2 , C3
real(8) :: CC1 , CC2 , CC3
INTEGER:: I

 !For uxx eq.
  CC1 =  3.0_8
  CC2 = -6.0_8
  CC3 =  3.0_8
 
 
  THIS%S(THIS%IS)  = (-3.5_8 * THIS%F(THIS%IS) + 4.0_8 * THIS%F(THIS%IS+1) - 0.5_8 * THIS%F(THIS%IS+2)) / THIS%Dx

  THIS%SS(THIS%IS) = ( 34.0_8 / 3.0_8 * THIS%F(THIS%IS) - 83.0_8 / 4.0_8 * THIS%F(THIS%IS+1) &
          &  +10.0_8 * THIS%F(THIS%IS+2) -7.0_8/12.0_8 * THIS%F(THIS%IS) ) / (THIS%Dx**2)

  !For I = 2 , N-1
  !$omp parallel do private(c3,c2,c1)
  DO I = THIS%IS+1 , THIS%IE-1
  
    C3 = THIS%C3(I)
    c2 = 15.0_8/8.0_8 - 2.0_8*C3
    c1 = C3 - 15.0_8 / 8.0_8
    
  	!U_X
    THIS%S(I)  = (C1*THIS%F(I-1) + C2*THIS%F(I) + C3*THIS%F(I+1)) / THIS%Dx
    !U_XX
    THIS%SS(I) = (CC1*THIS%F(I-1) + CC2*THIS%F(I) + CC3*THIS%F(I+1)) / (THIS%Dx**2)
  END DO
  !$omp end parallel do

  !For I = N
   !U_X
   THIS%S(THIS%IE)  = -(-3.5_8 * THIS%F(THIS%IE) + 4.0_8 * THIS%F(THIS%IE-1) - 0.5_8 * THIS%F(THIS%IE-2)) / THIS%Dx
   !U_XX
   THIS%SS(THIS%IE) = ( 34.0_8 / 3.0_8 * THIS%F(THIS%IE) - 83.0_8 / 4.0_8 * THIS%F(THIS%IE-1) &
         &  +10.0_8 * THIS%F(THIS%IE-2) -7.0_8/12.0_8 * THIS%F(THIS%IE) ) / (THIS%Dx**2)


END SUBROUTINE
   
SUBROUTINE SRKCCD_roots_ASSIGN_SOURCE_dn(THIS)
IMPLICIT NONE
CLASS(SRKCCD_roots) :: THIS
real(8) :: C1 , C2 , C3
real(8) :: CC1 , CC2 , CC3
INTEGER:: I

 !For uxx eq.
  CC1 =  3.0_8
  CC2 = -6.0_8
  CC3 =  3.0_8
 
 
  THIS%S(THIS%IS)  = (-3.5_8 * THIS%F(THIS%IS) + 4.0_8 * THIS%F(THIS%IS+1) - 0.5_8 * THIS%F(THIS%IS+2)) / THIS%Dx

  THIS%SS(THIS%IS) = ( 34.0_8 / 3.0_8 * THIS%F(THIS%IS) - 83.0_8 / 4.0_8 * THIS%F(THIS%IS+1) &
          &  +10.0_8 * THIS%F(THIS%IS+2) -7.0_8/12.0_8 * THIS%F(THIS%IS) ) / (THIS%Dx**2)

  !For I = 2 , N-1
  !$omp parallel do private(c3,c2,c1)
  DO I = THIS%IS+1 , THIS%IE-1
  
    C3 = THIS%C3(I)
    c2 = 15.0_8/8.0_8 - 2.0_8*C3
    c1 = C3 - 15.0_8 / 8.0_8
    
  	!U_X
    THIS%S(I)  = -(C3*THIS%F(I-1) + C2*THIS%F(I) + C1*THIS%F(I+1)) / THIS%Dx
    !U_XX
    THIS%SS(I) = (CC1*THIS%F(I-1) + CC2*THIS%F(I) + CC3*THIS%F(I+1)) / (THIS%Dx**2)
  END DO
  !$omp end parallel do

  !For I = N
   !U_X
   THIS%S(THIS%IE)  = -(-3.5_8 * THIS%F(THIS%IE) + 4.0_8 * THIS%F(THIS%IE-1) - 0.5_8 * THIS%F(THIS%IE-2)) / THIS%Dx
   !U_XX
   THIS%SS(THIS%IE) = ( 34.0_8 / 3.0_8 * THIS%F(THIS%IE) - 83.0_8 / 4.0_8 * THIS%F(THIS%IE-1) &
         &  +10.0_8 * THIS%F(THIS%IE-2) -7.0_8/12.0_8 * THIS%F(THIS%IE) ) / (THIS%Dx**2)


END SUBROUTINE

subroutine srkccd_solver_alloc(p,is,ie,js,je,dx,dy,dt)
implicit none
class(srkccd_solver) :: p
integer, intent(in) :: is, ie, js, je
real(8), intent(in) :: dx, dy, dt

 call p%x%alloc(is,ie,dx,dt)
 call p%y%alloc(js,je,dy,dt)

end subroutine
	 
END MODULE SRKCCD_solvers
