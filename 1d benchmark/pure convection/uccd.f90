MODULE SRKUCCD
use my_precision
USE MATRIX_SOLVER
IMPLICIT NONE

TYPE SRKUCCD_1D
INTEGER :: IS, IE
real(kind=ap) :: DX, DT
real(kind=ap),ALLOCATABLE,DIMENSION(:) :: F, C3
real(kind=ap),ALLOCATABLE,DIMENSION(:) :: S, SS 
real(kind=ap),ALLOCATABLE,DIMENSION(:) :: SU, SSU, SD, SSD
real(kind=ap),ALLOCATABLE,DIMENSION(:,:) :: A, AA, B, BB
real(kind=ap),ALLOCATABLE,DIMENSION(:,:) :: AU, AAU, BU, BBU
real(kind=ap),ALLOCATABLE,DIMENSION(:,:) :: AD, AAD, BD, BBD
CONTAINS
PROCEDURE INIT => uccd_INITIALIZED
PROCEDURE ASSIGN_u => uccd_ASSIGN_MATRIX_upwind
procedure assign_d => uccd_ASSIGN_MATRIX_downwind
PROCEDURE ASSIGN_SRC => uccd_ASSIGN_SOURCE
PROCEDURE SOLVE => uccd_SOLVE_1D
END TYPE SRKUCCD_1D

TYPE uccd_SOLVER
TYPE(SRKUCCD_1D) :: X, Y, Z
contains
procedure alloc => uccd_solver_alloc
END TYPE uccd_SOLVER

CONTAINS

subroutine uccd_solver_alloc(THIS,IS,IE,JS,JE,KS,KE,DX,DY,DZ,DT)
implicit none
class(uccd_SOLVER) :: this
integer :: IS,IE,JS,JE,KS,KE
real(kind=ap) :: dx,dy,dz,dt

CALL THIS%X%init(IS,IE,DX,DT)
CALL THIS%Y%init(JS,JE,DY,DT)
CALL THIS%Z%init(KS,KE,DZ,DT)
	
end subroutine

FUNCTION uccd_C3(X) RESULT(Y)
IMPLICIT NONE
real(kind=ap) :: a,b,c,d,e
real(kind=ap) :: X, Y

y=0.0_AP
 
END FUNCTION

SUBROUTINE uccd_INITIALIZED(THIS,IS,IE,DX,DT)
IMPLICIT NONE
CLASS(SRKUCCD_1D) :: THIS
INTEGER :: IS, IE
real(kind=ap) :: DX, dt

THIS%IS = IS; THIS%IE = IE; THIS%DX = DX; THIS%DT = DT

ALLOCATE( THIS%F(IS:IE), THIS%C3(IS:IE) )
ALLOCATE( THIS%SU(IS:IE), THIS%SSU(IS:IE), THIS%SD(IS:IE), THIS%SSD(IS:IE) )
ALLOCATE( THIS%S(IS:IE), THIS%SS(IS:IE) )
ALLOCATE( THIS%AU(3,IS:IE), THIS%AAU(3,IS:IE), THIS%BU(3,IS:IE), THIS%BBU(3,IS:IE) )
ALLOCATE( THIS%AD(3,IS:IE), THIS%AAD(3,IS:IE), THIS%BD(3,IS:IE), THIS%BBD(3,IS:IE) )
ALLOCATE( THIS%A(3,IS:IE), THIS%AA(3,IS:IE), THIS%B(3,IS:IE), THIS%BB(3,IS:IE) )

END SUBROUTINE

SUBROUTINE uccd_SOLVE_1D(THIS,U,F,FX,FXX)
IMPLICIT NONE
CLASS(SRKUCCD_1D) :: THIS
real(kind=ap),DIMENSION(THIS%IS:THIS%IE) :: U,F,FX
real(kind=ap),DIMENSION(THIS%IS:THIS%IE),optional :: FXX
INTEGER :: I

 DO I = THIS%IS, THIS%IE
   THIS%F(I) = F(I)
   THIS%C3(I) = -0.0609611900811_AP!uccd_C3(ABS(U(I))*THIS%DT/THIS%DX)
 ENDDO
 
 call this%ASSIGN_u
 call this%ASSIGN_d

 CALL THIS%ASSIGN_SRC
 
 CALL TWIN_DEC( THIS%AU, THIS%BU, THIS%AAU, THIS%BBU, THIS%IS, THIS%IE )
 CALL TWIN_BKS( THIS%AU, THIS%BU, THIS%AAU, THIS%BBU, THIS%SU, THIS%SSU, THIS%IS, THIS%IE )
 
 CALL TWIN_DEC( THIS%AD, THIS%BD, THIS%AAD, THIS%BBD, THIS%IS, THIS%IE )
 CALL TWIN_BKS( THIS%AD, THIS%BD, THIS%AAD, THIS%BBD, THIS%SD, THIS%SSD, THIS%IS, THIS%IE )
 
 DO I = THIS%IS, THIS%IE 
    IF( U(I)>=0.0_AP )THEN
	  FX(I) = THIS%SU(I)
	ELSE
	  FX(I) = THIS%SD(I)
	ENDIF
 ENDDO
 
 IF( PRESENT(FXX) )THEN
 	DO I = THIS%IS, THIS%IE 
		FXX(I) = 0.5_AP*(THIS%SSU(I)+THIS%SSD(I))
 	ENDDO 
 ENDIF

END SUBROUTINE

SUBROUTINE uccd_ASSIGN_MATRIX_upwind(THIS)
IMPLICIT NONE
CLASS(SRKUCCD_1D) :: THIS
real(kind=ap) :: A1 , B1 , B2 , B3, C3
real(kind=ap) :: A3D , B1D , B2D , B3D
real(kind=ap) :: AA1 , AA3 , BB1 , BB3
INTEGER :: I

  !For ux eq., a>0
  A1 =  0.875_AP
  B1 =  0.1251282341599089_AP
  B2 = -0.2487176584009104_AP
  B3 =  0.0001282341599089_AP

  !For uxx eq.
  AA1 = -9.0_AP / 8.0_AP
  AA3 =  9.0_AP / 8.0_AP
  BB1 = -1.0_AP / 8.0_AP
  BB3 = -1.0_AP / 8.0_AP
  
!For Upwind
  !For I = 1
  !I
   !For U_X
   this%AU(2,THIS%IS) = 1.0_AP
   this%BU(2,THIS%IS) = 0.0_AP
   !For U_XX
   this%AAU(2,THIS%IS) = 0.0_AP
   this%BBU(2,THIS%IS) = 1.0_AP
  !I+1
   !For U_X
   this%AU(3,THIS%IS) = 2.0_AP
   this%BU(3,THIS%IS) = -this%Dx
   !For U_XX
   this%AAU(3,THIS%IS) = -2.5_AP / this%Dx
   this%BBU(3,THIS%IS) =  8.5_AP

  !For I = 2 , N-1
  
  !$omp parallel do 
  DO I = THIS%IS+1 , THIS%IE-1
	
  !I-1
   !For U_X
    this%AU(1,I) = A1
    this%BU(1,I) = B1 * this%Dx
   !For U_XX
    this%AAU(1,I) = AA1 / this%Dx
    this%BBU(1,I) = BB1
  !I
    !For U_X
    this%AU(2,I) = 1.0_AP
    this%BU(2,I) = B2 * this%Dx
    !For U_XX
    this%AAU(2,I) = 0.0_AP
    this%BBU(2,I) = 1.0_AP
  !I+1
    !For U_X
    this%Au(3,I) = 0.0_AP
    this%Bu(3,I) = B3 * this%Dx
    !For U_XX
    this%AAu(3,I) = AA3 / this%Dx
    this%BBu(3,I) = BB3
	
  END DO
  !$omp end parallel do

  !For I = N
  !I-1
   !For U_X
   this%Au(1,THIS%IE) = 2.0_AP
   this%Bu(1,THIS%IE) = this%Dx
   !For U_XX
   this%AAu(1,THIS%IE) = 2.5_AP / THIS%Dx
   this%BBu(1,THIS%IE) = 8.5_AP
  !I
   !For U_X
   this%Au(2,THIS%IE) = 1.0_AP
   this%Bu(2,THIS%IE) = 0.0_AP
   !For U_XX
   this%AAu(2,THIS%IE) = 0.0_AP
   this%BBu(2,THIS%IE) = 1.0_AP  
   
END SUBROUTINE


SUBROUTINE uccd_ASSIGN_MATRIX_downwind(THIS)
IMPLICIT NONE
CLASS(SRKUCCD_1D) :: THIS
real(kind=ap) :: A1 , B1 , B2 , B3, C3
real(kind=ap) :: A3D , B1D , B2D , B3D
real(kind=ap) :: AA1 , AA3 , BB1 , BB3
INTEGER :: I

  !For ux eq., a>0
  A1 =  0.875_AP
  B1 =  0.1251282341599089_AP
  B2 = -0.2487176584009104_AP
  B3 =  0.0001282341599089_AP

  !For ux eq., a<0
  A3D =  0.875_AP
  B1D = -B3!-0.0001004428858241_AP
  B2D = -B2!0.248995571141759_AP
  B3D = -B1!-0.125100442885824_AP

  !For uxx eq.
  AA1 = -9.0_AP / 8.0_AP
  AA3 =  9.0_AP / 8.0_AP
  BB1 = -1.0_AP / 8.0_AP
  BB3 = -1.0_AP / 8.0_AP
  
!For Downwind
  !For I = 1
  !I
   !For U_X
   this%AD(2,THIS%IS) = 1.0_AP
   this%BD(2,THIS%IS) = 0.0_AP
   !For U_XX
   this%AAD(2,THIS%IS) = 0.0_AP
   this%BBD(2,THIS%IS) = 1.0_AP
  !I+1
   !For U_X
   this%AD(3,THIS%IS) = 2.0_AP
   this%BD(3,THIS%IS) = -this%Dx
   !For U_XX
   this%AAD(3,THIS%IS) = -2.5_AP / this%Dx
   this%BBD(3,THIS%IS) =  8.5_AP

  !For I = 2 , N-1
  
  !$omp parallel do 
  DO I = THIS%IS+1 , THIS%IE-1
	
  !I-1
   !For U_X
    this%AD(1,I) = 0.0_AP
    this%BD(1,I) = B1D * this%Dx
   !For U_XX
    this%AAD(1,I) = AA1 / this%Dx
    this%BBD(1,I) = BB1
	
  !I
    !For U_X
    this%AD(2,I) = 1.0_AP
    this%BD(2,I) = B2D * this%Dx
    !For U_XX
    this%AAD(2,I) = 0.0_AP
    this%BBD(2,I) = 1.0_AP
	
  !I+1
    !For U_X
    this%AD(3,I) = A3D
    this%BD(3,I) = B3D * this%Dx
    !For U_XX
    this%AAD(3,I) = AA3 / this%Dx
    this%BBD(3,I) = BB3
  END DO
  !$omp end parallel do

  !For I = N
  !I-1
   !For U_X
   this%AD(1,THIS%IE) = 2.0_AP
   this%BD(1,THIS%IE) = THIS%Dx
   !For U_XX
   this%AAD(1,THIS%IE) = 2.5_AP / this%Dx
   this%BBD(1,THIS%IE) = 8.5_AP
  !I
   !For U_X
   this%AD(2,THIS%IE) = 1.0_AP
   this%BD(2,THIS%IE) = 0.0_AP
   !For U_XX
   this%AAD(2,THIS%IE) = 0.0_AP
   this%BBD(2,THIS%IE) = 1.0_AP

END SUBROUTINE

SUBROUTINE uccd_ASSIGN_SOURCE(THIS)
IMPLICIT NONE
CLASS(SRKUCCD_1D) :: THIS
real(kind=ap) :: C1 , C2 , C3
real(kind=ap) :: C1D , C2D , C3D
real(kind=ap) :: CC1 , CC2 , CC3
INTEGER:: I

 !For ux eq., a>0
  C1 = -1.93596119008109_AP
  C2 =  1.99692238016218_AP
  C3 = -0.06096119008109_AP

 !For ux eq., a<0
  C1D =  -C3!0.061294685370110_AP
  C2D =  -C2!-1.99758937102742_AP
  C3D =  -C1!1.93629468537011_AP

 !For uxx eq.
  CC1 =  3.0_AP
  CC2 = -6.0_AP
  CC3 =  3.0_AP
 
!For upwind

  !For I = 1
   !U_X
    THIS%SU(THIS%IS)  = (-3.5_AP * THIS%F(THIS%IS) + 4.0_AP * THIS%F(THIS%IS+1) - 0.5_AP * THIS%F(THIS%IS+2)) / THIS%Dx
   !U_XX
    THIS%SSU(THIS%IS) = ( 34.0_AP / 3.0_AP * THIS%F(THIS%IS) - 83.0_AP / 4.0_AP * THIS%F(THIS%IS+1) &
          &  +10.0_AP * THIS%F(THIS%IS+2) -7.0_AP/12.0_AP * THIS%F(THIS%IS+3) ) / (THIS%Dx**2)

  !For I = 2 , N-1
  !$omp parallel do private(c3,c2,c1)
  DO I = THIS%IS+1 , THIS%IE-1
 
  	!U_X
    THIS%SU(I)  = (C1*THIS%F(I-1) + C2*THIS%F(I) + C3*THIS%F(I+1)) / THIS%Dx
    !U_XX
    THIS%SSU(I) = (CC1*THIS%F(I-1) + CC2*THIS%F(I) + CC3*THIS%F(I+1)) / (THIS%Dx**2)
  END DO
  !$omp end parallel do

  !For I = N
   !U_X
   THIS%SU(THIS%IE)  = -(-3.5_AP * THIS%F(THIS%IE) + 4.0_AP * THIS%F(THIS%IE-1) - 0.5_AP * THIS%F(THIS%IE-2)) / THIS%Dx
   !U_XX
   THIS%SSU(THIS%IE) = ( 34.0_AP / 3.0_AP * THIS%F(THIS%IE) - 83.0_AP / 4.0_AP * THIS%F(THIS%IE-1) &
         &  +10.0_AP * THIS%F(THIS%IE-2) -7.0_AP/12.0_AP * THIS%F(THIS%IE-3) ) / (THIS%Dx**2)


!For downwind

  !For I = 1
   !U_X
    THIS%SD(THIS%IS)  = (-3.5_AP * THIS%F(THIS%IS) + 4.0_AP * THIS%F(THIS%IS+1) - 0.5_AP * THIS%F(THIS%IS+2)) / THIS%Dx
   !U_XX
    THIS%SSD(THIS%IS) = ( 34.0_AP / 3.0_AP * THIS%F(THIS%IS) - 83.0_AP / 4.0_AP * THIS%F(THIS%IS+1) &
          &  +10.0_AP * THIS%F(THIS%IS+2) -7.0_AP/12.0_AP * THIS%F(THIS%IS+3) ) / (THIS%Dx**2)

  !For I = 2 , N-1
  !$omp parallel do private(c3,c2,c1,c1d,c2d,c3d)
  DO I = THIS%IS+1 , THIS%IE-1
  
  	!U_X
    THIS%SD(I)  = (C1D*THIS%F(I-1) + C2D*THIS%F(I) + C3D*THIS%F(I+1)) / THIS%Dx
    !U_XX
    THIS%SSD(I) = (CC1*THIS%F(I-1) + CC2*THIS%F(I) + CC3*THIS%F(I+1)) / (THIS%Dx**2)
  END DO
  !$omp end parallel do

  !For I = N
   !U_X
   THIS%SD(THIS%IE)  = -(-3.5_AP * THIS%F(THIS%IE) + 4.0_AP * THIS%F(THIS%IE-1) - 0.5_AP * THIS%F(THIS%IE-2)) / THIS%Dx
   !U_XX
   THIS%SSD(THIS%IE) = ( 34.0_AP / 3.0_AP * THIS%F(THIS%IE) - 83.0_AP / 4.0_AP * THIS%F(THIS%IE-1) &
   &  +10.0_AP * THIS%F(THIS%IE-2) -7.0_AP/12.0_AP * THIS%F(THIS%IE-3) ) / (THIS%Dx**2)		 
  
END SUBROUTINE

END MODULE SRKUCCD