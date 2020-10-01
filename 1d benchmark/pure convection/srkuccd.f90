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
PROCEDURE INIT => SRKuccd_INITIALIZED
PROCEDURE ASSIGN_u => SRKuccd_ASSIGN_MATRIX_upwind
procedure assign_d => SRKuccd_ASSIGN_MATRIX_downwind
PROCEDURE ASSIGN_SRC => SRKuccd_ASSIGN_SOURCE
PROCEDURE SOLVE => SRKuccd_SOLVE_1D
END TYPE SRKUCCD_1D

TYPE SRKUCCD_SOLVER
TYPE(SRKUCCD_1D) :: X, Y, Z
END TYPE SRKUCCD_SOLVER

CONTAINS

FUNCTION srkuccd_C3(X) RESULT(Y)
IMPLICIT NONE
real(kind=ap) :: a,b,c,d,e
real(kind=ap) :: X, Y

 A = -1.52104553243747110000E-01_ap
 B =  1.84486176189603670000E-03_ap
 C =  1.68201035107592190000E-02_ap
 D =  1.04876732091917260000E-02_ap

 Y = A + B*X + C*X**2.0_ap + D*X**3.0_ap
 
END FUNCTION

SUBROUTINE SRKuccd_INITIALIZED(THIS,IS,IE,DX,DT)
IMPLICIT NONE
CLASS(SRKUCCD_1D) :: THIS
INTEGER :: IS, IE
real(kind=ap) :: DX, dt

THIS%IS = IS; THIS%IE = IE; THIS%DX = DX; THIS%DT = DT

ALLOCATE( THIS%F(IS-3:IE+3), THIS%C3(IS-3:IE+3) )
ALLOCATE( THIS%SU(IS-3:IE+3), THIS%SSU(IS-3:IE+3), THIS%SD(IS-3:IE+3), THIS%SSD(IS-3:IE+3) )
ALLOCATE( THIS%S(IS-3:IE+3), THIS%SS(IS-3:IE+3) )
ALLOCATE( THIS%AU(3,IS-3:IE+3), THIS%AAU(3,IS-3:IE+3), THIS%BU(3,IS-3:IE+3), THIS%BBU(3,IS-3:IE+3) )
ALLOCATE( THIS%AD(3,IS-3:IE+3), THIS%AAD(3,IS-3:IE+3), THIS%BD(3,IS-3:IE+3), THIS%BBD(3,IS-3:IE+3) )
ALLOCATE( THIS%A(3,IS-3:IE+3), THIS%AA(3,IS-3:IE+3), THIS%B(3,IS-3:IE+3), THIS%BB(3,IS-3:IE+3) )

END SUBROUTINE

SUBROUTINE SRKuccd_SOLVE_1D(THIS,U,F,FX,FXX)
IMPLICIT NONE
CLASS(SRKUCCD_1D) :: THIS
real(kind=ap),DIMENSION(THIS%IS-3:THIS%IE+3) :: U,F,FX
real(kind=ap),DIMENSION(THIS%IS-3:THIS%IE+3),optional :: FXX
INTEGER :: I

 DO I = THIS%IS-3, THIS%IE+3
   THIS%F(I) = F(I)
   THIS%C3(I) = -0.0609611900811_ap!srkuccd_C3(ABS(U(I))*THIS%DT/THIS%DX)
 ENDDO
 
 call this%ASSIGN_u
 call this%ASSIGN_d

 CALL THIS%ASSIGN_SRC
 
 CALL TWIN_DEC( THIS%AU, THIS%BU, THIS%AAU, THIS%BBU, THIS%IS, THIS%IE )
 CALL TWIN_BKS( THIS%AU, THIS%BU, THIS%AAU, THIS%BBU, THIS%SU, THIS%SSU, THIS%IS, THIS%IE )
 
 CALL TWIN_DEC( THIS%AD, THIS%BD, THIS%AAD, THIS%BBD, THIS%IS, THIS%IE )
 CALL TWIN_BKS( THIS%AD, THIS%BD, THIS%AAD, THIS%BBD, THIS%SD, THIS%SSD, THIS%IS, THIS%IE )
 
 DO I = THIS%IS-3, THIS%IE+3 
    IF( U(I)>=0.0_ap )THEN
	  FX(I) = THIS%SU(I)
	ELSE
	  FX(I) = THIS%SD(I)
	ENDIF
 ENDDO
 
 IF( PRESENT(FXX) )THEN
 	DO I = THIS%IS-3, THIS%IE+3 
		FXX(I) = 0.5_ap*(THIS%SSU(I)+THIS%SSD(I))
 	ENDDO 
 ENDIF

END SUBROUTINE

SUBROUTINE SRKuccd_SOLVE_1D2(THIS,U,F,FX,FXX)
IMPLICIT NONE
CLASS(SRKUCCD_1D) :: THIS
real(kind=ap),DIMENSION(THIS%IS-3:THIS%IE+3) :: U,F,FX
real(kind=ap),DIMENSION(THIS%IS-3:THIS%IE+3),optional :: FXX
INTEGER :: I

 DO I = THIS%IS-3, THIS%IE+3
   THIS%F(I) = F(I)
   THIS%C3(I) = srkuccd_C3(ABS(U(I))*THIS%DT/THIS%DX)
 ENDDO
 
 call this%ASSIGN_u
 call this%ASSIGN_d

 CALL THIS%ASSIGN_SRC

 DO I = THIS%IS-3, THIS%IE+3
 	IF( U(I)>=0.0_ap)THEN
		THIS%A(:,I) = THIS%AU(:,I)
		THIS%AA(:,I) = THIS%AAU(:,I)
		THIS%B(:,I) = THIS%BU(:,I)
		THIS%BB(:,I) = THIS%BBU(:,I)
		THIS%S(I) = THIS%SU(i)
		THIS%SS(I) = THIS%SSU(i)
	ELSE
		THIS%A(:,I) = THIS%AD(:,I)
		THIS%AA(:,I) = THIS%AAD(:,I)
		THIS%B(:,I) = THIS%BD(:,I)
		THIS%BB(:,I) = THIS%BBD(:,I)
		THIS%S(I) = THIS%SD(i)
		THIS%SS(I) = THIS%SSD(i)
	ENDIF
 ENDDO
 
 CALL TWIN_DEC( THIS%A, THIS%B, THIS%AA, THIS%BB, THIS%IS, THIS%IE )
 CALL TWIN_BKS( THIS%A, THIS%B, THIS%AA, THIS%BB, THIS%S, THIS%SS, THIS%IS, THIS%IE )
 
 DO I = THIS%IS-3, THIS%IE+3 
 	FX(I) = THIS%S(I)
 ENDDO
 
 IF( PRESENT(FXX) )THEN
 	DO I = THIS%IS-3, THIS%IE+3 
		FXX(I) = THIS%SS(I)
 	ENDDO 
 ENDIF

END SUBROUTINE


SUBROUTINE SRKuccd_ASSIGN_MATRIX_upwind(THIS)
IMPLICIT NONE
CLASS(SRKUCCD_1D) :: THIS
real(kind=ap) :: A1 , B1 , B2 , B3, C3
real(kind=ap) :: A3D , B1D , B2D , B3D
real(kind=ap) :: AA1 , AA3 , BB1 , BB3
INTEGER :: I

  !For ux eq., a>0
  A1 =  0.875_ap
  B1 =  0.125100442885824_ap
  B2 = -0.248995571141759_ap
  B3 =  0.0001004428858241_ap

  !For ux eq., a<0
  A3D =  0.875_ap
  B1D = -0.0001004428858241_ap
  B2D =  0.248995571141759_ap
  B3D = -0.125100442885824_ap
  

  !For uxx eq.
  AA1 = -9.0_ap / 8.0_ap
  AA3 =  9.0_ap / 8.0_ap
  BB1 = -1.0_ap / 8.0_ap
  BB3 = -1.0_ap / 8.0_ap
  
!For Upwind
  !For I = 1
  !I
   !For U_X
   this%AU(2,THIS%IS-3) = 1.0_ap
   this%BU(2,THIS%IS-3) = 0.0_ap
   !For U_XX
   this%AAU(2,THIS%IS-3) = 0.0_ap
   this%BBU(2,THIS%IS-3) = 1.0_ap
  !I+1
   !For U_X
   this%AU(3,THIS%IS-3) = 2.0_ap
   this%BU(3,THIS%IS-3) = -this%Dx
   !For U_XX
   this%AAU(3,THIS%IS-3) = -2.5_ap / this%Dx
   this%BBU(3,THIS%IS-3) =  8.5_ap

  !For I = 2 , N-1
  
  !$omp parallel do private(c3,b1,b2,b3)
  DO I = this%IS-2 , this%IE+2
  
    C3 = THIS%C3(I)!F_C3(THIS%CR(I))
  
  	B1 = 25.0_ap / 192.0_ap + C3 / 12.0_ap
	B2 = -19.0_ap / 96.0_ap + 5.0_ap * C3 / 6.0_ap
	B3 = 1.0_ap / 192.0_ap + C3 / 12.0_ap	
	
  !I-1
   !For U_X
    this%AU(1,I) = A1
    this%BU(1,I) = B1 * this%Dx
   !For U_XX
    this%AAU(1,I) = AA1 / this%Dx
    this%BBU(1,I) = BB1
  !I
    !For U_X
    this%AU(2,I) = 1.0_ap
    this%BU(2,I) = B2 * this%Dx
    !For U_XX
    this%AAU(2,I) = 0.0_ap
    this%BBU(2,I) = 1.0_ap
  !I+1
    !For U_X
    this%Au(3,I) = 0.0_ap
    this%Bu(3,I) = B3 * this%Dx
    !For U_XX
    this%AAu(3,I) = AA3 / this%Dx
    this%BBu(3,I) = BB3
	
  END DO
  !$omp end parallel do

  !For I = N
  !I-1
   !For U_X
   this%Au(1,this%ie+3) = 2.0_ap
   this%Bu(1,this%ie+3) = this%Dx
   !For U_XX
   this%AAu(1,this%ie+3) = 2.5_ap / THIS%Dx
   this%BBu(1,this%ie+3) = 8.5_ap
  !I
   !For U_X
   this%Au(2,this%ie+3) = 1.0_ap
   this%Bu(2,this%ie+3) = 0.0_ap
   !For U_XX
   this%AAu(2,this%ie+3) = 0.0_ap
   this%BBu(2,this%ie+3) = 1.0_ap  
   
END SUBROUTINE


SUBROUTINE SRKuccd_ASSIGN_MATRIX_downwind(THIS)
IMPLICIT NONE
CLASS(SRKUCCD_1D) :: THIS
real(kind=ap) :: A1 , B1 , B2 , B3, C3
real(kind=ap) :: A3D , B1D , B2D , B3D
real(kind=ap) :: AA1 , AA3 , BB1 , BB3
INTEGER :: I

  !For ux eq., a>0
  A1 =  0.875_ap
  B1 =  0.125100442885824_ap
  B2 = -0.248995571141759_ap
  B3 =  0.0001004428858241_ap

  !For ux eq., a<0
  A3D =  0.875_ap
  B1D = -0.0001004428858241_ap
  B2D =  0.248995571141759_ap
  B3D = -0.125100442885824_ap

  !For uxx eq.
  AA1 = -9.0_ap / 8.0_ap
  AA3 =  9.0_ap / 8.0_ap
  BB1 = -1.0_ap / 8.0_ap
  BB3 = -1.0_ap / 8.0_ap
  
!For Downwind
  !For I = 1
  !I
   !For U_X
   this%AD(2,this%is-3) = 1.0_ap
   this%BD(2,this%is-3) = 0.0_ap
   !For U_XX
   this%AAD(2,this%is-3) = 0.0_ap
   this%BBD(2,this%is-3) = 1.0_ap
  !I+1
   !For U_X
   this%AD(3,this%is-3) = 2.0_ap
   this%BD(3,this%is-3) = -this%Dx
   !For U_XX
   this%AAD(3,this%is-3) = -2.5_ap / this%Dx
   this%BBD(3,this%is-3) =  8.5_ap

  !For I = 2 , N-1
  
  !$omp parallel do private(c3,b1,b2,b3,b1d,b2d,b3d)
  DO I = this%is-2 , this%ie+2

    C3 = THIS%C3(I)!F_C3(THIS%CR(I))
  
  	B1 = 25.0_ap / 192.0_ap + C3 / 12.0_ap
	B2 = -19.0_ap / 96.0_ap + 5.0_ap * C3 / 6.0_ap
	B3 = 1.0_ap / 192.0_ap + C3 / 12.0_ap	
	
	B1D = -B3
	B2D = -B2
	B3D = -B1
	
  !I-1
   !For U_X
    this%AD(1,I) = 0.0_ap
    this%BD(1,I) = B1D * this%Dx
   !For U_XX
    this%AAD(1,I) = AA1 / this%Dx
    this%BBD(1,I) = BB1
	
  !I
    !For U_X
    this%AD(2,I) = 1.0_ap
    this%BD(2,I) = B2D * this%Dx
    !For U_XX
    this%AAD(2,I) = 0.0_ap
    this%BBD(2,I) = 1.0_ap
	
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
   this%AD(1,this%ie+3) = 2.0_ap
   this%BD(1,this%ie+3) = THIS%Dx
   !For U_XX
   this%AAD(1,this%ie+3) = 2.5_ap / this%Dx
   this%BBD(1,this%ie+3) = 8.5_ap
  !I
   !For U_X
   this%AD(2,this%ie+3) = 1.0_ap
   this%BD(2,this%ie+3) = 0.0_ap
   !For U_XX
   this%AAD(2,this%ie+3) = 0.0_ap
   this%BBD(2,this%ie+3) = 1.0_ap

END SUBROUTINE

SUBROUTINE SRKuccd_ASSIGN_SOURCE(THIS)
IMPLICIT NONE
CLASS(SRKUCCD_1D) :: THIS
real(kind=ap) :: C1 , C2 , C3
real(kind=ap) :: C1D , C2D , C3D
real(kind=ap) :: CC1 , CC2 , CC3
INTEGER:: I

 !For ux eq., a>0
  C1 = -1.93629468537011_ap
  C2 =  1.99758937102742_ap
  C3 = -0.061294685370110_ap

 !For ux eq., a<0
  C1D =  0.061294685370110_ap
  C2D = -1.99758937102742_ap
  C3D =  1.93629468537011_ap

 !For uxx eq.
  CC1 =  3.0_ap
  CC2 = -6.0_ap
  CC3 =  3.0_ap
 
!For upwind

  !For I = 1
   !U_X
    THIS%SU(THIS%IS-3)  = (-3.5_ap * THIS%F(THIS%IS-3) + 4.0_ap * THIS%F(THIS%IS-2) - 0.5_ap * THIS%F(THIS%IS-1)) / THIS%Dx
   !U_XX
    THIS%SSU(THIS%IS-3) = ( 34.0_ap / 3.0_ap * THIS%F(THIS%IS-3) - 83.0_ap / 4.0_ap * THIS%F(THIS%IS-2) &
          &  +10.0_ap * THIS%F(THIS%IS-1) -7.0_ap/12.0_ap * THIS%F(THIS%IS) ) / (THIS%Dx**2)

  !For I = 2 , N-1
  !$omp parallel do private(c3,c2,c1)
  DO I = THIS%IS-2 , THIS%IE+2
  
    C3 = THIS%C3(I)!F_C3(THIS%CR(I))
    c2 = -2.0_ap*c3 + 15.0_ap / 8.0_ap
    c1 = c3 - 15.0_ap / 8.0_ap
    
  	!U_X
    THIS%SU(I)  = (C1*THIS%F(I-1) + C2*THIS%F(I) + C3*THIS%F(I+1)) / THIS%Dx
    !U_XX
    THIS%SSU(I) = (CC1*THIS%F(I-1) + CC2*THIS%F(I) + CC3*THIS%F(I+1)) / (THIS%Dx**2)
  END DO
  !$omp end parallel do

  !For I = N
   !U_X
   THIS%SU(THIS%IE+3)  = -(-3.5_ap * THIS%F(THIS%IE+3) + 4.0_ap * THIS%F(THIS%IE+2) - 0.5_ap * THIS%F(THIS%IE+1)) / THIS%Dx
   !U_XX
   THIS%SSU(THIS%IE+3) = ( 34.0_ap / 3.0_ap * THIS%F(THIS%IE+3) - 83.0_ap / 4.0_ap * THIS%F(THIS%IE+2) &
         &  +10.0_ap * THIS%F(THIS%IE+1) -7.0_ap/12.0_ap * THIS%F(THIS%IE) ) / (THIS%Dx**2)


!For downwind

  !For I = 1
   !U_X
    THIS%SD(THIS%IS-3)  = (-3.5_ap * THIS%F(THIS%IS-3) + 4.0_ap * THIS%F(THIS%IS-2) - 0.5_ap * THIS%F(THIS%IS-1)) / THIS%Dx
   !U_XX
    THIS%SSD(THIS%IS-3) = ( 34.0_ap / 3.0_ap * THIS%F(THIS%IS-3) - 83.0_ap / 4.0_ap * THIS%F(THIS%IS-2) &
          &  +10.0_ap * THIS%F(THIS%IS-1) -7.0_ap/12.0_ap * THIS%F(THIS%IS) ) / (THIS%Dx**2)

  !For I = 2 , N-1
  !$omp parallel do private(c3,c2,c1,c1d,c2d,c3d)
  DO I = THIS%IS-2 , THIS%IE+2
  
    C3 = THIS%C3(I)!F_C3(THIS%CR(I))
    c2 = -2.0_ap*c3 + 15.0_ap / 8.0_ap
    c1 = c3 - 15.0_ap / 8.0_ap
  
    c1d = -c3
    c2d = -c2
    c3d = -c1
  
  	!U_X
    THIS%SD(I)  = (C1D*THIS%F(I-1) + C2D*THIS%F(I) + C3D*THIS%F(I+1)) / THIS%Dx
    !U_XX
    THIS%SSD(I) = (CC1*THIS%F(I-1) + CC2*THIS%F(I) + CC3*THIS%F(I+1)) / (THIS%Dx**2)
  END DO
  !$omp end parallel do

  !For I = N
   !U_X
   THIS%SD(THIS%IE+3)  = -(-3.5_ap * THIS%F(THIS%IE+3) + 4.0_ap * THIS%F(THIS%IE+2) - 0.5_ap * THIS%F(THIS%IE+1)) / THIS%Dx
   !U_XX
   THIS%SSD(THIS%IE+3) = ( 34.0_ap / 3.0_ap * THIS%F(THIS%IE+3) - 83.0_ap / 4.0_ap * THIS%F(THIS%IE+2) &
   &  +10.0_ap * THIS%F(THIS%IE+1) -7.0_ap/12.0_ap * THIS%F(THIS%IE) ) / (THIS%Dx**2)		 
  
END SUBROUTINE

END MODULE SRKUCCD