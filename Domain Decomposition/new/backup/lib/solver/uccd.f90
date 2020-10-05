MODULE UCCD
USE MATRIX_SOLVER
IMPLICIT NONE

TYPE UCCD_root
INTEGER :: IS, IE
real(8) :: DX, DT
real(8),ALLOCATABLE,DIMENSION(:) :: F, C3
real(8),ALLOCATABLE,DIMENSION(:) :: S, SS 
real(8),ALLOCATABLE,DIMENSION(:) :: SU, SSU, SD, SSD
real(8),ALLOCATABLE,DIMENSION(:,:) :: A, AA, B, BB
real(8),ALLOCATABLE,DIMENSION(:,:) :: AU, AAU, BU, BBU
real(8),ALLOCATABLE,DIMENSION(:,:) :: AD, AAD, BD, BBD
CONTAINS
PROCEDURE INIT => SRKuccd_INITIALIZED
PROCEDURE ASSIGN_u => SRKuccd_ASSIGN_MATRIX_upwind
procedure assign_d => SRKuccd_ASSIGN_MATRIX_downwind
PROCEDURE ASSIGN_SRC => SRKuccd_ASSIGN_SOURCE
PROCEDURE SOLVE => SRKuccd_SOLVE_1D
END TYPE UCCD_root

TYPE UCCD_solver
TYPE(UCCD_root) :: X, Y, Z
contains
procedure alloc => uccd_solver_alloc
END TYPE UCCD_solver

CONTAINS

subroutine uccd_solver_alloc(THIS,IS,IE,JS,JE,KS,KE,DX,DY,DZ,DT)
implicit none
class(uccd_SOLVER) :: this
integer :: IS,IE,JS,JE,KS,KE
real(8) :: dx,dy,dz,dt

CALL THIS%X%init(IS,IE,DX,DT)
CALL THIS%Y%init(JS,JE,DY,DT)
CALL THIS%Z%init(KS,KE,DZ,DT)
	
end subroutine

SUBROUTINE SRKuccd_INITIALIZED(THIS,IS,IE,DX,DT)
IMPLICIT NONE
CLASS(UCCD_root) :: THIS
INTEGER :: IS, IE
real(8) :: DX, dt

THIS%IS = IS; THIS%IE = IE; THIS%DX = DX; THIS%DT = DT

ALLOCATE( THIS%F(IS:IE), THIS%C3(IS:IE) )
ALLOCATE( THIS%SU(IS:IE), THIS%SSU(IS:IE), THIS%SD(IS:IE), THIS%SSD(IS:IE) )
ALLOCATE( THIS%S(IS:IE), THIS%SS(IS:IE) )
ALLOCATE( THIS%AU(3,IS:IE), THIS%AAU(3,IS:IE), THIS%BU(3,IS:IE), THIS%BBU(3,IS:IE) )
ALLOCATE( THIS%AD(3,IS:IE), THIS%AAD(3,IS:IE), THIS%BD(3,IS:IE), THIS%BBD(3,IS:IE) )
ALLOCATE( THIS%A(3,IS:IE), THIS%AA(3,IS:IE), THIS%B(3,IS:IE), THIS%BB(3,IS:IE) )

END SUBROUTINE

SUBROUTINE SRKuccd_SOLVE_1D(THIS,U,F,FX,FXX)
IMPLICIT NONE
CLASS(UCCD_root) :: THIS
real(8),DIMENSION(this%IS:this%IE) :: U,F,FX
real(8),DIMENSION(this%IS:this%IE),optional :: FXX
INTEGER :: I

 DO I = THIS%IS, THIS%IE
   THIS%F(I) = F(I)
   THIS%C3(I) = -0.0609611900811d0
 ENDDO
 
 call this%ASSIGN_u
 call this%ASSIGN_d

 CALL THIS%ASSIGN_SRC
 
 CALL TWIN_DEC( THIS%AU, THIS%BU, THIS%AAU, THIS%BBU, THIS%IS, THIS%IE )
 CALL TWIN_BKS( THIS%AU, THIS%BU, THIS%AAU, THIS%BBU, THIS%SU, THIS%SSU, THIS%IS, THIS%IE )
 
 CALL TWIN_DEC( THIS%AD, THIS%BD, THIS%AAD, THIS%BBD, THIS%IS, THIS%IE )
 CALL TWIN_BKS( THIS%AD, THIS%BD, THIS%AAD, THIS%BBD, THIS%SD, THIS%SSD, THIS%IS, THIS%IE )
 
 DO I = THIS%IS, THIS%IE
    IF( U(I)>=0.0d0 )THEN
	  FX(I) = THIS%SU(I)
	ELSE
	  FX(I) = THIS%SD(I)
	ENDIF
 ENDDO
 
 IF( PRESENT(FXX) )THEN
 	DO I = THIS%IS, THIS%IE
		FXX(I) = 0.5d0*(THIS%SSU(I)+THIS%SSD(I))
 	ENDDO 
 ENDIF

END SUBROUTINE

SUBROUTINE SRKuccd_ASSIGN_MATRIX_upwind(THIS)
IMPLICIT NONE
CLASS(UCCD_root) :: THIS
real(8) :: A1 , B1 , B2 , B3, C3
real(8) :: A3D , B1D , B2D , B3D
real(8) :: AA1 , AA3 , BB1 , BB3
INTEGER :: I

  !For ux eq., a>0
  A1 =  0.875d0
  B1 =  0.125100442885824d0
  B2 = -0.248995571141759d0
  B3 =  0.0001004428858241d0

  !For ux eq., a<0
  A3D =  0.875d0
  B1D = -0.0001004428858241d0
  B2D =  0.248995571141759d0
  B3D = -0.125100442885824d0
  

  !For uxx eq.
  AA1 = -9.0d0 / 8.0d0
  AA3 =  9.0d0 / 8.0d0
  BB1 = -1.0d0 / 8.0d0
  BB3 = -1.0d0 / 8.0d0
  
!For Upwind
  !For I = 1
  !I
   !For U_X
   this%AU(2,THIS%IS) = 1.0d0
   this%BU(2,THIS%IS) = 0.0d0
   !For U_XX
   this%AAU(2,THIS%IS) = 0.0d0
   this%BBU(2,THIS%IS) = 1.0d0
  !I+1
   !For U_X
   this%AU(3,THIS%IS) = 2.0d0
   this%BU(3,THIS%IS) = -this%Dx
   !For U_XX
   this%AAU(3,THIS%IS) = -2.5d0 / this%Dx
   this%BBU(3,THIS%IS) =  8.5d0

  !For I = 2 , N-1
  
  !$omp parallel do private(c3,b1,b2,b3)
  DO I = this%IS+1 , this%IE-1
  
    C3 = THIS%C3(I)!F_C3(THIS%CR(I))
  
  	B1 = 25.0d0 / 192.0d0 + C3 / 12.0d0
	B2 = -19.0d0 / 96.0d0 + 5.0d0 * C3 / 6.0d0
	B3 = 1.0d0 / 192.0d0 + C3 / 12.0d0	
	
  !I-1
   !For U_X
    this%AU(1,I) = A1
    this%BU(1,I) = B1 * this%Dx
   !For U_XX
    this%AAU(1,I) = AA1 / this%Dx
    this%BBU(1,I) = BB1
  !I
    !For U_X
    this%AU(2,I) = 1.0d0
    this%BU(2,I) = B2 * this%Dx
    !For U_XX
    this%AAU(2,I) = 0.0d0
    this%BBU(2,I) = 1.0d0
  !I+1
    !For U_X
    this%Au(3,I) = 0.0d0
    this%Bu(3,I) = B3 * this%Dx
    !For U_XX
    this%AAu(3,I) = AA3 / this%Dx
    this%BBu(3,I) = BB3
	
  END DO
  !$omp end parallel do

  !For I = N
  !I-1
   !For U_X
   this%Au(1,this%ie) = 2.0d0
   this%Bu(1,this%ie) = this%Dx
   !For U_XX
   this%AAu(1,this%ie) = 2.5d0 / THIS%Dx
   this%BBu(1,this%ie) = 8.5d0
  !I
   !For U_X
   this%Au(2,this%ie) = 1.0d0
   this%Bu(2,this%ie) = 0.0d0
   !For U_XX
   this%AAu(2,this%ie) = 0.0d0
   this%BBu(2,this%ie) = 1.0d0  
   
END SUBROUTINE


SUBROUTINE SRKuccd_ASSIGN_MATRIX_downwind(THIS)
IMPLICIT NONE
CLASS(UCCD_root) :: THIS
real(8) :: A1 , B1 , B2 , B3, C3
real(8) :: A3D , B1D , B2D , B3D
real(8) :: AA1 , AA3 , BB1 , BB3
INTEGER :: I

  !For ux eq., a>0
  A1 =  0.875d0
  B1 =  0.125100442885824d0
  B2 = -0.248995571141759d0
  B3 =  0.0001004428858241d0

  !For ux eq., a<0
  A3D =  0.875d0
  B1D = -0.0001004428858241d0
  B2D =  0.248995571141759d0
  B3D = -0.125100442885824d0

  !For uxx eq.
  AA1 = -9.0d0 / 8.0d0
  AA3 =  9.0d0 / 8.0d0
  BB1 = -1.0d0 / 8.0d0
  BB3 = -1.0d0 / 8.0d0
  
!For Downwind
  !For I = 1
  !I
   !For U_X
   this%AD(2,THIS%IS) = 1.0d0
   this%BD(2,THIS%IS) = 0.0d0
   !For U_XX
   this%AAD(2,THIS%IS) = 0.0d0
   this%BBD(2,THIS%IS) = 1.0d0
  !I+1
   !For U_X
   this%AD(3,THIS%IS) = 2.0d0
   this%BD(3,THIS%IS) = -this%Dx
   !For U_XX
   this%AAD(3,THIS%IS) = -2.5d0 / this%Dx
   this%BBD(3,THIS%IS) =  8.5d0

  !For I = 2 , N-1
  
  !$omp parallel do private(c3,b1,b2,b3,b1d,b2d,b3d)
  DO I = this%is-2 , this%ie+2

    C3 = THIS%C3(I)!F_C3(THIS%CR(I))
  
  	B1 = 25.0d0 / 192.0d0 + C3 / 12.0d0
	B2 = -19.0d0 / 96.0d0 + 5.0d0 * C3 / 6.0d0
	B3 = 1.0d0 / 192.0d0 + C3 / 12.0d0	
	
	B1D = -B3
	B2D = -B2
	B3D = -B1
	
  !I-1
   !For U_X
    this%AD(1,I) = 0.0d0
    this%BD(1,I) = B1D * this%Dx
   !For U_XX
    this%AAD(1,I) = AA1 / this%Dx
    this%BBD(1,I) = BB1
	
  !I
    !For U_X
    this%AD(2,I) = 1.0d0
    this%BD(2,I) = B2D * this%Dx
    !For U_XX
    this%AAD(2,I) = 0.0d0
    this%BBD(2,I) = 1.0d0
	
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
   this%AD(1,this%ie) = 2.0d0
   this%BD(1,this%ie) = THIS%Dx
   !For U_XX
   this%AAD(1,this%ie) = 2.5d0 / this%Dx
   this%BBD(1,this%ie) = 8.5d0
  !I
   !For U_X
   this%AD(2,this%ie) = 1.0d0
   this%BD(2,this%ie) = 0.0d0
   !For U_XX
   this%AAD(2,this%ie) = 0.0d0
   this%BBD(2,this%ie) = 1.0d0

END SUBROUTINE

SUBROUTINE SRKuccd_ASSIGN_SOURCE(THIS)
IMPLICIT NONE
CLASS(UCCD_root) :: THIS
real(8) :: C1 , C2 , C3
real(8) :: C1D , C2D , C3D
real(8) :: CC1 , CC2 , CC3
INTEGER:: I

 !For ux eq., a>0
  C1 = -1.93629468537011d0
  C2 =  1.99758937102742d0
  C3 = -0.061294685370110d0

 !For ux eq., a<0
  C1D =  0.061294685370110d0
  C2D = -1.99758937102742d0
  C3D =  1.93629468537011d0

 !For uxx eq.
  CC1 =  3.0d0
  CC2 = -6.0d0
  CC3 =  3.0d0
 
!For upwind

  !For I = 1
   !U_X
    THIS%SU(THIS%IS)  = (-3.5d0 * THIS%F(THIS%IS) + 4.0d0 * THIS%F(THIS%IS+1) - 0.5d0 * THIS%F(THIS%IS+2)) / THIS%Dx
   !U_XX
    THIS%SSU(THIS%IS) = ( 34.0d0 / 3.0d0 * THIS%F(THIS%IS) - 83.0d0 / 4.0d0 * THIS%F(THIS%IS+1) &
          &  +10.0d0 * THIS%F(THIS%IS+2) -7.0d0/12.0d0 * THIS%F(THIS%IS+3) ) / (THIS%Dx**2)

  !For I = 2 , N-1
  !$omp parallel do private(c3,c2,c1)
  DO I = THIS%IS+1 , THIS%IE-1
  
    C3 = THIS%C3(I)
    c2 = -2.0d0*c3 + 15.0d0 / 8.0d0
    c1 = c3 - 15.0d0 / 8.0d0
    
  	!U_X
    THIS%SU(I)  = (C1*THIS%F(I-1) + C2*THIS%F(I) + C3*THIS%F(I+1)) / THIS%Dx
    !U_XX
    THIS%SSU(I) = (CC1*THIS%F(I-1) + CC2*THIS%F(I) + CC3*THIS%F(I+1)) / (THIS%Dx**2)
  END DO
  !$omp end parallel do

  !For I = N
   !U_X
   THIS%SU(THIS%IE)  = -(-3.5d0 * THIS%F(THIS%IE) + 4.0d0 * THIS%F(THIS%IE-1) - 0.5d0 * THIS%F(THIS%IE-2)) / THIS%Dx
   !U_XX
   THIS%SSU(THIS%IE) = ( 34.0d0 / 3.0d0 * THIS%F(THIS%IE) - 83.0d0 / 4.0d0 * THIS%F(THIS%IE-1) &
         &  +10.0d0 * THIS%F(THIS%IE-2) -7.0d0/12.0d0 * THIS%F(THIS%IE-3) ) / (THIS%Dx**2)


!For downwind

  !For I = 1
   !U_X
    THIS%SD(THIS%IS)  = (-3.5d0 * THIS%F(THIS%IS) + 4.0d0 * THIS%F(THIS%IS+1) - 0.5d0 * THIS%F(THIS%IS+2)) / THIS%Dx
   !U_XX
    THIS%SSD(THIS%IS) = ( 34.0d0 / 3.0d0 * THIS%F(THIS%IS) - 83.0d0 / 4.0d0 * THIS%F(THIS%IS+1) &
          &  +10.0d0 * THIS%F(THIS%IS+2) -7.0d0/12.0d0 * THIS%F(THIS%IS+3) ) / (THIS%Dx**2)

  !For I = 2 , N-1
  !$omp parallel do private(c3,c2,c1,c1d,c2d,c3d)
  DO I = THIS%IS+1 , THIS%IE-1
  
    C3 = THIS%C3(I)!F_C3(THIS%CR(I))
    c2 = -2.0d0*c3 + 15.0d0 / 8.0d0
    c1 = c3 - 15.0d0 / 8.0d0
  
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
   THIS%SD(THIS%IE)  = -(-3.5d0 * THIS%F(THIS%IE) + 4.0d0 * THIS%F(THIS%IE-1) - 0.5d0 * THIS%F(THIS%IE-2)) / THIS%Dx
   !U_XX
   THIS%SSD(THIS%IE) = ( 34.0d0 / 3.0d0 * THIS%F(THIS%IE) - 83.0d0 / 4.0d0 * THIS%F(THIS%IE-1) &
   &  +10.0d0 * THIS%F(THIS%IE-2) -7.0d0/12.0d0 * THIS%F(THIS%IE-3) ) / (THIS%Dx**2)		 
  
END SUBROUTINE

END MODULE UCCD