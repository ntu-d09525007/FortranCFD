subroutine weno_bc(f,fp,fm,Nx)
implicit none
integer :: Nx, i
real(8),dimension(-2:Nx+3) :: f, fp, fm
real(8) :: a1,a2,a3,b1,b2,b3,w1,w2,w3,e

e = 1.0d-12

!$omp parallel do private(b1,b2,b3,a1,a2,a3,w1,w2,w3)
do i = 1, Nx, Nx-1
	
	b1 = 13.0*(f(i-2)-2.0*f(i-1)+f(i))**2 + 3.0*(f(i-2)-4*f(i-1)+3.0*f(i))**2
	b2 = 13.0*(f(i-1)-2.0*f(i)+f(i+1))**2 + 3.0*(f(i-1)-f(i+1))**2
	b3 = 13.0*(f(i)-2.0*f(i+1)+f(i+2))**2 + 3.0*(3.0*f(i)-4.0*f(i+1)+f(i+2))**2
	
	a1 = 1.0*(1.0+abs(b3-b1)/(e+b1))
	a2 = 6.0*(1.0+abs(b3-b1)/(e+b2))
	a3 = 3.0*(1.0+abs(b3-b1)/(e+b3))
	
	w1 = a1/(a1+a2+a3)
	w2 = a2/(a1+a2+a3)
	w3 = a3/(a1+a2+a3)
	
	fm(i) = w1/3.0*f(i-2) - (7.0*w1+w2)/6.0*f(i-1) + (11.0*w1+5.0*w2+2.0*w3)/6.0*f(i) &
			+ (2.0*w2+5.0*w3)/6.0*f(i+1) - w3/6.0*f(i+2)
			
end do
!$omp end parallel do 

!$omp parallel do private(b1,b2,b3,a1,a2,a3,w1,w2,w3)	
do i = 1, Nx, Nx-1
	
	b3 = 13.0*(f(i-1)-2.0*f(i)  +f(i+1))**2 + 3.0*(f(i-1)-4.0*f(i)+3.0*f(i+1))**2
	b2 = 13.0*(f(i)  -2.0*f(i+1)+f(i+2))**2 + 3.0*(f(i)-f(i+2))**2
	b1 = 13.0*(f(i+1)-2.0*f(i+2)+f(i+3))**2 + 3.0*(3.0*f(i+1)-4.0*f(i+2)+f(i+3))**2
	
  a1 = 1.0 * ( 1.0 + abs(b1-b3)/(e+b1) )
  a2 = 6.0 * ( 1.0 + abs(b1-b3)/(e+b2) )
  a3 = 3.0 * ( 1.0 + abs(b1-b3)/(e+b3) )	
	
  w1 = a1 / ( a1+a2+a3)
  w2 = a2 / ( a1+a2+a3)
  w3 = a3 / ( a1+a2+a3)
 	
 	fp(i) = w3*(-f(i-1)+5.0*f(i)+2.0*f(i+1))/6.0 &
 				 +w2*(2.0*f(i)+5.0*f(i+1)-f(i+2))/6.0 &
 				 +w1*(11.0*f(i+1)-7.0*f(i+2)+2.0*f(i+3))/6.0	
	
	
end do
!$omp end parallel do 

end subroutine

subroutine crweno(f,fp,fm,Nx)
implicit none
integer :: Nx, i
real(8),dimension(-2:Nx+3) :: f, fp, fm
real(8) :: a1,a2,a3,b1,b2,b3,w1,w2,w3,e,c1,c2,c3
real(8),dimension(2:Nx-1) :: A,B,C,S

e = 1.0e-12

	c1 = 2.0
	c2 = 5.0
	c3 = 3.0

 call weno_bc(f,fp,fm,Nx)

!$omp parallel do private(b1,b2,b3,a1,a2,a3,w1,w2,w3)
do i = 2, Nx-1
	
	b1 = 13.0*(f(i-2)-2.0*f(i-1)+f(i))**2   + 3.0*(    f(i-2)-4*f(i-1)+3.0*f(i))**2
	b2 = 13.0*(f(i-1)-2.0*f(i)  +f(i+1))**2 + 3.0*(    f(i-1)  -f(i+1))**2
	b3 = 13.0*(f(i)  -2.0*f(i+1)+f(i+2))**2 + 3.0*(3.0*f(i)-4.0*f(i+1)+f(i+2))**2
	
	a1 = c1*(  1.0 +  abs(b3-b1) / (e+b1)  )
	a2 = c2*(  1.0 +  abs(b3-b1) / (e+b2)  )
	a3 = c3*(  1.0 +  abs(b3-b1) / (e+b3)  )
	
	w1 = a1/(a1+a2+a3)
	w2 = a2/(a1+a2+a3)
	w3 = a3/(a1+a2+a3)	
	
	A(i) = (2.0*w1+w2)/3.0
	B(i) = (w1+2.0*(w2+w3))/3.0
	C(i) =  w3/3.0
	S(i) = w1/6.0*f(i-1) + (5.0*(w1+w2)+w3)/6.0*f(i) + (w2+5.0*w3)/6.0*f(i+1)
	
end do
!$omp end parallel do 

  S(2) = S(2) - A(2)*fm(1)
  S(Nx-1) = S(Nx-1) - C(Nx-1)*fm(Nx)
  
  
  call solve_tridiagonal(A,B,C,S,fm(2:Nx-1),2,Nx-1)

!$omp parallel do private(b1,b2,b3,a1,a2,a3,w1,w2,w3)  
do i = 2, Nx-1
	
	b3 = 13.0*(f(i-1)-2.0*f(i)  +f(i+1))**2 + 3.0*(f(i-1)-4.0*f(i)+3.0*f(i+1))**2
	b2 = 13.0*(f(i)  -2.0*f(i+1)+f(i+2))**2 + 3.0*(f(i)-f(i+2))**2
	b1 = 13.0*(f(i+1)-2.0*f(i+2)+f(i+3))**2 + 3.0*(3.0*f(i+1)-4.0*f(i+2)+f(i+3))**2
	
	a1 = c1*(1.0+abs(b3-b1)/(e+b1))
	a2 = c2*(1.0+abs(b3-b1)/(e+b2))
	a3 = c3*(1.0+abs(b3-b1)/(e+b3))
	
	w1 = a1/(a1+a2+a3)
	w2 = a2/(a1+a2+a3)
	w3 = a3/(a1+a2+a3)	
	
	A(i) = (w3)/3.0
	B(i) = (w1+2.0*(w2+w3))/3.0
	C(i) = (w2+2.0*w1)/3.0
	S(i) = (5.0*w3+w2)/6.0*f(i) + (w3+5.0*(w2+w1))/6.0*f(i+1) + w1/6.0*f(i+2)
	
end do
!$omp end parallel do 

  S(2) = S(2) - A(2)*fp(1)
  S(Nx-1) = S(Nx-1) - C(Nx-1)*fp(Nx)
  
  call solve_tridiagonal(A,B,C,S,fp(2:Nx-1),2,Nx-1)
  
end subroutine

subroutine ocrweno(f,fp,fm,Nx)
implicit none
integer :: Nx, i
real(8),dimension(-2:Nx+3) :: f, fp, fm
real(8) :: a1,a2,a3,b1,b2,b3,w1,w2,w3,e,c1,c2,c3
real(8),dimension(2:Nx-1) :: A,B,C,S

e = 1.0d-12

	c1 = 0.2089141306
	c2 = 0.4999999998
	c3 = 0.2910858692

 call weno_bc(f,fp,fm,Nx)

!$omp parallel do private(b1,b2,b3,a1,a2,a3,w1,w2,w3)
do i = 2, Nx-1
	
	b1 = 13.0*(f(i-2)-2.0*f(i-1)+f(i))**2   + 3.0*(    f(i-2)-4*f(i-1)+3.0*f(i))**2
	b2 = 13.0*(f(i-1)-2.0*f(i)  +f(i+1))**2 + 3.0*(    f(i-1)  -f(i+1))**2
	b3 = 13.0*(f(i)  -2.0*f(i+1)+f(i+2))**2 + 3.0*(3.0*f(i)-4.0*f(i+1)+f(i+2))**2
	
	a1 = c1*(  1.0 +  abs(b3-b1) / (e+b1)  )
	a2 = c2*(  1.0 +  abs(b3-b1) / (e+b2)  )
	a3 = c3*(  1.0 +  abs(b3-b1) / (e+b3)  )
	
	w1 = a1/(a1+a2+a3)
	w2 = a2/(a1+a2+a3)
	w3 = a3/(a1+a2+a3)	
	
	A(i) = (2.0*w1+w2)/3.0
	B(i) = (w1+2.0*(w2+w3))/3.0
	C(i) =  w3/3.0
	S(i) = w1/6.0*f(i-1) + (5.0*(w1+w2)+w3)/6.0*f(i) + (w2+5.0*w3)/6.0*f(i+1)
	
end do
!$omp end parallel do 

  S(2) = S(2) - A(2)*fm(1)
  S(Nx-1) = S(Nx-1) - C(Nx-1)*fm(Nx)
  
  call solve_tridiagonal(A,B,C,S,fm(2:Nx-1),2,Nx-1)
  
!$omp parallel do private(b1,b2,b3,a1,a2,a3,w1,w2,w3)
do i = 2, Nx-1
	
	b3 = 13.0*(f(i-1)-2.0*f(i)  +f(i+1))**2 + 3.0*(f(i-1)-4.0*f(i)+3.0*f(i+1))**2
	b2 = 13.0*(f(i)  -2.0*f(i+1)+f(i+2))**2 + 3.0*(f(i)-f(i+2))**2
	b1 = 13.0*(f(i+1)-2.0*f(i+2)+f(i+3))**2 + 3.0*(3.0*f(i+1)-4.0*f(i+2)+f(i+3))**2
	
	a1 = c1*(1.0+abs(b3-b1)/(e+b1))
	a2 = c2*(1.0+abs(b3-b1)/(e+b2))
	a3 = c3*(1.0+abs(b3-b1)/(e+b3))
	
	w1 = a1/(a1+a2+a3)
	w2 = a2/(a1+a2+a3)
	w3 = a3/(a1+a2+a3)	
	
	A(i) = (w3)/3.0
	B(i) = (w1+2.0*(w2+w3))/3.0
	C(i) = (w2+2.0*w1)/3.0
	S(i) = (5.0*w3+w2)/6.0*f(i) + (w3+5.0*(w2+w1))/6.0*f(i+1) + w1/6.0*f(i+2)
	
end do
!$omp end parallel do 

  S(2) = S(2) - A(2)*fp(1)
  S(Nx-1) = S(Nx-1) - C(Nx-1)*fp(Nx)
  
  call solve_tridiagonal(A,B,C,S,fp(2:Nx-1),2,Nx-1)
  
end subroutine

subroutine ocrweno_LD(f,fp,fm,Nx)
implicit none
integer :: Nx, i
real(8),dimension(-2:Nx+3) :: f, fp, fm
real(8) :: a1,a2,a3,a4,b1,b2,b3,b4,w1,w2,w3,w4,e,c1,c2,c3,c4
real(8),dimension(2:Nx-1) :: A,B,C,S

e = 1.0d-12

	c1 = 0.09918522298
	c2 = 0.3991852223
	c3 = 0.4008147914
	c4 = 0.1008147633


 call weno_bc(f,fp,fm,Nx)

!$omp parallel do private(b1,b2,b3,b4,a1,a2,a3,a4,w1,w2,w3,w4)
do i = 2, Nx-1
	
	b1 = 13.0*(f(i-2)-2.0*f(i-1)+f(i))**2   + 3.0*(    f(i-2)-4*f(i-1)+3.0*f(i))**2
	b2 = 13.0*(f(i-1)-2.0*f(i)  +f(i+1))**2 + 3.0*(    f(i-1)  -f(i+1))**2
	b3 = 13.0*(f(i)  -2.0*f(i+1)+f(i+2))**2 + 3.0*(3.0*f(i)-4.0*f(i+1)+f(i+2))**2
	b4 = 13.0*(f(i+1)-2.0*f(i+2)+f(i+3))**2 + 3.0*(-5.0*f(i+1)+8.0*f(i+2)-3.0*f(i+3))**2
	b4 = max(b3,b4)
	
	a1 = c1*(  1.0 +  abs(b4-b1) / (e+b1)  )
	a2 = c2*(  1.0 +  abs(b4-b1) / (e+b2)  )
	a3 = c3*(  1.0 +  abs(b4-b1) / (e+b3)  )
	a4 = c4*(  1.0 +  abs(b4-b1) / (e+b4)  )
	
	w1 = a1/(a1+a2+a3+a4)
	w2 = a2/(a1+a2+a3+a4)
	w3 = a3/(a1+a2+a3+a4)	
	w4 = a4/(a1+a2+a3+a4)
	
	A(i) = (2.0*w1+w2)/3.0
	B(i) = (w1+2.0*(w2+w3)+w4)/3.0
	C(i) =  (w3+2.0*w4)/3.0
	S(i) =   w1/6.0*f(i-1) + (5.0*(w1+w2)+w3)/6.0*f(i) + (w2+5.0*(w3+w4))/6.0*f(i+1) + w4/6.0*f(i+2)
	
end do
!$omp end parallel do 

  S(2) = S(2) - A(2)*fm(1)
  S(Nx-1) = S(Nx-1) - C(Nx-1)*fm(Nx)
  
  
  call solve_tridiagonal(A,B,C,S,fm(2:Nx-1),2,Nx-1)

!$omp parallel do private(b1,b2,b3,b4,a1,a2,a3,a4,w1,w2,w3,w4) 
do i = 2, Nx-1
	
	b4 = 13.0*(f(i-2)-2.0*f(i-1)+f(i)  )**2 + 3.0*(-3.0*f(i-2)+8.0*f(i-1)-5.0*f(i))**2
	b3 = 13.0*(f(i-1)-2.0*f(i)  +f(i+1))**2 + 3.0*(f(i-1)-4.0*f(i)+3.0*f(i+1))**2
	b2 = 13.0*(f(i)  -2.0*f(i+1)+f(i+2))**2 + 3.0*(f(i)-f(i+2))**2
	b1 = 13.0*(f(i+1)-2.0*f(i+2)+f(i+3))**2 + 3.0*(3.0*f(i+1)-4.0*f(i+2)+f(i+3))**2
	b4 = max(b3,b4)
	
	a1 = c1*(  1.0 +  abs(b4-b1) / (e+b1)  )
	a2 = c2*(  1.0 +  abs(b4-b1) / (e+b2)  )
	a3 = c3*(  1.0 +  abs(b4-b1) / (e+b3)  )
	a4 = c4*(  1.0 +  abs(b4-b1) / (e+b4)  )
	
	w1 = a1/(a1+a2+a3+a4)
	w2 = a2/(a1+a2+a3+a4)
	w3 = a3/(a1+a2+a3+a4)	
	w4 = a4/(a1+a2+a3+a4)
	
	A(i) = (w3+2.0*w4)/3.0
	B(i) = (w1+2.0*(w2+w3)+w4)/3.0
	C(i) = (w2+2.0*w1)/3.0
	S(i) = w4/6.0*f(i-1) + (5.0*(w3+w4)+w2)/6.0*f(i) + (w3+5.0*(w2+w1))/6.0*f(i+1) + w1/6.0*f(i+2)
	
end do
!$omp end parallel do 

  S(2) = S(2) - A(2)*fp(1)
  S(Nx-1) = S(Nx-1) - C(Nx-1)*fp(Nx)
  
  call solve_tridiagonal(A,B,C,S,fp(2:Nx-1),2,Nx-1)
  
end subroutine

subroutine crweno_LD(f,fp,fm,Nx)
implicit none
integer :: Nx, i
real(8),dimension(-2:Nx+3) :: f, fp, fm
real(8) :: a1,a2,a3,a4,b1,b2,b3,b4,w1,w2,w3,w4,e,c1,c2,c3,c4
real(8),dimension(2:Nx-1) :: A,B,C,S

e = 1.0d-12
	
	c1 = 3.0/20.0
	c2 = 9.0/20.0
	c3 = 7.0/20.0
	c4 = 1.0/20.0

 call weno_bc(f,fp,fm,Nx)

!$omp parallel do private(b1,b2,b3,b4,a1,a2,a3,a4,w1,w2,w3,w4)
do i = 2, Nx-1
	
	b1 = 13.0*(f(i-2)-2.0*f(i-1)+f(i))**2   + 3.0*(    f(i-2)-4*f(i-1)+3.0*f(i))**2
	b2 = 13.0*(f(i-1)-2.0*f(i)  +f(i+1))**2 + 3.0*(    f(i-1)  -f(i+1))**2
	b3 = 13.0*(f(i)  -2.0*f(i+1)+f(i+2))**2 + 3.0*(3.0*f(i)-4.0*f(i+1)+f(i+2))**2
	b4 = 13.0*(f(i+1)-2.0*f(i+2)+f(i+3))**2 + 3.0*(-5.0*f(i+1)+8.0*f(i+2)-3.0*f(i+3))**2
	b4 = max(b3,b4)
	
	a1 = c1*(  1.0 +  abs(b4-b1) / (e+b1)  )
	a2 = c2*(  1.0 +  abs(b4-b1) / (e+b2)  )
	a3 = c3*(  1.0 +  abs(b4-b1) / (e+b3)  )
	a4 = c4*(  1.0 +  abs(b4-b1) / (e+b4)  )
	
	w1 = a1/(a1+a2+a3+a4)
	w2 = a2/(a1+a2+a3+a4)
	w3 = a3/(a1+a2+a3+a4)	
	w4 = a4/(a1+a2+a3+a4)
	
	A(i) = (2.0*w1+w2)/3.0
	B(i) = (w1+2.0*(w2+w3)+w4)/3.0
	C(i) =  (w3+2.0*w4)/3.0
	S(i) =   w1/6.0*f(i-1) + (5.0*(w1+w2)+w3)/6.0*f(i) + (w2+5.0*(w3+w4))/6.0*f(i+1) + w4/6.0*f(i+2)
	
end do
!$omp end parallel do 

  S(2) = S(2) - A(2)*fm(1)
  S(Nx-1) = S(Nx-1) - C(Nx-1)*fm(Nx)
  
  
  call solve_tridiagonal(A,B,C,S,fm(2:Nx-1),2,Nx-1)

!$omp parallel do private(b1,b2,b3,b4,a1,a2,a3,a4,w1,w2,w3,w4) 
do i = 2, Nx-1
	
	b4 = 13.0*(f(i-2)-2.0*f(i-1)+f(i)  )**2 + 3.0*(-3.0*f(i-2)+8.0*f(i-1)-5.0*f(i))**2
	b3 = 13.0*(f(i-1)-2.0*f(i)  +f(i+1))**2 + 3.0*(f(i-1)-4.0*f(i)+3.0*f(i+1))**2
	b2 = 13.0*(f(i)  -2.0*f(i+1)+f(i+2))**2 + 3.0*(f(i)-f(i+2))**2
	b1 = 13.0*(f(i+1)-2.0*f(i+2)+f(i+3))**2 + 3.0*(3.0*f(i+1)-4.0*f(i+2)+f(i+3))**2
	b4 = max(b3,b4)
	
	a1 = c1*(  1.0 +  abs(b4-b1) / (e+b1)  )
	a2 = c2*(  1.0 +  abs(b4-b1) / (e+b2)  )
	a3 = c3*(  1.0 +  abs(b4-b1) / (e+b3)  )
	a4 = c4*(  1.0 +  abs(b4-b1) / (e+b4)  )
	
	w1 = a1/(a1+a2+a3+a4)
	w2 = a2/(a1+a2+a3+a4)
	w3 = a3/(a1+a2+a3+a4)	
	w4 = a4/(a1+a2+a3+a4)
	
	A(i) = (w3+2.0*w4)/3.0
	B(i) = (w1+2.0*(w2+w3)+w4)/3.0
	C(i) = (w2+2.0*w1)/3.0
	S(i) = w4/6.0*f(i-1) + (5.0*(w3+w4)+w2)/6.0*f(i) + (w3+5.0*(w2+w1))/6.0*f(i+1) + w1/6.0*f(i+2)
	
end do
!$omp end parallel do 

  S(2) = S(2) - A(2)*fp(1)
  S(Nx-1) = S(Nx-1) - C(Nx-1)*fp(Nx)
  
  call solve_tridiagonal(A,B,C,S,fp(2:Nx-1),2,Nx-1)
  
end subroutine
