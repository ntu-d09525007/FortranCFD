program main
implicit none
integer :: i,n,j
real(16) :: dx, xstart, fx, fx_old, factor
real(16),dimension(:),allocatable :: f,x

n = 40
n = 2*n+1

dx=1.0q-3
factor=0.5q0

allocate( x(n), f(n) )

do j = 1, 5

dx=dx*factor

xstart = -(n-1)*dx/2.0q0

do i = 1, n
    x(i) = xstart + (i-1)*dx
    f(i) = x(i)**3.0q0 + qcos(x(i))
end do

!call ocrweno_ld(f,fx,N,dx)
!call crweno_ld(f,fx,N,dx)
call crweno(f,fx,N,dx)

if(j==1)then
	write(*,'(ES20.5,ES20.7)')dx,abs(fx)
else
	write(*,'(ES20.5,ES20.7,F8.4)')dx,abs(fx),qlog(abs(fx)/abs(fx_old))/qlog(factor)
endif

fx_old=fx

end do

contains

subroutine weno_js(f,fp,fm,Nx)
implicit none
integer :: Nx, i, N
real(16),dimension(Nx) :: f, fp, fm
real(16) :: a1,a2,a3,b1,b2,b3,w1,w2,w3,fx,dx,eps


eps = 1.0q-40

N = (Nx-1)/2+1

!$OMP PARALLEL dO PRIVATE(B1,B2,B3,A1,A2,A3,W1,W2,W3)
do i = 3, Nx-2
	
	b1 = 13.0q0*(f(i-2)-2.0q0*f(i-1)+f(i))**2.0q0 + 3.0q0*(f(i-2)-4*f(i-1)+3.0*f(i))**2.0q0
	b2 = 13.0q0*(f(i-1)-2.0q0*f(i)+f(i+1))**2.0q0 + 3.0q0*(f(i-1)-f(i+1))**2.0q0
	b3 = 13.0q0*(f(i)-2.0q0*f(i+1)+f(i+2))**2.0q0 + 3.0q0*(3.0q0*f(i)-4.0*f(i+1)+f(i+2))**2.0q0
	
	a1 = 1.0q0/(EPS+b1)**2.0q0
	a2 = 6.0q0/(EPS+b2)**2.0q0
	a3 = 3.0q0/(EPS+b3)**2.0q0
	
	w1 = a1/(a1+a2+a3)
	w2 = a2/(a1+a2+a3)
	w3 = a3/(a1+a2+a3)
	
	fm(i) = w1/3.0q0*f(i-2) - (7.0q0*w1+w2)/6.0q0*f(i-1) + (11.0q0*w1+5.0q0*w2+2.0q0*w3)/6.0q0*f(i) &
			+ (2.0q0*w2+5.0q0*w3)/6.0q0*f(i+1) - w3/6.0q0*f(i+2)
			
end do
!$OMP ENd PARALLEL dO

!$OMP PARALLEL dO PRIVATE(B1,B2,B3,A1,A2,A3,W1,W2,W3)	
do i = 2, Nx-3
	
	b3 = 13.0q0*(f(i-1)-2.0q0*f(i)  +f(i+1))**2.0q0 + 3.0q0*(f(i-1)-4.0q0*f(i)+3.0q0*f(i+1))**2
	b2 = 13.0q0*(f(i)  -2.0q0*f(i+1)+f(i+2))**2.0q0 + 3.0q0*(f(i)-f(i+2))**2
	b1 = 13.0q0*(f(i+1)-2.0q0*f(i+2)+f(i+3))**2.0q0 + 3.0q0*(3.0*f(i+1)-4.0q0*f(i+2)+f(i+3))**2
	
    a1 = 1.0q0/(EPS+b1)**2
    a2 = 6.0q0/(EPS+b2)**2
    a3 = 3.0q0/(EPS+b3)**2	
	
    w1 = a1 / (a1+a2+a3)
    w2 = a2 / (a1+a2+a3)
    w3 = a3 / (a1+a2+a3)
 	
 	fp(i) = w3*(-f(i-1)+5.0q0*f(i)+2.0q0*f(i+1))/6.0q0 &
 				 +w2*(2.0q0*f(i)+5.0q0*f(i+1)-f(i+2))/6.0q0 &
 				 +w1*(11.0q0*f(i+1)-7.0q0*f(i+2)+2.0q0*f(i+3))/6.0q0	
	
	
end do
!$OMP ENd PARALLEL dO


fx = 0.5q0*(fm(N)+fp(N)-fm(N-1)-fp(N-1))/dx
!fx = (fp(n)-fp(n-1))/dx
!fx = (fm(n)-fm(n-1))/dx
!fx = 0.5q0*( abs((fp(n)-fp(n-1))/dx)+abs((fm(n)-fm(n-1))/dx) )


end subroutine

subroutine crweno(f,fx,Nx,dx)
implicit none
integer :: Nx, i, N
real(16),dimension(Nx) :: f, fp, fm
real(16) :: a1,a2,a3,a4,b1,b2,b3,b4,w1,w2,w3,w4,e,c1,c2,c3,c4,fx,dx
real(16),dimension(2:Nx-2) :: A,B,C,S

	e = 1.0q-40

	N = (Nx-1)/2+1

	!c1 = 0.09918522298q0
	!c2 = 0.3991852223q0
	!c3 = 0.4008147914q0
	!c4 = 0.1008147633q0
	
	c1 = 2.0q0
	c2 = 5.0q0
	c3 = 3.0q0

 call weno_js(f,fp,fm,Nx)
 
 ! p : 3, n-2
 ! m : 2, n-3

!$omp parallel do private(b1,b2,b3,b4,a1,a2,a3,a4,w1,w2,w3,w4)
do i = 4, Nx-3
	
	b1 = 13.0q0*(f(i-2)-2.0q0*f(i-1)+f(i))**2.0q0   + 3.0q0*(    f(i-2)-4.0q0*f(i-1)+3.0q0*f(i))**2.0q0
	b2 = 13.0q0*(f(i-1)-2.0q0*f(i)  +f(i+1))**2.0q0 + 3.0q0*(    f(i-1)  -f(i+1))**2.0q0
	b3 = 13.0q0*(f(i)  -2.0q0*f(i+1)+f(i+2))**2.0q0 + 3.0q0*(3.0q0*f(i)-4.0q0*f(i+1)+f(i+2))**2.0q0
	!b4 = 13.0q0*(f(i+1)-2.0q0*f(i+2)+f(i+3))**2.0q0 + 3.0q0*(-5.0q0*f(i+1)+8.0q0*f(i+2)-3.0q0*f(i+3))**2.0q0
	!b4 = max(b3,b4)
	
	a1 = c1*(  1.0q0 +  abs(b3-b1) / (e+b1)  )
	a2 = c2*(  1.0q0 +  abs(b3-b1) / (e+b2)  )
	a3 = c3*(  1.0q0 +  abs(b3-b1) / (e+b3)  )
	!a4 = c4*(  1.0q0 +  abs(b4-b1) / (e+b4)  )
	
	w1 = a1/(a1+a2+a3+a4)
	w2 = a2/(a1+a2+a3+a4)
	w3 = a3/(a1+a2+a3+a4)	
	w4 = 0.0q0!a4/(a1+a2+a3+a4)
	
	A(i) = (2.0q0*w1+w2)/3.0q0
	B(i) = (w1+2.0q0*(w2+w3)+w4)/3.0q0
	C(i) =  (w3+2.0q0*w4)/3.0q0
	S(i) =   w1/6.0q0*f(i-1) + (5.0q0*(w1+w2)+w3)/6.0q0*f(i) + (w2+5.0q0*(w3+w4))/6.0q0*f(i+1) + w4/6.0q0*f(i+2)
	
end do
!$omp end parallel do 

  S(4) = S(4) - A(4)*fm(3)
  S(Nx-3) = S(Nx-3) - C(Nx-3)*fm(Nx-2)
   
  call solve_tridiagonal(A(4:Nx-3),B(4:Nx-3),C(4:Nx-3),S(4:Nx-3),fm(4:Nx-3),4,Nx-3)

!$omp parallel do private(b1,b2,b3,b4,a1,a2,a3,a4,w1,w2,w3,w4) 
do i = 3, Nx-4
	
	!b4 = 13.0q0*(f(i-2)-2.0q0*f(i-1)+f(i)  )**2.0q0 + 3.0q0*(-3.0q0*f(i-2)+8.0q0*f(i-1)-5.0q0*f(i))**2.0q0
	b3 = 13.0q0*(f(i-1)-2.0q0*f(i)  +f(i+1))**2.0q0 + 3.0q0*(f(i-1)-4.0q0*f(i)+3.0q0*f(i+1))**2.0q0
	b2 = 13.0q0*(f(i)  -2.0q0*f(i+1)+f(i+2))**2.0q0 + 3.0q0*(f(i)-f(i+2))**2.0q0
	b1 = 13.0q0*(f(i+1)-2.0q0*f(i+2)+f(i+3))**2.0q0 + 3.0q0*(3.0q0*f(i+1)-4.0q0*f(i+2)+f(i+3))**2.0q0
	!b4 = max(b3,b4)
	
	a1 = c1*(  1.0q0 +  abs(b4-b1) / (e+b1)  )
	a2 = c2*(  1.0q0 +  abs(b4-b1) / (e+b2)  )
	a3 = c3*(  1.0q0 +  abs(b4-b1) / (e+b3)  )
	!a4 = c4*(  1.0q0 +  abs(b4-b1) / (e+b4)  )
	
	w1 = a1/(a1+a2+a3+a4)
	w2 = a2/(a1+a2+a3+a4)
	w3 = a3/(a1+a2+a3+a4)	
	w4 = 0.0q0!a4/(a1+a2+a3+a4)
	
	A(i) = (w3+2.0q0*w4)/3.0q0
	B(i) = (w1+2.0q0*(w2+w3)+w4)/3.0q0
	C(i) = (w2+2.0q0*w1)/3.0q0
	S(i) = w4/6.0q0*f(i-1) + (5.0q0*(w3+w4)+w2)/6.0q0*f(i) + (w3+5.0q0*(w2+w1))/6.0q0*f(i+1) + w1/6.0q0*f(i+2)
	
end do
!$omp end parallel do 

  S(3) = S(3) - A(3)*fp(2)
  S(Nx-4) = S(Nx-4) - C(Nx-4)*fp(Nx-3)
  
  call solve_tridiagonal(A(3:Nx-4),B(3:Nx-4),C(3:Nx-4),S(3:Nx-4),fp(3:Nx-4),3,Nx-4)

	fx = 0.5q0*(fm(N)+fp(N)-fm(N-1)-fp(N-1))/dx

  
end subroutine

subroutine crweno_ld(f,fx,Nx,dx)
implicit none
integer :: Nx, i, N
real(16),dimension(Nx) :: f, fp, fm
real(16) :: a1,a2,a3,a4,b1,b2,b3,b4,w1,w2,w3,w4,e,c1,c2,c3,c4,fx,dx
real(16),dimension(2:Nx-2) :: A,B,C,S

	e = 1.0q-40

	N = (Nx-1)/2+1

	!c1 = 0.09918522298q0
	!c2 = 0.3991852223q0
	!c3 = 0.4008147914q0
	!c4 = 0.1008147633q0
	
	c1 = 3.0q0/20.0q0
	c2 = 9.0q0/20.0q0
	c3 = 7.0q0/20.0q0
	c4 = 1.0q0/20.0q0

 call weno_js(f,fp,fm,Nx)
 
 ! p : 3, n-2
 ! m : 2, n-3

!$omp parallel do private(b1,b2,b3,b4,a1,a2,a3,a4,w1,w2,w3,w4)
do i = 4, Nx-3
	
	b1 = 13.0q0*(f(i-2)-2.0q0*f(i-1)+f(i))**2.0q0   + 3.0q0*(    f(i-2)-4.0q0*f(i-1)+3.0q0*f(i))**2.0q0
	b2 = 13.0q0*(f(i-1)-2.0q0*f(i)  +f(i+1))**2.0q0 + 3.0q0*(    f(i-1)  -f(i+1))**2.0q0
	b3 = 13.0q0*(f(i)  -2.0q0*f(i+1)+f(i+2))**2.0q0 + 3.0q0*(3.0q0*f(i)-4.0q0*f(i+1)+f(i+2))**2.0q0
	b4 = 13.0q0*(f(i+1)-2.0q0*f(i+2)+f(i+3))**2.0q0 + 3.0q0*(-5.0q0*f(i+1)+8.0q0*f(i+2)-3.0q0*f(i+3))**2.0q0
	b4 = max(b3,b4)
	
	a1 = c1*(  1.0q0 +  abs(b4-b1) / (e+b1)  )
	a2 = c2*(  1.0q0 +  abs(b4-b1) / (e+b2)  )
	a3 = c3*(  1.0q0 +  abs(b4-b1) / (e+b3)  )
	a4 = c4*(  1.0q0 +  abs(b4-b1) / (e+b4)  )
	
	w1 = a1/(a1+a2+a3+a4)
	w2 = a2/(a1+a2+a3+a4)
	w3 = a3/(a1+a2+a3+a4)	
	w4 = a4/(a1+a2+a3+a4)
	
	A(i) = (2.0q0*w1+w2)/3.0q0
	B(i) = (w1+2.0q0*(w2+w3)+w4)/3.0q0
	C(i) =  (w3+2.0q0*w4)/3.0q0
	S(i) =   w1/6.0q0*f(i-1) + (5.0q0*(w1+w2)+w3)/6.0q0*f(i) + (w2+5.0q0*(w3+w4))/6.0q0*f(i+1) + w4/6.0q0*f(i+2)
	
end do
!$omp end parallel do 

  S(4) = S(4) - A(4)*fm(3)
  S(Nx-3) = S(Nx-3) - C(Nx-3)*fm(Nx-2)
   
  call solve_tridiagonal(A(4:Nx-3),B(4:Nx-3),C(4:Nx-3),S(4:Nx-3),fm(4:Nx-3),4,Nx-3)

!$omp parallel do private(b1,b2,b3,b4,a1,a2,a3,a4,w1,w2,w3,w4) 
do i = 3, Nx-4
	
	b4 = 13.0q0*(f(i-2)-2.0q0*f(i-1)+f(i)  )**2.0q0 + 3.0q0*(-3.0q0*f(i-2)+8.0q0*f(i-1)-5.0q0*f(i))**2.0q0
	b3 = 13.0q0*(f(i-1)-2.0q0*f(i)  +f(i+1))**2.0q0 + 3.0q0*(f(i-1)-4.0q0*f(i)+3.0q0*f(i+1))**2.0q0
	b2 = 13.0q0*(f(i)  -2.0q0*f(i+1)+f(i+2))**2.0q0 + 3.0q0*(f(i)-f(i+2))**2.0q0
	b1 = 13.0q0*(f(i+1)-2.0q0*f(i+2)+f(i+3))**2.0q0 + 3.0q0*(3.0q0*f(i+1)-4.0q0*f(i+2)+f(i+3))**2.0q0
	b4 = max(b3,b4)
	
	a1 = c1*(  1.0q0 +  abs(b4-b1) / (e+b1)  )
	a2 = c2*(  1.0q0 +  abs(b4-b1) / (e+b2)  )
	a3 = c3*(  1.0q0 +  abs(b4-b1) / (e+b3)  )
	a4 = c4*(  1.0q0 +  abs(b4-b1) / (e+b4)  )
	
	w1 = a1/(a1+a2+a3+a4)
	w2 = a2/(a1+a2+a3+a4)
	w3 = a3/(a1+a2+a3+a4)	
	w4 = a4/(a1+a2+a3+a4)
	
	A(i) = (w3+2.0q0*w4)/3.0q0
	B(i) = (w1+2.0q0*(w2+w3)+w4)/3.0q0
	C(i) = (w2+2.0q0*w1)/3.0q0
	S(i) = w4/6.0q0*f(i-1) + (5.0q0*(w3+w4)+w2)/6.0q0*f(i) + (w3+5.0q0*(w2+w1))/6.0q0*f(i+1) + w1/6.0q0*f(i+2)
	
end do
!$omp end parallel do 

  S(3) = S(3) - A(3)*fp(2)
  S(Nx-4) = S(Nx-4) - C(Nx-4)*fp(Nx-3)
  
  call solve_tridiagonal(A(3:Nx-4),B(3:Nx-4),C(3:Nx-4),S(3:Nx-4),fp(3:Nx-4),3,Nx-4)

	fx = 0.5q0*(fm(N)+fp(N)-fm(N-1)-fp(N-1))/dx

  
end subroutine

subroutine ocrweno_ld(f,fx,Nx,dx)
implicit none
integer :: Nx, i, N
real(16),dimension(Nx) :: f, fp, fm
real(16) :: a1,a2,a3,a4,b1,b2,b3,b4,w1,w2,w3,w4,e,c1,c2,c3,c4,fx,dx
real(16),dimension(2:Nx-2) :: A,B,C,S

	e = 1.0q-40

	N = (Nx-1)/2+1

	c1 = 0.09918522298q0
	c2 = 0.3991852223q0
	c3 = 0.4008147914q0
	c4 = 0.1008147633q0

 call weno_js(f,fp,fm,Nx)
 
 ! p : 3, n-2
 ! m : 2, n-3

!$omp parallel do private(b1,b2,b3,b4,a1,a2,a3,a4,w1,w2,w3,w4)
do i = 4, Nx-3
	
	b1 = 13.0q0*(f(i-2)-2.0q0*f(i-1)+f(i))**2.0q0   + 3.0q0*(    f(i-2)-4.0q0*f(i-1)+3.0q0*f(i))**2.0q0
	b2 = 13.0q0*(f(i-1)-2.0q0*f(i)  +f(i+1))**2.0q0 + 3.0q0*(    f(i-1)  -f(i+1))**2.0q0
	b3 = 13.0q0*(f(i)  -2.0q0*f(i+1)+f(i+2))**2.0q0 + 3.0q0*(3.0q0*f(i)-4.0q0*f(i+1)+f(i+2))**2.0q0
	b4 = 13.0q0*(f(i+1)-2.0q0*f(i+2)+f(i+3))**2.0q0 + 3.0q0*(-5.0q0*f(i+1)+8.0q0*f(i+2)-3.0q0*f(i+3))**2.0q0
	b4 = max(b3,b4)
	
	a1 = c1*(  1.0q0 +  abs(b4-b1) / (e+b1)  )
	a2 = c2*(  1.0q0 +  abs(b4-b1) / (e+b2)  )
	a3 = c3*(  1.0q0 +  abs(b4-b1) / (e+b3)  )
	a4 = c4*(  1.0q0 +  abs(b4-b1) / (e+b4)  )
	
	w1 = a1/(a1+a2+a3+a4)
	w2 = a2/(a1+a2+a3+a4)
	w3 = a3/(a1+a2+a3+a4)	
	w4 = a4/(a1+a2+a3+a4)
	
	A(i) = (2.0q0*w1+w2)/3.0q0
	B(i) = (w1+2.0q0*(w2+w3)+w4)/3.0q0
	C(i) =  (w3+2.0q0*w4)/3.0q0
	S(i) =   w1/6.0q0*f(i-1) + (5.0q0*(w1+w2)+w3)/6.0q0*f(i) + (w2+5.0q0*(w3+w4))/6.0q0*f(i+1) + w4/6.0q0*f(i+2)
	
end do
!$omp end parallel do 

  S(4) = S(4) - A(4)*fm(3)
  S(Nx-3) = S(Nx-3) - C(Nx-3)*fm(Nx-2)
   
  call solve_tridiagonal(A(4:Nx-3),B(4:Nx-3),C(4:Nx-3),S(4:Nx-3),fm(4:Nx-3),4,Nx-3)

!$omp parallel do private(b1,b2,b3,b4,a1,a2,a3,a4,w1,w2,w3,w4) 
do i = 3, Nx-4
	
	b4 = 13.0q0*(f(i-2)-2.0q0*f(i-1)+f(i)  )**2.0q0 + 3.0q0*(-3.0q0*f(i-2)+8.0q0*f(i-1)-5.0q0*f(i))**2.0q0
	b3 = 13.0q0*(f(i-1)-2.0q0*f(i)  +f(i+1))**2.0q0 + 3.0q0*(f(i-1)-4.0q0*f(i)+3.0q0*f(i+1))**2.0q0
	b2 = 13.0q0*(f(i)  -2.0q0*f(i+1)+f(i+2))**2.0q0 + 3.0q0*(f(i)-f(i+2))**2.0q0
	b1 = 13.0q0*(f(i+1)-2.0q0*f(i+2)+f(i+3))**2.0q0 + 3.0q0*(3.0q0*f(i+1)-4.0q0*f(i+2)+f(i+3))**2.0q0
	b4 = max(b3,b4)
	
	a1 = c1*(  1.0q0 +  abs(b4-b1) / (e+b1)  )
	a2 = c2*(  1.0q0 +  abs(b4-b1) / (e+b2)  )
	a3 = c3*(  1.0q0 +  abs(b4-b1) / (e+b3)  )
	a4 = c4*(  1.0q0 +  abs(b4-b1) / (e+b4)  )
	
	w1 = a1/(a1+a2+a3+a4)
	w2 = a2/(a1+a2+a3+a4)
	w3 = a3/(a1+a2+a3+a4)	
	w4 = a4/(a1+a2+a3+a4)
	
	A(i) = (w3+2.0q0*w4)/3.0q0
	B(i) = (w1+2.0q0*(w2+w3)+w4)/3.0q0
	C(i) = (w2+2.0q0*w1)/3.0q0
	S(i) = w4/6.0q0*f(i-1) + (5.0q0*(w3+w4)+w2)/6.0q0*f(i) + (w3+5.0q0*(w2+w1))/6.0q0*f(i+1) + w1/6.0q0*f(i+2)
	
end do
!$omp end parallel do 

  S(3) = S(3) - A(3)*fp(2)
  S(Nx-4) = S(Nx-4) - C(Nx-4)*fp(Nx-3)
  
  call solve_tridiagonal(A(3:Nx-4),B(3:Nx-4),C(3:Nx-4),S(3:Nx-4),fp(3:Nx-4),3,Nx-4)

	fx = 0.5q0*(fm(N)+fp(N)-fm(N-1)-fp(N-1))/dx

  
end subroutine

 subroutine solve_tridiagonal(AA,BB,CC,SS,X,M,N)
!cccccccccccccccccccccccccccccccccc
!
! A, first  coefficient matrix
! B, second coefficient matrix
! C, third  coefficient matrix
! S, sorce matrix
! X, solution
! [M,N], index domain
!
!cccccccccccccccccccccccccccccccccc
implicit none
integer :: M, N, i
real(kind=16),dimension(m:n) :: AA,BB,CC,SS
real(kind=16),dimension(m:n) :: A,B,C,S,X
 
 A = AA
 B = BB
 C = CC
 S = SS
 
 C(M) = C(M)/ B(M)
 S(M) = S(M)/ B(M)
 
 do i = M+1, N, 1
	 
	 C(i) = C(i) / ( B(I)-A(I)*C(I-1) )
	 S(I) = (S(I) - A(i)*S(I-1))/(B(i)-A(I)*C(i-1))

 end do
 
 X(N) = S(N)
 
 do i = N-1, M, -1
	 
	 X(i) = S(i) - C(i)*X(i+1)	 
	 
 end do


end subroutine


end program
