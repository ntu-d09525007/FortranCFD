module srk6_solver
implicit none

type srk6_roots
integer :: is, ie, js, je, ghc
real(8) :: dt, w
real(8),dimension(:,:),allocatable :: target
real(8),dimension(:,:),allocatable :: s1,s2,s3,ss1,ss2,ss3,l1,l2,l3
contains
procedure alloc => srk6_roots_alloc
procedure init => srk6_roots_init
procedure solve => srk6_roots_solve
procedure final => srk6_roots_final
end type srk6_roots

type srk6_data
integer :: is, ie, js, je
logical :: is_vector_solver
type(srk6_roots) :: x, y
contains
procedure alloc => srk6_data_alloc
procedure init => srk6_data_init
procedure solve => srk6_data_solve
procedure final => srk6_data_final
end type srk6_data

contains

subroutine srk6_roots_alloc(p,is,ie,js,je,dt,w,ghc)
implicit none
class(srk6_roots) :: p
integer, intent(in) :: is, ie, js, je, ghc
real(8), intent(in) :: dt, w

p%is = is; p%ie = ie
p%js = js; p%je = je
p%dt = dt; p%w = w
p%ghc = ghc

 allocate( p%s1(is:ie,js:je), p%s2(is:ie,js:je), p%s3(is:ie,js:je))
 allocate( p%ss1(is:ie,js:je), p%ss2(is:ie,js:je), p%ss3(is:ie,js:je))
 allocate( p%l1(is:ie,js:je), p%l2(is:ie,js:je), p%l3(is:ie,js:je))
 allocate( p%target(is:ie,js:je) )

end subroutine

subroutine srk6_data_alloc(p,is,ie,js,je,dt,w,ghc)
implicit none
class(srk6_data) :: p
integer, intent(in) :: is, ie, js, je, ghc
real(8), intent(in) :: dt, w

p%is = is; p%ie = ie
p%js = js; p%je = je

call p%x%alloc(is,ie,js,je,dt,w,ghc)
call p%y%alloc(is,ie,js,je,dt,w,ghc)

end subroutine

subroutine srk6_roots_init(p,u)
implicit none
class(srk6_roots) :: p
real(8),dimension(p%is:p%ie,p%js:p%je) :: u
integer :: i,j

 do j = p%js, p%je
 do i = p%is, p%ie
	p%target(i,j) = u(i,j)
	p%s1(i,j) = u(i,j)
	p%s2(i,j) = u(i,j)
	p%s3(i,j) = u(i,j)
 enddo
 enddo

end subroutine

subroutine srk6_data_init(p,btn,x,y)
implicit none
class(srk6_data) :: p
logical :: btn
real(8),dimension(p%is:p%ie,p%js:p%je) :: x
real(8),dimension(p%is:p%ie,p%js:p%je),optional :: y
integer :: i,j
 
p%is_vector_solver = btn

call p%x%init(x)

if( p%is_vector_solver ) then
	call p%y%init(y)
endif

end subroutine

subroutine srk6_roots_solve(p,err)
implicit none
class(srk6_roots) :: p
integer :: i,j
real(8),intent(out) :: err
real(8) :: coef1_1, coef1_2, coef1_3
real(8) :: coef2_1, coef2_2, coef2_3
real(8) :: coef3_1, coef3_2, coef3_3
real(8) :: cc

cc = dsqrt(0.6_8)*0.5_8

coef1_1 = 5.0_8/36.0_8; coef1_2 = 2.0_8/9.0_8 + 2.0_8*cc/3.0_8; coef1_3 = 5.0_8/36.0_8 + cc/3.0_8
coef2_1 = 5.0_8/36.0_8 - 5.0_8*cc/12.0_8; coef2_2 = 2.0_8/9.0_8; coef2_3 = 5.0_8/36.0_8 + 5.0_8*cc/12.0_8 
coef3_1 = 5.0_8/36.0_8 - cc/3.0_8; coef3_2 = 2.0_8/9.0_8 - 2.0_8*cc/3.0_8; coef3_3 = 5.0_8/36.0_8

 err = 0.0_8

 do j = p%js, p%je
 do i = p%is, p%ie
 
	p%ss1(i,j) = p%s1(i,j)
	p%ss2(i,j) = p%s2(i,j)
	p%ss3(i,j) = p%s3(i,j)
	
	p%s1(i,j) = p%target(i,j) + p%dt*( coef1_1*p%l1(i,j) + coef1_2*p%l2(i,j) + coef1_3*p%l3(i,j) )
	p%s2(i,j) = p%target(i,j) + p%dt*( coef2_1*p%l1(i,j) + coef2_2*p%l2(i,j) + coef2_3*p%l3(i,j) )
	p%s3(i,j) = p%target(i,j) + p%dt*( coef3_1*p%l1(i,j) + coef3_2*p%l2(i,j) + coef3_3*p%l3(i,j) )
	
	p%s1(i,j) = p%w * p%s1(i,j) + (1.0_8-p%w) * p%ss1(i,j) 
	p%s2(i,j) = p%w * p%s2(i,j) + (1.0_8-p%w) * p%ss2(i,j) 
	p%s3(i,j) = p%w * p%s3(i,j) + (1.0_8-p%w) * p%ss3(i,j) 
	
 enddo
 enddo

 do j = p%js+p%ghc, p%je-p%ghc
 do i = p%is+p%ghc, p%ie-p%ghc
	
	err = max( err, abs( p%s1(i,j)-p%ss1(i,j) ), abs( p%s2(i,j)-p%ss2(i,j) ), abs( p%s3(i,j)-p%ss3(i,j) ) )
	
 enddo
 enddo
 
 
end subroutine 

subroutine srk6_data_solve(p,err)
implicit none
class(srk6_data) :: p
real(8), intent(out) :: err
real(8) :: err1, err2

 call p%x%solve(err1)

 if( p%is_vector_solver )then
	call p%y%solve(err2)
	err = max( err1, err2 )
 else
	err = err1
 endif

end subroutine 

subroutine srk6_roots_final(p)
implicit none
class(srk6_roots) :: p
integer :: i,j

 do j = p%js, p%je
 do i = p%is, p%ie
	 p%target(i,j) = p%target(i,j) + p%dt/18.0_8 * ( 5.0_8*p%l1(i,j) + 8.0_8*p%l2(i,j) + 5.0_8*p%l3(i,j) )
 enddo
 enddo

end subroutine

subroutine srk6_data_final(p)
implicit none
class(srk6_data) :: p
integer :: i,j

 call p%x%final

 if( p%is_vector_solver )then
	 call p%y%final
 endif
 
 p%is_vector_solver = .false.
 
 end subroutine

end module srk6_solver