module time_solver
implicit none

type tsolver_roots
integer :: is, ie, js, je, ghc
real(8) :: dt, w
real(8),dimension(:,:),allocatable :: target
real(8),dimension(:,:),allocatable :: s1,s2,s3,ss1,ss2,ss3,l1,l2,l3
contains
procedure alloc => tsolver_roots_alloc
procedure init => tsolver_roots_init
procedure solve_srk6 => tsolver_roots_solve_srk6
procedure final_srk6 => tsolver_roots_final_srk6
procedure solve_srk4 => tsolver_roots_solve_srk4
procedure final_srk4 => tsolver_roots_final_srk4
end type tsolver_roots

type tsolver_data
integer :: is, ie, js, je
logical :: is_vector_solver
type(tsolver_roots) :: x, y
contains
procedure alloc => tsolver_data_alloc
procedure init => tsolver_data_init
procedure solve_srk6 => tsolver_data_solve_srk6
procedure final_srk6 => tsolver_data_final_srk6
procedure solve_srk4 => tsolver_data_solve_srk4
procedure final_srk4 => tsolver_data_final_srk4
end type tsolver_data

contains

subroutine tsolver_roots_alloc(p,is,ie,js,je,dt,w,ghc)
implicit none
class(tsolver_roots) :: p
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

subroutine tsolver_data_alloc(p,is,ie,js,je,dt,w,ghc)
implicit none
class(tsolver_data) :: p
integer, intent(in) :: is, ie, js, je, ghc
real(8), intent(in) :: dt, w

p%is = is; p%ie = ie
p%js = js; p%je = je

call p%x%alloc(is,ie,js,je,dt,w,ghc)
call p%y%alloc(is,ie,js,je,dt,w,ghc)

end subroutine

subroutine tsolver_roots_init(p,u)
implicit none
class(tsolver_roots) :: p
real(8),dimension(p%is:p%ie,p%js:p%je) :: u
integer :: i,j

 do j = p%js, p%je
 do i = p%is, p%ie
    p%target(i,j) = u(i,j)
    p%s1(i,j) = u(i,j)
    p%s2(i,j) = u(i,j)
    p%s3(i,j) = u(i,j)
 enddo
 end do 

end subroutine

subroutine tsolver_data_init(p,btn,x,y)
implicit none
class(tsolver_data) :: p
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

include 'srk6.f90'
include 'srk4.f90'

end module time_solver
