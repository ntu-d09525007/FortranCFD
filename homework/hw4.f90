module data
implicit none
save
type :: work_data
    integer :: n
    real(8),dimension(:),allocatable :: now, old, x, xx
    contains
    procedure init => data_init
    procedure switch => data_switch
    procedure get_fx  => find_first_dervaitve
    procedure get_fxx => find_second_dervaitve
    procedure bc => boundary_condition
end type work_data
integer :: n, pltid
real(8),dimension(:,:),allocatable :: A
real(8),dimension(:),allocatable :: x
type(work_data) :: h,u,phi,eta,psi,hu
real(8) :: g,alpha
real(8) :: xstart, xend, dx
real(8) :: dt,t,t2s,t2p

contains

subroutine data_init(p,n)
implicit none
class(work_data) :: p
integer :: n

p%n = n

allocate(p%now(-1:n+2),p%old(-1:n+2),p%x(-1:n+2),p%xx(-1:n+2))

end subroutine

subroutine data_switch(p)
implicit none
class(work_data) :: p
integer :: i

!$omp parallel do 
do i = -1, p%n+2
    p%old(i) = p%now(i)
enddo
!$omp end parallel do 

end subroutine

subroutine find_first_dervaitve(p,dx)
implicit none
class(work_data) :: p
real(8) :: dx
integer :: i

call p%bc()

!$omp parallel do
do i = 1, p%n
    p%x(i) = (-p%now(i-2)+8.0*p%now(i-1)-8.0*p%now(i+1)+p%now(i+2))/(12.0*dx)
enddo
!$omp end parallel do

end subroutine

subroutine find_second_dervaitve(p,dx)
implicit none
class(work_data) :: p
real(8) :: dx
integer :: i

call p%bc()

!$omp parallel do
do i = 1, p%n
    p%xx(i) = (-p%now(i-2)+16.0*p%now(i-1)-30.0*p%now(i)+16.0*p%now(i+1)-p%now(i+2))/(12.0*dx**2)
enddo
!$omp end parallel do

end subroutine

subroutine boundary_condition(p)
implicit none
class(work_data) :: p

p%now(0)  = p%now(2)
p%now(-1) = p%now(3)
p%now(p%n+1) = p%now(p%n-1)
p%now(p%n+2) = p%now(p%n-2)

end subroutine

end module data

program main
use data
implicit none
integer :: i

xstart=-6
xend=6
dx=0.06
g=9.81
alpha=-0.531
n=(xend-xstart)/dx + 1

allocate(A(n,n),x(n))
call h%init(n)
call u%init(n)
call phi%init(n)
call eta%init(n)
call psi%init(n)
call hu%init(n)

!$omp parallel do 
do i =  1, n
    x(i) = xstart + (i-1)*dx
    h%now(i) = 0.3
    eta%now(i) = 0.1*dexp(-18.0*x(i)**2)
    phi%now(i) = 0.0
enddo
!$omp end parallel do 
call h%bc
call eta%bc
call phi%bc

pltid=0
dt=0.9*dx/dsqrt(g*0.4)
t=0.0
t2s=2.0
t2p=0.5

call plot

do 
    t=t+dt
    call ssp_rk3_solve
    call plot
    if(t>t2s)exit

enddo

contains

subroutine find_u()
use data, only : h, alpha, dx, n, A, phi, u
implicit none
integer :: i,j,info
real(8) :: a1,a2,a3,z
real(8),dimension(n) :: ipiv, work

!$omp parallel do
do i = 1, n
do j = 1, n
    A(i,j) = 0.0
enddo
enddo
!$omp end parallel do

!$omp parallel do private(a1,a2,a3,z)
do i = 1, n

    z = -alpha*h%now(i)

    a1 = (z*h(i-1)+0.5*z**2) / dx**2
    a2 = 1.0 - (2.0*z*h%now(i)+z**2) / dx**2
    a3 = (z*h(i+1)+0.5*z**2) / dx**2

    A(i,i) = a2

    if (i>1 .and. i<N) then
        A(i,i-1) = a1
        A(i,i+1) = a3
    else if (i==1)then
        A(i,i) = a2+a3
        A(i,i+1) = a3
    else 
        A(i,i) = a1+a2
        A(i,i-1) = a1
    endif

enddo
!$omp end parallel do

call SGETRF(n,n,A,n,ipiv,info)
call SGETRI(n,A,n,ipiv,work,n,info)

u%now(1:n) = matmul(A,phi%now(1:n))

call u%bc()

end subroutine

subroutine find_src()
use data, only : eta,phi,psi,u,h,hu,n,dx,alpha
implicit none
integer :: i
real(8) :: z, u1

call find_u()

!$omp parallel do
do i = 1, n
    hu%now(i) = u%now(i) * h%now(i)
enddo
!$omp end parallel do

call hu%get_fxx(dx)
call u%get_fxx(dx)

!$omp parallel do private(z,u1)
do i = 1, n
    z =  - alpha * h%now(i)
    u1 = (z+0.5)*hu%xx(i)+(z**2/2.0-h%now(I)**2/6.0)*u%xx(i)
    psi%now(i) = h%now(i)*(u%now(i)+u1)
enddo
!$omp end parallel do

call psi%get_fx(dx)
call eta%get_fx(dx)

end subroutine

subroutine ssp_rk3_solve()
use data
implicit none
integer :: i

call eta%switch
call phi%switch

call find_src()
!$omp parallel do 
do i = 1, n
    eta%now(i) = eta%old(i) - dt * psi%x(i)
    phi%now(i) = phi%old(i) - dt * g * eta%x(i)
enddo
!$omp end parallel do

call find_src()
!$omp parallel do 
do i = 1, n
    eta%now(i) = ( ( eta%now(i) - dt * psi%x(i) ) + 3.0 * eta%old(i) )/4.0
    phi%now(i) = ( ( phi%now(i) - dt * g * eta%x(i) ) + 3.0 * phi%old(i) )/4.0
enddo
!$omp end parallel do

call find_src()
!$omp parallel do 
do i = 1, n
    eta%now(i) = ( ( eta%now(i) - dt * psi%x(i) ) * 2.0 + eta%old(i) )/4.0
    phi%now(i) = ( ( phi%now(i) - dt * g * eta%x(i) ) * 2.0 + phi%old(i) )/4.0
enddo
!$omp end parallel do
call eta%bc; call phi%bc

end subroutine

subroutine plot()
use data
implicit none
integer :: i
character(2) :: name

if(abs(t-pltid*t2p)>dt)return

write(name,'(i2.2)')pltid
open(unit=66,file='hw4_'//name//'.dat')

do i = 1, n
    write(66,*)x(i),eta%now(i)
enddo
close(unit=66)
pltid=pltid+1

end subroutine

end program