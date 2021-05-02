module data
implicit none
save

type :: work_data
! ==============================
! Define structure for physical properties
! ==============================
    integer :: n
    real(8),dimension(:),allocatable :: now, old, x, xx
    contains
    procedure init => data_init
    procedure switch => data_switch
    procedure get_fx  => find_first_dervaitve
    procedure get_fxx => find_second_dervaitve
    procedure bc => boundary_condition
end type work_data

type :: task
! ==============================
! Define all data needed to accomplish a simulation
! ==============================
    integer :: n, pltid
    real(8),dimension(:),allocatable :: x,A,B,C,S,diff
    type(work_data) :: h,u,phi,eta,psi,hu
    real(8) :: g,alpha
    real(8) :: xstart, xend, dx
    real(8) :: dt,t,t2s,t2p
    real(8) :: error
    character(20) :: name
    logical :: sponge
    contains
    procedure run => task_run
    procedure solve => ssp_rk3_solve
    procedure find_u => task_find_u
    procedure find_src => task_find_src
    procedure plot => task_plot
end type task

contains

subroutine data_init(p,n)
! ==============================
! Allocate the data
! ==============================
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
! ==============================
! Get first derivative with fourth-order central difference scheme
! ==============================
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
! ==============================
! Get first derivative with fourth-order central difference scheme
! ==============================
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
! ==============================
! Open boundary condition, ( d{\phi}/dx = 0 at boundary )
! ==============================
implicit none
class(work_data) :: p

p%now(0)  = p%now(1)
p%now(-1) = p%now(1)
p%now(p%n+1) = p%now(p%n)
p%now(p%n+2) = p%now(p%n)

end subroutine

subroutine task_run(p,dx,name,sponge)
implicit none
class(task) :: p
real(8) :: dx
character(*) :: name
integer :: i
logical,optional :: sponge


p%sponge = .false.
if(present(sponge)) p%sponge = sponge

p%xstart=-6
p%xend=6

! ==============================
! Create virutal doamin for sponge layer
! ==============================
if(p%sponge)then
    p%xstart = -8
    p%xend = 8
endif

p%dx=dx
p%g=9.81
p%alpha=0.531
p%n=(p%xend-p%xstart)/p%dx + 1
p%name = trim(name)


allocate(p%A(p%n),p%B(p%n),p%C(p%n),p%S(p%n),p%x(p%n),p%diff(p%n))

call p%h%init(p%n)
call p%u%init(p%n)
call p%phi%init(p%n)
call p%eta%init(p%n)
call p%psi%init(p%n)
call p%hu%init(p%n)

!$omp parallel do 
do i =  1, p%n
    p%x(i) = p%xstart + (i-1)*dx
    p%h%now(i) = 0.3
    p%eta%now(i) = 0.1*dexp(-18.0*p%x(i)**2)
    p%u%now(i)=0.0
    p%phi%now(i) = 0.0
enddo
!$omp end parallel do 
call p%h%bc
call p%eta%bc
call p%phi%bc
call p%u%bc

! ==============================
! Define dacay coefficient for sponge layer
! ==============================
if(.not. p%sponge)then
    !$omp parallel do
    do i = 1, p%n
        p%diff(i) = 1.0d0
    enddo
    !$omp end parallel do
else
    !$omp parallel do
    do i = 1, p%n
        if( p%x(i) > 6.0 )then
            p%diff(i) = 0.5 * (1.0 - dcos(dacos(-1.0d0)*(p%x(i)-6.0)/(p%xend-6.0)) )
        else if ( p%x(i) < -6.0 )then
            p%diff(i) = 0.5 * (1.0 - dcos(dacos(-1.0d0)*(p%x(i)+6.0)/(p%xstart+6.0)) )
        else
            p%diff(i) = 1.0
        endif
    enddo
    !$omp end parallel do    
endif

p%pltid=0
p%t=0.0
if(present(sponge))then
    ! ==============================
    ! Make the solutions to reach boundary to see
    ! how sponge layer would affect the solutions.
    ! ==============================
    p%dt=0.9*p%dx/dsqrt(p%g*0.4)
    p%t2s=4.0
    p%t2p=1.0
else
    p%dt=0.9*0.001/dsqrt(p%g*0.4)
    p%t2s=2.0
    p%t2p=0.5
endif

call p%plot

do 
    p%t=p%t+p%dt
    call p%solve
    call p%plot
    if(p%t>p%t2s)exit

enddo

end subroutine

subroutine ssp_rk3_solve(p)
! ==============================
! Solve the equations with SSP-RK3
! ==============================
implicit none
class(task) :: p
integer :: i

call p%eta%switch
call p%phi%switch

call p%find_src()
!$omp parallel do 
do i = 1, p%n
    p%eta%now(i) = p%eta%old(i) - p%dt * p%psi%x(i)
    p%phi%now(i) = p%phi%old(i) - p%dt * p%g * p%eta%x(i)
    ! ==============================
    !  Sponge layer
    ! ==============================
    p%eta%now(i) = p%diff(i) * p%eta%now(i)
enddo
!$omp end parallel do

call p%find_src()
!$omp parallel do 
do i = 1, p%n
    p%eta%now(i) = ( ( p%eta%now(i) - p%dt * p%psi%x(i) ) + 3.0 * p%eta%old(i) )/4.0
    p%phi%now(i) = ( ( p%phi%now(i) - p%dt * p%g * p%eta%x(i) ) + 3.0 * p%phi%old(i) )/4.0
    ! ==============================
    !  Sponge layer
    ! ==============================
    p%eta%now(i) = p%diff(i) * p%eta%now(i)
enddo
!$omp end parallel do

call p%find_src()
!$omp parallel do 
do i = 1, p%n
    p%eta%now(i) = ( ( p%eta%now(i) - p%dt * p%psi%x(i) ) * 2.0 + p%eta%old(i) )/3.0
    p%phi%now(i) = ( ( p%phi%now(i) - p%dt * p%g * p%eta%x(i) ) * 2.0 + p%phi%old(i) )/3.0
    ! ==============================
    !  Sponge layer
    ! ==============================
    p%eta%now(i) = p%diff(i) * p%eta%now(i)
enddo
!$omp end parallel do

end subroutine

subroutine task_find_u(p)
! ==============================
! Get the horizontal velocit by second-order difference
! scheme and solve the matrix with Thomas tridiagonal algorithm.
! ==============================
implicit none
class(task) :: p
integer :: i,j
real(8) :: a1,a2,a3,z

!$omp parallel do private(a1,a2,a3,z)
do i = 1, p%n

    z = - p%alpha * p%h%now(i)

    a1 = (z*p%h%now(i-1)+0.5*z**2) / p%dx**2
    a2 = 1.0 - (2.0*z*p%h%now(i)+z**2) / p%dx**2
    a3 = (z*p%h%now(i+1)+0.5*z**2) / p%dx**2

    p%B(i) = a2
    p%A(i) = a1
    p%C(i) = a3

    if (i==1)then
        p%B(i) = a2+a1
        p%C(i) = a3
    else if (i==p%N)then
        p%B(i) = a2+a3
        p%A(i) = a1
    endif

    p%S(i) = p%phi%now(i)

enddo
!$omp end parallel do

call solve_tridiagonal(p%A,p%B,p%C,p%S,p%u%now(1:p%n),1,p%n)

call p%u%bc()

end subroutine

subroutine solve_tridiagonal(A,B,C,S,X,M,N)
!cccccccccccccccccccccccccccccccccc
!
!  Solve of Thomas tridiagonal algorithm
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
real(kind=8),dimension(m:n) :: A,B,C,S,X
 
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

subroutine task_find_src(p)
! ==============================
! Find the source to march [eta, (u+u_2)]
! ==============================
implicit none
class(task) :: p
integer :: i
real(8) :: z, u1

call p%find_u()

!$omp parallel do
do i = 1, p%n
    p%hu%now(i) = p%u%now(i) * p%h%now(i)
enddo
!$omp end parallel do

call p%hu%get_fxx(p%dx)
call p%u%get_fxx(p%dx)

!$omp parallel do private(z,u1)
do i = 1, p%n
    z =  - p%alpha * p%h%now(i)
    u1 = (z+0.5*p%h%now(i))*p%hu%xx(i)+(z**2/2.0-p%h%now(I)**2/6.0)*p%u%xx(i)
    p%psi%now(i) = p%h%now(i)*(p%u%now(i)+u1)
enddo
!$omp end parallel do

call p%psi%get_fx(p%dx)
call p%eta%get_fx(p%dx)

end subroutine

subroutine task_plot(p)
! ==============================
! Plot the data with user-defined period
! ==============================
implicit none
class(task) :: p
integer :: i
character(2) :: name

if(abs(p%t-p%pltid*p%t2p)>p%dt)return

write(name,'(i2.2)')p%pltid
open(unit=66,file='hw4_'//trim(p%name)//'_'//name//'.plt')

do i = 1, p%n
    write(66,*)p%x(i),p%eta%now(i)
enddo
close(unit=66)
p%pltid=p%pltid+1

end subroutine

end module data