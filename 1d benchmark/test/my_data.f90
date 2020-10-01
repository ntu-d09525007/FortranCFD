include 'my_precision.f90'
include 'matrix_solver.f90'
include 'ccd_solvers.f90'

module mydata
use CCD_SOLVERS
implicit none
save
real(kind=ap) :: tol, w
logical :: IsPloting, IsShowing

type my_solver
integer :: id, is, ie
real(kind=ap) :: dx, dt
type(ccd_root) :: ccd
contains
procedure init => solver_init
procedure solve => solver_solve
end type my_solver

type pure_adv
integer :: nx, turn
integer :: id, pid, sid
real(kind=ap) :: xstart, xend, xl
real(kind=ap) :: c, dx, dt, cfl, norm(3), cpu_time
real(kind=ap) :: t, t2s, t2p, t2show, t_shock
real(kind=ap),allocatable,dimension(:) :: x, cc, u, u_exact
real(kind=ap),allocatable,dimension(:) :: u1, u2, u3
real(kind=ap),allocatable,dimension(:) :: s1, s2, s3
real(kind=ap),allocatable,dimension(:) :: uu1, uu2, uu3
type(my_solver) :: solver
contains 
procedure init => initialize
procedure set => set_value
procedure alloc => alloc_data
procedure bc => boundary_condition
procedure find_error => calculate_error
procedure plot => plot_data
procedure srk4 => solve_srk4
procedure srk6 => solve_srk6
procedure find_dt => get_time_step
end type pure_adv

contains

subroutine solver_init(this,is,ie,dx,dt,id)
implicit none
class(my_solver) :: this
integer :: is, ie, id
real(kind=ap) :: dx,dt
 
 this%id=id; this%is=is; this%ie=ie
 this%dx=dx; this%dt=dt
 
 this%ccd%is = is-3; this%ccd%ie = ie+3
 this%ccd%dx = dx; this%ccd%dt = dt;

end subroutine

subroutine solver_solve(this,u,f,s)
implicit none
class(my_solver) :: this
real(kind=ap), dimension(this%is:this%ie) :: u, f, s, tmp
integer :: i

call this%ccd%solve("srkccd",f,s,tmp,u)
 
 !$omp parallel do
 do i = this%is, this%ie
    s(i) = s(i)*u(i)
 enddo
 !$omp end parallel do
 
end subroutine

subroutine initialize(this,id,num)
implicit none
class(pure_adv) :: this
integer :: id, i, num
real(kind=ap) :: c3

this%id = id
this%pid = -1;if(IsPloting)this%pid=0;
this%sid = -1;if(IsShowing)this%sid=0;

this%xstart = -1.0_ap
this%xend   = 1.0_ap
this%xl = this%xend - this%xstart

this%nx = 16*2**num
this%dx = this%xl / real(this%nx-1,kind=ap)
this%cfl = 0.1_ap !+ 0.2_ap*(id-1)
this%c = 1.0_ap

this%turn = 1
this%t2s = this%xl / this%c
this%t2p = 0.2_ap*this%t2s
this%t = 0.0_ap
this%t_shock = 2.0_ap

call this%alloc

!$omp parallel do
do i = -2, this%nx+3
  this%x(i) = this%xstart + real(i-1,kind=ap)*this%dx
  this%cc(i) = this%c
enddo
!$omp end parallel do

call this%set(this%u,this%c,0.0_ap)

call this%solver%init(1, this%nx, this%dx, this%dt, this%id )

call this%plot

end subroutine

subroutine alloc_data(this)
implicit none
class(pure_adv) :: this

 allocate(this%x(-2:this%nx+3), this%cc(-2:this%nx+3), this%u(-2:this%nx+3), this%u_exact(-2:this%nx+3) )
 allocate(this%u1(-2:this%nx+3), this%u2(-2:this%nx+3), this%u3(-2:this%nx+3))
 allocate(this%s1(-2:this%nx+3), this%s2(-2:this%nx+3), this%s3(-2:this%nx+3))
 allocate(this%uu1(-2:this%nx+3), this%uu2(-2:this%nx+3), this%uu3(-2:this%nx+3))
 
end subroutine

subroutine boundary_condition(this,f)
implicit none
class(pure_adv) :: this
integer :: i
real(kind=ap),dimension(-2:this%nx+3) :: f

 f(this%nx+1) = f(2)
 f(this%nx+2) = f(3)
 f(this%nx+3) = f(4)
 
 f(0) = f(this%nx-1)
 f(-1) = f(this%nx-2)
 f(-2) = f(this%nx-3)
 
 if(this%c>=0.0_ap)then
   f(1) = f(this%nx)
 else
   f(this%nx) = f(1)
 endif
 
end subroutine

subroutine set_value(this,u,c,t)
use my_precision
implicit none
class(pure_adv) :: this
real(kind=ap),dimension(-2:this%nx+3) :: u
integer :: i
real(kind=ap) :: x,t,c

 !$omp parallel do private(x)
 do i = 1, this%nx
    x = this%x(i)-c*t
    u(i) =  sin(acos(-1.0_ap)*x)
 enddo
 !$omp end parallel do

 
end subroutine

subroutine calculate_error(this,a,b,tmp1,tmp2,tmp3)
use my_precision
implicit none
class(pure_adv) :: this
integer :: i
real(kind=ap), dimension(-2:this%nx+3) :: a,b
real(kind=ap) :: tmp1, tmp2, tmp3, err

  tmp1 = 0.0_ap
  tmp2 = 0.0_ap
  tmp3 = 0.0_ap

  !$omp parallel do reduction(+:tmp1,tmp2), reduction(max:tmp3), private(err)
  do i = 2, this%nx-1
      err = a(i) - b(i)
      tmp1 = tmp1 + abs(err)
      tmp2 = tmp2 + err**2.0_ap
      tmp3 = max(tmp3, abs(err))
  enddo
  !$omp end parallel do
  
 tmp1 = tmp1 / real(this%nx-2,kind=ap)
 tmp2 = sqrt( tmp2 / real(this%nx-2,kind=ap) )
 tmp3 = tmp3
 
end subroutine

subroutine plot_data(this)
implicit none
class(pure_adv) :: this
integer :: i, uid
character(10) :: str

 if(abs(this%t-this%pid*this%t2p)>this%dt)return

 write(str,'(I2,A,I2,A)')this%id,'_',this%pid,'.dat'
 
 uid = this%id + 11
 
 open(unit=uid,file=trim(str))
 
 do i = 1, this%nx
   write(uid,'(2es15.4)')this%x(I),this%u(i)
 enddo
 
 close(unit=uid)
 
 this%pid = this%pid + 1

end subroutine

subroutine solve_srk6(this)
use my_precision
implicit none
class(pure_adv) :: this
integer :: i, iter
real(kind=ap)  :: CC , COE1_1 , COE1_2 , COE1_3  , COE2_1 , COE2_2 , COE2_3 , COE3_1 , COE3_2 , COE3_3
real(kind=ap), dimension(3) :: err1, err2, err3

  CC = 0.5_ap * SQRT(3.0_ap / 5.0_ap)
     
  Coe1_1 = 5.0_ap / 36.0_ap
  Coe1_2 = 2.0_ap /  9.0_ap + 2.0_ap * CC / 3.0_ap
  Coe1_3 = 5.0_ap / 36.0_ap +          CC / 3.0_ap
    
  Coe2_1 = 5.0_ap / 36.0_ap - 5.0_ap * CC / 12.0_ap
  Coe2_2 = 2.0_ap /  9.0_ap
  Coe2_3 = 5.0_ap / 36.0_ap + 5.0_ap * CC / 12.0_ap
    
  Coe3_1 = 5.0_ap / 36.0_ap - CC / 3.0_ap          
  Coe3_2 = 2.0_ap /  9.0_ap - 2.0_ap * CC / 3.0_ap
  Coe3_3 = 5.0_ap / 36.0_ap   

 !$omp parallel do
 do i = 1, this%nx
    this%u1(I) = this%u(i)
    this%u2(i) = this%u(i)
    this%u3(i) = this%u(i)
 enddo
 !$omp end parallel do
 
 call this%bc(this%u1)
 call this%bc(this%u2)
 call this%bc(this%u3)
 
 iter = 0
 
 do
     iter=iter+1
     
    !$omp parallel do
    do i = 1, this%nx
        this%uu1(I) = this%u1(i)
        this%uu2(i) = this%u2(i)
        this%uu3(i) = this%u3(i)
    enddo
    !$omp end parallel do
    
    call this%solver%solve(this%cc,this%u1,this%s1)
    call this%solver%solve(this%cc,this%u2,this%s2)
    call this%solver%solve(this%cc,this%u3,this%s3)
    
    !$omp parallel do 
    do i = 1, this%nx
        
        this%u1(i) = this%u(i) - this%dt*(coe1_1*this%s1(i)+&
                                          coe1_2*this%s2(i)+&
                                          coe1_3*this%s3(i))
        
        this%u2(i) = this%u(i) - this%dt*(coe2_1*this%s1(i)+&
                                          coe2_2*this%s2(i)+&
                                          coe2_3*this%s3(i))
        
        this%u3(i) = this%u(i) - this%dt*(coe3_1*this%s1(i)+&
                                          coe3_2*this%s2(i)+&
                                          coe3_3*this%s3(i))
        
        this%u1(i) = w*this%u1(i) + (1.0_ap-w)*this%uu1(i)
        this%u2(i) = w*this%u2(i) + (1.0_ap-w)*this%uu2(i)
        this%u3(i) = w*this%u3(i) + (1.0_ap-w)*this%uu3(i)
    enddo
    !$omp end parallel do
    
    call this%bc(this%u1)
    call this%bc(this%u2)
    call this%bc(this%u3)
     
    call this%find_error(this%u1,this%uu1,err1(1),err1(2),err1(3))
    call this%find_error(this%u2,this%uu2,err2(1),err2(2),err2(3))
    call this%find_error(this%u3,this%uu3,err3(1),err3(2),err3(3))
    
    if(max(err1(3),err2(3),err3(3))<tol)exit
    
    !if(mod(iter,100).eq.0)write(*,*)iter,max(err1(3),err2(3),err3(3))
     
 enddo
 
 call this%solver%solve(this%cc,this%u1,this%s1)
 call this%solver%solve(this%cc,this%u2,this%s2)
 call this%solver%solve(this%cc,this%u3,this%s3)
 
 !$omp parallel do 
 do i = 1, this%nx
     this%u(i) = this%u(i) - this%dt/18.0_ap*(5.0_ap*this%s1(i)+&
                                              8.0_ap*this%s2(i)+&
                                              5.0_ap*this%s3(i))
 enddo
 !$omp end parallel do
 
 call this%bc(this%u)
 
end subroutine

subroutine solve_srk4(this)
use my_precision
implicit none
class(pure_adv) :: this
integer :: i, iter
real(kind=ap)  :: CC , COE1_1 , COE1_2 ,  COE2_1 , COE2_2
real(kind=ap), dimension(3) :: err1, err2, err3

  CC = 0.5_ap / SQRT(3.0_AP)
     
  Coe1_1 = 0.25_ap
  Coe1_2 = 0.25_ap - cc
 
  Coe2_1 = 0.25_ap + cc
  Coe2_2 = 0.25_ap

 !$omp parallel do
 do i = 1, this%nx
    this%u1(I) = this%u(i)
    this%u2(i) = this%u(i)
 enddo
 !$omp end parallel do
 
 call this%bc(this%u1)
 call this%bc(this%u2)
 
 iter = 0
 
 do
     iter=iter+1
     
    !$omp parallel do
    do i = 1, this%nx
        this%uu1(I) = this%u1(i)
        this%uu2(i) = this%u2(i)
    enddo
    !$omp end parallel do
    
    call this%solver%solve(this%cc,this%u1,this%s1)
    call this%solver%solve(this%cc,this%u2,this%s2)
        
    !$omp parallel do 
    do i = 1, this%nx
        
        this%u1(i) = this%u(i) - this%dt*(coe1_1*this%s1(i)+&
                                          coe1_2*this%s2(i))
        
        this%u2(i) = this%u(i) - this%dt*(coe2_1*this%s1(i)+&
                                          coe2_2*this%s2(i))
        
        this%u1(i) = w*this%u1(i) + (1.0_ap-w)*this%uu1(i)
        this%u2(i) = w*this%u2(i) + (1.0_ap-w)*this%uu2(i)
    enddo
    !$omp end parallel do
    
    call this%bc(this%u1)
    call this%bc(this%u2)
     
    call this%find_error(this%u1,this%uu1,err1(1),err1(2),err1(3))
    call this%find_error(this%u2,this%uu2,err2(1),err2(2),err2(3))
    
    if(max(err1(3),err2(3))<tol)exit
    
    !if(mod(iter,100).eq.0)write(*,*)iter,max(err1(3),err2(3))
     
 enddo
 
 call this%solver%solve(this%cc,this%u1,this%s1)
 call this%solver%solve(this%cc,this%u2,this%s2)
 
 !$omp parallel do 
 do i = 1, this%nx
     this%u(i) = this%u(i) - this%dt/2.0_ap*(this%s1(i)+this%s2(i))
 enddo
 !$omp end parallel do
 
 call this%bc(this%u)
 
end subroutine

subroutine get_time_step(this) 
use my_precision
implicit none
class(pure_adv) :: this
integer :: i
real(kind=ap) :: maxc
 
 this%dt = this%cfl * this%dx / abs(this%c)

end subroutine


end module mydata




