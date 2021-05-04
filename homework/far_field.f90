subroutine task_run2(p,dx,name)
implicit none
class(task) :: p
real(8) :: dx
character(*) :: name
integer :: i

p%xstart=-20.0
p%xend=20.0

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
    p%diff(i) = 1.0d0
enddo
!$omp end parallel do 
call p%h%bc
call p%eta%bc
call p%phi%bc
call p%u%bc

p%pltid=0
p%dt=0.9*p%dx/dsqrt(p%g*0.4)
p%t2s=10.0
p%t2p=2.5d0

call p%plot

do 
    p%t=p%t+p%dt
    call p%solve
    call p%plot
    if(p%t>p%t2s)exit

enddo

end subroutine

subroutine task_plot2(p)
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