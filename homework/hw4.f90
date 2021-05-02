include 'data.f90'

program main
use data
implicit none
type(task) :: pr
type(task),allocatable :: p(:)
integer :: id, num
character(1) :: name

num=5

call pr%run(0.001d0,"ref")

allocate(p(num))

do id = 1, num
    write(name,'(i1.1)')id
    call p(id)%run(0.06d0/2.0d0**(id-1),name)
    write(*,'(2ES15.4)')p(id)%dx,l2_error(pr,p(id)) 
enddo

contains

function l2_error(pr,p)  result(l2)
use data
implicit none
type(task) :: pr, p
integer :: i,j
real(8) :: int, l2

!!$omp parallel do private(j,int), reduction(+:l2)
do i = 2, p%n-1
    do j = 1, pr%n
        if( pr%x(j)>p%x(i) .and. pr%x(j-1)<p%x(i) )then
            int = pr%eta%now(j-1) + (pr%eta%now(j)-pr%eta%now(j-1))/pr%dx*(p%x(i)-pr%x(j-1))
            l2 = l2 + (int-p%eta%now(i))**2
            exit
        endif
    enddo
enddo
!!$omp end parallel do

l2 = dsqrt(l2/p%n)


end function

end program
