include 'sorts.f90'

program main
use omp_lib
implicit none
integer :: num_of_test, data_size
real(8) :: data_range

data_size = 12000
num_of_test = 10000

call omp_set_dynamic(.false.)
call omp_set_num_threads(omp_get_max_threads())

call dec_stump(10000,12000,2.0d0,0.0d0)
call dec_stump(100000,12000,2.0d0,0.0d0)

contains

function sign(x)
implicit none
integer :: sign
real(8) :: x

    if(x>0)then
        sign = 1
    else
        sign = -1
    endif

end function

subroutine dec_stump(num_of_test,data_size,data_range,tau)
use omp_lib
use sorts
implicit none
integer :: num_of_test,data_size
integer :: tid, i, j, id, s
real(8) :: data_range, tau, sum, randnum
real(8) :: E, E_total, Emin, mtheta, ms
real(8),allocatable :: x(:,:),theta(:,:)

allocate(x(omp_get_max_threads(),data_size),&
         theta(omp_get_max_threads(),data_size))

call random_seed()

sum=0.0d0
!$omp parallel do reduction(+:sum,E_total), private(id,i,j,randnum,s,E,Emin)
do tid = 1, num_of_test

    id = omp_get_thread_num()
    do i = 1, data_size
        call random_number(randnum)
        x(id,i) = data_range*randnum-data_range*0.5d0
    enddo

    call heap_sort(x(id,:))

    theta(tid,1) = -1.0d0
    do i = 2, data_size
        theta(id,i) = 0.5d0*( x(id,i)+x(id,i-1) )
    enddo

    Emin=data_size
    do s = 1, -1, 2
    do j = 1, data_size
        E = 0.0d0
        do i = 1, data_size
            if( sign(x(id,i)) .ne. s*sign(x(id,i)-theta(id,j)) )then
                E = E + 1.0d0
            endif
        enddo
        if( E<Emin )then
            Emin=E
            ms=s
            mtheta = theta(id,j)
        endif
    enddo
    enddo
    E_total = E_total + E / real(data_size,8)

enddo
!$omp end parallel do

write(*,'("Error:",ES11.4,I2,ES11.4)')E_total / real(num_of_test,8), ms, mtheta


end subroutine

end program