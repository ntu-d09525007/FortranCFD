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

call dec_stump(num_of_test,data_size,2.0d0,0.0d0)
call dec_stump(num_of_test,data_size,20.0d0,0.0d0)

call dec_stump(num_of_test,data_size,20.0d0,0.1d0)
call dec_stump(num_of_test,data_size,200.0d0,0.1d0)

contains

function sign2(x,tau)
implicit none
integer :: sign2
real(8) :: x,tau,rand

call random_number(rand)
    
if( rand > tau )then

    if(x>0)then
        sign2 = -1
    else
        sign2 = 1
    endif

else

    if(x>0)then
        sign2 = 1
    else
        sign2 = -1
    endif

endif

end function

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
integer :: num_of_test,data_size,out_size
integer :: tid, i, j, id, s
real(8) :: data_range, tau, sum, randnum
real(8) :: E, E_total, Emin, mtheta, ms
real(8),allocatable :: x(:,:),theta(:,:),Ein(:),xx(:,:)

out_size = 10

allocate(x(omp_get_max_threads(),data_size),xx(out_size,data_size),&
         theta(omp_get_max_threads(),data_size),Ein(num_of_test))

call random_seed()

!$omp parallel do private(randnum)
do j = 1, out_size
do i = 1, data_size
    call random_number(randnum)
    xx(j,i) = data_range*randnum-data_range*0.5d0
enddo
enddo
!$omp end parallel do

!$omp parallel do private(id,i,j,randnum,s,E,Emin)
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

    ! Training

    Emin=data_size
    do s = 1, -1, 2
    do j = 1, data_size
        E = 0.0d0
        do i = 1, data_size
            if( sign2(x(id,i),tau) .ne. s*sign(x(id,i)-theta(id,j)) )then
                E = E + 1.0d0
            endif
        enddo
        E = E / real( data_size, 8)
        if( E<Emin )then
            Emin=E
            ms=s
            mtheta = theta(id,j)
        endif
    enddo
    enddo

    ! Calculate Eout

    E=0.0d0
    do i = 1, out_size
    do j = 1, data_size
        if( sign2(xx(i,j),tau) .ne. ms*sign(xx(i,j)-mtheta) )then
            E=E+1.0d0
        endif
    enddo
    enddo
    E=E/real(out_size*data_size,8)

    Ein(tid) = E-Emin

enddo
!$omp end parallel do

!$omp parallel do reduction(+:E)
do i = 1, num_of_test
    E = E + Ein(tid) / real(num_of_test,8)
enddo
!$omp end parallel do

write(*,*)E

deallocate(x,xx,theta,Ein)

end subroutine

end program