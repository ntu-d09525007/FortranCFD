include 'sorts.f90'

program main
use omp_lib
implicit none
integer :: num_of_test, data_size
real(8) :: data_range

data_size = 12000
num_of_test = 10000

call dec_stump(num_of_test,data_size,2.0d0,0.0d0)
call dec_stump(num_of_test,data_size,20.0d0,0.0d0)

call dec_stump(num_of_test,data_size,20.0d0,0.1d0)
call dec_stump(num_of_test,data_size,200.0d0,0.1d0)

contains

function sign2(x,tau)
implicit none
integer :: sign2
real(8) :: x,tau,rand

if( tau < 1.0d-10 )then
    sign2=sign(x)
    return
endif

call random_number(rand)
    
if( rand > tau )then

    if(x>0)then
        sign2 = 1
    else
        sign2 = -1
    endif

else

    if(x>0)then
        sign2 = -1
    else
        sign2 = 1
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
integer :: num_of_test,data_size,out_size,threads
integer :: tid, i, j, id, s, ms
real(8) :: data_range, tau, randnum
real(8) :: E, Emin, mtheta
real(8),allocatable :: x(:),theta(:),Ein(:)

out_size = 10

threads = omp_get_max_threads()

call omp_set_dynamic(.false.)
call omp_set_num_threads(threads)

allocate(x(data_size),theta(data_size),Ein(num_of_test))

call random_seed()

do tid = 1, num_of_test

    !$omp parallel do private(randnum)
    do i = 1, data_size
        call random_number(randnum)
        x(i) = data_range*randnum-data_range*0.5d0
    enddo
    !$omp end parallel do

    call heap_sort(x(:))

    theta(1) = -1.0d0
    !$omp parallel do 
    do i = 2, data_size
        theta(i) = 0.5d0*( x(i)+x(i-1) )
    enddo
    !$omp end parallel do

    ! Training

    Emin=data_size
    do s = 1, -1, -2
    do j = 1, data_size
        E = 0.0d0
        !$omp parallel do reduction(+:E)
        do i = 1, data_size
            if( sign2(x(i),tau) .ne. s*sign(x(i)-theta(j)) )then
                E = E + 1.0d0
            endif
        enddo
        !$omp end parallel do
        E = E / real( data_size, 8)                                                             
        if( E<Emin )then
            Emin=E
            ms=s
            mtheta = theta(j)
        endif
    enddo
    enddo

    !write(*,'("Emin:",ES11.4,"  s:",I5,"  theta:",ES11.4)')Emin,ms,mtheta

    ! Calculate Eout
    E=0.0d0
    !$omp parallel do reduction(+:E), private(randnum)
    do i = 1, out_size*data_size
        call random_number(randnum)
        randnum = data_range*randnum-data_range*0.5d0
        if( sign(randnum) .ne. ms*sign(randnum-mtheta) )then
            E=E+1.0d0
        endif
    enddo
    !$omp end parallel do

    E = E / real(out_size*data_size,8)

    Ein(tid) = E-Emin

    !write(*,*)E, Ein(tid)

enddo

E=0.0d0
!$omp parallel do reduction(+:E)
do i = 1, num_of_test
    E = E + Ein(tid) / real(num_of_test,8)
enddo
!$omp end parallel do

write(*,*)"Ein-Eout",E

end subroutine

end program