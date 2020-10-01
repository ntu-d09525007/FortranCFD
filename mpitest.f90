program main
use omp_lib
implicit none
include 'mpif.h'
integer :: id, ierror,istat,threads,maxthreads(2),maxt
integer,dimension(3,3) :: dat
integer :: buffer(9),num
integer :: i,j


call MPI_INIT(ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierror)
threads = omp_get_max_threads()

if( id==0 )then

do i = 1, 3
do j = 1, 3
    dat(i,j) = i+j-1
enddo
enddo

else

dat=0
endif

write(*,*)">>",id
do i = 1, 3
    write(*,*)dat(i,:)
enddo

buffer = reshape(dat,(/9/))

CALL MPI_GATHER(threads,1,MPI_INTEGER, &
             maxthreads,1,MPI_INTEGER, &
             0, MPI_COMM_WORLD, ierror)
buffer = sum(maxthreads)

call mpi_sendrecv(buffer,1,MPI_INT,1,0, maxt,1,MPI_INT,0,0,  MPI_COMM_WORLD,istat,ierror)

if( id==0 )then
    CALL MPI_SEND(buffer,9,MPI_INT,1,0,MPI_COMM_WORLD,ierror)
else
    CALL MPI_RECV(buffer,9,MPI_INT,0,0,MPI_COMM_WORLD,istat,ierror)
endif

if( id==1 )then
dat=reshape(buffer,(/3,3/))
write(*,*)">>",id
do i = 1, 3
    write(*,*)dat(i,:)
enddo
endif

CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

write(*,*)id,threads,maxt,"Hello World"

CALL MPI_FINALIZE()

end program
