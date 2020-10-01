subroutine ptrpart_mpi_assign(p)
implicit none
class(pointer_parent) :: p
integer :: ix,iy,iz,i,j,k

!              o-------o
!             /|      /|
! send_m     / |     / |    send_p
! ------ <= o--|----o  | => -------
! recv_p    |  o----|--o    recv_m
!           | /     | /
!           |/      |/
!           o-------o
!       i-1/2   i+1/2     

if(p%split==1)then

    ix=p%grids(1,1)
    do iy = p%grids(2,1), p%grids(2,2)
    do iz = p%grids(3,1), p%grids(3,2)
        do i = 1, p%ghc
        do j = p%of(ix,iy,iz)%js, p%of(ix,iy,iz)%je
        do k = p%of(ix,iy,iz)%ks, p%of(ix,iy,iz)%ke
            p%mpi_sendm(i,j,k) = p%of(ix,iy,iz)%dat(i,j,k)
        enddo
        enddo
        enddo
    enddo
    enddo

    ix=p%grids(1,2)
    do iy = p%grids(2,1), p%grids(2,2)
    do iz = p%grids(3,1), p%grids(3,2)
        do i = 1, p%ghc
        do j = p%of(ix,iy,iz)%js, p%of(ix,iy,iz)%je
        do k = p%of(ix,iy,iz)%ks, p%of(ix,iy,iz)%ke
            p%mpi_sendp(i,j,k) = p%of(ix,iy,iz)%dat(p%node_x-p%ghc+i,j,k)
        enddo
        enddo
        enddo
    enddo
    enddo

else if (p%split==2)then

    iy=p%grids(2,1)
    do ix = p%grids(1,1), p%grids(1,2)
    do iz = p%grids(3,1), p%grids(3,2)
        do j = 1, p%ghc
        do i = p%of(ix,iy,iz)%is, p%of(ix,iy,iz)%ie
        do k = p%of(ix,iy,iz)%ks, p%of(ix,iy,iz)%ke
            p%mpi_sendm(i,j,k) = p%of(ix,iy,iz)%dat(i,j,k)
        enddo
        enddo
        enddo
    enddo
    enddo

    iy=p%grids(2,2)
    do ix = p%grids(1,1), p%grids(1,2)
    do iz = p%grids(3,1), p%grids(3,2)
        do j = 1, p%ghc
        do i = p%of(ix,iy,iz)%is, p%of(ix,iy,iz)%ie
        do k = p%of(ix,iy,iz)%ks, p%of(ix,iy,iz)%ke
            p%mpi_sendp(i,j,k) = p%of(ix,iy,iz)%dat(i,p%node_y-p%ghc+j,k)
        enddo
        enddo
        enddo
    enddo
    enddo

else

    iz=p%grids(3,1)
    do ix = p%grids(1,1), p%grids(1,2)
    do iy = p%grids(2,1), p%grids(2,2)
        do k = 1, p%ghc
        do i = p%of(ix,iy,iz)%is, p%of(ix,iy,iz)%ie
        do j = p%of(ix,iy,iz)%js, p%of(ix,iy,iz)%je
            p%mpi_sendm(i,j,k) = p%of(ix,iy,iz)%dat(i,j,k)
        enddo
        enddo
        enddo
    enddo
    enddo

    iz=p%grids(3,2)
    do ix = p%grids(1,1), p%grids(1,2)
    do iy = p%grids(2,1), p%grids(2,2)
        do k = 1, p%ghc
        do i = p%of(ix,iy,iz)%is, p%of(ix,iy,iz)%ie
        do j = p%of(ix,iy,iz)%js, p%of(ix,iy,iz)%je
            p%mpi_sendm(i,j,k) = p%of(ix,iy,iz)%dat(i,j,p%node_z-p%ghc+k)
        enddo
        enddo
        enddo
    enddo
    enddo

endif

end subroutine

subroutine ptrpart_mpi_sync(p)
implicit none
include 'mpif.h'
class(pointer_parent) :: p
integer :: n,istat,ierr

!              o-------o
!             /|      /|
! send_m     / |     / |    send_p
! ------ <= o--|----o  | => -------
! recv_p    |  o----|--o    recv_m
!           | /     | /
!           |/      |/
!           o-------o
!       i-1/2   i+1/2     


n = p%ghc*p%mpidim(1)*p%mpidim(2)

!---------------------------------------------------------

if( p%mpirank < p%mpisize-1 )then
    p%mpi_buffer = reshape(p%mpi_sendp,(/n/))
    call mpi_send(p%mpi_buffer, n, mpi_real8, p%mpirank+1, 100*p%mpirank+1, mpi_comm_world, ierr)
endif

if( p%mpirank > 0 )then
    call mpi_recv(p%mpi_buffer2, n, mpi_real8, p%mpirank-1, 100*(p%mpirank-1)+1, mpi_comm_world, istat, ierr)
    p%mpi_recvp = reshape( p%mpi_buffer2, (/p%ghc,p%mpidim(1),p%mpidim(2)/) )
endif

!---------------------------------------------------------

if( p%mpirank > 0 )then
    p%mpi_buffer = reshape(p%mpi_sendm,(/n/))
    call mpi_send(p%mpi_buffer, n, mpi_real8, p%mpirank-1, 200*p%mpirank+2, mpi_comm_world, ierr)
endif

if( p%mpirank < p%mpisize-1 )then
    call mpi_recv(p%mpi_buffer2, n, mpi_real8, p%mpirank+1, 200*(p%mpirank+1)+2, mpi_comm_world, istat, ierr)
    p%mpi_recvm = reshape( p%mpi_buffer2, (/p%ghc,p%mpidim(1),p%mpidim(2)/) )
endif

!---------------------------------------------------------

call mpi_barrier(mpi_comm_world, ierr)

end subroutine

subroutine ptrpart_mpi_finalize(p)
implicit none
include 'mpif.h'
class(pointer_parent) :: p
integer :: ix,iy,iz,i,j,k
integer :: ierr

!              o-------o
!             /|      /|
! send_m     / |     / |    send_p
! ------ <= o--|----o  | => -------
! recv_p    |  o----|--o    recv_m
!           | /     | /
!           |/      |/
!           o-------o
!       i-1/2   i+1/2     


if(p%split==1)then

    if( p%mpirank>0)then

        ix=p%grids(1,1)
        do iy = p%grids(2,1), p%grids(2,2)
        do iz = p%grids(3,1), p%grids(3,2)
            do i = 1, p%ghc
            do j = p%of(ix,iy,iz)%js, p%of(ix,iy,iz)%je
            do k = p%of(ix,iy,iz)%ks, p%of(ix,iy,iz)%ke
                p%of(ix,iy,iz)%dat(p%of(ix,iy,iz)%ie+i,j,k) = p%mpi_recvm(i,j,k)
            enddo
            enddo
            enddo
        enddo
        enddo

    endif

    if( p%mpirank < p%mpisize-1 )then

        ix=p%grids(1,2)
        do iy = p%grids(2,1), p%grids(2,2)
        do iz = p%grids(3,1), p%grids(3,2)
            do i = 1, p%ghc
            do j = p%of(ix,iy,iz)%js, p%of(ix,iy,iz)%je
            do k = p%of(ix,iy,iz)%ks, p%of(ix,iy,iz)%ke
                p%of(ix,iy,iz)%dat(-p%ghc+i,j,k) = p%mpi_recvp(i,j,k)
            enddo
            enddo
            enddo
        enddo
        enddo

    endif

else if (p%split==2)then

    if( p%mpirank > 0 )then

        iy=p%grids(2,1)
        do ix = p%grids(1,1), p%grids(1,2)
        do iz = p%grids(3,1), p%grids(3,2)
            do j = 1, p%ghc
            do i = p%of(ix,iy,iz)%is, p%of(ix,iy,iz)%ie
            do k = p%of(ix,iy,iz)%ks, p%of(ix,iy,iz)%ke
                p%of(ix,iy,iz)%dat(i,p%of(ix,iy,iz)%je+j,k) = p%mpi_recvm(j,i,k)
            enddo
            enddo
            enddo
        enddo
        enddo

    endif

    if( p%mpirank < p%mpisize-1 )then

        iy=p%grids(2,2)
        do ix = p%grids(1,1), p%grids(1,2)
        do iz = p%grids(3,1), p%grids(3,2)
            do j = 1, p%ghc
            do i = p%of(ix,iy,iz)%is, p%of(ix,iy,iz)%ie
            do k = p%of(ix,iy,iz)%ks, p%of(ix,iy,iz)%ke
                p%of(ix,iy,iz)%dat(i,-p%ghc+j,k) = p%mpi_recvp(j,i,k)
            enddo
            enddo
            enddo
        enddo
        enddo

    endif

else

    if( p%mpirank > 0 )then

        iz=p%grids(3,1)
        do ix = p%grids(1,1), p%grids(1,2)
        do iy = p%grids(2,1), p%grids(2,2)
            do k = 1, p%ghc
            do i = p%of(ix,iy,iz)%is, p%of(ix,iy,iz)%ie
            do j = p%of(ix,iy,iz)%js, p%of(ix,iy,iz)%je
                p%of(ix,iy,iz)%dat(i,j,p%of(ix,iy,iz)%ke+k) = p%mpi_recvm(k,i,j)
            enddo
            enddo
            enddo
        enddo
        enddo

    endif

    if( p%mpirank < p%mpisize-1 )then

        iz=p%grids(3,2)
        do ix = p%grids(1,1), p%grids(1,2)
        do iy = p%grids(2,1), p%grids(2,2)
            do k = 1, p%ghc
            do i = p%of(ix,iy,iz)%is, p%of(ix,iy,iz)%ie
            do j = p%of(ix,iy,iz)%js, p%of(ix,iy,iz)%je
                p%of(ix,iy,iz)%dat(i,j,-p%ghc+k) = p%mpi_recvp(k,i,j)
            enddo
            enddo
            enddo
        enddo
        enddo

    endif

endif

call mpi_barrier(mpi_comm_world, ierr)

end subroutine
