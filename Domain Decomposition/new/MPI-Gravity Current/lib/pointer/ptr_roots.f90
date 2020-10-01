module ptr_roots
!$ use omp_lib
use tree
implicit none

type pointer_child
integer :: is, ie, js, je, ks, ke, ghc
real(8), dimension(:,:,:), pointer :: dat
contains
end type pointer_child

type pointer_parent
integer :: grids(3,2)
integer(8) :: cpurate
real(8) :: cputime
type(pointer_child),allocatable :: of(:,:,:)
real(8),allocatable :: mpi_sendp(:,:,:), mpi_recvp(:,:,:), mpi_sendm(:,:,:), mpi_recvm(:,:,:), mpi_buffer(:), mpi_buffer2(:)
integer :: mpirank, mpisize, mpidim(2), split
integer :: node_x, node_y, node_z, ghc
contains
procedure init => ptrpart_init
procedure sync => ptrpart_sync
procedure check => ptrpart_chk
procedure release => ptrpart_release
procedure reset => ptrpart_reset
procedure mpi_assign => ptrpart_mpi_assign
procedure mpi_sync => ptrpart_mpi_sync
procedure mpi_finalize => ptrpart_mpi_finalize
end type pointer_parent

type pointer_vector_parent
type(pointer_parent),allocatable :: nodes(:)
integer :: num
integer(8) :: cpurate
real(8) :: cputime
contains
procedure init => ptrvecpart_init
procedure sync => ptrvecpart_sync
procedure check => ptrvecpart_chk
procedure release => ptrvecpart_release
procedure reset => ptrvecpart_reset
end type pointer_vector_parent

contains 

include 'ptr_mpi.f90'

subroutine ptrpart_init(p,src)
implicit none
class(pointer_parent) :: p
type(manager) :: src
integer :: ix,iy,iz,id,m,n,ghc

    p%cputime = 0.0d0
    p%cpurate = src%glb%cpurate
    
    p%grids(1,1) = src%glb%grids(1,1)
    p%grids(1,2) = src%glb%grids(1,2)
    p%grids(2,1) = src%glb%grids(2,1)
    p%grids(2,2) = src%glb%grids(2,2)
    p%grids(3,1) = src%glb%grids(3,1)
    p%grids(3,2) = src%glb%grids(3,2)

    p%ghc = src%glb%ghc

    p%mpidim(1) = src%glb%mpidim(1)
    p%mpidim(2) = src%glb%mpidim(2)

    p%node_x = src%glb%node_x
    p%node_y = src%glb%node_y
    p%node_z = src%glb%node_z
    
    allocate( p%of(p%grids(1,1):p%grids(1,2),p%grids(2,1):p%grids(2,2),p%grids(3,1):p%grids(3,2)) )
    
    !$omp parallel do collapse(3), private(id)
    do ix = p%grids(1,1), p%grids(1,2)
    do iy = p%grids(2,1), p%grids(2,2)
    do iz = p%grids(3,1), p%grids(3,2)
        
        id = src%glb%tid(ix,iy,iz)

        p%of(ix,iy,iz)%is = src%of(id)%loc%is
        p%of(ix,iy,iz)%ie = src%of(id)%loc%ie 
        p%of(ix,iy,iz)%js = src%of(id)%loc%js
        p%of(ix,iy,iz)%je = src%of(id)%loc%je
        p%of(ix,iy,iz)%ks = src%of(id)%loc%ks
        p%of(ix,iy,iz)%ke = src%of(id)%loc%ke
        
        p%of(ix,iy,iz)%ghc = src%glb%ghc
        
    enddo
    enddo
    enddo
    !$omp end parallel do

    p%mpirank = src%glb%mpirank
    p%mpisize = src%glb%mpisize
    p%split = src%glb%split

    m=p%mpidim(1); n=p%mpidim(2); ghc=p%ghc

    allocate(p%mpi_sendp(ghc,m,n),p%mpi_recvp(ghc,m,n))
    allocate(p%mpi_sendm(ghc,m,n),p%mpi_recvm(ghc,m,n))
    allocate(p%mpi_buffer(ghc*m*n), p%mpi_buffer2(ghc*m*n) )

end subroutine

subroutine ptrvecpart_init(p,src,num)
implicit none
class(pointer_vector_parent) :: p 
type(manager) :: src
integer :: num, comp

    p%cputime = 0.0d0
    p%cpurate = src%glb%cpurate
    p%num = num
    
    allocate( p%nodes(1:p%num) )
    
    do comp = 1, p%num
        call p%nodes(comp)%init(src)
    enddo
    
end subroutine


subroutine ptrpart_chk(p) 
implicit none
class(pointer_parent) :: p
integer :: ix,iy,iz
    
    !$omp parallel do collapse(3)
    do iz = p%grids(3,1), p%grids(3,2)
    do iy = p%grids(2,1), p%grids(2,2)
    do ix = p%grids(1,1), p%grids(1,2)

        if( .not. associated(p%of(ix,iy,iz)%dat) )then
            write(*,*)" The pointer is not associated. Stop the program. "
            stop
        endif
    
    enddo
    enddo
    enddo    
    !$omp end parallel do

end subroutine

subroutine ptrvecpart_chk(p)
implicit none
class(pointer_vector_parent) :: p
integer :: comp

    do comp = 1, p%num
        call p%nodes(comp)%check
    enddo
    
end subroutine

subroutine ptrpart_release(p)
implicit none
class(pointer_parent) :: p
integer :: ix,iy,iz
    
    !$omp parallel do collapse(3)
    do iz = p%grids(3,1), p%grids(3,2)
    do iy = p%grids(2,1), p%grids(2,2)
    do ix = p%grids(1,1), p%grids(1,2)

        nullify(p%of(ix,iy,iz)%dat)
    
    enddo
    enddo
    enddo    
    !$omp end parallel do

end subroutine

subroutine ptrvecpart_release(p)
implicit none
class(pointer_vector_parent) :: p
integer :: comp

    do comp = 1, p%num
        call p%nodes(comp)%release
    enddo
    
end subroutine

subroutine ptrpart_sync(p)
implicit none
class(pointer_parent) :: p
integer :: ix,iy,iz,i,j,k
integer(8) :: cpustart, cpuend
integer :: n

    call system_clock(cpustart)

    call p%check

    !$omp parallel do collapse(3), private(i,j,k)
    do ix = p%grids(1,1), p%grids(1,2)
    do iy = p%grids(2,1), p%grids(2,2)
    do iz = p%grids(3,1), p%grids(3,2)
        
        ! x direction
        if( ix<p%grids(1,2) )then
            
            do k = p%of(ix,iy,iz)%ks-p%of(ix,iy,iz)%ghc, p%of(ix,iy,iz)%ke+p%of(ix,iy,iz)%ghc
            do j = p%of(ix,iy,iz)%js-p%of(ix,iy,iz)%ghc, p%of(ix,iy,iz)%je+p%of(ix,iy,iz)%ghc
            do i = p%of(ix,iy,iz)%ie+1, p%of(ix,iy,iz)%ie+p%of(ix,iy,iz)%ghc
                p%of(ix,iy,iz)%dat(i,j,k) = p%of(ix+1,iy,iz)%dat(i,j,k)
            end do
            end do
            end do
        
        endif

        if( ix>p%grids(1,1) )then
            
            do k = p%of(ix,iy,iz)%ks-p%of(ix,iy,iz)%ghc, p%of(ix,iy,iz)%ke+p%of(ix,iy,iz)%ghc
            do j = p%of(ix,iy,iz)%js-p%of(ix,iy,iz)%ghc, p%of(ix,iy,iz)%je+p%of(ix,iy,iz)%ghc
            do i = p%of(ix,iy,iz)%is-1, p%of(ix,iy,iz)%is-p%of(ix,iy,iz)%ghc, -1
                p%of(ix,iy,iz)%dat(i,j,k) = p%of(ix-1,iy,iz)%dat(i,j,k)
            end do
            end do
            end do
        
        endif
        
        ! y direction
        if( iy<p%grids(2,2) )then
            
            do k = p%of(ix,iy,iz)%ks-p%of(ix,iy,iz)%ghc, p%of(ix,iy,iz)%ke+p%of(ix,iy,iz)%ghc
            do i = p%of(ix,iy,iz)%is-p%of(ix,iy,iz)%ghc, p%of(ix,iy,iz)%ie+p%of(ix,iy,iz)%ghc
            do j = p%of(ix,iy,iz)%je+1, p%of(ix,iy,iz)%je+p%of(ix,iy,iz)%ghc
                p%of(ix,iy,iz)%dat(i,j,k) = p%of(ix,iy+1,iz)%dat(i,j,k)
            end do
            end do
            end do
        
        endif
        
        if( iy>p%grids(2,1) )then
            
            do k = p%of(ix,iy,iz)%ks-p%of(ix,iy,iz)%ghc, p%of(ix,iy,iz)%ke+p%of(ix,iy,iz)%ghc
            do i = p%of(ix,iy,iz)%is-p%of(ix,iy,iz)%ghc, p%of(ix,iy,iz)%ie+p%of(ix,iy,iz)%ghc
            do j = p%of(ix,iy,iz)%js-1, p%of(ix,iy,iz)%js-p%of(ix,iy,iz)%ghc, -1
                p%of(ix,iy,iz)%dat(i,j,k) = p%of(ix,iy-1,iz)%dat(i,j,k)
            end do
            end do
            end do
        
        endif
        
        ! z direction
        if( iz<p%grids(3,2) )then
        
            do j = p%of(ix,iy,iz)%js-p%of(ix,iy,iz)%ghc, p%of(ix,iy,iz)%je+p%of(ix,iy,iz)%ghc
            do i = p%of(ix,iy,iz)%is-p%of(ix,iy,iz)%ghc, p%of(ix,iy,iz)%ie+p%of(ix,iy,iz)%ghc       
            do k = p%of(ix,iy,iz)%ke+1, p%of(ix,iy,iz)%ke+p%of(ix,iy,iz)%ghc
                p%of(ix,iy,iz)%dat(i,j,k) = p%of(ix,iy,iz+1)%dat(i,j,k)
            end do
            end do
            end do
        
        endif

        if( iz>p%grids(3,1) )then
            
            do j = p%of(ix,iy,iz)%js-p%of(ix,iy,iz)%ghc, p%of(ix,iy,iz)%je+p%of(ix,iy,iz)%ghc
            do i = p%of(ix,iy,iz)%is-p%of(ix,iy,iz)%ghc, p%of(ix,iy,iz)%ie+p%of(ix,iy,iz)%ghc
            do k = p%of(ix,iy,iz)%ks-1, p%of(ix,iy,iz)%ks-p%of(ix,iy,iz)%ghc, -1
                p%of(ix,iy,iz)%dat(i,j,k) = p%of(ix,iy,iz-1)%dat(i,j,k)
            end do
            end do
            end do
        
        endif

    enddo
    enddo
    enddo
    !$omp end parallel do
    
    if( p%mpisize > 1)then

        call p%mpi_assign
        call p%mpi_sync
        call p%mpi_finalize

    endif

    call system_clock(cpuend)
    p%cputime = p%cputime + real(cpuend-cpustart,kind=8)/real(p%cpurate,kind=8)

end subroutine

subroutine ptrvecpart_sync(p)
implicit none
class(pointer_vector_parent) :: p
integer(8) :: cpustart, cpuend
integer :: comp
    
    call system_clock(cpustart)

    call p%check

    do comp = 1, p%num
        call p%nodes(comp)%sync
    enddo
    
    call system_clock(cpuend)
    p%cputime = p%cputime + real(cpuend-cpustart,kind=8)/real(p%cpurate,kind=8)
    
end subroutine 

subroutine ptrpart_reset(p)
implicit none
class(pointer_parent) :: p

    p%cputime=0.0d0

end subroutine

subroutine ptrvecpart_reset(p)
implicit none
class(pointer_vector_parent) :: p
integer :: comp

    p%cputime=0.0d0
    do comp = 1, p%num
        call p%nodes(comp)%reset
    enddo

end subroutine

end module ptr_roots



