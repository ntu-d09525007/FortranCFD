module mutligrid_root

type multigrid_root
integer :: nx, ny, nz, n
integer :: idx, idy, idz, gx, gy, gz
real(8)  :: dx, dy, dz
integer, dimension(:,:,:), allocatable :: node
integer, dimension(:), allocatable :: i, j, k, ipiv
real(8), dimension(:,:,:), allocatable :: sol, res, pol, src
real(8),dimension(:,:),allocatable :: A, AA
real(8),dimension(:),allocatable :: B, work
real(8) :: L2norm, L2norm0
contains
procedure :: init => mg_init
procedure :: solve => mg_solve
end type multigrid_root

contains

subroutine mg_init(p)
implicit none
class(multigrid_root) :: p
integer :: i,j,k
real(8) :: dx,dy,dz

p%n = p%nx * p%ny * p%nz

dx = p%dx;dy = p%dy;dz = p%dz

allocate( p%sol(0:p%nx+1,0:p%ny+1,0:p%nz+1), p%src(0:p%nx+1,0:p%ny+1,0:p%nz+1) )
allocate( p%res(0:p%nx+1,0:p%ny+1,0:p%nz+1), p%pol(0:p%nx+1,0:p%ny+1,0:p%nz+1) )
allocate( p%B(p%n), p%work(p%n))
allocate( p%A(p%n,p%n), p%AA(p%n,p%n) )
allocate( p%node(p%nx,p%ny,p%nz), p%i(p%n), p%j(p%n), p%k(p%n), p%ipiv(p%n))

do k = 1, p%nz
do j = 1, p%ny
do i = 1, p%nx
    p%node(i,j,k) = (k-1)*p%nx*p%ny + (j-1)*p%nx + i
    p%i(p%node(i,j,k)) = i
    p%j(p%node(i,j,k)) = j
    p%k(p%node(i,j,k)) = k
enddo
enddo
enddo

do i = 1, p%n

    p%A(i,i) = -2.0d0/dx**2.0d0 -2.0d0/dy**2.0d0 -2.0d0/dz**2.0d0

    if( p%i(i)>1 )then
        p%A(i,p%node(p%i(i)-1,p%j(i),p%k(i))) = 1.0d0 / dx**2.0d0
    else !if ( p%idx==0 )then
        p%A(i,i) = p%A(i,i) + 1.0d0 / dx**2.0d0
    endif

    if( p%i(i)<p%nx )then
        p%A(i,p%node(p%i(i)+1,p%j(i),p%k(i))) = 1.0d0 / dx**2.0d0
    else !if ( p%idx==p%gx-1 )then
        p%A(i,i) = p%A(i,i) + 1.0d0 / dx**2.0d0
    endif

    if( p%j(i)>1 )then
        p%A(i,p%node(p%i(i),p%j(i)-1,p%k(i))) = 1.0d0 / dy**2.0d0
    else !if ( p%idy==0 )then
        p%A(i,i) = p%A(i,i) + 1.0d0 / dx**2.0d0
    endif

    if( p%j(i)<p%ny )then
        p%A(i,p%node(p%i(i),p%j(i)+1,p%k(i))) = 1.0d0 / dy**2.0d0
    else !if ( p%idy==p%gy-1 )then
        p%A(i,i) = p%A(i,i) + 1.0d0 / dx**2.0d0
    endif

    if( p%k(i)>1 )then
        p%A(i,p%node(p%i(i),p%j(i),p%k(i)-1)) = 1.0d0 / dz**2.0d0
    else !if ( p%idz==0 )then
        p%A(i,i) = p%A(i,i) + 1.0d0 / dx**2.0d0
    endif

    if( p%k(i)<p%nz )then
        p%A(i,p%node(p%i(i),p%j(i),p%k(i)+1)) = 1.0d0 / dz**2.0d0
    else !if ( p%idz==p%gz-1 )then
        p%A(i,i) = p%A(i,i) + 1.0d0 / dx**2.0d0
    endif

enddo

end subroutine

subroutine mg_solve(p)
implicit none
class(multigrid_root) :: p
integer :: i,j,info
real(8) :: dx,dy,dz

dx = p%dx
dy = p%dy
dz = p%dz

do i = 1, p%n

    p%B(i) = p%src(p%i(i),p%j(i),p%k(i))

    do j = 1, p%n
        p%AA(i,j) = p%A(i,j)
    enddo

    ! if( p%i(i) == 1 .and. p%idx.ne.0 )then
    !     p%B(i) = p%B(i) - p%sol(p%i(i)-1,p%j(i),p%k(i))/dx**2.0d0
    ! endif

    ! if( p%i(i) == p%nx .and. p%idx.ne.p%gx-1 )then
    !     p%B(i) = p%B(i) - p%sol(p%i(i)+1,p%j(i),p%k(i))/dx**2.0d0
    ! endif

    ! if( p%j(i) == 1 .and. p%idy.ne.0 )then
    !     p%B(i) = p%B(i) - p%sol(p%i(i),p%j(i)-1,p%k(i))/dy**2.0d0
    ! endif

    ! if( p%j(i) == p%ny .and. p%idy.ne.p%gy-1 )then
    !     p%B(i) = p%B(i) - p%sol(p%i(i),p%j(i)+1,p%k(i))/dy**2.0d0
    ! endif

    ! if( p%k(i) == 1 .and. p%idz.ne.0 )then
    !     p%B(i) = p%B(i) - p%sol(p%i(i),p%j(i),p%k(i)-1)/dz**2.0d0
    ! endif

    ! if( p%k(i) == p%nz .and. p%idz.ne.p%gz-1 )then
    !     p%B(i) = p%B(i) - p%sol(p%i(i),p%j(i),p%k(i)+1)/dz**2.0d0
    ! endif

enddo

call dgesv(p%n,1,p%AA,p%n,p%ipiv,p%B,p%n,info)
if(info.ne.0)stop "Something Wrong with Multigrid Solver!!"

do i = 1, p%n
    p%sol(p%i(i),p%j(i),p%k(i)) = p%B(i)
enddo

end subroutine

end module mutligrid_root
