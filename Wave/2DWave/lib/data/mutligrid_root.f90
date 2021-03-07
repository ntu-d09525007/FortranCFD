module mutligrid_root

type multigrid_root
integer :: nx, ny, n
integer :: idx, idy, gx, gy
real(8)  :: dx, dy
real(8), dimension(:,:), allocatable :: sol, res, pol, src
real(8), dimension(:,:), allocatable :: node, node2
real(8), dimension(:), allocatable :: i, j
real(8) :: L2norm, L2norm0
contains
procedure :: init => mg_init
end type multigrid_root

contains

subroutine mg_init(p)
implicit none
class(multigrid_root) :: p
integer :: i,j

allocate(p%sol(0:p%nx+1,0:p%ny+1), p%src(0:p%nx+1,0:p%ny+1), &
        &p%res(0:p%nx+1,0:p%ny+1), p%pol(0:p%nx+1,0:p%ny+1) )
            
p%n = p%nx * p%ny
allocate( p%node(p%nx,p%ny), p%i(p%n), p%j(p%n) )

do j = 1, p%ny
do i = 1, p%nx
    p%node(i,j) = (j-1)*p%nx + i
    p%i(p%node(i,j)) = i
    p%j(p%node(i,j)) = j
enddo
enddo

end subroutine

end module mutligrid_root
