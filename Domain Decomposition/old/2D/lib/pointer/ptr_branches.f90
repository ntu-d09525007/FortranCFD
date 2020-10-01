module ptr_branches
!$ use omp_lib
implicit none

type pointer_vector_child
integer :: id, is, ie, js, je, ghc
real(8), dimension(:,:), pointer :: x, y, div
contains
end type pointer_vector_child 

type pointer_vector_parent
integer :: threads
type(pointer_vector_child),allocatable :: of(:)
contains
procedure init => ptrvecpart_init
procedure sync => ptrvecpart_sync
procedure check => ptrvecpart_chk
procedure release => ptrvecpart_release
end type pointer_vector_parent

contains

subroutine ptrvecpart_init(p,src)
use tree
implicit none
class(pointer_vector_parent) :: p 
type(manager) :: src
integer :: id

	p%threads = src%glb%threads
	
	allocate( p%of(0:p%threads-1) )
	
	!$omp parallel private(id), num_threads(p%threads)
	
		id = 0
		!$ id = omp_get_thread_num()
		
		p%of(id)%is = src%of(id)%loc%is
		p%of(id)%ie = src%of(id)%loc%ie
		p%of(id)%js = src%of(id)%loc%js
		p%of(id)%je = src%of(id)%loc%je
		
		p%of(id)%ghc = src%of(id)%glb%ghc
		
	!$omp end parallel 
	
end subroutine

subroutine ptrvecpart_chk(p)
implicit none
class(pointer_vector_parent) :: p
integer :: id

	!$omp parallel private(id), num_threads(p%threads)
	
		id = 0
		!$ id = omp_get_thread_num()
		
		if( .not. associated(p%of(id)%x) .or. .not. associated(p%of(id)%y) .or. .not. associated(p%of(id)%div) )then
			write(*,*)" The pointer is not associated. Stop the program. "
			stop
		endif
		
	!$omp end parallel

end subroutine

subroutine ptrvecpart_sync(p)
implicit none
class(pointer_vector_parent) :: p
integer :: id,i,j

	call p%check

	do id = p%threads-2, 0, -1
		do j = p%of(id)%js-p%of(id)%ghc, p%of(id)%je+p%of(id)%ghc
		do i = p%of(id)%ie+1, p%of(id)%ie+p%of(id)%ghc
			p%of(id)%x(i,j) = p%of(id+1)%x(i,j)
			p%of(id)%y(i,j) = p%of(id+1)%y(i,j)
			p%of(id)%div(i,j) = p%of(id+1)%div(i,j)
		end do
		end do
	end do
	
	do id = 1, p%threads-1
		do j = p%of(id)%js-p%of(id)%ghc, p%of(id)%je+p%of(id)%ghc
		do i = p%of(id)%is-1, p%of(id)%is-p%of(id)%ghc, -1
			p%of(id)%x(i,j) = p%of(id-1)%x(i,j)
			p%of(id)%y(i,j) = p%of(id-1)%y(i,j)
			p%of(id)%div(i,j) = p%of(id-1)%div(i,j)
		end do
		end do
	end do
	
end subroutine 


subroutine ptrvecpart_release(p)
implicit none
class(pointer_vector_parent) :: p
integer :: id
	
	!$omp parallel private(id), num_threads(p%threads)
		
		id=0
		!$ id = omp_get_thread_num()
		nullify( p%of(id)%x, p%of(id)%y, p%of(id)%div )
	
	!$omp end parallel 
	
end subroutine

end module ptr_branches