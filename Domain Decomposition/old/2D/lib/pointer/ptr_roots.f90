module ptr_roots
!$ use omp_lib
implicit none

type pointer_child
integer :: id, is, ie, js, je, ghc
real(8), dimension(:,:), pointer :: dat
contains
end type pointer_child

type pointer_parent
integer :: threads
type(pointer_child),allocatable :: of(:)
contains
procedure init => ptrpart_init
procedure sync => ptrpart_sync
procedure check => ptrpart_chk
procedure release => ptrpart_release
end type pointer_parent

contains 

subroutine ptrpart_init(p,src)
use tree
implicit none
class(pointer_parent) :: p
type(manager) :: src
integer :: i, id
integer :: is, ie, js, je, ghc

	p%threads = src%glb%threads
	
	allocate( p%of(0:p%threads-1) )
	
	!$omp parallel private(id,is,ie,js,je,ghc), num_threads(p%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		p%of(id)%is = src%of(id)%loc%is
		p%of(id)%ie = src%of(id)%loc%ie
		
		p%of(id)%js = src%of(id)%loc%js
		p%of(id)%je = src%of(id)%loc%je
		
		p%of(id)%ghc = src%glb%ghc
	
	!$omp end parallel

end subroutine

subroutine ptrpart_chk(p) 
implicit none
class(pointer_parent) :: p
integer :: id

	!$omp parallel private(id), num_threads(p%threads) 
		
		id=0
		!$ id = omp_get_thread_num()
		
		if( .not. associated(p%of(id)%dat) )then
			write(*,*)" The pointer is not associated. Stop the program. "
			stop
		endif

	!$omp end parallel

end subroutine

subroutine ptrpart_release(p)
implicit none
class(pointer_parent) :: p
integer :: id
	
	!$omp parallel private(id), num_threads(p%threads) 
		
		id = 0
		!$ id = omp_get_thread_num()
		
		nullify(p%of(id)%dat)
		
	!$omp end parallel

end subroutine

subroutine ptrpart_sync(p)
implicit none
class(pointer_parent) :: p
integer :: id,i,j

	call p%check

	do id = p%threads-2, 0, -1
		do j = p%of(id)%js-p%of(id)%ghc, p%of(id)%je+p%of(id)%ghc
		do i = p%of(id)%ie+1, p%of(id)%ie+p%of(id)%ghc
			p%of(id)%dat(i,j) = p%of(id+1)%dat(i,j)
		end do
		end do
	end do
	
	do id = 1, p%threads-1
		do j = p%of(id)%js-p%of(id)%ghc, p%of(id)%je+p%of(id)%ghc
		do i = p%of(id)%is-1, p%of(id)%is-p%of(id)%ghc, -1
			p%of(id)%dat(i,j) = p%of(id-1)%dat(i,j)
		end do
		end do
	end do

end subroutine

end module ptr_roots



