module ptr_roots
!$ use omp_lib
use tree
implicit none

type pointer_child
integer :: idx, idy, idz
integer :: is, ie, js, je, ks, ke, ghc
real(8), dimension(:,:,:), pointer :: dat
contains
end type pointer_child

type pointer_parent
integer :: threads, gx, gy, gz
integer(8) :: cpurate
real(8) :: cputime
type(pointer_child),allocatable :: of(:)
contains
procedure init => ptrpart_init
procedure sync => ptrpart_sync
procedure check => ptrpart_chk
procedure release => ptrpart_release
procedure reset => ptrpart_reset
end type pointer_parent

type pointer_vector_parent
type(pointer_parent) :: x,y,z
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

subroutine ptrpart_init(p,src)
implicit none
class(pointer_parent) :: p
type(manager) :: src
integer :: id

	p%threads = src%glb%threads
	p%cputime = 0.0d0
	p%cpurate = src%glb%cpurate
	
	p%gx = src%glb%grid_x
	p%gy = src%glb%grid_y
	p%gz = src%glb%grid_z
	
	allocate( p%of(0:p%threads-1) )
	
	!$omp parallel private(id), num_threads(p%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		p%of(id)%is = src%of(id)%loc%is
		p%of(id)%ie = src%of(id)%loc%ie	
		p%of(id)%js = src%of(id)%loc%js
		p%of(id)%je = src%of(id)%loc%je
		p%of(id)%ks = src%of(id)%loc%ks
		p%of(id)%ke = src%of(id)%loc%ke
		
		p%of(id)%ghc = src%glb%ghc
		
		p%of(id)%idx = src%of(id)%loc%idx
		p%of(id)%idy = src%of(id)%loc%idy
		p%of(id)%idz = src%of(id)%loc%idz
	
	!$omp end parallel

end subroutine

subroutine ptrvecpart_init(p,src)
implicit none
class(pointer_vector_parent) :: p 
type(manager) :: src
integer :: id

	p%cputime = 0.0d0
	p%cpurate = src%glb%cpurate

	call p%x%init(src)
	call p%y%init(src)
	call p%z%init(src)
	
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

subroutine ptrvecpart_chk(p)
implicit none
class(pointer_vector_parent) :: p
integer :: id

	call p%x%check
	call p%y%check
	call p%z%check

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

subroutine ptrvecpart_release(p)
implicit none
class(pointer_vector_parent) :: p
integer :: id
	
	call p%x%release
	call p%y%release
	call p%z%release
	
end subroutine

subroutine ptrpart_sync(p)
implicit none
class(pointer_parent) :: p
integer :: id,i,j,k
integer(8) :: cpustart, cpuend

	call system_clock(cpustart)

	!call p%check

	!$omp parallel private(id,i,j,k), num_threads(p%threads)	
	
		id=0
		!$ id = omp_get_thread_num()
		
		! x direction
		if(p%of(id)%idx<p%gx-1)then
			
			do k = p%of(id)%ks-p%of(id)%ghc, p%of(id)%ke+p%of(id)%ghc
			do j = p%of(id)%js-p%of(id)%ghc, p%of(id)%je+p%of(id)%ghc
			do i = p%of(id)%ie+1, p%of(id)%ie+p%of(id)%ghc
				p%of(id)%dat(i,j,k) = p%of(id+1)%dat(i,j,k)
			end do
			end do
			end do
		
		endif

		if(p%of(id)%idx>0)then
			
			do k = p%of(id)%ks-p%of(id)%ghc, p%of(id)%ke+p%of(id)%ghc
			do j = p%of(id)%js-p%of(id)%ghc, p%of(id)%je+p%of(id)%ghc
			do i = p%of(id)%is-1, p%of(id)%is-p%of(id)%ghc, -1
				p%of(id)%dat(i,j,k) = p%of(id-1)%dat(i,j,k)
			end do
			end do
			end do
		
		endif
		
		! y direction
		if(p%of(id)%idy<p%gy-1)then
			
			do k = p%of(id)%ks-p%of(id)%ghc, p%of(id)%ke+p%of(id)%ghc
			do i = p%of(id)%is-p%of(id)%ghc, p%of(id)%ie+p%of(id)%ghc
			do j = p%of(id)%je+1, p%of(id)%je+p%of(id)%ghc
				p%of(id)%dat(i,j,k) = p%of(id+p%gx)%dat(i,j,k)
			end do
			end do
			end do
		
		endif
		
		if(p%of(id)%idy>0)then
			
			do k = p%of(id)%ks-p%of(id)%ghc, p%of(id)%ke+p%of(id)%ghc
			do i = p%of(id)%is-p%of(id)%ghc, p%of(id)%ie+p%of(id)%ghc
			do j = p%of(id)%js-1, p%of(id)%js-p%of(id)%ghc, -1
				p%of(id)%dat(i,j,k) = p%of(id-p%gx)%dat(i,j,k)
			end do
			end do
			end do
		
		endif
		
		! z direction
		if(p%of(id)%idz<p%gz-1)then
		
			do j = p%of(id)%js-p%of(id)%ghc, p%of(id)%je+p%of(id)%ghc
			do i = p%of(id)%is-p%of(id)%ghc, p%of(id)%ie+p%of(id)%ghc		
			do k = p%of(id)%ke+1, p%of(id)%ke+p%of(id)%ghc
				p%of(id)%dat(i,j,k) = p%of(id+p%gx*p%gy)%dat(i,j,k)
			end do
			end do
			end do
		
		endif

		if(p%of(id)%idz>0)then
			
			do j = p%of(id)%js-p%of(id)%ghc, p%of(id)%je+p%of(id)%ghc
			do i = p%of(id)%is-p%of(id)%ghc, p%of(id)%ie+p%of(id)%ghc
			do k = p%of(id)%ks-1, p%of(id)%ks-p%of(id)%ghc, -1
				p%of(id)%dat(i,j,k) = p%of(id-p%gx*p%gy)%dat(i,j,k)
			end do
			end do
			end do
		
		endif

	!$omp end parallel
	
	call system_clock(cpuend)
	p%cputime = p%cputime + real(cpuend-cpustart,kind=8)/real(p%cpurate,kind=8)

end subroutine

subroutine ptrvecpart_sync(p)
implicit none
class(pointer_vector_parent) :: p
integer :: id,i,j,k
integer(8) :: cpustart, cpuend

	call system_clock(cpustart)

	!call p%check

	call p%x%sync
	call p%y%sync
	call p%z%sync
	
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

	p%cputime=0.0d0
	call p%x%reset
	call p%y%reset
	call p%z%reset

end subroutine

end module ptr_roots



