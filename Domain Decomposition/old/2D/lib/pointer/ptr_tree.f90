module ptr_tree
use ptr_roots
use ptr_branches
!$ use omp_lib
implicit none

type ptr_family
integer :: threads
type(pointer_parent) :: phi, p
type(pointer_vector_parent) :: vel, nvel, nveltmp
type(pointer_vector_parent) :: normals_1, normals_2
type(pointer_vector_parent) :: srk6_x, srk6_y
contains 
procedure init => ptr_family_init
end type ptr_family

contains

subroutine ptr_family_init(this,src)
use tree
implicit none
class(ptr_family) :: this
type(manager),target :: src
integer :: id

	this%threads = src%glb%threads
	
	call this%phi%init(src)

	call this%p%init(src)
	
	call this%vel%init(src)
	
	call this%nvel%init(src)
	call this%nveltmp%init(src)
	
	call this%normals_1%init(src)
	call this%normals_2%init(src)
		
	call this%srk6_x%init(src)
	call this%srk6_y%init(src)
	
	!$omp parallel private(id), num_threads(this%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		this%phi%of(id)%dat => src%of(id)%loc%phi%now
		
		this%p%of(id)%dat => src%of(id)%loc%p%now
		
		this%vel%of(id)%x => src%of(id)%loc%vel%x%now
		this%vel%of(id)%y => src%of(id)%loc%vel%y%now
		this%vel%of(id)%div => src%of(id)%loc%vel%div%now
	
		this%nvel%of(id)%x => src%of(id)%loc%nvel%x%now
		this%nvel%of(id)%y => src%of(id)%loc%nvel%y%now
		this%nvel%of(id)%div => src%of(id)%loc%nvel%div%now
		
		this%nveltmp%of(id)%x => src%of(id)%loc%nvel%x%tmp
		this%nveltmp%of(id)%y => src%of(id)%loc%nvel%y%tmp
		this%nveltmp%of(id)%div => src%of(id)%loc%nvel%div%tmp
		
		this%normals_1%of(id)%x => src%of(id)%loc%normals%x%now
		this%normals_1%of(id)%y => src%of(id)%loc%normals%y%now
		this%normals_1%of(id)%div => src%of(id)%loc%normals%xy%now
		
		this%normals_2%of(id)%x => src%of(id)%loc%normals%xx%now
		this%normals_2%of(id)%y => src%of(id)%loc%normals%yy%now
		this%normals_2%of(id)%div => src%of(id)%loc%normals%curv%now
		
		this%srk6_x%of(id)%x => src%of(id)%loc%srk6%x%s1
		this%srk6_x%of(id)%y => src%of(id)%loc%srk6%x%s2
		this%srk6_x%of(id)%div => src%of(id)%loc%srk6%x%s3
		
		this%srk6_y%of(id)%x => src%of(id)%loc%srk6%y%s1
		this%srk6_y%of(id)%y => src%of(id)%loc%srk6%y%s2
		this%srk6_y%of(id)%div => src%of(id)%loc%srk6%y%s3
		
	!$omp end parallel 

end subroutine

end module 

