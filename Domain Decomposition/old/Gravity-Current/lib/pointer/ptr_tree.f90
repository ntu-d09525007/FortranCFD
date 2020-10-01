module ptr_tree
use ptr_roots
!$ use omp_lib
implicit none

type ptr_family
integer :: threads
type(pointer_parent) :: phi, p, curv, vof
type(pointer_vector_parent) :: vel, nvel, velsrc_old, velsrc, nvel_old
type(pointer_vector_parent) :: normals_1, normals_2, normals_3
type(pointer_vector_parent) :: tdatax, tdatay, tdataz
contains 
procedure init => ptr_family_init
procedure cputime => ptr_family_cputime
procedure reset => ptr_family_reset
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
	call this%vof%init(src)
	call this%p%init(src)
	call this%curv%init(src)
	
	call this%vel%init(src)
	call this%velsrc_old%init(src)
	call this%velsrc%init(src)
	call this%nvel%init(src)
	call this%nvel_old%init(src)
	
	call this%normals_1%init(src)
	call this%normals_2%init(src)
	call this%normals_3%init(src)
		
	call this%tdatax%init(src)
	call this%tdatay%init(src)
	call this%tdataz%init(src)

	!$omp parallel private(id), num_threads(this%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		this%phi%of(id)%dat => src%of(id)%loc%phi%now
		this%vof%of(id)%dat => src%of(id)%loc%vof%now
		this%curv%of(id)%dat => src%of(id)%loc%normals%curv%now
		this%p%of(id)%dat => src%of(id)%loc%p%now
		
		this%vel%x%of(id)%dat => src%of(id)%loc%vel%x%now
		this%vel%y%of(id)%dat => src%of(id)%loc%vel%y%now
		this%vel%z%of(id)%dat => src%of(id)%loc%vel%z%now

		this%velsrc%x%of(id)%dat => src%of(id)%loc%velsrc%x%now
		this%velsrc%y%of(id)%dat => src%of(id)%loc%velsrc%y%now
		this%velsrc%z%of(id)%dat => src%of(id)%loc%velsrc%z%now

		this%velsrc_old%x%of(id)%dat => src%of(id)%loc%velsrc%x%old
		this%velsrc_old%y%of(id)%dat => src%of(id)%loc%velsrc%y%old
		this%velsrc_old%z%of(id)%dat => src%of(id)%loc%velsrc%z%old
		
		this%nvel%x%of(id)%dat => src%of(id)%loc%nvel%x%now
		this%nvel%y%of(id)%dat => src%of(id)%loc%nvel%y%now
		this%nvel%z%of(id)%dat => src%of(id)%loc%nvel%z%now
		
		this%nvel_old%x%of(id)%dat => src%of(id)%loc%nvel%x%old
		this%nvel_old%y%of(id)%dat => src%of(id)%loc%nvel%y%old
		this%nvel_old%z%of(id)%dat => src%of(id)%loc%nvel%z%old
		
		this%normals_1%x%of(id)%dat => src%of(id)%loc%normals%x%now
		this%normals_1%y%of(id)%dat => src%of(id)%loc%normals%y%now
		this%normals_1%z%of(id)%dat => src%of(id)%loc%normals%z%now
		
		this%normals_2%x%of(id)%dat => src%of(id)%loc%normals%xx%now
		this%normals_2%y%of(id)%dat => src%of(id)%loc%normals%yy%now
		this%normals_2%z%of(id)%dat => src%of(id)%loc%normals%zz%now
		
		this%normals_3%x%of(id)%dat => src%of(id)%loc%normals%xy%now
		this%normals_3%y%of(id)%dat => src%of(id)%loc%normals%xz%now
		this%normals_3%z%of(id)%dat => src%of(id)%loc%normals%yz%now
				
		this%tdatax%x%of(id)%dat => src%of(id)%loc%tdata%x%s1
		this%tdatax%y%of(id)%dat => src%of(id)%loc%tdata%x%s2
		this%tdatax%z%of(id)%dat => src%of(id)%loc%tdata%x%s3
		
		this%tdatay%x%of(id)%dat => src%of(id)%loc%tdata%y%s1
		this%tdatay%y%of(id)%dat => src%of(id)%loc%tdata%y%s2
		this%tdatay%z%of(id)%dat => src%of(id)%loc%tdata%y%s3

		this%tdataz%x%of(id)%dat => src%of(id)%loc%tdata%z%s1
		this%tdataz%y%of(id)%dat => src%of(id)%loc%tdata%z%s2
		this%tdataz%z%of(id)%dat => src%of(id)%loc%tdata%z%s3
				
	!$omp end parallel 

end subroutine

subroutine ptr_family_cputime(p,y) 
implicit none
class(ptr_family) :: p
real(8)  :: y

	y = p%phi%cputime + p%p%cputime + p%curv%cputime + p%vof%cputime
	y = y + p%vel%cputime + p%nvel%cputime + p%velsrc%cputime + p%velsrc_old%cputime
	y = y + p%normals_1%cputime + p%normals_2%cputime + p%normals_3%cputime
	y = y + p%tdatax%cputime + p%tdatay%cputime + p%tdataz%cputime

end subroutine

subroutine ptr_family_reset(p)
implicit none
class(ptr_family) :: p

	call p%phi%reset
	call p%p%reset
	call p%curv%reset
	call p%vof%reset
	call p%vel%reset
	call p%nvel%reset
	call p%velsrc%reset
	call p%velsrc_old%reset
	call p%normals_1%reset
	call p%normals_2%reset
	call p%normals_3%reset
	call p%tdatax%reset
	call p%tdatay%reset
	call p%tdataz%reset
	
end subroutine

end module 

