module ptr_tree
use ptr_roots
!$ use omp_lib
implicit none

type ptr_family
integer :: threads
type(pointer_parent) :: phi, p, vof, c
type(pointer_vector_parent) :: vel, nvel, nvel_old, velsrc_old, velsrc
type(pointer_vector_parent) :: normals, vort
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
    call this%c%init(src)
    
    call this%vel%init(src,3)
    call this%velsrc_old%init(src,3)
    call this%velsrc%init(src,3)
    call this%nvel%init(src,3)
    call this%nvel_old%init(src,3)
    
    call this%normals%init(src,10)
    call this%vort%init(src,3)
        
    call this%tdatax%init(src,3)
    call this%tdatay%init(src,3)
    call this%tdataz%init(src,3)

    !$omp parallel do 
    do id = 0, this%threads-1

        this%phi%of(id)%dat => src%of(id)%loc%phi%now
        this%vof%of(id)%dat => src%of(id)%loc%vof%now
        this%p%of(id)%dat => src%of(id)%loc%p%now
        this%c%of(id)%dat => src%of(id)%loc%c%now
        
        this%vel%nodes(1)%of(id)%dat => src%of(id)%loc%vel%x%now
        this%vel%nodes(2)%of(id)%dat => src%of(id)%loc%vel%y%now
        this%vel%nodes(3)%of(id)%dat => src%of(id)%loc%vel%z%now
        
        this%vort%nodes(1)%of(id)%dat => src%of(id)%loc%vort%x%now
        this%vort%nodes(2)%of(id)%dat => src%of(id)%loc%vort%y%now
        this%vort%nodes(3)%of(id)%dat => src%of(id)%loc%vort%z%now

        this%velsrc%nodes(1)%of(id)%dat => src%of(id)%loc%velsrc%x%now
        this%velsrc%nodes(2)%of(id)%dat => src%of(id)%loc%velsrc%y%now
        this%velsrc%nodes(3)%of(id)%dat => src%of(id)%loc%velsrc%z%now

        this%velsrc_old%nodes(1)%of(id)%dat => src%of(id)%loc%velsrc%x%old
        this%velsrc_old%nodes(2)%of(id)%dat => src%of(id)%loc%velsrc%y%old
        this%velsrc_old%nodes(3)%of(id)%dat => src%of(id)%loc%velsrc%z%old
        
        this%nvel%nodes(1)%of(id)%dat => src%of(id)%loc%nvel%x%now
        this%nvel%nodes(2)%of(id)%dat => src%of(id)%loc%nvel%y%now
        this%nvel%nodes(3)%of(id)%dat => src%of(id)%loc%nvel%z%now
        
        this%nvel_old%nodes(1)%of(id)%dat => src%of(id)%loc%nvel%x%old
        this%nvel_old%nodes(2)%of(id)%dat => src%of(id)%loc%nvel%y%old
        this%nvel_old%nodes(3)%of(id)%dat => src%of(id)%loc%nvel%z%old
        
        this%normals%nodes(1)%of(id)%dat => src%of(id)%loc%normals%x%now
        this%normals%nodes(2)%of(id)%dat => src%of(id)%loc%normals%y%now
        this%normals%nodes(3)%of(id)%dat => src%of(id)%loc%normals%z%now
        
        this%normals%nodes(4)%of(id)%dat => src%of(id)%loc%normals%xx%now
        this%normals%nodes(5)%of(id)%dat => src%of(id)%loc%normals%yy%now
        this%normals%nodes(6)%of(id)%dat => src%of(id)%loc%normals%zz%now
        
        this%normals%nodes(7)%of(id)%dat => src%of(id)%loc%normals%xy%now
        this%normals%nodes(8)%of(id)%dat => src%of(id)%loc%normals%xz%now
        this%normals%nodes(9)%of(id)%dat => src%of(id)%loc%normals%yz%now
        
        this%normals%nodes(10)%of(id)%dat => src%of(id)%loc%normals%curv%now
                
        this%tdatax%nodes(1)%of(id)%dat => src%of(id)%loc%tdata%x%s1
        this%tdatax%nodes(2)%of(id)%dat => src%of(id)%loc%tdata%x%s2
        this%tdatax%nodes(3)%of(id)%dat => src%of(id)%loc%tdata%x%s3
        
        this%tdatay%nodes(1)%of(id)%dat => src%of(id)%loc%tdata%y%s1
        this%tdatay%nodes(2)%of(id)%dat => src%of(id)%loc%tdata%y%s2
        this%tdatay%nodes(3)%of(id)%dat => src%of(id)%loc%tdata%y%s3

        this%tdataz%nodes(1)%of(id)%dat => src%of(id)%loc%tdata%z%s1
        this%tdataz%nodes(2)%of(id)%dat => src%of(id)%loc%tdata%z%s2
        this%tdataz%nodes(3)%of(id)%dat => src%of(id)%loc%tdata%z%s3
                
    end do
    !$omp end parallel do

end subroutine


subroutine ptr_family_cputime(p,y) 
implicit none
class(ptr_family) :: p
real(8)  :: y

    y = p%phi%cputime + p%p%cputime + p%vof%cputime
    y = y + p%vel%cputime + p%nvel%cputime + p%velsrc%cputime + p%velsrc_old%cputime + p%vort%cputime
    y = y + p%normals%cputime 
    y = y + p%tdatax%cputime + p%tdatay%cputime + p%tdataz%cputime

end subroutine

subroutine ptr_family_reset(p)
implicit none
class(ptr_family) :: p

    call p%phi%reset
    call p%p%reset
    call p%vof%reset
    call p%vel%reset
    call p%vort%reset
    call p%nvel%reset
    call p%velsrc%reset
    call p%velsrc_old%reset
    call p%normals%reset
    call p%tdatax%reset
    call p%tdatay%reset
    call p%tdataz%reset
    

end subroutine

end module 

