module ptr_tree
use ptr_roots
implicit none

type ptr_family
type(pointer_parent) :: phi, p, vof
type(pointer_vector_parent) :: vel, nvel, nvel_old, velsrc_old, velsrc
type(pointer_vector_parent) :: normals
type(pointer_vector_parent) :: tdatax, tdatay, tdataz
type(pointer_mg_parent) :: mg
contains 
procedure init => ptr_family_init
procedure cputime => ptr_family_cputime
procedure reset => ptr_family_reset
end type ptr_family

contains

subroutine ptr_family_init(p,src)
use tree
implicit none
class(ptr_family) :: p
type(manager),target :: src
integer :: idx,idy,id,level
    
    call p%phi%init(src)
    call p%vof%init(src)
    call p%p%init(src)
    
    call p%vel%init(src,2)
    call p%velsrc_old%init(src,2)
    call p%velsrc%init(src,2)
    call p%nvel%init(src,2)
    call p%nvel_old%init(src,2)
    
    call p%normals%init(src,6)
        
    call p%tdatax%init(src,3)
    call p%tdatay%init(src,3)
    call p%tdataz%init(src,3)

    call p%mg%init(src)

    !$omp parallel do private(id,level), collapse(2)
    do idy = 0, src%glb%grid_y-1
    do idx = 0, src%glb%grid_x-1

        id = src%glb%id(idx,idy)

        do level = 1, src%glb%level
            p%mg%of(idx,idy)%at(level)%dat => src%of(id)%loc%mg(level)%sol
        enddo

        p%phi%of(idx,idy)%dat => src%of(id)%loc%phi%now
        p%vof%of(idx,idy)%dat => src%of(id)%loc%vof%now
        p%p%of(idx,idy)%dat => src%of(id)%loc%p%now

        p%vel%nodes(1)%of(idx,idy)%dat => src%of(id)%loc%vel%x%now
        p%vel%nodes(2)%of(idx,idy)%dat => src%of(id)%loc%vel%y%now

        p%velsrc%nodes(1)%of(idx,idy)%dat => src%of(id)%loc%velsrc%x%now
        p%velsrc%nodes(2)%of(idx,idy)%dat => src%of(id)%loc%velsrc%y%now

        p%velsrc_old%nodes(1)%of(idx,idy)%dat => src%of(id)%loc%velsrc%x%old
        p%velsrc_old%nodes(2)%of(idx,idy)%dat => src%of(id)%loc%velsrc%y%old
        
        p%nvel%nodes(1)%of(idx,idy)%dat => src%of(id)%loc%nvel%x%now
        p%nvel%nodes(2)%of(idx,idy)%dat => src%of(id)%loc%nvel%y%now
        
        p%nvel_old%nodes(1)%of(idx,idy)%dat => src%of(id)%loc%nvel%x%old
        p%nvel_old%nodes(2)%of(idx,idy)%dat => src%of(id)%loc%nvel%y%old
        
        p%normals%nodes(1)%of(idx,idy)%dat => src%of(id)%loc%normals%x%now
        p%normals%nodes(2)%of(idx,idy)%dat => src%of(id)%loc%normals%y%now
        p%normals%nodes(3)%of(idx,idy)%dat => src%of(id)%loc%normals%xx%now
        p%normals%nodes(4)%of(idx,idy)%dat => src%of(id)%loc%normals%yy%now
        p%normals%nodes(5)%of(idx,idy)%dat => src%of(id)%loc%normals%xy%now
        p%normals%nodes(6)%of(idx,idy)%dat => src%of(id)%loc%normals%curv%now
                
        p%tdatax%nodes(1)%of(idx,idy)%dat => src%of(id)%loc%tdata%x%s1
        p%tdatax%nodes(2)%of(idx,idy)%dat => src%of(id)%loc%tdata%x%s2
        p%tdatax%nodes(3)%of(idx,idy)%dat => src%of(id)%loc%tdata%x%s3
        
        p%tdatay%nodes(1)%of(idx,idy)%dat => src%of(id)%loc%tdata%y%s1
        p%tdatay%nodes(2)%of(idx,idy)%dat => src%of(id)%loc%tdata%y%s2
        p%tdatay%nodes(3)%of(idx,idy)%dat => src%of(id)%loc%tdata%y%s3

        p%tdataz%nodes(1)%of(idx,idy)%dat => src%of(id)%loc%tdata%z%s1
        p%tdataz%nodes(2)%of(idx,idy)%dat => src%of(id)%loc%tdata%z%s2
        p%tdataz%nodes(3)%of(idx,idy)%dat => src%of(id)%loc%tdata%z%s3
        
    end do
    end do
    !$omp end parallel do

end subroutine


subroutine ptr_family_cputime(p,y) 
implicit none
class(ptr_family) :: p
real(8)  :: y

    y = p%phi%cputime + p%p%cputime + p%vof%cputime
    y = y + p%vel%cputime + p%nvel%cputime + p%velsrc%cputime + p%velsrc_old%cputime
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
    call p%nvel%reset
    call p%velsrc%reset
    call p%velsrc_old%reset
    call p%normals%reset
    call p%tdatax%reset
    call p%tdatay%reset
    call p%tdataz%reset
    

end subroutine

end module 

