module ptr_roots
use tree
implicit none

type pointer_child
integer :: is, ie, js, je, ghc
real(8), dimension(:,:), pointer :: dat
end type pointer_child

type pointer_parent
integer :: gx, gy
integer(8) :: cpurate
real(8) :: cputime
logical :: xper, yper
type(pointer_child),allocatable :: of(:,:)
contains
procedure init => ptrpart_init
procedure sync => ptrpart_sync
procedure reset => ptrpart_reset
end type pointer_parent

type pointer_vector_parent
type(pointer_parent),allocatable :: nodes(:)
integer :: num
integer(8) :: cpurate
real(8) :: cputime
contains
procedure init => ptrvecpart_init
procedure sync => ptrvecpart_sync
procedure reset => ptrvecpart_reset
end type pointer_vector_parent

type pointer_mg_child
integer :: n
real(8), dimension(:,:), pointer :: dat
end type pointer_mg_child

type pointer_mg_child2
type(pointer_mg_child),allocatable :: at(:)
end type pointer_mg_child2

type pointer_mg_parent
integer :: gx, gy
integer(8) :: cpurate
real(8) :: cputime
logical :: xper, yper, zper
type(pointer_mg_child2),allocatable :: of(:,:)
contains
procedure init => ptrmgpart_init
procedure sync => ptrmgpart_sync
procedure reset => ptrmgpart_reset
end type pointer_mg_parent

contains 

subroutine ptrpart_init(p,src)
implicit none
class(pointer_parent) :: p
type(manager) :: src
integer :: idx,idy,id

    p%cputime = 0.0d0
    p%cpurate = src%glb%cpurate
    
    p%gx = src%glb%grid_x
    p%gy = src%glb%grid_y

    p%xper = src%glb%xper
    p%yper = src%glb%yper
    
    allocate( p%of(0:p%gx-1,0:p%gy-1) )
    
    !$omp parallel do private(id)
    do idy = 0, p%gy-1
    do idx = 0, p%gx-1

        id = src%glb%id(idx,idy)
        
        p%of(idx,idy)%is = src%of(id)%loc%is
        p%of(idx,idy)%ie = src%of(id)%loc%ie 
        p%of(idx,idy)%js = src%of(id)%loc%js
        p%of(idx,idy)%je = src%of(id)%loc%je
        
        p%of(idx,idy)%ghc = src%glb%ghc

    enddo
    enddo
    !$omp end parallel do

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

subroutine ptrmgpart_init(p,src)
implicit none
class(pointer_mg_parent) :: p
type(manager) :: src
integer :: idx,idy,id,level

p%cputime = 0.0d0
p%cpurate = src%glb%cpurate

p%gx = src%glb%grid_x
p%gy = src%glb%grid_y

p%xper = src%glb%xper
p%yper = src%glb%yper

allocate( p%of(0:p%gx-1,0:p%gy-1) )

!$omp parallel do private(id,level)
do idy = 0, p%gy-1
do idx = 0, p%gx-1

    id = src%glb%id(idx,idy)

    allocate( p%of(idx,idy)%at(src%glb%level) )

    do level = 1, src%glb%level
        p%of(idx,idy)%at(level)%n = src%of(id)%loc%mg(level)%nx
    enddo

enddo
enddo
!$omp end parallel do

end subroutine

subroutine ptrpart_sync(p)
implicit none
class(pointer_parent) :: p
integer :: idx,idy,i,j,ghc
integer(8) :: cpustart, cpuend

    call system_clock(cpustart)

    ghc = p%of(0,0)%ghc

    !$omp parallel do private(i,j), collapse(2)
    do idy = 0, p%gy-1
    do idx = 0, p%gx-1
        
        ! x direction
        if(idx<p%gx-1)then
            
            do j = p%of(idx,idy)%js-p%of(idx,idy)%ghc, p%of(idx,idy)%je+p%of(idx,idy)%ghc
            do i = p%of(idx,idy)%ie+1, p%of(idx,idy)%ie+p%of(idx,idy)%ghc
                p%of(idx,idy)%dat(i,j) = p%of(idx+1,idy)%dat(i,j)
            end do
            end do
        
        endif

        if(idx>0)then
            
            do j = p%of(idx,idy)%js-p%of(idx,idy)%ghc, p%of(idx,idy)%je+p%of(idx,idy)%ghc
            do i = p%of(idx,idy)%is-1, p%of(idx,idy)%is-p%of(idx,idy)%ghc, -1
                p%of(idx,idy)%dat(i,j) = p%of(idx-1,idy)%dat(i,j)
            end do
            end do
        
        endif
        
        ! y direction
        if(idy<p%gy-1)then
            
            do i = p%of(idx,idy)%is-p%of(idx,idy)%ghc, p%of(idx,idy)%ie+p%of(idx,idy)%ghc
            do j = p%of(idx,idy)%je+1, p%of(idx,idy)%je+p%of(idx,idy)%ghc
                p%of(idx,idy)%dat(i,j) = p%of(idx,idy+1)%dat(i,j)
            end do
            end do
        
        endif
        
        if(idy>0)then
            
            do i = p%of(idx,idy)%is-p%of(idx,idy)%ghc, p%of(idx,idy)%ie+p%of(idx,idy)%ghc
            do j = p%of(idx,idy)%js-1, p%of(idx,idy)%js-p%of(idx,idy)%ghc, -1
                p%of(idx,idy)%dat(i,j) = p%of(idx,idy-1)%dat(i,j)
            end do
            end do
        
        endif
        
    end do
    end do
    !$omp end parallel do

    if(p%xper)then
        
        !$omp parallel do private(i,j)
        do idy = 0, p%gy-1

            do j = p%of(0,idy)%js - ghc, p%of(0,idy)%je + ghc
            do i = 1, ghc
                p%of(0,idy)%dat(1-i,j) = p%of(p%gx-1,idy)%dat(p%of(p%gx-1,idy)%ie+1-i,j)
                p%of(p%gx-1,idy)%dat(p%of(p%gx-1,idy)%ie+i,j) = p%of(0,idy)%dat(i,j)
            enddo
            enddo

        enddo
        !$omp end parallel do

    endif

    if(p%yper)then

        !$omp parallel do private(i,j)
        do idx = 0, p%gx-1

            do i = p%of(idx,0)%is - ghc, p%of(idx,0)%ie + ghc
            do j = 1, ghc
                p%of(idx,0)%dat(i,1-j) = p%of(idx,p%gy-1)%dat(i,p%of(idx,p%gy-1)%je+1-j)
                p%of(idx,p%gy-1)%dat(i,p%of(idx,p%gy-1)%je+j) = p%of(idx,0)%dat(i,j)
            enddo
            enddo

        enddo
        !$omp end parallel do


    endif
    
    call system_clock(cpuend)
    p%cputime = p%cputime + real(cpuend-cpustart,kind=8)/real(p%cpurate,kind=8)

end subroutine

subroutine ptrvecpart_sync(p)
implicit none
class(pointer_vector_parent) :: p
integer :: id,i,j,k
integer(8) :: cpustart, cpuend
integer :: comp
    
    call system_clock(cpustart)

    do comp = 1, p%num
        call p%nodes(comp)%sync
    enddo
    
    call system_clock(cpuend)
    p%cputime = p%cputime + real(cpuend-cpustart,kind=8)/real(p%cpurate,kind=8)
    
end subroutine

subroutine ptrmgpart_sync(p,level)
implicit none
class(pointer_mg_parent) :: p
integer,intent(in) :: level
integer :: idx,idy,i,j,n
integer(8) :: cpustart, cpuend

call system_clock(cpustart)

n = p%of(0,0)%at(level)%n

!$omp parallel do collapse(2), private(i,j)
do idx = 0, p%gx-1
do idy = 0, p%gy-1

    if( idx>0 )then
        do j = 1, n
            p%of(idx,idy)%at(level)%dat(0,j) = p%of(idx-1,idy)%at(level)%dat(n,j)
        enddo
    else
        do j = 1, n
            p%of(idx,idy)%at(level)%dat(0,j) = p%of(idx,idy)%at(level)%dat(1,j)
        enddo
    endif

    if( idx<p%gx-1 )then
        do j = 1, n
            p%of(idx,idy)%at(level)%dat(n+1,j) = p%of(idx+1,idy)%at(level)%dat(1,j)
        enddo
    else
        do j = 1, n
            p%of(idx,idy)%at(level)%dat(n+1,j) = p%of(idx,idy)%at(level)%dat(n,j)
        enddo
    endif

    !--------------------------

    if( idy>0 )then
        do i = 1, n
            p%of(idx,idy)%at(level)%dat(i,0) = p%of(idx,idy-1)%at(level)%dat(i,n)
        enddo
    else
        do i = 1, n
            p%of(idx,idy)%at(level)%dat(i,0) = p%of(idx,idy)%at(level)%dat(i,1)
        enddo
    endif

    if( idy<p%gy-1 )then
        do i = 1, n
            p%of(idx,idy)%at(level)%dat(i,n+1) = p%of(idx,idy+1)%at(level)%dat(i,1)
        enddo
    else
        do i = 1, n
            p%of(idx,idy)%at(level)%dat(i,n+1) = p%of(idx,idy)%at(level)%dat(i,n)
        enddo
    endif

enddo
enddo
!$omp end parallel do

if(p%xper)then
        
    !$omp parallel do private(i,j)
    do idy = 0, p%gy-1

        do j = 1, n
            p%of(0,idy)%at(level)%dat(0,j) = p%of(p%gx-1,idy)%at(level)%dat(n,j)
            p%of(p%gx-1,idy)%at(level)%dat(n+1,j) = p%of(0,idy)%at(level)%dat(1,j)
        enddo

    enddo
    !$omp end parallel do

endif

if(p%yper)then
        
    !$omp parallel do private(i,j)
    do idx = 0, p%gx-1

        do i = 1, n
            p%of(idx,0)%at(level)%dat(i,0) = p%of(idx,p%gy-1)%at(level)%dat(i,n)
            p%of(idx,p%gy-1)%at(level)%dat(i,n+1) = p%of(idx,0)%at(level)%dat(i,1)
        enddo

    enddo
    !$omp end parallel do

endif

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

subroutine ptrmgpart_reset(p)
implicit none
class(pointer_mg_parent) :: p

p%cputime = 0.0d0

end subroutine

end module ptr_roots



