subroutine manager_mgsync_setup(p)
implicit none
class(manager) :: p
integer :: level

allocate( p%mg%syncer(p%glb%level) )

if( p%glb%split == 1)then

endif

do level = 1, p%glb%level
    allocate( p%mg%syncer(level)%sendp()
enddo

end subroutine

subroutine manager_mg_setup(p)
implicit none
class(manager) :: p
integer :: i,j,k,id,level,ind,ii,jj,kk
integer :: nx,ny,nz,n,info
integer :: imax,imin,jmax,jmin,kmax,kmin
real(8) :: dx,dy,dz
integer :: ix,iy,iz

    level = p%glb%level

    if( p%of(0)%loc%mg(level)%n .ne. 1)stop "Coarest grid in MG solver should be 1."
    
    p%mg%n  = p%glb%num_domains*p%of(0)%loc%mg(level)%n
    p%mg%nx = p%glb%grid_x*p%of(0)%loc%mg(level)%nx
    p%mg%ny = p%glb%grid_y*p%of(0)%loc%mg(level)%ny
    p%mg%nz = p%glb%grid_z*p%of(0)%loc%mg(level)%nz
    
    write(*,'(A,4I6)')"Local A size",p%of(0)%loc%mg(level)%n,p%of(0)%loc%mg(level)%nx,p%of(0)%loc%mg(level)%ny,p%of(0)%loc%mg(level)%nz
    write(*,'(A,4I6)')"Global A size",p%mg%n,p%mg%nx,p%mg%ny,p%mg%nz
    
    allocate(p%mg%A(p%mg%n,p%mg%n),p%mg%AA(p%mg%n,p%mg%n),p%mg%B(p%mg%n),p%mg%BB(p%mg%n),p%mg%ipiv(p%mg%n)) 
    allocate(p%mg%i(p%mg%n),p%mg%j(p%mg%n),p%mg%k(p%mg%n),p%mg%node(p%mg%nx,p%mg%ny,p%mg%nz))
    allocate(p%mg%error(p%mg%n),p%mg%sol(p%mg%n))

    allocate(p%mg%ix(0:p%glb%num_domains-1), p%mg%iy(0:p%glb%num_domains-1), p%mg%iz(0:p%glb%num_domains-1) )
    
    dx = p%of(0)%loc%mg(level)%dx
    dy = p%of(0)%loc%mg(level)%dy
    dz = p%of(0)%loc%mg(level)%dz
    
    do i = 1, p%mg%n
    do j = 1, p%mg%n
        p%mg%A(i,j)=0.0d0
    enddo
    enddo

    do id = 0, p%glb%num_domains-1

        if( p%glb%split == 1 )then
            ix = id / p%glb%grid_A
            iy = ( id - ix*p%glb%grid_A ) / p%glb%grid_z
            iz = id - ix*p%glb%grid_A - iy*p%glb%grid_z
        else if ( p%glb%split == 2)then
            iy = id / p%glb%grid_A
            iz = ( id - iy*p%glb%grid_A ) / p%glb%grid_x
            ix = id - ix*p%glb%grid_A - iz*p%glb%grid_x
        else
            iz = id / p%glb%grid_A
            ix = ( id - iz*p%glb%grid_A ) / p%glb%grid_y
            iy = id - ix*p%glb%grid_A - ix*p%glb%grid_y
        endif  

        p%mg%ix(id) = ix
        p%mg%iy(id) = iy
        p%mg%iz(id) = iz

    enddo
    
    imax=0; imin=10000
    jmax=0; jmin=10000
    kmax=0; kmin=10000
    do i = 1, p%mg%n    
        id = (i-1) / p%of(0)%loc%mg(level)%n
        ! global index of local A
        ind = i - id*p%of(0)%loc%mg(level)%n   
        
        ! global (i,j,k) of global A
        ii = p%of(0)%loc%mg(level)%i(ind) + p%of(0)%loc%mg(level)%nx * p%mg%ix(id) 
        jj = p%of(0)%loc%mg(level)%j(ind) + p%of(0)%loc%mg(level)%ny * p%mg%iy(id) 
        kk = p%of(0)%loc%mg(level)%k(ind) + p%of(0)%loc%mg(level)%nz * p%mg%iz(id) 
        
        p%mg%i(i) = ii
        p%mg%j(i) = jj
        p%mg%k(i) = kk
        p%mg%node(ii,jj,kk) = i
        
        p%mg%A(i,i) = -2.0d0/dx**2.0d0-2.0d0/dy**2.0d0-2.0d0/dz**2.0d0
        
        imax=max(imax,ii); imin=min(imin,ii)
        jmax=max(jmax,jj); jmin=min(jmin,jj)
        kmax=max(kmax,kk); kmin=min(kmin,kk)
    end do
    
    write(*,'(A15,2I10)')"I matrix",imin,imax
    write(*,'(A15,2I10)')"J matrix",jmin,jmax
    write(*,'(A15,2I10)')"K matrix",kmin,kmax
    
    do i = 1, p%mg%n
    
        if(p%mg%i(i)>1)then
            p%mg%A(i, p%mg%node(p%mg%i(i)-1,p%mg%j(i),p%mg%k(i)) ) = 1.0d0/dx**2.0d0
        else
            p%mg%A(i,i) = p%mg%A(i,i) + 1.0d0/dx**2.0d0
        endif
        
        if(p%mg%i(i)<p%mg%nx)then
            p%mg%A(i, p%mg%node(p%mg%i(i)+1,p%mg%j(i),p%mg%k(i)) ) = 1.0d0/dx**2.0d0
        else
            p%mg%A(i,i) = p%mg%A(i,i) + 1.0d0/dx**2.0d0
        endif
        
        if(p%mg%j(i)>1)then
            p%mg%A(i, p%mg%node(p%mg%i(i),p%mg%j(i)-1,p%mg%k(i)) ) = 1.0d0/dy**2.0d0
        else
            p%mg%A(i,i) = p%mg%A(i,i) + 1.0d0/dy**2.0d0
        endif
        
        if(p%mg%j(i)<p%mg%ny)then
            p%mg%A(i, p%mg%node(p%mg%i(i),p%mg%j(i)+1,p%mg%k(i)) ) = 1.0d0/dy**2.0d0
        else
            p%mg%A(i,i) = p%mg%A(i,i) + 1.0d0/dy**2.0d0
        endif   

        if(p%mg%k(i)>1)then
            p%mg%A(i, p%mg%node(p%mg%i(i),p%mg%j(i),p%mg%k(i)-1) ) = 1.0d0/dz**2.0d0
        else
            p%mg%A(i,i) = p%mg%A(i,i) + 1.0d0/dz**2.0d0
        endif
        
        if(p%mg%k(i)<p%mg%nz)then
            p%mg%A(i, p%mg%node(p%mg%i(i),p%mg%j(i),p%mg%k(i)+1) ) = 1.0d0/dz**2.0d0
        else
            p%mg%A(i,i) = p%mg%A(i,i) + 1.0d0/dz**2.0d0
        endif   
        
    end do
    
    end subroutine

subroutine manager_mg_solve_exact(p,show)
! Only for n=1
implicit none
include 'mpif.h'
class(manager) :: p
integer :: rank, ierr, istat
integer :: i,j,id,level,ind,info
real(8) :: buffer
logical :: show

level = p%glb%level

if( p%glb%mpisize > 1)then

    if( p%glb%mpirank == 0 )then
        do rank = 1, p%glb%mpisize-1
            do id = 0, p%glb%mpithreads(rank)-1
                call mpi_recv(buffer,1,mpi_real8,rank,1000*rank+id, mpi_comm_world, istat, ierr)
                ind = sum(p%glb%mpithreads(0:rank-1))+id
                p%mg%B(ind+1) = buffer
            enddo
        enddo
    else
        do id = 0, p%glb%threads-1
            call mpi_send(p%of(id)%loc%mg(level)%src(1,1,1),1,MPI_REAL8,0, 1000*p%glb%mpirank+id, mpi_comm_world, ierr)
        enddo
    endif

endif

if( p%glb%mpirank == 0)then

    do i = 1, p%glb%threads
        p%mg%B(i) = p%of(i-1)%loc%mg(level)%src(1,1,1)
    enddo

    do i = 1, p%mg%n
        p%mg%BB(i) = p%mg%B(i)
        do j = 1, p%mg%n
            p%mg%AA(i,j) = p%mg%A(i,j)
        end do
    end do

    call dgesv(p%mg%n,1,p%mg%AA,p%mg%n,p%mg%ipiv,p%mg%B,p%mg%n,info)
    if(info.ne.0)stop "Something Wrong with Multigrid Solver!!"

    do i = 1, p%mg%n
        p%mg%Sol(i) = p%mg%B(i)
    enddo

    call p%mg_find_error(show)

    do id = 0, p%glb%threads-1
        p%of(id)%loc%mg(level)%sol(1,1,1) = p%mg%sol(id+1)
    enddo

endif

if( p%glb%mpisize > 1)then

    if( p%glb%mpirank == 0 )then
        do rank = 1, p%glb%mpisize-1
            do id = 0, p%glb%mpithreads(rank)-1
                ind = sum(p%glb%mpithreads(0:rank-1))+id
                call mpi_send(p%mg%sol(1+ind),1,MPI_REAL8,rank, 2000*rank+id, mpi_comm_world, ierr)
            enddo
        enddo
    else
        do id = 0, p%glb%threads-1
            call mpi_recv(buffer,1,mpi_real8,0,2000*p%glb%mpirank+id, mpi_comm_world, istat, ierr)   
            p%of(id)%loc%mg(level)%sol(1,1,1) = buffer
        enddo
    endif

    call mpi_barrier(mpi_comm_world,ierr)

endif

end subroutine

subroutine manager_mg_find_error(p,show)
class(manager) :: p
integer :: i,j,k,id,level
integer :: ii,jj,kk,ind,info
real(8) :: error, tmp
logical :: show

level = p%glb%level

error=0.0d0
do i = 1, p%mg%n
    tmp= p%mg%BB(i)
    do j = 1, p%mg%n
        tmp = tmp - p%mg%A(i,j)*p%mg%sol(j)
    enddo
    p%mg%error(i) = tmp
    error = max(error,abs(tmp))
enddo

if(show)write(*,*)"Exact solver (LAPACK) error:",error
if(error>1.0d-10)call p%mg_solve_correct

end subroutine

subroutine manager_mg_solve_correct(p)
implicit none
class(manager) :: p
integer :: i,j,k,id,level
integer :: ii,jj,kk,ind,info

level = p%glb%level

do i = 1, p%mg%n
do j = 1, p%mg%n
    p%mg%AA(i,j) = p%mg%A(i,j)
end do
end do

call dgesv(p%mg%n,1,p%mg%AA,p%mg%n,p%mg%ipiv,p%mg%error,p%mg%n,info)
if(info.ne.0)stop "Something Wrong with Multigrid Solver!!"

do i = 1, p%mg%n
    p%mg%sol(i) = p%mg%sol(i) + p%mg%error(i) 
enddo

end subroutine