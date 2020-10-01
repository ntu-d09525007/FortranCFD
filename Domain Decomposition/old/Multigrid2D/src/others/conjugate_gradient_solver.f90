subroutine ConjugateGradient_sync()
use all
implicit none
integer :: id,i,j,k,level
integer :: nx,ny,nz

!$omp parallel private(id,i,j,k,nx,ny,nz), num_threads(p%glb%threads)

	id=0
	!$ id = omp_get_thread_num()
	
	nz = p%of(id)%loc%cg%nz
	ny = p%of(id)%loc%cg%ny
	nx = p%of(id)%loc%cg%nx
	
	!==================================================================
	
	if( p%of(id)%loc%idx>0 )then
		do k = 1, nz
		do j = 1, ny
			p%of(id)%loc%cg%x(0,j,k)=p%of(id-1)%loc%cg%x(p%of(id-1)%loc%cg%nx,j,k)
			p%of(id)%loc%cg%a(0,j,k)=p%of(id-1)%loc%cg%a(p%of(id-1)%loc%cg%nx,j,k)
			p%of(id)%loc%cg%b(0,j,k)=p%of(id-1)%loc%cg%b(p%of(id-1)%loc%cg%nx,j,k)
			p%of(id)%loc%cg%f(0,j,k)=p%of(id-1)%loc%cg%f(p%of(id-1)%loc%cg%nx,j,k)
		enddo
		enddo
	else
		do k = 1, nz
		do j = 1, ny
			p%of(id)%loc%cg%x(0,j,k)=p%of(id)%loc%cg%x(1,j,k)
			p%of(id)%loc%cg%a(0,j,k)=p%of(id)%loc%cg%a(1,j,k)
			p%of(id)%loc%cg%b(0,j,k)=p%of(id)%loc%cg%b(1,j,k)
			p%of(id)%loc%cg%f(0,j,k)=p%of(id)%loc%cg%f(1,j,k)
		enddo
		enddo
	endif
	
	if( p%of(id)%loc%idx<p%glb%grid_x-1 )then
		do k = 1, nz
		do j = 1, ny
			p%of(id)%loc%cg%x(nx+1,j,k)=p%of(id+1)%loc%cg%x(1,j,k)
			p%of(id)%loc%cg%a(nx+1,j,k)=p%of(id+1)%loc%cg%a(1,j,k)
			p%of(id)%loc%cg%b(nx+1,j,k)=p%of(id+1)%loc%cg%b(1,j,k)
			p%of(id)%loc%cg%f(nx+1,j,k)=p%of(id+1)%loc%cg%f(1,j,k)
		enddo
		enddo
	else
		do k = 1, nz
		do j = 1, ny
			p%of(id)%loc%cg%x(nx+1,j,k)=p%of(id)%loc%cg%x(nx,j,k)
			p%of(id)%loc%cg%a(nx+1,j,k)=p%of(id)%loc%cg%a(nx,j,k)
			p%of(id)%loc%cg%b(nx+1,j,k)=p%of(id)%loc%cg%b(nx,j,k)
			p%of(id)%loc%cg%f(nx+1,j,k)=p%of(id)%loc%cg%f(nx,j,k)
		enddo
		enddo
	endif
	
	!==============================================================
	
	if( p%of(id)%loc%idy>0 )then
		do k = 1, nz
		do i = 1, nx
			p%of(id)%loc%cg%x(i,0,k)=p%of(id-p%glb%grid_x)%loc%cg%x(i,p%of(id-p%glb%grid_x)%loc%cg%ny,k)
			p%of(id)%loc%cg%a(i,0,k)=p%of(id-p%glb%grid_x)%loc%cg%a(i,p%of(id-p%glb%grid_x)%loc%cg%ny,k)
			p%of(id)%loc%cg%b(i,0,k)=p%of(id-p%glb%grid_x)%loc%cg%b(i,p%of(id-p%glb%grid_x)%loc%cg%ny,k)
			p%of(id)%loc%cg%f(i,0,k)=p%of(id-p%glb%grid_x)%loc%cg%f(i,p%of(id-p%glb%grid_x)%loc%cg%ny,k)
		enddo
		enddo
	else
		do k = 1, nz
		do i = 1, nx
			p%of(id)%loc%cg%x(i,0,k)=p%of(id)%loc%cg%x(i,1,k)
			p%of(id)%loc%cg%a(i,0,k)=p%of(id)%loc%cg%a(i,1,k)
			p%of(id)%loc%cg%b(i,0,k)=p%of(id)%loc%cg%b(i,1,k)
			p%of(id)%loc%cg%f(i,0,k)=p%of(id)%loc%cg%f(i,1,k)
		enddo
		enddo
	endif
	
	if( p%of(id)%loc%idy<p%glb%grid_y-1 )then
		do k = 1, nz
		do i = 1, nx
			p%of(id)%loc%cg%x(i,ny+1,k)=p%of(id+p%glb%grid_x)%loc%cg%x(i,1,k)
			p%of(id)%loc%cg%a(i,ny+1,k)=p%of(id+p%glb%grid_x)%loc%cg%a(i,1,k)
			p%of(id)%loc%cg%b(i,ny+1,k)=p%of(id+p%glb%grid_x)%loc%cg%b(i,1,k)
			p%of(id)%loc%cg%f(i,ny+1,k)=p%of(id+p%glb%grid_x)%loc%cg%f(i,1,k)
		enddo
		enddo
	else
		do k = 1, nz
		do i = 1, nx
			p%of(id)%loc%cg%x(i,ny+1,k)=p%of(id)%loc%cg%x(i,ny,k)
			p%of(id)%loc%cg%a(i,ny+1,k)=p%of(id)%loc%cg%a(i,ny,k)
			p%of(id)%loc%cg%b(i,ny+1,k)=p%of(id)%loc%cg%b(i,ny,k)
			p%of(id)%loc%cg%f(i,ny+1,k)=p%of(id)%loc%cg%f(i,ny,k)
		enddo
		enddo
	endif	
	
	!==============================================================
	
	if( p%of(id)%loc%idz>0 )then
		do j = 1, ny
		do i = 1, nx
			p%of(id)%loc%cg%x(i,j,0)=&
			&p%of(id-p%glb%grid_x*p%glb%grid_y)%loc%cg%x(i,j,p%of(id-p%glb%grid_x*p%glb%grid_y)%loc%cg%nz)
			p%of(id)%loc%cg%a(i,j,0)=&
			&p%of(id-p%glb%grid_x*p%glb%grid_y)%loc%cg%a(i,j,p%of(id-p%glb%grid_x*p%glb%grid_y)%loc%cg%nz)
			p%of(id)%loc%cg%b(i,j,0)=&
			&p%of(id-p%glb%grid_x*p%glb%grid_y)%loc%cg%b(i,j,p%of(id-p%glb%grid_x*p%glb%grid_y)%loc%cg%nz)
			p%of(id)%loc%cg%f(i,j,0)=&
			&p%of(id-p%glb%grid_x*p%glb%grid_y)%loc%cg%f(i,j,p%of(id-p%glb%grid_x*p%glb%grid_y)%loc%cg%nz)
		enddo
		enddo
	else
		do j = 1, ny
		do i = 1, nx
			p%of(id)%loc%cg%x(i,j,0)=p%of(id)%loc%cg%x(i,j,1)
			p%of(id)%loc%cg%a(i,j,0)=p%of(id)%loc%cg%a(i,j,1)
			p%of(id)%loc%cg%b(i,j,0)=p%of(id)%loc%cg%b(i,j,1)
			p%of(id)%loc%cg%f(i,j,0)=p%of(id)%loc%cg%f(i,j,1)
		enddo
		enddo
	endif
	
	if( p%of(id)%loc%idz<p%glb%grid_z-1 )then
		do j = 1, ny
		do i = 1, nx
			p%of(id)%loc%cg%x(i,j,nz+1)=p%of(id+p%glb%grid_x*p%glb%grid_y)%loc%cg%x(i,j,1)
			p%of(id)%loc%cg%a(i,j,nz+1)=p%of(id+p%glb%grid_x*p%glb%grid_y)%loc%cg%a(i,j,1)
			p%of(id)%loc%cg%b(i,j,nz+1)=p%of(id+p%glb%grid_x*p%glb%grid_y)%loc%cg%b(i,j,1)
			p%of(id)%loc%cg%f(i,j,nz+1)=p%of(id+p%glb%grid_x*p%glb%grid_y)%loc%cg%f(i,j,1)
		enddo
		enddo
	else
		do j = 1, ny
		do i = 1, nx
			p%of(id)%loc%cg%x(i,j,nz+1)=p%of(id)%loc%cg%x(i,j,nz)
			p%of(id)%loc%cg%a(i,j,nz+1)=p%of(id)%loc%cg%a(i,j,nz)
			p%of(id)%loc%cg%b(i,j,nz+1)=p%of(id)%loc%cg%b(i,j,nz)
			p%of(id)%loc%cg%f(i,j,nz+1)=p%of(id)%loc%cg%f(i,j,nz)
		enddo
		enddo
	endif
	
!$omp end parallel

end subroutine