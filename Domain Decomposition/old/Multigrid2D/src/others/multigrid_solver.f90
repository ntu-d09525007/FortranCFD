subroutine multigrid_relax(level,num)
use all
implicit none
integer :: relax_iter, num, id, level
integer :: i,j,k
real(8) :: dx,dy,dz,a

do relax_iter = 1, num

	call multigrid_sync(level) 
		
	!$omp parallel private(id,i,j,k,dx,dy,dz,a), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		dx = p%of(id)%loc%mg(level)%dx
		dy = p%of(id)%loc%mg(level)%dy
		dz = p%of(id)%loc%mg(level)%dz

		!a = -2.0d0/dx**2.0d0 -2.0d0/dy**2.0d0 -2.0d0/dz**2.0d0
		a = -2.0d0/dx**2.0d0 -2.0d0/dy**2.0d0
		
		do k = 1, p%of(id)%loc%mg(level)%nz
		do j = 1, p%of(id)%loc%mg(level)%ny
		do i = 1, p%of(id)%loc%mg(level)%nx
	
			! 3D
			! p%of(id)%loc%mg(level)%sol(i,j,k) = ( p%of(id)%loc%mg(level)%src(i,j,k) &
			! & - (p%of(id)%loc%mg(level)%sol(i+1,j,k)+p%of(id)%loc%mg(level)%sol(i-1,j,k))/dx**2.0d0 &
			! & - (p%of(id)%loc%mg(level)%sol(i,j+1,k)+p%of(id)%loc%mg(level)%sol(i,j-1,k))/dy**2.0d0 &
			! & - (p%of(id)%loc%mg(level)%sol(i,j,k+1)+p%of(id)%loc%mg(level)%sol(i,j,k-1))/dz**2.0d0 ) / a
				
			! 2D
			p%of(id)%loc%mg(level)%sol(i,j,k) = ( p%of(id)%loc%mg(level)%src(i,j,k) &
			& - (p%of(id)%loc%mg(level)%sol(i+1,j,k)+p%of(id)%loc%mg(level)%sol(i-1,j,k))/dx**2.0d0 &
			& - (p%of(id)%loc%mg(level)%sol(i,j+1,k)+p%of(id)%loc%mg(level)%sol(i,j-1,k))/dy**2.0d0 )/a
		enddo
		enddo
		enddo

	!$omp end parallel
			
	call multigrid_compatibility_sol(level)
			
enddo

call multigrid_sync(level) 
		
end subroutine

subroutine multigrid_restriction(level)
use all
implicit none
integer :: id,level
integer :: i,j,k

call multigrid_residual(level,.false.)

!$omp parallel private(id,i,j,k), num_threads(p%glb%threads)

	id=0
	!$ id = omp_get_thread_num()
	
	do k = 1, p%of(id)%loc%mg(level+1)%nz
	do j = 1, p%of(id)%loc%mg(level+1)%ny
	do i = 1, p%of(id)%loc%mg(level+1)%nx
	
		! 3D
		! p%of(id)%loc%mg(level+1)%src(i,j,k) = &
		! &( p%of(id)%loc%mg(level)%res(2*i-1,2*j-1,2*k-1)+p%of(id)%loc%mg(level)%res(2*i,2*j-1,2*k-1)+ &
		! &  p%of(id)%loc%mg(level)%res(2*i-1,2*j  ,2*k-1)+p%of(id)%loc%mg(level)%res(2*i,2*j  ,2*k-1)+ &
		! &  p%of(id)%loc%mg(level)%res(2*i-1,2*j-1,2*k)  +p%of(id)%loc%mg(level)%res(2*i,2*j-1,2*k  )+ &
		! &  p%of(id)%loc%mg(level)%res(2*i-1,2*j  ,2*k)  +p%of(id)%loc%mg(level)%res(2*i,2*j  ,2*k  )) / 8.0d0
		
		! 2D
		p%of(id)%loc%mg(level+1)%src(i,j,k) = &
		&( p%of(id)%loc%mg(level)%res(2*i-1,2*j-1,2*k-1)+p%of(id)%loc%mg(level)%res(2*i,2*j-1,2*k-1)+ &
		&  p%of(id)%loc%mg(level)%res(2*i-1,2*j  ,2*k-1)+p%of(id)%loc%mg(level)%res(2*i,2*j  ,2*k-1)) / 4.0d0
		 
		p%of(id)%loc%mg(level+1)%sol(i,j,k) = 0.0d0
		
	enddo
	enddo
	enddo
	
!$omp end parallel	
		
call multigrid_compatibility_src(level+1)

end subroutine

subroutine multigrid_final_exact()
use all
implicit none
integer :: id,level,relax_iter
integer :: i,j,k
real(8) :: dx,dy,dz,a,w,pnew,rms0

level=p%glb%level
   
call p%mg_solve_exact   
call multigrid_compatibility_sol(level)
call multigrid_residual(level,.false.)
call multigrid_sync(level)
write(*,*)"Exact solver,",relax_iter,p%of(0)%loc%mg(level)%l2norm

end subroutine

subroutine multigrid_iterative()
use all
implicit none
integer :: id,level,relax_iter
integer :: i,j,k
real(8) :: dx,dy,dz,a,w,pnew

level=p%glb%level
	
w=1.0d0

do relax_iter = 1, 100000
		
    if(relax_iter>2)then
        w=0.5d0
    endif

	call multigrid_sync(level) 
			
	!$omp parallel private(id,i,j,k,dx,dy,dz,a,pnew), num_threads(p%glb%threads)
		
		id=0
		!$ id = omp_get_thread_num()
			
		dx = p%of(id)%loc%mg(level)%dx
		dy = p%of(id)%loc%mg(level)%dy
		dz = p%of(id)%loc%mg(level)%dz

		!a = -2.0d0/dx**2.0d0 -2.0d0/dy**2.0d0 -2.0d0/dz**2.0d0
		a = -2.0d0/dx**2.0d0 -2.0d0/dy**2.0d0
			
		do k = 1, p%of(id)%loc%mg(level)%nz
		do j = 1, p%of(id)%loc%mg(level)%ny
		do i = 1, p%of(id)%loc%mg(level)%nx
		
			! 3D
			! pnew = ( p%of(id)%loc%mg(level)%src(i,j,k) &
			! & - (p%of(id)%loc%mg(level)%sol(i+1,j,k)+p%of(id)%loc%mg(level)%sol(i-1,j,k))/dx**2.0d0 &
			! & - (p%of(id)%loc%mg(level)%sol(i,j+1,k)+p%of(id)%loc%mg(level)%sol(i,j-1,k))/dy**2.0d0 &
			! & - (p%of(id)%loc%mg(level)%sol(i,j,k+1)+p%of(id)%loc%mg(level)%sol(i,j,k-1))/dz**2.0d0 ) / a
		
			! 2D
			pnew = ( p%of(id)%loc%mg(level)%src(i,j,k) &
			& - (p%of(id)%loc%mg(level)%sol(i+1,j,k)+p%of(id)%loc%mg(level)%sol(i-1,j,k))/dx**2.0d0 &
			& - (p%of(id)%loc%mg(level)%sol(i,j+1,k)+p%of(id)%loc%mg(level)%sol(i,j-1,k))/dy**2.0d0 )/a
			
			p%of(id)%loc%mg(level)%sol(i,j,k) = w*p%of(id)%loc%mg(level)%sol(i,j,k) + (1.0d0-w)*pnew
			
		enddo
		enddo
		enddo
			
	!$omp end parallel
				
	call multigrid_compatibility_sol(level)
	call multigrid_residual(level,.false.)
			
	if( p%of(0)%loc%mg(level)%l2norm .lt. p%glb%p_tol*0.001d0 )exit
	if(mod(relax_iter,5000).eq.0)write(*,*)"MG error",relax_iter,p%of(0)%loc%mg(level)%l2norm
	
enddo

!write(*,*)"Finding exact solution, error=",p%of(0)%loc%mg(level)%l2norm
	
end subroutine

subroutine multigrid_prolongation(level)
use all 
implicit none
integer :: id,level
integer :: i,j,k
real(8) :: mx,my,mz

call multigrid_sync(level+1)
		
!$omp parallel private(id,i,j,k,mx,my,mz), num_threads(p%glb%threads)

	id=0
	!$ id = omp_get_thread_num()

	do k = 1, p%of(id)%loc%mg(level+1)%nz
	do j = 1, p%of(id)%loc%mg(level+1)%ny
	do i = 1, p%of(id)%loc%mg(level+1)%nx
	
		! 3D
		! mx = 0.5d0*( p%of(id)%loc%mg(level+1)%sol(i+1,j,k)-p%of(id)%loc%mg(level+1)%sol(i-1,j,k) )
		! my = 0.5d0*( p%of(id)%loc%mg(level+1)%sol(i,j+1,k)-p%of(id)%loc%mg(level+1)%sol(i,j-1,k) )
		! mz = 0.5d0*( p%of(id)%loc%mg(level+1)%sol(i,j,k+1)-p%of(id)%loc%mg(level+1)%sol(i,j,k-1) )
		! p%of(id)%loc%mg(level)%pol(2*i-1,2*j-1,2*k-1) = p%of(id)%loc%mg(level+1)%sol(i,j,k) - 0.25d0*mx - 0.25d0*my - 0.25d0*mz
		! p%of(id)%loc%mg(level)%pol(2*i-1,2*j,2*k-1)   = p%of(id)%loc%mg(level+1)%sol(i,j,k) - 0.25d0*mx + 0.25d0*my - 0.25d0*mz
		! p%of(id)%loc%mg(level)%pol(2*i,2*j-1,2*k-1)   = p%of(id)%loc%mg(level+1)%sol(i,j,k) + 0.25d0*mx - 0.25d0*my - 0.25d0*mz
		! p%of(id)%loc%mg(level)%pol(2*i,2*j,2*k-1)     = p%of(id)%loc%mg(level+1)%sol(i,j,k) + 0.25d0*mx + 0.25d0*my - 0.25d0*mz
		! p%of(id)%loc%mg(level)%pol(2*i-1,2*j-1,2*k)   = p%of(id)%loc%mg(level+1)%sol(i,j,k) - 0.25d0*mx - 0.25d0*my + 0.25d0*mz
		! p%of(id)%loc%mg(level)%pol(2*i-1,2*j,2*k)     = p%of(id)%loc%mg(level+1)%sol(i,j,k) - 0.25d0*mx + 0.25d0*my + 0.25d0*mz
		! p%of(id)%loc%mg(level)%pol(2*i,2*j-1,2*k)     = p%of(id)%loc%mg(level+1)%sol(i,j,k) + 0.25d0*mx - 0.25d0*my + 0.25d0*mz
		! p%of(id)%loc%mg(level)%pol(2*i,2*j,2*k)       = p%of(id)%loc%mg(level+1)%sol(i,j,k) + 0.25d0*mx + 0.25d0*my + 0.25d0*mz

		! 2D
		mx = 0.5d0*( p%of(id)%loc%mg(level+1)%sol(i+1,j,k)-p%of(id)%loc%mg(level+1)%sol(i-1,j,k) )
		my = 0.5d0*( p%of(id)%loc%mg(level+1)%sol(i,j+1,k)-p%of(id)%loc%mg(level+1)%sol(i,j-1,k) )		
		p%of(id)%loc%mg(level)%pol(2*i-1,2*j-1,2*k-1) = p%of(id)%loc%mg(level+1)%sol(i,j,k) - 0.25d0*mx - 0.25d0*my 
		p%of(id)%loc%mg(level)%pol(2*i-1,2*j  ,2*k-1) = p%of(id)%loc%mg(level+1)%sol(i,j,k) - 0.25d0*mx + 0.25d0*my 
		p%of(id)%loc%mg(level)%pol(2*i  ,2*j-1,2*k-1) = p%of(id)%loc%mg(level+1)%sol(i,j,k) + 0.25d0*mx - 0.25d0*my 
		p%of(id)%loc%mg(level)%pol(2*i  ,2*j  ,2*k-1) = p%of(id)%loc%mg(level+1)%sol(i,j,k) + 0.25d0*mx + 0.25d0*my 
		
	enddo
	enddo
	enddo

	do k = 1, p%of(id)%loc%mg(level)%nz
	do j = 1, p%of(id)%loc%mg(level)%ny
	do i = 1, p%of(id)%loc%mg(level)%nx	
		p%of(id)%loc%mg(level)%sol(i,j,k) = p%of(id)%loc%mg(level)%sol(i,j,k) + p%of(id)%loc%mg(level)%pol(i,j,k)
	enddo
	enddo
	enddo
	
!$omp end parallel	
		

end subroutine

subroutine multigrid_residual(level,reset)
use all
implicit none
integer :: level,id,i,j,k,nx,ny,nz
real(8) :: dx,dy,dz,l2norm
logical :: reset

call multigrid_sync(level)

nx = 0
ny = 0
nz = 0
l2norm = 0

!$omp parallel private(id,i,j,k,dx,dy,dz), num_threads(p%glb%threads), reduction(+:nx,ny,nz,l2norm), shared(level)
	
	id=0
	!$ id = omp_get_thread_num()
	
	dx = p%of(id)%loc%mg(level)%dx
	dy = p%of(id)%loc%mg(level)%dy
	dz = p%of(id)%loc%mg(level)%dz
		
	nx = nx + p%of(id)%loc%mg(level)%nx
	ny = ny + p%of(id)%loc%mg(level)%ny
	nz = nz + p%of(id)%loc%mg(level)%nz
	
	L2norm = 0.0d0
	do k = 1, p%of(id)%loc%mg(level)%nz
	do j = 1, p%of(id)%loc%mg(level)%ny
	do i = 1, p%of(id)%loc%mg(level)%nx
		! 3D
		! p%of(id)%loc%mg(level)%res(i,j,k) = p%of(id)%loc%mg(level)%src(i,j,k) &
		! &   - (p%of(id)%loc%mg(level)%sol(i+1,j,k)-2.0d0*p%of(id)%loc%mg(level)%sol(i,j,k)+p%of(id)%loc%mg(level)%sol(i-1,j,k))/dx**2.0d0 &
		! &   - (p%of(id)%loc%mg(level)%sol(i,j+1,k)-2.0d0*p%of(id)%loc%mg(level)%sol(i,j,k)+p%of(id)%loc%mg(level)%sol(i,j-1,k))/dy**2.0d0 &
		! &   - (p%of(id)%loc%mg(level)%sol(i,j,k+1)-2.0d0*p%of(id)%loc%mg(level)%sol(i,j,k)+p%of(id)%loc%mg(level)%sol(i,j,k-1))/dz**2.0d0

		! 2D
		p%of(id)%loc%mg(level)%res(i,j,k) = p%of(id)%loc%mg(level)%src(i,j,k) &
		&   - (p%of(id)%loc%mg(level)%sol(i+1,j,k)-2.0d0*p%of(id)%loc%mg(level)%sol(i,j,k)+p%of(id)%loc%mg(level)%sol(i-1,j,k))/dx**2.0d0 &
		&   - (p%of(id)%loc%mg(level)%sol(i,j+1,k)-2.0d0*p%of(id)%loc%mg(level)%sol(i,j,k)+p%of(id)%loc%mg(level)%sol(i,j-1,k))/dy**2.0d0
									
		L2norm = L2norm + p%of(id)%loc%mg(level)%res(i,j,k)**2.0d0
	enddo
	enddo
	enddo

!$omp end parallel

!l2norm = dsqrt( l2norm / (nx*ny*nz) )
l2norm = dsqrt( l2norm / real(nx*ny,8) )

if(reset)then
	!$omp parallel do num_threads(p%glb%threads)
	do id = 0, p%glb%threads-1
		p%of(id)%loc%mg(level)%l2norm0 = l2norm
	enddo
	!$omp end parallel do
else
	!$omp parallel do num_threads(p%glb%threads)
	do id = 0, p%glb%threads-1
		p%of(id)%loc%mg(level)%l2norm = l2norm
	enddo
	!$omp end parallel do 
endif

end subroutine

subroutine multigrid_sync(level)
use all
implicit none
integer :: id,i,j,k,level
integer :: nx,ny,nz

!$omp parallel private(id,i,j,k,nx,ny,nz), num_threads(p%glb%threads)

	id=0
	!$ id = omp_get_thread_num()
	
	nz = p%of(id)%loc%mg(level)%nz
	ny = p%of(id)%loc%mg(level)%ny
	nx = p%of(id)%loc%mg(level)%nx
	
	!==================================================================
	
	if( p%of(id)%loc%idx>0 )then
		do k = 1, nz
		do j = 1, ny
			p%of(id)%loc%mg(level)%sol(0,j,k)=p%of(id-1)%loc%mg(level)%sol(p%of(id-1)%loc%mg(level)%nx,j,k)
		enddo
		enddo
	else
		do k = 1, nz
		do j = 1, ny
			p%of(id)%loc%mg(level)%sol(0,j,k)=p%of(id)%loc%mg(level)%sol(1,j,k)
		enddo
		enddo
	endif
	
	if( p%of(id)%loc%idx<p%glb%grid_x-1 )then
		do k = 1, nz
		do j = 1, ny
			p%of(id)%loc%mg(level)%sol(nx+1,j,k)=p%of(id+1)%loc%mg(level)%sol(1,j,k)
		enddo
		enddo
	else
		do k = 1, nz
		do j = 1, ny
			p%of(id)%loc%mg(level)%sol(nx+1,j,k)=p%of(id)%loc%mg(level)%sol(nx,j,k)
		enddo
		enddo
	endif
	
	!==============================================================
	
	if( p%of(id)%loc%idy>0 )then
		do k = 1, nz
		do i = 1, nx
			p%of(id)%loc%mg(level)%sol(i,0,k)=p%of(id-p%glb%grid_x)%loc%mg(level)%sol(i,p%of(id-p%glb%grid_x)%loc%mg(level)%ny,k)
		enddo
		enddo
	else
		do k = 1, nz
		do i = 1, nx
			p%of(id)%loc%mg(level)%sol(i,0,k)=p%of(id)%loc%mg(level)%sol(i,1,k)
		enddo
		enddo
	endif
	
	if( p%of(id)%loc%idy<p%glb%grid_y-1 )then
		do k = 1, nz
		do i = 1, nx
			p%of(id)%loc%mg(level)%sol(i,ny+1,k)=p%of(id+p%glb%grid_x)%loc%mg(level)%sol(i,1,k)
		enddo
		enddo
	else
		do k = 1, nz
		do i = 1, nx
			p%of(id)%loc%mg(level)%sol(i,ny+1,k)=p%of(id)%loc%mg(level)%sol(i,ny,k)
		enddo
		enddo
	endif	
	
	!==============================================================
	
	if( p%of(id)%loc%idz>0 )then
		do j = 1, ny
		do i = 1, nx
			p%of(id)%loc%mg(level)%sol(i,j,0)=&
			&p%of(id-p%glb%grid_x*p%glb%grid_y)%loc%mg(level)%sol(i,j,p%of(id-p%glb%grid_x*p%glb%grid_y)%loc%mg(level)%nz)
		enddo
		enddo
	else
		do j = 1, ny
		do i = 1, nx
			p%of(id)%loc%mg(level)%sol(i,j,0)=p%of(id)%loc%mg(level)%sol(i,j,1)
		enddo
		enddo
	endif
	
	if( p%of(id)%loc%idz<p%glb%grid_z-1 )then
		do j = 1, ny
		do i = 1, nx
			p%of(id)%loc%mg(level)%sol(i,j,nz+1)=p%of(id+p%glb%grid_x*p%glb%grid_y)%loc%mg(level)%sol(i,j,1)
		enddo
		enddo
	else
		do j = 1, ny
		do i = 1, nx
			p%of(id)%loc%mg(level)%sol(i,j,nz+1)=p%of(id)%loc%mg(level)%sol(i,j,nz)
		enddo
		enddo
	endif
	
!$omp end parallel

end subroutine

subroutine multigrid_compatibility_src(level)
use all
implicit none
integer :: id,i,j,k,level
integer :: nx,ny,nz
real(8) :: s

nx = 0
ny = 0
nz = 0

!$omp parallel private(id,i,j,k), num_threads(p%glb%threads), reduction(+:s,nx,ny,nz)	
	
	id=0
	!$ id = omp_get_thread_num()
	
	nx = nx + p%of(id)%loc%mg(level)%nx
	ny = ny + p%of(id)%loc%mg(level)%ny
	nz = nz + p%of(id)%loc%mg(level)%nz
	
	s=0.0d0
	
	do k = 1, p%of(id)%loc%mg(level)%nz
	do j = 1, p%of(id)%loc%mg(level)%ny
	do i = 1, p%of(id)%loc%mg(level)%nx
		s = s + p%of(id)%loc%mg(level)%src(i,j,k)
	enddo
	enddo
	enddo
	
!$omp end parallel

!s=s/(nx*ny*nz)
s=s/(nx*ny)

!$omp parallel private(id,i,j,k), num_threads(p%glb%threads)
		
	id=0
	!$ id = omp_get_thread_num()
			
	do k = 1,  p%of(id)%loc%mg(level)%nz
	do j = 1,  p%of(id)%loc%mg(level)%ny
	do i = 1,  p%of(id)%loc%mg(level)%nx	
		p%of(id)%loc%mg(level)%src(i,j,k) = p%of(id)%loc%mg(level)%src(i,j,k) - s
	enddo
	enddo
	enddo

!$omp end parallel

end subroutine

subroutine multigrid_compatibility_sol(level)
use all
implicit none
integer :: id,i,j,k,level
integer :: nx,ny,nz
real(8) :: s

nx = 0
ny = 0
nz = 0

!$omp parallel private(id,i,j,k), num_threads(p%glb%threads), reduction(+:s,nx,ny,nz)	
	
	id=0
	!$ id = omp_get_thread_num()
	
	nx = nx + p%of(id)%loc%mg(level)%nx
	ny = ny + p%of(id)%loc%mg(level)%ny
	nz = nz + p%of(id)%loc%mg(level)%nz
	
	s=0.0d0
	
	do k = 1, p%of(id)%loc%mg(level)%nz
	do j = 1, p%of(id)%loc%mg(level)%ny
	do i = 1, p%of(id)%loc%mg(level)%nx
		s = s + p%of(id)%loc%mg(level)%sol(i,j,k)
	enddo
	enddo
	enddo
	
!$omp end parallel

!s=s/(nx*ny*nz)
s=s/(nx*ny)

!$omp parallel private(id,i,j,k), num_threads(p%glb%threads)
		
	id=0
	!$ id = omp_get_thread_num()
			
	do k = 1,  p%of(id)%loc%mg(level)%nz
	do j = 1,  p%of(id)%loc%mg(level)%ny
	do i = 1,  p%of(id)%loc%mg(level)%nx	
		p%of(id)%loc%mg(level)%sol(i,j,k) = p%of(id)%loc%mg(level)%sol(i,j,k) - s
	enddo
	enddo
	enddo

!$omp end parallel

end subroutine
