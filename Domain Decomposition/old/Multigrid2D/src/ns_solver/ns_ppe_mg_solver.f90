subroutine ppe_mg_solver
use all
!$ use omp_lib
implicit none
integer :: iter, level, id
integer :: nx,ny,nz
real(8) :: s, rms0, x, y, exact, pi
real(8) :: px, py, rms1
integer :: i,j,k
integer(8) :: cpustart, cpuend

pi = dacos(-1.0d0)


call ppe_mg_solver_src
	
iter = 0
call system_clock(cpustart)	
do iter = 1, 200

	p%glb%piter=p%glb%piter+1
	
	do level = 1, p%glb%level-1
		call multigrid_relax(level,2)
		call multigrid_restriction(level)
	enddo
	
	!call multigrid_iterative
	call multigrid_final_exact()
	
	do level = p%glb%level-1, 1, -1
		call multigrid_prolongation(level)
		call multigrid_relax(level,2)
	enddo

	call multigrid_residual(1,.false.)
	
	write(*,*)iter,p%of(0)%loc%mg(1)%l2norm
	if( p%of(0)%loc%mg(1)%l2norm .lt. p%glb%p_tol )exit

enddo
call system_clock(cpuend)

p%glb%ppe_linf = p%of(0)%loc%mg(1)%l2norm

!$omp parallel private(id,i,j,k,x,y,exact), num_threads(p%glb%threads), reduction(+:rms0)

	id=0
	!$ id = omp_get_thread_num()
	
	k = 1	
	!do k = p%of(id)%loc%ks, p%of(id)%loc%ke
	do j = p%of(id)%loc%js, p%of(id)%loc%je
	do i = p%of(id)%loc%is, p%of(id)%loc%ie 
		p%of(id)%loc%p%now(i,j,k) = p%of(id)%loc%mg(1)%sol(i-p%of(id)%loc%is+1,&
														&  j-p%of(id)%loc%js+1,&
														&  k-p%of(id)%loc%ks+1)
															
	!enddo
	enddo
	enddo
	
	call p%of(id)%bc(0,p%of(id)%loc%p%now)

!$omp end parallel

call pt%p%sync

rms0=0.0d0
rms1=0.0d0

!$omp parallel private(id,i,j,k,x,y,exact,px,py), num_threads(p%glb%threads), reduction(max:rms0,rms1)

	id=0
	!$ id = omp_get_thread_num()
	
	do j = p%of(id)%loc%js, p%of(id)%loc%je
	do i = p%of(id)%loc%is, p%of(id)%loc%ie 
		p%of(id)%loc%p%now(i,j,k) = p%of(id)%loc%mg(1)%sol(i-p%of(id)%loc%is+1,&
														&  j-p%of(id)%loc%js+1,&
														&  k-p%of(id)%loc%ks+1)
		x = p%glb%x(i)
		y = p%glb%y(j)
		
		exact = -pi*dsin(pi*x)*dcos(pi*y)
		px = 0.5d0*( p%of(id)%loc%p%now(i+1,j,k)-p%of(id)%loc%p%now(i-1,j,k))/p%glb%dx
		
		rms0 = max(	rms0, abs(px - exact ) )
		!rms0 = rms0 + (	px - exact )**2.0d0
		
		exact = -pi*dsin(pi*y)*dcos(pi*x)
		py = 0.5d0*( p%of(id)%loc%p%now(i,j+1,k)-p%of(id)%loc%p%now(i,j-1,k))/p%glb%dx
		
		rms0 = max(	rms0, abs(py - exact ) )
		!rms0 = rms0 + (	py - exact )**2.0d0
		
	enddo
	enddo		 
	
!$omp end parallel

rms0 = dsqrt( rms0 / (p%glb%node_x*p%glb%node_y) )

write(*,'(A,ES11.4)')"Error respect to exact solution:",rms0
write(*,'(A,ES11.4)')"CPU time: ",real(cpuend-cpustart,kind=8)/real(p%glb%cpurate,kind=8)
!call pt%p%sync
!call ppe_mg_correction


p%glb%ppe = p%glb%ppe + real(cpuend-cpustart,kind=8)/real(p%glb%cpurate,kind=8)
	
end subroutine

subroutine ppe_mg_solver_init
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k

	!$omp parallel private(id,i,j,k), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		do k = p%of(id)%loc%ks-p%of(id)%glb%ghc, p%of(id)%loc%ke+p%of(id)%glb%ghc
		do j = p%of(id)%loc%js-p%of(id)%glb%ghc, p%of(id)%loc%je+p%of(id)%glb%ghc
		do i = p%of(id)%loc%is-p%of(id)%glb%ghc, p%of(id)%loc%ie+p%of(id)%glb%ghc
	
			p%of(id)%loc%p%tmp(i,j,k) = 2.0d0*p%of(id)%loc%p%old(i,j,k) - p%of(id)%loc%p%old2(i,j,k)
		
		enddo
		enddo
		enddo
		
		do k = p%of(id)%loc%ks, p%of(id)%loc%ke
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			p%of(id)%loc%coe%r(i,j,k) = 1.0d0 - 2.0d0*p%of(id)%glb%rho_12 &
			& / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i+1,j,k))			
			p%of(id)%loc%coe%l(i,j,k) = 1.0d0 - 2.0d0*p%of(id)%glb%rho_12 &
			& / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i-1,j,k))
			p%of(id)%loc%coe%f(i,j,k) = 1.0d0 - 2.0d0*p%of(id)%glb%rho_12 &
			& / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j+1,k))
			p%of(id)%loc%coe%b(i,j,k) = 1.0d0 - 2.0d0*p%of(id)%glb%rho_12 &
			& / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j-1,k))
			p%of(id)%loc%coe%u(i,j,k) = 1.0d0 - 2.0d0*p%of(id)%glb%rho_12 &
			& / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j,k+1))
			p%of(id)%loc%coe%d(i,j,k) = 1.0d0 - 2.0d0*p%of(id)%glb%rho_12 &
			& / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j,k-1))
						
			p%of(id)%loc%coe%r(i,j,k) = p%of(id)%loc%coe%r(i,j,k) / p%glb%dx**2.0d0
			p%of(id)%loc%coe%l(i,j,k) = p%of(id)%loc%coe%l(i,j,k) / p%glb%dx**2.0d0
			p%of(id)%loc%coe%f(i,j,k) = p%of(id)%loc%coe%f(i,j,k) / p%glb%dy**2.0d0
			p%of(id)%loc%coe%b(i,j,k) = p%of(id)%loc%coe%b(i,j,k) / p%glb%dy**2.0d0
			p%of(id)%loc%coe%u(i,j,k) = p%of(id)%loc%coe%u(i,j,k) / p%glb%dz**2.0d0
			p%of(id)%loc%coe%d(i,j,k) = p%of(id)%loc%coe%d(i,j,k) / p%glb%dz**2.0d0
			
			p%of(id)%loc%coe%c(i,j,k) = - ( p%of(id)%loc%coe%r(i,j,k) + p%of(id)%loc%coe%l(i,j,k) + &
									&	    p%of(id)%loc%coe%f(i,j,k) + p%of(id)%loc%coe%b(i,j,k) + &
									& 	    p%of(id)%loc%coe%u(i,j,k) + p%of(id)%loc%coe%d(i,j,k) )
							
			if( i==1 )then
				p%of(id)%loc%coe%c(i,j,k)=p%of(id)%loc%coe%c(i,j,k)+p%of(id)%loc%coe%l(i,j,k)
				p%of(id)%loc%coe%l(i,j,k)=0.0d0
			endif
			
			if( i==p%glb%node_x )then
				p%of(id)%loc%coe%c(i,j,k)=p%of(id)%loc%coe%c(i,j,k)+p%of(id)%loc%coe%r(i,j,k)
				p%of(id)%loc%coe%r(i,j,k)=0.0d0
			endif
			
			if( j==1 )then
				p%of(id)%loc%coe%c(i,j,k)=p%of(id)%loc%coe%c(i,j,k)+p%of(id)%loc%coe%b(i,j,k)
				p%of(id)%loc%coe%b(i,j,k)=0.0d0
			endif
			
			if( j==p%glb%node_y )then
				p%of(id)%loc%coe%c(i,j,k)=p%of(id)%loc%coe%c(i,j,k)+p%of(id)%loc%coe%f(i,j,k)
				p%of(id)%loc%coe%f(i,j,k)=0.0d0
			endif
			
			if( k==1 )then
				p%of(id)%loc%coe%c(i,j,k)=p%of(id)%loc%coe%c(i,j,k)+p%of(id)%loc%coe%d(i,j,k)
				p%of(id)%loc%coe%d(i,j,k)=0.0d0
			endif
			
			if( k==p%glb%node_z )then
				p%of(id)%loc%coe%c(i,j,k)=p%of(id)%loc%coe%c(i,j,k)+p%of(id)%loc%coe%u(i,j,k)
				p%of(id)%loc%coe%u(i,j,k)=0.0d0
			endif
			
		end do
		end do
		end do
		
	!$omp end parallel 
	
end subroutine

subroutine ppe_mg_solver_src
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: pi,x,y

pi = dacos(-1.0d0)

!$omp parallel private(id,i,j,k,x,y), num_threads(p%glb%threads)

	id=0
	!$ id = omp_get_thread_num()
	
	! do k = p%of(id)%loc%ks, p%of(id)%loc%ke
	! do j = p%of(id)%loc%js, p%of(id)%loc%je
	! do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
		! p%of(id)%loc%coe%src(i,j,k) = ( ( p%of(id)%loc%vel%x%now(i,j,k) - p%of(id)%loc%vel%x%now(i-1,j,k) ) / p%glb%dx + &
									! &	( p%of(id)%loc%vel%y%now(i,j,k) - p%of(id)%loc%vel%y%now(i,j-1,k) ) / p%glb%dy + & 
									! &	( p%of(id)%loc%vel%z%now(i,j,k) - p%of(id)%loc%vel%z%now(i,j,k-1) ) / p%glb%dz ) / p%glb%dt
										
		! p%of(id)%loc%coe%src(i,j,k) = p%of(id)%loc%coe%src(i,j,k) * p%of(id)%glb%rho_12 + &
									! & p%of(id)%loc%coe%c(i,j,k) * p%of(id)%loc%p%tmp(i,j,k)   + &
									! & p%of(id)%loc%coe%r(i,j,k) * p%of(id)%loc%p%tmp(i+1,j,k) + p%of(id)%loc%coe%l(i,j,k) * p%of(id)%loc%p%tmp(i-1,j,k) + &
									! & p%of(id)%loc%coe%f(i,j,k) * p%of(id)%loc%p%tmp(i,j+1,k) + p%of(id)%loc%coe%b(i,j,k) * p%of(id)%loc%p%tmp(i,j-1,k) + &
									! & p%of(id)%loc%coe%u(i,j,k) * p%of(id)%loc%p%tmp(i,j,k+1) + p%of(id)%loc%coe%d(i,j,k) * p%of(id)%loc%p%tmp(i,j,k-1) 
		
	! enddo
	! enddo
	! enddo
	
	do k = 1, p%of(id)%loc%mg(1)%nz
	do j = 1, p%of(id)%loc%mg(1)%ny
	do i = 1, p%of(id)%loc%mg(1)%nx
		
		!p%of(id)%loc%mg(1)%src(i,j,k) = p%of(id)%loc%coe%src(i+p%of(id)%loc%is-1,j+p%of(id)%loc%js-1,k+p%of(id)%loc%ks-1)
		!p%of(id)%loc%mg(1)%sol(i,j,k) = p%of(id)%loc%p%tmp(i+p%of(id)%loc%is-1,j+p%of(id)%loc%js-1,k+p%of(id)%loc%ks-1)
		
		x = p%glb%x(i+p%of(id)%loc%is-1)
		y = p%glb%y(j+p%of(id)%loc%js-1)
		
		p%of(id)%loc%mg(1)%src(i,j,k) = -2.0d0*pi**2.0d0*dcos(pi*x)*dcos(pi*y)
		p%of(id)%loc%mg(1)%sol(i,j,k) = 0.0d0
		
	enddo
	enddo
	enddo
	
!$omp end parallel

	call multigrid_sync(1)

end subroutine

subroutine ppe_mg_correction
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: src

	!$omp parallel private(id,i,j,k,src), num_threads(p%glb%threads)
		
		id=0
		!$ id = omp_get_thread_num()
		
		do k = p%of(id)%loc%ks, p%of(id)%loc%ke
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
			
			src = 1.0d0 / p%glb%rho_12 * ( p%of(id)%loc%p%now(i+1,j,k) - p%of(id)%loc%p%now(i,j,k) )/p%glb%dx
			src = src + ( 2.0d0/(p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i+1,j,k)) - 1.0d0/p%glb%rho_12 ) * &
			& ( p%of(id)%loc%p%tmp(i+1,j,k) - p%of(id)%loc%p%tmp(i,j,k) )/p%glb%dx
			
			p%of(id)%loc%vel%x%now(i,j,k) = p%of(id)%loc%vel%x%now(i,j,k) - src * p%glb%dt
			
			src = 1.0d0 / p%glb%rho_12 * ( p%of(id)%loc%p%now(i,j+1,k) - p%of(id)%loc%p%now(i,j,k) )/p%glb%dy
			src = src + ( 2.0d0/(p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j+1,k)) - 1.0d0/p%glb%rho_12 ) * &
			& ( p%of(id)%loc%p%tmp(i,j+1,k) - p%of(id)%loc%p%tmp(i,j,k) )/p%glb%dy
			
			p%of(id)%loc%vel%y%now(i,j,k) = p%of(id)%loc%vel%y%now(i,j,k) - src * p%glb%dt
			
			src = 1.0d0 / p%glb%rho_12 * ( p%of(id)%loc%p%now(i,j,k+1) - p%of(id)%loc%p%now(i,j,k) )/p%glb%dz
			src = src + ( 2.0d0/(p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j,k+1)) - 1.0d0/p%glb%rho_12 ) * &
			& ( p%of(id)%loc%p%tmp(i,j,k+1) - p%of(id)%loc%p%tmp(i,j,k) )/p%glb%dz
			
			p%of(id)%loc%vel%z%now(i,j,k) = p%of(id)%loc%vel%z%now(i,j,k) - src * p%glb%dt
			
		enddo
		enddo
		enddo
		
		call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now,p%of(id)%loc%vel%z%now)
	
	!$omp end parallel
	
	call pt%vel%sync
	
end subroutine
