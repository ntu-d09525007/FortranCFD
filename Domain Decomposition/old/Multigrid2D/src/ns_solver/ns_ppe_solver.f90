subroutine ppe_sor_init
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k

	!$omp parallel private(id,i,j,k), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		do k = p%of(id)%loc%ks, p%of(id)%loc%ke
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			p%of(id)%loc%coe%r(i,j,k) = 2.0d0 / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i+1,j,k))
			p%of(id)%loc%coe%l(i,j,k) = 2.0d0 / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i-1,j,k))
			p%of(id)%loc%coe%f(i,j,k) = 2.0d0 / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j+1,k))
			p%of(id)%loc%coe%b(i,j,k) = 2.0d0 / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j-1,k))
			p%of(id)%loc%coe%u(i,j,k) = 2.0d0 / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j,k+1))
			p%of(id)%loc%coe%d(i,j,k) = 2.0d0 / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j,k-1))
			
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

subroutine ppe_sor_solver
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k,iter
integer(8) :: cpustart, cpuend
real(8) :: sump, err, w, pcal
real(8) :: CC,CR,CL,CF,CB,CU,CD,RR,RL,RF,RB,RU,RD

	call system_clock(cpustart)

	!$omp parallel private(id,i,j,k), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		do k = p%of(id)%loc%ks, p%of(id)%loc%ke
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			p%of(id)%loc%coe%src(i,j,k) = ( ( p%of(id)%loc%vel%x%now(i,j,k) - p%of(id)%loc%vel%x%now(i-1,j,k) ) / p%glb%dx + &
										    ( p%of(id)%loc%vel%y%now(i,j,k) - p%of(id)%loc%vel%y%now(i,j-1,k) ) / p%glb%dy + & 
											( p%of(id)%loc%vel%z%now(i,j,k) - p%of(id)%loc%vel%z%now(i,j,k-1) ) / p%glb%dz ) / p%glb%dt
										  			
		end do
		end do
		end do
		
	!$omp end parallel
	
	w = p%glb%p_w1
	
do
	
	p%glb%piter=p%glb%piter+1

	!$omp parallel private(id,i,j,k,CC,CR,CL,CF,CB,CU,CD,RR,RL,RF,RB,RU,RD), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()

		do k = p%of(id)%loc%ks, p%of(id)%loc%ke
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			p%of(id)%loc%p%tmp(i,j,k) = p%of(id)%loc%p%now(i,j,k)
			
			p%of(id)%loc%p%now(i,j,k) = p%of(id)%loc%coe%src(i,j,k) - p%of(id)%loc%coe%r(i,j,k)*p%of(id)%loc%p%now(i+1,j,k) &
																&   - p%of(id)%loc%coe%l(i,j,k)*p%of(id)%loc%p%now(i-1,j,k) &
																&   - p%of(id)%loc%coe%f(i,j,k)*p%of(id)%loc%p%now(i,j+1,k) &
																&   - p%of(id)%loc%coe%b(i,j,k)*p%of(id)%loc%p%now(i,j-1,k) &
																&   - p%of(id)%loc%coe%u(i,j,k)*p%of(id)%loc%p%now(i,j,k+1) &
																&   - p%of(id)%loc%coe%d(i,j,k)*p%of(id)%loc%p%now(i,j,k-1)
															
			p%of(id)%loc%p%now(i,j,k) = p%of(id)%loc%p%now(i,j,k) / p%of(id)%loc%coe%c(i,j,k)	
  
		end do
		end do
		end do
	
	!$omp end parallel 
	
	sump=0.0d0
		
	!$omp parallel private(id,i,j,k), num_threads(p%glb%threads), reduction(+:sump)
	
		id=0
		!$ id = omp_get_thread_num()
		sump=0.0d0
		
		do k = p%of(id)%loc%ks, p%of(id)%loc%ke
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			p%of(id)%loc%p%now(i,j,k) = w * p%of(id)%loc%p%now(i,j,k) + (1.0d0-w)*p%of(id)%loc%p%tmp(i,j,k)
			
			sump = sump + p%of(id)%loc%p%now(i,j,k)

		end do
		end do
		end do
	
	!$omp end parallel 
	
	sump = sump / ( p%glb%node_x * p%glb%node_y * p%glb%node_z )
	err=0.0d0
	
	!$omp parallel private(id,i,j,k), num_threads(p%glb%threads), reduction(+:err), shared(sump)
	
		id=0
		!$ id = omp_get_thread_num()
		
		err=0.0d0
		
		do k = p%of(id)%loc%ks, p%of(id)%loc%ke
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
			
			p%of(id)%loc%p%now(i,j,k) = p%of(id)%loc%p%now(i,j,k) - sump
			
			p%of(id)%loc%p%tmp(i,j,k) = p%of(id)%loc%coe%src(i,j,k) - p%of(id)%loc%coe%r(i,j,k)*p%of(id)%loc%p%now(i+1,j,k) &
																&   - p%of(id)%loc%coe%l(i,j,k)*p%of(id)%loc%p%now(i-1,j,k) &
																&   - p%of(id)%loc%coe%f(i,j,k)*p%of(id)%loc%p%now(i,j+1,k) &
																&   - p%of(id)%loc%coe%b(i,j,k)*p%of(id)%loc%p%now(i,j-1,k) &
																&   - p%of(id)%loc%coe%u(i,j,k)*p%of(id)%loc%p%now(i,j,k+1) &
																&   - p%of(id)%loc%coe%d(i,j,k)*p%of(id)%loc%p%now(i,j,k-1) &
																&   - p%of(id)%loc%coe%c(i,j,k)*p%of(id)%loc%p%now(i,j,k) 
																
			err = err + p%of(id)%loc%p%tmp(i,j,k)**2.0d0
			
		end do
		end do
		end do
		
		call p%of(id)%bc(0,p%of(id)%loc%p%now)
		
	!$omp end parallel
	
	err = dsqrt( err / (p%glb%node_x*p%glb%node_y*p%glb%node_z) )
	
	call pt%p%sync
	
	if( err < p%glb%p_tol ) exit
	if( err < p%glb%p_b ) w = p%glb%p_w2
	if( err > 10 .and. p%glb%piter > 10000 )then
		write(*,*)"The solution can not converge in PPE :",err
		stop
	end if
	
	if( mod(p%glb%piter,5000) .eq. 0 )then
		write(*,'("PPE iter:",I8,",error:",ES15.7)')p%glb%piter,err
	endif
	
end do

	p%glb%ppe_linf = err
	
	call system_clock(cpuend)
	p%glb%ppe = p%glb%ppe + real(cpuend-cpustart,kind=8)/real(p%glb%cpurate,kind=8)
	
	call ns_momentum_correction
	
end subroutine

subroutine ns_momentum_correction
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: px,py,pz,rho,ux,vy,wz

	!$omp parallel private(id,i,j,k,px,py,pz,rho), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		do k = p%of(id)%loc%ks, p%of(id)%loc%ke
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
			
			rho = (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i+1,j,k))/2.0d0
			px  = (p%of(id)%loc%p%now(i+1,j,k)-p%of(id)%loc%p%now(i,j,k)) / p%glb%dx  
			p%of(id)%loc%vel%x%now(i,j,k) = p%of(id)%loc%vel%x%now(i,j,k) - p%glb%dt*px/rho
			
			rho = (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j+1,k))/2.0d0
			py  = (p%of(id)%loc%p%now(i,j+1,k)-p%of(id)%loc%p%now(i,j,k)) / p%glb%dy 		
			p%of(id)%loc%vel%y%now(i,j,k) = p%of(id)%loc%vel%y%now(i,j,k) - p%glb%dt*py/rho
			
			rho = (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j,k+1))/2.0d0
			pz  = (p%of(id)%loc%p%now(i,j,k+1)-p%of(id)%loc%p%now(i,j,k)) / p%glb%dz 		
			p%of(id)%loc%vel%z%now(i,j,k) = p%of(id)%loc%vel%z%now(i,j,k) - p%glb%dt*pz/rho
			
		end do
		end do
		end do
	
		call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now,p%of(id)%loc%vel%z%now)
		
	!$omp end parallel
	
	call pt%vel%sync
	
end subroutine