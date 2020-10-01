subroutine deforming_field_velocity()
use all 
!$ use omp_lib
implicit none
integer :: i,j,id
real(8) :: x, y, xx, yy, pi, ct

	pi = dacos(-1.0_8)
	ct = dcos( pi*p%glb%time / p%glb%t2s )

	!$omp parallel private(id,i,j,x,y,xx,yy), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			x = p%glb%x(i)
			y = p%glb%y(j)
			
			xx = 0.5_8 * ( p%glb%x(i) + p%glb%x(i+1) )
			yy = 0.5_8 * ( p%glb%y(j) + p%glb%y(j+1) )
			
			p%of(id)%loc%nvel%x%now(i,j) =   dsin(pi*x)**2.0_8 * dsin(2.0_8*pi*y) * ct
			p%of(id)%loc%nvel%y%now(i,j) = - dsin(2.0_8*pi*x) * dsin(pi*y)**2.0_8 * ct
			
			p%of(id)%loc%vel%x%now(i,j) =   dsin(pi*xx)**2.0_8 * dsin(2.0_8*pi*y) * ct
			p%of(id)%loc%vel%y%now(i,j) = - dsin(2.0_8*pi*x) * dsin(pi*yy)**2.0_8 * ct
			
		enddo
		enddo
		
		call p%of(id)%bc(2,p%of(id)%loc%nvel%x%now)
		call p%of(id)%bc(2,p%of(id)%loc%nvel%y%now)
		
		call p%of(id)%bc(2,p%of(id)%loc%vel%x%now)
		call p%of(id)%bc(2,p%of(id)%loc%vel%y%now)
	
	!$omp end parallel 
	
	call pt%vel%sync
	
	call pt%nvel%sync	

end subroutine

