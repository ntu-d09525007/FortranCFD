subroutine mass_preserving_level_set
use all
!$ use omp_lib
implicit none
integer :: i,j,id,iter
real(8) :: lam, plam

do iter = 1, 2
    
	call p%ls_mv()
	call p%surface_norms()
	call pt%normals_1%sync
	
	lam = 0.0_8
	
	!$omp parallel private(id,i,j,plam), num_threads(p%glb%threads), reduction(+:lam)
	
		id=0
		!$ id = omp_get_thread_num()
		
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			p%of(id)%loc%grad%now(i,j) = dsqrt( p%of(id)%loc%normals%x%now(i,j)**2.0_8 +  p%of(id)%loc%normals%y%now(i,j)**2.0_8 )
			
			plam =  p%of(id)%loc%delta%now(i,j)**2.0_8 * p%of(id)%loc%grad%now(i,j) 
			plam = plam * ( 2.0_8*(1.0_8-p%glb%rho_12)*p%of(id)%loc%heavy%now(i,j) + p%glb%rho_12 )*p%glb%dx*p%glb%dy
			
			lam = lam + plam
			
		end do
		end do
		
	!$omp end parallel
	
	lam = ( p%glb%imass - p%glb%mass ) / lam 

	!$omp parallel private(id,i,j), num_threads(p%glb%threads), shared(lam)
	
		id=0
		!$ id = omp_get_thread_num()
			
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%phi%now(i,j) + lam * p%of(id)%loc%delta%now(i,j) * p%of(id)%loc%grad%now(i,j)
			
		end do
		end do
		
		call p%of(id)%bc(0,p%of(id)%loc%phi%now)
		
	!$omp end parallel
	
	call pt%phi%sync
	
end do	

end subroutine

