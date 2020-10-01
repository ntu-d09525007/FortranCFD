subroutine problem_init()
use all
!$ use omp_lib
implicit none
integer :: id, i, j
real(8) :: x, y, err
CHARACTER(100) :: NAME_OF_FILE
	
	NAME_OF_FILE='d'
	
	WRITE(*,*)"============================================"
	WRITE(*,*)'Returning the files in directory "/input" '
	WRITE(*,*)
	CALL SYSTEM("ls ./input")
	WRITE(*,*)
	WRITE(*,*)"============================================"
	WRITE(*,*)'Selet an input file (<d> : default.txt), with full name'
	WRITE(*,*)
	!READ(*,*)NAME_OF_FILE
	if( name_of_file == 'd' )name_of_file = "default.txt"
    
	
	call p%init( "./input/"//trim(name_of_file) )
	call plt%init(p)
	call pt%init(p)
	
	!---------------------------------------------------
	
	!$omp parallel private(id,i,j,x,y), num_threads(p%glb%threads)
		
		id=0
		!$ id = omp_get_thread_num()
		
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			x = p%glb%x(i)
			y = p%glb%y(j)
			
			!p%of(id)%loc%phi%now(i,j) = - dsqrt( (x-0.5_8)**2 + (y-0.5_8)**2 ) + 0.25d0
			
			if( x<=1.0d0 .and. y<=1.0d0 )then
				p%of(id)%loc%phi%now(i,j) = 1.0_8
			else
				p%of(id)%loc%phi%now(i,j) = -1.0_8
			end if
			
			!p%of(id)%loc%phi%now(i,j) =   dsqrt( x**2.0d0 + (y-1.5d0)**2.0d0 ) - 1.0d0
			
			x = 0.5_8*( p%glb%x(i)+p%glb%x(i+1) )
			y = 0.5_8*( p%glb%y(j)+p%glb%y(j+1) )
			
			p%of(id)%loc%vel%x%now(i,j) = 0.0_8
			p%of(id)%loc%vel%y%now(i,j) = 0.0_8
			
		end do
		end do
		
		call p%of(id)%bc(0,p%of(id)%loc%phi%now)
		call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now)
		
	!$omp end parallel
	
	call pt%phi%sync
	call level_set_rk3_redis(0)
	
	call pt%vel%sync
	call p%node_vel
	call pt%nvel%sync
	call ns_init
	
	!---------------------------------------------------

	!$omp parallel private(id), num_threads(p%glb%threads)	
		id=0
		!$ id = omp_get_thread_num()
		call p%of(id)%loc%init	
	!$omp end parallel 
	
	call ns_adv_source
	call ns_diff_source
	call ns_final
	
	call p%ls_mv
	p%glb%ivol = p%glb%vol
	p%glb%imass = p%glb%mass
	call p%sync

	call plt%plot

		
end subroutine
