subroutine level_set_srk6_solver()
use all
!$ use omp_lib
implicit none
integer :: i,j,id,iter
integer(8) :: cpustart, cpuend
real(8) :: err,perr

	call system_clock(cpustart)

	!$omp parallel private(id,i,j), num_threads(p%glb%threads)
	
		id = 0
		!$ id = omp_get_thread_num()
		
		call p%of(id)%loc%srk6%init(.false.,p%of(id)%loc%phi%old)

	!$omp end parallel 
	
	iter=0
	
	do 
	
		iter = iter + 1
		
		call level_set_srk6_source()
		
		err=0.0_8
		
		!$omp parallel private(id), num_threads(p%glb%threads), reduction(max:err)
			
			id = 0
			!$ id = omp_get_thread_num()
		
			call p%of(id)%loc%srk6%solve(err)

		!$omp end parallel 
	
		if( err < p%glb%srk6_tol )exit
		
		if( mod(iter,500) == 0)write(*,*)"LS solver,",iter,err
		
	end do
	
	!$omp parallel private(id,i,j), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		call p%of(id)%loc%srk6%final()
		
		do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
		do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
			p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%srk6%x%target(i,j)
		end do
		end do
		
		call p%of(id)%bc(0,p%of(id)%loc%phi%now)

	!$omp end parallel

	call pt%phi%sync
	
	call system_clock(cpuend)
	p%glb%ls_adv = p%glb%ls_adv + real(cpuend-cpustart,kind=8)/real(p%glb%cpurate,kind=8)

end subroutine

subroutine level_set_srk6_source()
use all 
!$ use omp_lib
implicit none
integer :: i,j,id

	!$omp parallel private(id,i,j), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		call p%of(id)%bc(0,p%of(id)%loc%srk6%x%s1)
		call p%of(id)%bc(0,p%of(id)%loc%srk6%x%s2)
		call p%of(id)%bc(0,p%of(id)%loc%srk6%x%s3)
		
	!$omp end parallel 
	
	call pt%srk6_x%sync
	
	!$omp parallel private(id,i,j), num_threads(p%glb%threads)
		
		id=0
		!$ id = omp_get_thread_num()
		
		! calculate x derivatives
		do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
		
			call p%of(id)%loc%ccd%x%SOLVE_SRKCCD(p%of(id)%loc%nvel%x%old(:,j),&
												&p%of(id)%loc%srk6%x%s1(:,j),&
												&p%of(id)%loc%srk6%x%ss1(:,j))
									 	
			call p%of(id)%loc%ccd%x%SOLVE_SRKCCD(p%of(id)%loc%nvel%x%old(:,j),&
												&p%of(id)%loc%srk6%x%s2(:,j),&
												&p%of(id)%loc%srk6%x%ss2(:,j))
									 	
			call p%of(id)%loc%ccd%x%SOLVE_SRKCCD(p%of(id)%loc%nvel%x%old(:,j),&
												&p%of(id)%loc%srk6%x%s3(:,j),&
												&p%of(id)%loc%srk6%x%ss3(:,j))
			
		enddo
		
		do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
		do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
			
			p%of(id)%loc%srk6%x%l1(i,j) = - p%of(id)%loc%srk6%x%ss1(i,j) * p%of(id)%loc%nvel%x%old(i,j)
			p%of(id)%loc%srk6%x%l2(i,j) = - p%of(id)%loc%srk6%x%ss2(i,j) * p%of(id)%loc%nvel%x%old(i,j)
			p%of(id)%loc%srk6%x%l3(i,j) = - p%of(id)%loc%srk6%x%ss3(i,j) * p%of(id)%loc%nvel%x%old(i,j)
			
		end do
		end do
		
		! calculate y derivatives
		do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
			
			call p%of(id)%loc%ccd%y%SOLVE_SRKCCD(p%of(id)%loc%nvel%y%old(i,:),&
												&p%of(id)%loc%srk6%x%s1(i,:),&
												&p%of(id)%loc%srk6%x%ss1(i,:))
										
			call p%of(id)%loc%ccd%y%SOLVE_SRKCCD(p%of(id)%loc%nvel%y%old(i,:),&
												&p%of(id)%loc%srk6%x%s2(i,:),&
												&p%of(id)%loc%srk6%x%ss2(i,:))
										
			call p%of(id)%loc%ccd%y%SOLVE_SRKCCD(p%of(id)%loc%nvel%y%old(i,:),&
												&p%of(id)%loc%srk6%x%s3(i,:),&
												&p%of(id)%loc%srk6%x%ss3(i,:))
		end do
		
		do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
		do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
		
			p%of(id)%loc%srk6%x%l1(i,j) = p%of(id)%loc%srk6%x%l1(i,j) - p%of(id)%loc%srk6%x%ss1(i,j) * p%of(id)%loc%nvel%y%old(i,j)
			p%of(id)%loc%srk6%x%l2(i,j) = p%of(id)%loc%srk6%x%l2(i,j) - p%of(id)%loc%srk6%x%ss2(i,j) * p%of(id)%loc%nvel%y%old(i,j)
			p%of(id)%loc%srk6%x%l3(i,j) = p%of(id)%loc%srk6%x%l3(i,j) - p%of(id)%loc%srk6%x%ss3(i,j) * p%of(id)%loc%nvel%y%old(i,j)
			
		end do
		end do
		
	!$omp end parallel
	
end subroutine