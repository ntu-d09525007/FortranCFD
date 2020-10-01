subroutine ns_solver
use all
implicit none
integer :: iter
integer(8) :: cpustart, cpuend

	call system_clock(cpustart)

	call ns_init
	call ns_diff_source
	call ppe_init
	
	iter=0
	p%glb%piter=0
	
	do 
	
		iter=iter+1
		
		call ns_adv_source
		call ns_predictor
		call ppe_solver
		call ns_check_convergence

		if(iter>2)exit
		
	end do

	call ns_final
	
	call system_clock(cpuend)
	p%glb%ns = p%glb%ns + real(cpuend-cpustart,kind=8)/real(p%glb%cpurate,kind=8)
	
end subroutine

subroutine ns_init
use all
!$ use omp_lib
implicit none
integer :: id,i,j

	call p%rho_mu
	
	call p%surface_norms
	call pt%normals_1%sync
	call pt%normals_2%sync
	
	call p%curv
	call pt%normals_2%sync

end subroutine

subroutine ns_adv_source
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: u,v,xp,xm,yp,ym

	call ns_adv_source_quick
	return

	!$omp parallel private(id,i,j,u,v,xp,xm,yp,ym), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
	
		do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
		do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
			p%of(id)%loc%vel%x%tmp(i,j) = p%of(id)%loc%vel%x%now(i,j)
			p%of(id)%loc%vel%y%tmp(i,j) = p%of(id)%loc%vel%y%now(i,j)
		end do
		end do
		
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			u = p%of(id)%loc%vel%x%tmp(i,j)
			v = (p%of(id)%loc%vel%y%tmp(i+1,j-1)+p%of(id)%loc%vel%y%tmp(i+1,j)+p%of(id)%loc%vel%y%tmp(i,j)+p%of(id)%loc%vel%y%tmp(i,j-1))/4.0_8
			
			xp = 0.5_8*(-p%of(id)%loc%vel%x%old(i+2,j)+4.0_8*p%of(id)%loc%vel%x%old(i+1,j)-3.0_8*p%of(id)%loc%vel%x%old(i,j))/p%glb%dx 
			xm = 0.5_8*( p%of(id)%loc%vel%x%old(i-2,j)-4.0_8*p%of(id)%loc%vel%x%old(i-1,j)+3.0_8*p%of(id)%loc%vel%x%old(i,j))/p%glb%dx
			
			yp = 0.5_8*(-p%of(id)%loc%vel%x%old(i,j+2)+4.0_8*p%of(id)%loc%vel%x%old(i,j+1)-3.0_8*p%of(id)%loc%vel%x%old(i,j))/p%glb%dy
			ym = 0.5_8*( p%of(id)%loc%vel%x%old(i,j-2)-4.0_8*p%of(id)%loc%vel%x%old(i,j-1)+3.0_8*p%of(id)%loc%vel%x%old(i,j))/p%glb%dy
			
			p%of(id)%loc%velsrc%x%tmp(i,j) = - ((u+abs(u))*xm+(u-abs(u))*xp)/2.0_8 - ((v+abs(v))*ym+(v-abs(v))*yp)/2.0_8
					
			!-----------------------------------------------------------
			
			u = (p%of(id)%loc%vel%x%tmp(i-1,j+1)+p%of(id)%loc%vel%x%tmp(i,j+1)+p%of(id)%loc%vel%x%tmp(i,j)+p%of(id)%loc%vel%x%tmp(i-1,j))/4.0_8
			v = p%of(id)%loc%vel%y%tmp(i,j)
			
			xp = 0.5_8*(-p%of(id)%loc%vel%y%old(i+2,j)+4.0_8*p%of(id)%loc%vel%y%old(i+1,j)-3.0_8*p%of(id)%loc%vel%y%old(i,j))/p%glb%dx 
			xm = 0.5_8*( p%of(id)%loc%vel%y%old(i-2,j)-4.0_8*p%of(id)%loc%vel%y%old(i-1,j)+3.0_8*p%of(id)%loc%vel%y%old(i,j))/p%glb%dx
			
			yp = 0.5_8*(-p%of(id)%loc%vel%y%old(i,j+2)+4.0_8*p%of(id)%loc%vel%y%old(i,j+1)-3.0_8*p%of(id)%loc%vel%y%old(i,j))/p%glb%dy
			ym = 0.5_8*( p%of(id)%loc%vel%y%old(i,j-2)-4.0_8*p%of(id)%loc%vel%y%old(i,j-1)+3.0_8*p%of(id)%loc%vel%y%old(i,j))/p%glb%dy
			
			p%of(id)%loc%velsrc%y%tmp(i,j) = - ((u+abs(u))*xm+(u-abs(u))*xp)/2.0_8 - ((v+abs(v))*ym+(v-abs(v))*yp)/2.0_8
			
		end do
		end do 
		
	!$omp end parallel
	
end subroutine

subroutine ns_adv_source_quick
!--------------------------------------
! u*Px = [ u_{i}+u_{i+1} ] * [ P_{i+1/2}-P_{i-1/2} ] / dx
!--------------------------------------
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: uh, vh

	!$omp parallel private(id,i,j), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
	
		do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
		do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
			p%of(id)%loc%vel%x%tmp(i,j) = p%of(id)%loc%vel%x%now(i,j)
			p%of(id)%loc%vel%y%tmp(i,j) = p%of(id)%loc%vel%y%now(i,j)
		end do
		end do
		
		call p%of(id)%find_stag_vel( p%of(id)%loc%nvel%x%tmp,p%of(id)%loc%nvel%y%tmp,&
									&p%of(id)%loc%vel%x%tmp ,p%of(id)%loc%vel%y%tmp )
		
	!$omp end parallel
	
	call pt%nveltmp%sync
	
	!$omp parallel private(id,i,j,uh,vh), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		do j = p%of(id)%loc%js-1, p%of(id)%loc%je
		do i = p%of(id)%loc%is-1, p%of(id)%loc%ie
		
			uh = 0.5d0*( p%of(id)%loc%vel%x%tmp(i,j) + p%of(id)%loc%vel%x%tmp(i+1,j) )
			vh = 0.5d0*( p%of(id)%loc%vel%y%tmp(i,j) + p%of(id)%loc%vel%y%tmp(i+1,j) )
			
			if( uh>=0.0d0 )then
				p%of(id)%loc%srk6%x%s1(i,j) = ( -p%of(id)%loc%vel%x%old(i-1,j)+6.0d0*p%of(id)%loc%vel%x%old(i,j)+3.0d0*p%of(id)%loc%vel%x%old(i+1,j) )/8.0d0
			else
				p%of(id)%loc%srk6%x%s1(i,j) = ( -p%of(id)%loc%vel%x%old(i+2,j)+6.0d0*p%of(id)%loc%vel%x%old(i+1,j)+3.0d0*p%of(id)%loc%vel%x%old(i,j) )/8.0d0
			endif
			
			if( vh>=0.0d0 )then
				p%of(id)%loc%srk6%x%s2(i,j) = ( -p%of(id)%loc%vel%x%old(i,j-1)+6.0d0*p%of(id)%loc%vel%x%old(i,j)+3.0d0*p%of(id)%loc%vel%x%old(i,j+1) )/8.0d0
			else
				p%of(id)%loc%srk6%x%s2(i,j) = ( -p%of(id)%loc%vel%x%old(i,j+2)+6.0d0*p%of(id)%loc%vel%x%old(i,j+1)+3.0d0*p%of(id)%loc%vel%x%old(i,j) )/8.0d0
			endif
			
			!-----------------------------------------------------------
			
			uh = 0.5d0*( p%of(id)%loc%vel%x%tmp(i,j) + p%of(id)%loc%vel%x%tmp(i,j+1) )
			vh = 0.5d0*( p%of(id)%loc%vel%y%tmp(i,j) + p%of(id)%loc%vel%y%tmp(i,j+1) )
			
			if( uh>=0.0d0 )then
				p%of(id)%loc%srk6%y%s1(i,j) = ( -p%of(id)%loc%vel%y%old(i-1,j)+6.0d0*p%of(id)%loc%vel%y%old(i,j)+3.0d0*p%of(id)%loc%vel%y%old(i+1,j) )/8.0d0
			else
				p%of(id)%loc%srk6%y%s1(i,j) = ( -p%of(id)%loc%vel%y%old(i+2,j)+6.0d0*p%of(id)%loc%vel%y%old(i+1,j)+3.0d0*p%of(id)%loc%vel%y%old(i,j) )/8.0d0
			endif
			
			if( vh>=0.0d0 )then
				p%of(id)%loc%srk6%y%s2(i,j) = ( -p%of(id)%loc%vel%y%old(i,j-1)+6.0d0*p%of(id)%loc%vel%y%old(i,j)+3.0d0*p%of(id)%loc%vel%y%old(i,j+1) )/8.0d0
			else
				p%of(id)%loc%srk6%y%s2(i,j) = ( -p%of(id)%loc%vel%y%old(i,j+2)+6.0d0*p%of(id)%loc%vel%y%old(i,j+1)+3.0d0*p%of(id)%loc%vel%y%old(i,j) )/8.0d0				!p%of(id)%loc%srk6%y%s2(i,j) = ( -p%of(id)%loc%srk6%y%l2(i,j+2)+6.0d0*p%of(id)%loc%srk6%y%l2(i,j+1)+3.0d0*p%of(id)%loc%srk6%y%l2(i,j) )/8.0d0
			endif
			
			!-----------------------------------------------------------
			
		end do
		end do
		
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			uh = 0.5d0*( p%of(id)%loc%vel%x%tmp(i,j) + p%of(id)%loc%vel%x%tmp(i+1,j) )
			vh = 0.5d0*( p%of(id)%loc%vel%y%tmp(i,j) + p%of(id)%loc%vel%y%tmp(i+1,j) )
		
			p%of(id)%loc%velsrc%x%tmp(i,j) = - uh*( p%of(id)%loc%srk6%x%s1(i,j)-p%of(id)%loc%srk6%x%s1(i-1,j) )/p%glb%dx &
											&- vh*( p%of(id)%loc%srk6%x%s2(i,j)-p%of(id)%loc%srk6%x%s2(i,j-1) )/p%glb%dy
			
			uh = 0.5d0*( p%of(id)%loc%vel%x%tmp(i,j) + p%of(id)%loc%vel%x%tmp(i,j+1) )
			vh = 0.5d0*( p%of(id)%loc%vel%y%tmp(i,j) + p%of(id)%loc%vel%y%tmp(i,j+1) )
			
			p%of(id)%loc%velsrc%y%tmp(i,j) = - uh*( p%of(id)%loc%srk6%y%s1(i,j)-p%of(id)%loc%srk6%y%s1(i-1,j) )/p%glb%dx &
											&- vh*( p%of(id)%loc%srk6%y%s2(i,j)-p%of(id)%loc%srk6%y%s2(i,j-1) )/p%glb%dy
																			
		end do
		end do
	
	!$omp end parallel
	
end subroutine

subroutine ns_adv_source_quick2
!--------------------------------------
! u*Px = [ (uP)_{i+1/2}-(uP)_{i-1/2} ] / dx
!--------------------------------------
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: uh, vh

	!$omp parallel private(id,i,j), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
	
		do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
		do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
			p%of(id)%loc%vel%x%tmp(i,j) = p%of(id)%loc%vel%x%now(i,j)
			p%of(id)%loc%vel%y%tmp(i,j) = p%of(id)%loc%vel%y%now(i,j)
		end do
		end do
		
		call p%of(id)%find_stag_vel( p%of(id)%loc%nvel%x%tmp,p%of(id)%loc%nvel%y%tmp,&
									&p%of(id)%loc%vel%x%tmp ,p%of(id)%loc%vel%y%tmp )
		
	!$omp end parallel
	
	call pt%nveltmp%sync
	
	!$omp parallel private(id,i,j,uh,vh), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
		do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
			p%of(id)%loc%srk6%x%l1(i,j) = p%of(id)%loc%vel%x%tmp(i,j)*p%of(id)%loc%vel%x%old(i,j)
			p%of(id)%loc%srk6%x%l2(i,j) = p%of(id)%loc%nvel%y%tmp(i,j)*p%of(id)%loc%vel%x%old(i,j)
			
			p%of(id)%loc%srk6%y%l1(i,j) = p%of(id)%loc%nvel%x%tmp(i,j)*p%of(id)%loc%vel%y%old(i,j)
			p%of(id)%loc%srk6%y%l2(i,j) = p%of(id)%loc%vel%y%tmp(i,j)*p%of(id)%loc%vel%y%old(i,j)
		end do
		end do
		
		do j = p%of(id)%loc%js-1, p%of(id)%loc%je
		do i = p%of(id)%loc%is-1, p%of(id)%loc%ie
		
			uh = 0.5d0*( p%of(id)%loc%vel%x%tmp(i,j) + p%of(id)%loc%vel%x%tmp(i+1,j) )
			vh = 0.5d0*( p%of(id)%loc%vel%y%tmp(i,j) + p%of(id)%loc%vel%y%tmp(i+1,j) )
			
			if( uh>=0.0d0 )then
				p%of(id)%loc%srk6%x%s1(i,j) = ( -p%of(id)%loc%srk6%x%l1(i-1,j)+6.0d0*p%of(id)%loc%srk6%x%l1(i,j)+3.0d0*p%of(id)%loc%srk6%x%l1(i+1,j) )/8.0d0
			else
				p%of(id)%loc%srk6%x%s1(i,j) = ( -p%of(id)%loc%srk6%x%l1(i+2,j)+6.0d0*p%of(id)%loc%srk6%x%l1(i+1,j)+3.0d0*p%of(id)%loc%srk6%x%l1(i,j) )/8.0d0
			endif
			
			if( vh>=0.0d0 )then
				p%of(id)%loc%srk6%x%s2(i,j) = ( -p%of(id)%loc%srk6%x%l2(i,j-1)+6.0d0*p%of(id)%loc%srk6%x%l2(i,j)+3.0d0*p%of(id)%loc%srk6%x%l2(i,j+1) )/8.0d0
			else
				p%of(id)%loc%srk6%x%s2(i,j) = ( -p%of(id)%loc%srk6%x%l2(i,j+2)+6.0d0*p%of(id)%loc%srk6%x%l2(i,j+1)+3.0d0*p%of(id)%loc%srk6%x%l2(i,j) )/8.0d0
			endif
			
			!-----------------------------------------------------------
			
			uh = 0.5d0*( p%of(id)%loc%vel%x%tmp(i,j) + p%of(id)%loc%vel%x%tmp(i,j+1) )
			vh = 0.5d0*( p%of(id)%loc%vel%y%tmp(i,j) + p%of(id)%loc%vel%y%tmp(i,j+1) )
			
			if( uh>=0.0d0 )then
				p%of(id)%loc%srk6%y%s1(i,j) = ( -p%of(id)%loc%srk6%y%l1(i-1,j)+6.0d0*p%of(id)%loc%srk6%y%l1(i,j)+3.0d0*p%of(id)%loc%srk6%y%l1(i+1,j) )/8.0d0
			else
				p%of(id)%loc%srk6%y%s1(i,j) = ( -p%of(id)%loc%srk6%y%l1(i+2,j)+6.0d0*p%of(id)%loc%srk6%y%l1(i+1,j)+3.0d0*p%of(id)%loc%srk6%y%l1(i,j) )/8.0d0
			endif
			
			if( vh>=0.0d0 )then
				p%of(id)%loc%srk6%y%s2(i,j) = ( -p%of(id)%loc%srk6%y%l2(i,j-1)+6.0d0*p%of(id)%loc%srk6%y%l2(i,j)+3.0d0*p%of(id)%loc%srk6%y%l2(i,j+1) )/8.0d0
			else
				p%of(id)%loc%srk6%y%s2(i,j) = ( -p%of(id)%loc%srk6%y%l2(i,j+2)+6.0d0*p%of(id)%loc%srk6%y%l2(i,j+1)+3.0d0*p%of(id)%loc%srk6%y%l2(i,j) )/8.0d0
			endif
			
			!-----------------------------------------------------------
			
		end do
		end do
		
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			p%of(id)%loc%velsrc%x%tmp(i,j) = - ( p%of(id)%loc%srk6%x%s1(i,j)-p%of(id)%loc%srk6%x%s1(i-1,j) )/p%glb%dx &
											&- ( p%of(id)%loc%srk6%x%s2(i,j)-p%of(id)%loc%srk6%x%s2(i,j-1) )/p%glb%dy
			
			p%of(id)%loc%velsrc%y%tmp(i,j) = - ( p%of(id)%loc%srk6%y%s1(i,j)-p%of(id)%loc%srk6%y%s1(i-1,j) )/p%glb%dx &
											&- ( p%of(id)%loc%srk6%y%s2(i,j)-p%of(id)%loc%srk6%y%s2(i,j-1) )/p%glb%dy
																			
		end do
		end do
	
	!$omp end parallel
	
end subroutine

subroutine ns_diff_source
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: xx,yy,rho,mu
real(8) :: curv,delta,phix,phiy

	!$omp parallel private(id,i,j,xx,yy,rho,mu,curv,delta,phix,phiy), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
	
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			rho   = (p%of(id)%loc%rho%old(i,j)+p%of(id)%loc%rho%old(i+1,j))/2.0_8
			mu    = (p%of(id)%loc%mu%old(i,j)+p%of(id)%loc%mu%old(i+1,j))/2.0_8
			
			curv  = (p%of(id)%loc%normals%curv%old(i,j)+p%of(id)%loc%normals%curv%old(i+1,j))/2.0_8
			delta = (p%of(id)%loc%delta%old(i+1,j)+p%of(id)%loc%delta%old(i,j))/2.0_8
			phix  = (p%of(id)%loc%phi%old(i+1,j)-p%of(id)%loc%phi%old(i,j))/p%glb%dx
			
			xx = (p%of(id)%loc%vel%x%old(i+1,j)-2.0_8*p%of(id)%loc%vel%x%old(i,j)+p%of(id)%loc%vel%x%old(i-1,j))/p%glb%dx**2
			yy = (p%of(id)%loc%vel%x%old(i,j+1)-2.0_8*p%of(id)%loc%vel%x%old(i,j)+p%of(id)%loc%vel%x%old(i,j-1))/p%glb%dy**2

			p%of(id)%loc%velsrc%x%now(i,j) = mu*(xx+yy)/(rho*p%glb%re)
			p%of(id)%loc%velsrc%x%now(i,j) = p%of(id)%loc%velsrc%x%now(i,j) - p%glb%btn_sf*curv*delta*phix/(p%glb%we*rho)
			
			!=============================================================
			
			rho   = (p%of(id)%loc%rho%old(i,j)+p%of(id)%loc%rho%old(i,j+1))/2.0_8		
			mu    = (p%of(id)%loc%mu%old(i,j)+p%of(id)%loc%mu%old(i,j+1))/2.0_8
			
			curv  = (p%of(id)%loc%normals%curv%old(i,j)+p%of(id)%loc%normals%curv%old(i,j+1))/2.0_8
			delta = (p%of(id)%loc%delta%old(i,j)+p%of(id)%loc%delta%old(i,j+1))/2.0_8
			phix  = (p%of(id)%loc%phi%old(i,j+1)-p%of(id)%loc%phi%old(i,j))/p%glb%dy
			
			xx = (p%of(id)%loc%vel%y%old(i+1,j)-2.0_8*p%of(id)%loc%vel%y%old(i,j)+p%of(id)%loc%vel%y%old(i-1,j))/p%glb%dx**2
			yy = (p%of(id)%loc%vel%y%old(i,j+1)-2.0_8*p%of(id)%loc%vel%y%old(i,j)+p%of(id)%loc%vel%y%old(i,j-1))/p%glb%dy**2
			
			p%of(id)%loc%velsrc%y%now(i,j) = mu*(xx+yy)/(rho*p%glb%re)
			p%of(id)%loc%velsrc%y%now(i,j) = p%of(id)%loc%velsrc%y%now(i,j) - p%glb%btn_sf*curv*delta*phiy/(p%glb%we*rho)
			p%of(id)%loc%velsrc%y%now(i,j) = p%of(id)%loc%velsrc%y%now(i,j) - p%glb%btn_g / p%glb%fr
			
			!=============================================================
		
		end do
		end do
		
	!$omp end parallel
	
	call ns_diff_source2

end subroutine

subroutine ns_diff_source2()
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: rho,delta,phix,phiy
real(8) :: ux,uy,vx,vy

	!$omp parallel private(id,i,j,rho,delta,phix,phiy,ux,uy,vx,vy), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
	
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
			
			rho   = (p%of(id)%loc%rho%old(i,j)+p%of(id)%loc%rho%old(i+1,j))/2.0d0
			delta = (p%of(id)%loc%delta%old(i,j)+p%of(id)%loc%delta%old(i+1,j))/2.0d0
			
			phix  = (p%of(id)%loc%phi%old(i+1,j)-p%of(id)%loc%phi%old(i,j))/p%glb%dx
			phiy  = (p%of(id)%loc%phi%old(i+1,j+1)-p%of(id)%loc%phi%old(i+1,j-1)+&
					&p%of(id)%loc%phi%old(i  ,j+1)-p%of(id)%loc%phi%old(i  ,j-1))/4.0d0/p%glb%dy
			
			ux = 0.5d0*( p%of(id)%loc%vel%x%old(i+1,j)-p%of(id)%loc%vel%x%old(i-1,j))/p%glb%dx
			uy = 0.5d0*( p%of(id)%loc%vel%x%old(i,j+1)-p%of(id)%loc%vel%x%old(i,j-1))/p%glb%dy
			
			vx = ( p%of(id)%loc%vel%y%old(i+1,j)+p%of(id)%loc%vel%y%old(i+1,j-1)&
				 &-p%of(id)%loc%vel%y%old(i  ,j)-p%of(id)%loc%vel%y%old(i  ,j-1) )/2.0d0/p%glb%dx
			vy = ( p%of(id)%loc%vel%y%old(i+1,j  )+p%of(id)%loc%vel%y%old(i,j  )&
				 &-p%of(id)%loc%vel%y%old(i+1,j-1)-p%of(id)%loc%vel%y%old(i,j-1) )/2.0d0/p%glb%dy
				  
			p%of(id)%loc%velsrc%x%now(i,j) = p%of(id)%loc%velsrc%x%now(i,j) + &
											&(1.0d0-p%glb%rho_12)*delta*(phix*(2.0d0*ux)+phiy*(uy+vx))/(p%glb%re*rho)
											
			!=============================================================
			
			rho   = (p%of(id)%loc%rho%old(i,j)+p%of(id)%loc%rho%old(i,j+1))/2.0d0
			delta = (p%of(id)%loc%delta%old(i,j)+p%of(id)%loc%delta%old(i,j+1))/2.0d0
			
			phix  = ( p%of(id)%loc%phi%old(i+1,j+1)-p%of(id)%loc%phi%old(i-1,j+1) &
					&+p%of(id)%loc%phi%old(i+1,j  )-p%of(id)%loc%phi%old(i-1,j  ) )/4.0d0/p%glb%dx
			phiy  = (p%of(id)%loc%phi%old(i,j+1)-p%of(id)%loc%phi%old(i,j))/p%glb%dy
			
			ux = ( p%of(id)%loc%vel%x%old(i  ,j+1)+p%of(id)%loc%vel%x%old(i  ,j)&
				  -p%of(id)%loc%vel%x%old(i-1,j+1)-p%of(id)%loc%vel%x%old(i-1,j) )/2.0d0/p%glb%dx
			uy = ( p%of(id)%loc%vel%x%old(i,j+1)+p%of(id)%loc%vel%x%old(i-1,j+1)&
				  -p%of(id)%loc%vel%x%old(i,j  )-p%of(id)%loc%vel%x%old(i-1,j  ) )/2.0d0/p%glb%dy
			
			vx = 0.5d0*(p%of(id)%loc%vel%y%old(i+1,j)-p%of(id)%loc%vel%y%old(i-1,j))/p%glb%dx
			vy = 0.5d0*(p%of(id)%loc%vel%y%old(i,j+1)-p%of(id)%loc%vel%y%old(i,j-1))/p%glb%dy
				  
			p%of(id)%loc%velsrc%y%now(i,j) = p%of(id)%loc%velsrc%y%now(i,j) + &
											(1.0d0-p%glb%rho_12)*delta*(phiy*(2.0d0*vy)+phix*(uy+vx))/(p%glb%re*rho)
											
			!=============================================================
			
		end do
		end do
		
	!$omp end parallel

end subroutine

subroutine ns_predictor
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: src

	!$omp parallel private(id,i,j,src), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			src = 1.5_8*( p%of(id)%loc%velsrc%x%now(i,j)+p%of(id)%loc%velsrc%x%tmp(i,j) ) - 0.5_8*p%of(id)%loc%velsrc%x%old(i,j)			
			p%of(id)%loc%vel%x%now(i,j) = p%of(id)%loc%vel%x%old(i,j) + p%glb%dt * src
			
			src = 1.5_8*( p%of(id)%loc%velsrc%y%now(i,j)+p%of(id)%loc%velsrc%y%tmp(i,j) ) - 0.5_8*p%of(id)%loc%velsrc%y%old(i,j)
			p%of(id)%loc%vel%y%now(i,j) = p%of(id)%loc%vel%y%old(i,j) + p%glb%dt * src
			
		end do
		end do 
		
		call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now)
		
	!$omp end parallel
	
	call pt%vel%sync
	
end subroutine

subroutine ppe_init
use all
!$ use omp_lib
implicit none
integer :: id,i,j

	!$omp parallel private(id,i,j), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			p%of(id)%loc%coe%r(i,j) = 2.0_8 / (p%glb%dx**2*(p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i+1,j)))
			p%of(id)%loc%coe%l(i,j) = 2.0_8 / (p%glb%dx**2*(p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i-1,j)))
			p%of(id)%loc%coe%f(i,j) = 2.0_8 / (p%glb%dy**2*(p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i,j+1)))
			p%of(id)%loc%coe%b(i,j) = 2.0_8 / (p%glb%dy**2*(p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i,j-1)))
			
			p%of(id)%loc%coe%c(i,j) = - ( p%of(id)%loc%coe%r(i,j) + p%of(id)%loc%coe%l(i,j) + &
									&	  p%of(id)%loc%coe%f(i,j) + p%of(id)%loc%coe%b(i,j) )
							
			if( i==1 )then
				p%of(id)%loc%coe%c(i,j)=p%of(id)%loc%coe%c(i,j)+p%of(id)%loc%coe%l(i,j)
				p%of(id)%loc%coe%l(i,j)=0.0_8
			endif
			
			if( i==p%glb%node_x )then
				p%of(id)%loc%coe%c(i,j)=p%of(id)%loc%coe%c(i,j)+p%of(id)%loc%coe%r(i,j)
				p%of(id)%loc%coe%r(i,j)=0.0_8
			endif
			
			if( j==1 )then
				p%of(id)%loc%coe%c(i,j)=p%of(id)%loc%coe%c(i,j)+p%of(id)%loc%coe%b(i,j)
				p%of(id)%loc%coe%b(i,j)=0.0_8
			endif
			
			if( j==p%glb%node_y )then
				p%of(id)%loc%coe%c(i,j)=p%of(id)%loc%coe%c(i,j)+p%of(id)%loc%coe%f(i,j)
				p%of(id)%loc%coe%f(i,j)=0.0_8
			endif
			
		end do
		end do
		
	!$omp end parallel 
	
end subroutine

subroutine ppe_solver
use all
!$ use omp_lib
implicit none
integer :: id,i,j,iter
integer(8) :: cpustart, cpuend
real(8) :: sump, err, w

	call system_clock(cpustart)

	!$omp parallel private(id,i,j), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			p%of(id)%loc%coe%src(i,j) = ( ( p%of(id)%loc%vel%x%now(i,j) - p%of(id)%loc%vel%x%now(i-1,j) ) / p%glb%dx + &
										  ( p%of(id)%loc%vel%y%now(i,j) - p%of(id)%loc%vel%y%now(i,j-1) ) / p%glb%dy ) / p%glb%dt
										  			
		end do
		end do
		
	!$omp end parallel
	
	w = p%glb%p_w1
	
do
	
	p%glb%piter=p%glb%piter+1
	sump=0.0_8

	!$omp parallel private(id,i,j), num_threads(p%glb%threads), reduction(+:sump)
	
		id=0
		!$ id = omp_get_thread_num()
		sump=0.0_8
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			p%of(id)%loc%p%tmp(i,j) = p%of(id)%loc%p%now(i,j)
			
			p%of(id)%loc%p%now(i,j) = p%of(id)%loc%coe%src(i,j) - p%of(id)%loc%coe%r(i,j)*p%of(id)%loc%p%now(i+1,j) &
															&   - p%of(id)%loc%coe%l(i,j)*p%of(id)%loc%p%now(i-1,j) &
															&   - p%of(id)%loc%coe%f(i,j)*p%of(id)%loc%p%now(i,j+1) &
															&   - p%of(id)%loc%coe%b(i,j)*p%of(id)%loc%p%now(i,j-1) 
															
			p%of(id)%loc%p%now(i,j) = p%of(id)%loc%p%now(i,j) / p%of(id)%loc%coe%c(i,j)									
			
			sump = sump + p%of(id)%loc%p%now(i,j)
			
		end do
		end do
	
	!$omp end parallel 
	
	sump=sump/(p%glb%node_x*p%glb%node_y)
	err=0.0_8
	
	!$omp parallel private(id,i,j), num_threads(p%glb%threads), reduction(max:err), shared(sump)
	
		id=0
		!$ id = omp_get_thread_num()
		
		err=0.0_8
		
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
			
			p%of(id)%loc%p%now(i,j) = p%of(id)%loc%p%now(i,j) - sump
			p%of(id)%loc%p%now(i,j) = w * p%of(id)%loc%p%now(i,j) + (1.0_8-w)*p%of(id)%loc%p%tmp(i,j)
			err = max( err, abs(p%of(id)%loc%p%now(i,j)-p%of(id)%loc%p%tmp(i,j)) )
			
		end do
		end do
		
		call p%of(id)%bc(0,p%of(id)%loc%p%now)
		
	!$omp end parallel
	
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
	
end subroutine

subroutine ns_check_convergence
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: linf, l2f, div
real(8) :: px,py,rho

	!$omp parallel private(id,i,j,px,py,rho), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
			
			rho = (p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i+1,j))/2.0_8
			px  = (p%of(id)%loc%p%now(i+1,j)-p%of(id)%loc%p%now(i,j)) / p%glb%dx  
			p%of(id)%loc%vel%x%now(i,j) = p%of(id)%loc%vel%x%now(i,j) - p%glb%dt*px/rho
			
			rho = (p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i,j+1))/2.0_8
			py  = (p%of(id)%loc%p%now(i,j+1)-p%of(id)%loc%p%now(i,j)) / p%glb%dx 		
			p%of(id)%loc%vel%y%now(i,j) = p%of(id)%loc%vel%y%now(i,j) - p%glb%dt*py/rho
			
		end do
		end do
	
		call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now)
		
	!$omp end parallel
	
	call pt%vel%sync
	
	linf=0.0_8
	l2f=0.0_8
	div=0.0_8
	
	!$omp parallel private(id,i,j), num_threads(p%glb%threads), reduction(max:linf,div), reduction(+:l2f)
	
		id=0
		!$ id = omp_get_thread_num()
		div=0.0d0
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			linf = max(linf,abs(p%of(id)%loc%vel%x%now(i,j)-p%of(id)%loc%vel%x%tmp(i,j)))
			linf = max(linf,abs(p%of(id)%loc%vel%y%now(i,j)-p%of(id)%loc%vel%y%tmp(i,j)))
			
			l2f = l2f + (p%of(id)%loc%vel%x%now(i,j)-p%of(id)%loc%vel%x%tmp(i,j))**2
			l2f = l2f + (p%of(id)%loc%vel%y%now(i,j)-p%of(id)%loc%vel%y%tmp(i,j))**2
			
			div = max( div, abs( (p%of(id)%loc%vel%x%now(i,j)-p%of(id)%loc%vel%x%now(i-1,j))/p%glb%dx + &
							   & (p%of(id)%loc%vel%y%now(i,j)-p%of(id)%loc%vel%y%now(i,j-1))/p%glb%dy ) )						
											
		end do
		end do
		
	!$omp end parallel
	
	l2f = dsqrt( l2f / (2.0_8*p%glb%node_x*p%glb%node_y) )
	
	p%glb%vel_div = div
	p%glb%ns_linf = linf
	p%glb%ns_l2f = l2f
	
end subroutine

subroutine ns_final
use all
!$ use omp_lib
implicit none
integer :: id,i,j

	call p%node_vel

	!$omp parallel private(id,i,j), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		call p%of(id)%nvelbc(p%of(id)%loc%nvel%x%now,p%of(id)%loc%nvel%y%now)
		
		do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
		do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
	
			!p%of(id)%loc%p%now(i,j) = p%of(id)%loc%p%now(i,j) !+  p%of(id)%loc%p%old(i,j)
			
			p%of(id)%loc%velsrc%x%now(i,j) = p%of(id)%loc%velsrc%x%now(i,j) + p%of(id)%loc%velsrc%x%tmp(i,j)
			p%of(id)%loc%velsrc%y%now(i,j) = p%of(id)%loc%velsrc%y%now(i,j) + p%of(id)%loc%velsrc%y%tmp(i,j)
			
		end do
		end do
		
		call p%of(id)%bc(0,p%of(id)%loc%p%now)
		
	!$omp end parallel
	
	call pt%nvel%sync
	call pt%p%sync

end subroutine

