module tree
!$ use omp_lib
use branches
implicit none

type filemanager
integer :: ls_mv
end type filemanager

type manager
type(global) :: glb
type(filemanager) :: fil
type(job),allocatable :: of(:)
contains
procedure show => manager_show
procedure read => manager_read
procedure init => manager_init
procedure sync => manager_sync
procedure ls_mv => manager_ls_mv
procedure ls_funs => manager_ls_funs
procedure rho_mu => manager_rho_mu
procedure surface_norms => manager_surface_norms
procedure node_vel => manager_node_vel
procedure curv => manager_curv
procedure switch => manager_switch
end type manager

contains 

subroutine manager_show(p)
implicit none
class(manager) :: p
integer :: id

    write(*,*)" --- Problem information --- "
    write(*,'("Number of computing threads:",I5)')p%glb%threads
    write(*,'("Number of ploting threads:",I5)')p%glb%pthreads
    write(*,'("dx:",ES15.4)')p%glb%dx
    write(*,'("dy:",ES15.4)')p%glb%dy
    write(*,'("dt:",ES15.4)')p%glb%dt
    write(*,'("Re:",F8.3)')p%glb%re
    write(*,'("We:",F8.3)')p%glb%we
    write(*,'("Fr:",F8.3)')p%glb%fr
    write(*,'("Grids",I5,"x",I5)')p%glb%node_x,p%glb%node_y
    do id = 0, p%glb%threads-1
        write(*,'(I3,4I5)')id,p%of(id)%loc%is,p%of(id)%loc%ie,p%of(id)%loc%js,p%of(id)%loc%je
    end do

end subroutine

subroutine manager_read(p,path)
implicit none
class(manager) :: p
character(*) :: path

 open(unit=526,file=trim(path),status='old')
 
 read(526,*)
 read(526,*)p%glb%name
 read(526,*)
 read(526,*)p%glb%threads
 read(526,*)
 read(526,*)p%glb%pthreads
 read(526,*)
 read(526,*)p%glb%ug
 read(526,*)
 read(526,*)p%glb%ghc
 read(526,*)
 read(526,*)p%glb%xstart, p%glb%xend
 read(526,*)
 read(526,*)p%glb%ystart, p%glb%yend
 read(526,*)
 read(526,*)p%glb%t2s, p%glb%t2p
 read(526,*)
 read(526,*)p%glb%dt, p%glb%rdt
 read(526,*)
 read(526,*)p%glb%p_tol, p%glb%p_w1, p%glb%p_w2, p%glb%p_b
 read(526,*)
 read(526,*)p%glb%srk6_tol, p%glb%srk6_w
 read(526,*)
 read(526,*)p%glb%ls_wid
 read(526,*)
 read(526,*)p%glb%how_to_paras
 read(526,*)
 read(526,*)p%glb%rho_1, p%glb%mu_1
 read(526,*)
 read(526,*)p%glb%rho_2, p%glb%mu_2
 read(526,*)
 read(526,*)p%glb%sigma, p%glb%btn_sf
 read(526,*)
 read(526,*)p%glb%g, p%glb%btn_g
 read(526,*)
 read(526,*)p%glb%L, p%glb%U, p%glb%T
 read(526,*)
 read(526,*)p%glb%Re, p%glb%Fr, p%glb%We
 read(526,*)
 read(526,*)p%glb%rho_12, p%glb%mu_12
 read(526,*)
 read(526,*)p%glb%ubc(1), p%glb%ubc(2)
 read(526,*)
 read(526,*)p%glb%vbc(1), p%glb%vbc(2)
 
 close(unit=526)
 
end subroutine

subroutine manager_init(p,path)
implicit none
class(manager) :: p
character(*) :: path
integer :: max_threads
integer :: i, j, id

 	call p%read(path)
	
	p%fil%ls_mv = 15
	open(unit=p%fil%ls_mv,file="./out/"//trim(p%glb%name)//"_MVloss.plt")
	write(p%fil%ls_mv,*)'variables = "T" "Loss of mass" "Loss of Volume" '
	
	max_threads = 1
	!$ max_threads = omp_get_max_threads()

  	if( p%glb%threads > max_threads )then
		write(*,*)"Warning >> number of CPU threads are not the same as you requested"
		p%glb%threads = max_threads
	endif
	
	if( p%glb%pthreads > p%glb%threads )then
		write(*,*)" Warning >> number of pieces of single plot are too many. "
		p%glb%pthreads = p%glb%threads
	endif
	
	allocate( p%of(0:p%glb%threads-1) )
	
	p%glb%node_x = p%glb%ug * ( p%glb%xend - p%glb%xstart )
	p%glb%node_y = p%glb%ug * ( p%glb%yend - p%glb%ystart )
	
	allocate( p%glb%x(0:p%glb%node_x+1), p%glb%y(0:p%glb%node_y+1) )
	
	!$omp parallel private(id), num_threads(p%glb%threads)
 		id = 0
		!$ id = omp_get_thread_num()
		allocate( p%of(id)%glb%x(0:p%glb%node_x+1), p%of(id)%glb%y(0:p%glb%node_y+1) )
	!$omp end parallel
		
	p%glb%dx = ( p%glb%xend - p%glb%xstart ) / p%glb%node_x
	p%glb%dy = ( p%glb%yend - p%glb%ystart ) / p%glb%node_y
		
	!$omp parallel do
	do i = 1, p%glb%node_x
		p%glb%x(i) = p%glb%xstart + (i-0.5)*p%glb%dx
	enddo
	!$omp end parallel do
	
	!$omp parallel do
	do j = 1, p%glb%node_y
		p%glb%y(j) = p%glb%ystart + (j-0.5)*p%glb%dy
	enddo
	!$omp end parallel do
	
	p%glb%x(0)=p%glb%xstart; p%glb%x(p%glb%node_x+1)=p%glb%xend
	p%glb%y(0)=p%glb%ystart; p%glb%y(p%glb%node_y+1)=p%glb%yend
	
	p%glb%dt = p%glb%dt * p%glb%dx
	p%glb%rdt = p%glb%rdt * p%glb%dx
	p%glb%ls_wid = p%glb%ls_wid * p%glb%dx
	
	p%glb%p_tol = 0.1_8 ** p%glb%p_tol
	p%glb%srk6_tol = 0.1_8 ** p%glb%srk6_tol

	p%glb%p_b = 0.1_8 ** p%glb%p_b

	
	select case ( p%glb%how_to_paras )
	
		case (1)
			continue
		case (2) ! L
			p%glb%U = dsqrt( p%glb%L * p%glb%g )
			p%glb%T = p%glb%L / p%glb%U
		case (3) ! L+U
			p%glb%T = p%glb%L / p%glb%U
		case (4) ! L+T
			p%glb%U = p%glb%L / p%glb%T
		case (5) ! U+T
			p%glb%L = p%glb%U * p%glb%T
		case default
			write(*,*)"Error >> Wrong parameter selector "
			stop
			
 	end select 
	
	if( p%glb%how_to_paras > 1 )then
		p%glb%mu_12 = p%glb%mu_2 / p%glb%mu_1
		p%glb%rho_12 = p%glb%rho_2 / p%glb%rho_1
		p%glb%re = p%glb%rho_1 * p%glb%u * p%glb%L / p%glb%mu_1
		p%glb%we = p%glb%rho_1 * p%glb%u**2 * p%glb%L / p%glb%sigma
		p%glb%fr = p%glb%u**2 / ( p%glb%g * p%glb%L ) 
	endif
	
	call p%sync
	
	!$omp parallel private(id), num_threads(p%glb%threads)
		id=0
		!$ id = omp_get_thread_num()
		call p%of(id)%init(id)
	!$omp end parallel
	
	p%glb%time = 0.0_8
	p%glb%iter = 0
	p%glb%pid = 0
	
	p%glb%ls_adv = 0.0d0
	p%glb%ls_red = 0.0d0
	p%glb%ppe    = 0.0d0
	p%glb%ns     = 0.0d0
	
	call system_clock( count_rate=p%glb%cpurate )

end subroutine

subroutine manager_sync(p)
implicit none
class(manager) :: p
integer :: i, id

 !$omp parallel private(id), num_threads(p%glb%threads)
 	id = 0
	!$ id = omp_get_thread_num()
	p%of(id)%glb = p%glb
 !$omp end parallel
 
end subroutine

subroutine manager_ls_funs(p)
implicit none
class(manager) :: p
integer :: id, i, j
real(8) :: x, heavy, hp, pi, eps


    eps = 1.0d-12
    pi = dacos(-1.0_8)

    !$omp parallel private(id,i,j,x,heavy,hp), num_threads(p%glb%threads)
    
        id=0
        !$ id = omp_get_thread_num()
        
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        
            x = p%of(id)%loc%phi%now(i,j) / p%glb%ls_wid
            
            if( x > 1.0_8-eps )then
                heavy = 1.0_8
                hp = 0.0_8
            else if ( x < -1.0_8+eps )then
                heavy = 0.0_8
                hp = 0.0_8
            else
                heavy = 0.5_8 * (1.0_8 + x + dsin(pi*x) / pi )
                hp = 0.5_8 * ( 1.0_8 + dcos(pi*x) ) / p%glb%ls_wid 
            endif
            
            p%of(id)%loc%heavy%now(i,j) = heavy 
            p%of(id)%loc%delta%now(i,j) = hp
            p%of(id)%loc%sign%now(i,j) = 2.0_8*heavy-1.0_8
            
        end do
        end do
    
    !$omp end parallel

end subroutine

subroutine manager_rho_mu(p)
implicit none
class(manager) :: p
integer :: id,i,j

    call p%ls_funs
    
    !$omp parallel private(id,i,j), num_threads(p%glb%threads)
    
        id=0
        !$ id = omp_get_thread_num()
        
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
            
            p%of(id)%loc%rho%now(i,j) = p%of(id)%loc%heavy%now(i,j) + p%glb%rho_12 * (1.0_8 - p%of(id)%loc%heavy%now(i,j) )
            p%of(id)%loc%mu%now(i,j) = p%of(id)%loc%heavy%now(i,j) + p%glb%mu_12 * (1.0_8 - p%of(id)%loc%heavy%now(i,j) )
			
			!p%of(id)%loc%rho%now(i,j) = p%glb%rho_12 * p%of(id)%loc%heavy%now(i,j) +  (1.0_8 - p%of(id)%loc%heavy%now(i,j) )
            !p%of(id)%loc%mu%now(i,j) = p%glb%mu_12 * p%of(id)%loc%heavy%now(i,j) +  (1.0_8 - p%of(id)%loc%heavy%now(i,j) )
            
        end do
        end do
        
    !$omp end parallel
    
    
end subroutine

subroutine manager_ls_mv(p)
implicit none
class(manager) :: p
integer :: id, i, j
real(8) :: mass, vol

    call p%rho_mu

    mass = 0.0_8; vol=0.0_8
    
    !$omp parallel private(id,i,j), num_threads(p%glb%threads), reduction(+:mass,vol)
        
        id=0
        !$ id = omp_get_thread_num()
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            
            mass = mass + p%of(id)%loc%rho%now(i,j)*p%of(id)%loc%heavy%now(i,j)*p%glb%dx*p%glb%dy
            vol = vol + p%of(id)%loc%heavy%now(i,j)*p%glb%dx*p%glb%dy
        
        enddo
        enddo
    
    !$omp end parallel
    
    p%glb%mass = mass
    p%glb%vol = vol
    
    call p%sync


end subroutine

subroutine manager_surface_norms(p)
implicit none
class(manager) :: p
integer :: id,I,J   
    
    !$omp parallel private(id,i,j), num_threads(p%glb%threads)
    
        id=0
        !$ id = omp_get_thread_num()
        
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        
            call p%of(id)%loc%ccd%x%solve_fixed_central(15.0_8/16.0_8, p%of(id)%loc%phi%now(:,j),&
                                                                      &p%of(id)%loc%normals%x%now(:,j),&
                                                                      &p%of(id)%loc%normals%xx%now(:,j) )
                                                                    
        enddo
        
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        
            call p%of(id)%loc%ccd%y%solve_fixed_central(15.0_8/16.0_8, p%of(id)%loc%phi%now(i,:),&
                                                                      &p%of(id)%loc%normals%y%now(i,:),&
                                                                      &p%of(id)%loc%normals%yy%now(i,:) )
                                                                    
            call p%of(id)%loc%ccd%y%solve_fixed_central(15.0_8/16.0_8, p%of(id)%loc%normals%x%now(i,:),&
                                                                      &p%of(id)%loc%normals%xy%now(i,:))
                                                                    
        enddo
        
        call p%of(id)%bc(0,p%of(id)%loc%normals%x%now)
        call p%of(id)%bc(0,p%of(id)%loc%normals%xx%now)
        call p%of(id)%bc(0,p%of(id)%loc%normals%y%now)
        call p%of(id)%bc(0,p%of(id)%loc%normals%yy%now)
        call p%of(id)%bc(0,p%of(id)%loc%normals%xy%now)
        
    !$omp end parallel
    
end subroutine

subroutine manager_curv(p)
implicit none
class(manager) :: p
integer :: id,i,j
real(8) :: fx,fxx,fy,fyy,fxy

	!$omp parallel private(id,i,j,fx,fxx,fy,fyy,fxy), num_threads(p%glb%threads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
		
			fx  = p%of(id)%loc%normals%x%now(i,j)
			fy  = p%of(id)%loc%normals%y%now(i,j)
			fxx = p%of(id)%loc%normals%xx%now(i,j)
			fyy = p%of(id)%loc%normals%yy%now(i,j)
			fxy = p%of(id)%loc%normals%xy%now(i,j)
			
			p%of(id)%loc%normals%curv%now(i,j) = (fx*fx*fyy-2.0_8*fx*fy*fxy+fy*fy*fxx) / (fx*fx+fy*fy+1.0d-10)**1.5_8
			
		end do
		end do
		
		call p%of(id)%bc(0,p%of(id)%loc%normals%curv%now)
	
	!$omp end parallel 

end subroutine

subroutine manager_node_vel(p)
implicit none
class(manager) :: p
integer :: id,i,j

    !$omp parallel private(id,i,j), num_threads(p%glb%threads)
    
        id=0
        !$ id = omp_get_thread_num()
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            p%of(id)%loc%nvel%x%now(i,j) = 0.5_8 * ( p%of(id)%loc%vel%x%now(i-1,j) + p%of(id)%loc%vel%x%now(i,j) )
            p%of(id)%loc%nvel%y%now(i,j) = 0.5_8 * ( p%of(id)%loc%vel%y%now(i,j-1) + p%of(id)%loc%vel%y%now(i,j) )
        end do
        end do
    
        call p%of(id)%nvelbc(p%of(id)%loc%nvel%x%now,p%of(id)%loc%nvel%y%now)
		
    !$omp end parallel

end subroutine

subroutine manager_switch(p)
implicit none
class(manager) :: p
integer :: id

    !$omp parallel private(id), num_threads(p%glb%threads)
    
        id=0
        !$ id = omp_get_thread_num()
        
        call p%of(id)%loc%phi%switch
        
        call p%of(id)%loc%rho%switch
		call p%of(id)%loc%mu%switch
		
		call p%of(id)%loc%delta%switch
		call p%of(id)%loc%heavy%switch
		call p%of(id)%loc%sign%switch
		
		call p%of(id)%loc%normals%switch
		
		call p%of(id)%loc%vel%switch
		call p%of(id)%loc%nvel%switch
		call p%of(id)%loc%velsrc%switch
		
		call p%of(id)%loc%p%switch
		    
    !$omp end parallel

end subroutine

end module tree