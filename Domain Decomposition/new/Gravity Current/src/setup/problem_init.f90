subroutine problem_init()
use all
!$ use omp_lib
implicit none
integer :: id, i, j, k, ug, ii,jj,kk
real(8) :: x, y, z, err, heavy, pi
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
    write(*,*)"data init finish"
    call pt%init(p)
    write(*,*)"pointer init finish"
        
    call p%show
    !---------------------------------------------------
    
    ug=30
    pi=dacos(-1.0d0)
    !$omp parallel do private(i,j,k,ii,jj,kk,x,y,z,heavy)
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            p%of(id)%loc%vof%now(i,j,k) = 0.0d0
        
            do ii = 1, ug
            do jj = 1, ug
            do kk = 1, ug
                
                x = 0.5d0*( p%glb%x(i,j,k)+p%glb%x(i-1,j,k) ) + real(ii,8)*p%glb%dx/real(ug,8)
                y = 0.5d0*( p%glb%y(i,j,k)+p%glb%y(i,j-1,k) ) + real(jj,8)*p%glb%dy/real(ug,8)
                z = 0.5d0*( p%glb%z(i,j,k)+p%glb%z(i,j,k-1) ) + real(kk,8)*p%glb%dz/real(ug,8)
                
                if( - dsqrt( (x-0.35d0)**2.0d0 + (y-0.35d0)**2.0d0 + (z-0.35d0)**2.0d0 ) + 0.15d0 >= 0.0d0 )then
                    p%of(id)%loc%vof%now(i,j,k) = p%of(id)%loc%vof%now(i,j,k) + 1.0d0/real(ug,8)**3.0d0
                end if
                
            end do
            end do
            end do
        
            x = p%glb%x(i,j,k)
            y = p%glb%y(i,j,k)
            z = p%glb%z(i,j,k)
            
            ! vortex
            !=========================================
            !p%of(id)%loc%phi%now(i,j,k) = - dsqrt( (x-0.5d0)**2.0d0 + (y-0.5d0)**2.0d0 + (z-0.5d0)**2.0d0 ) + 0.15d0
            
            ! dambreak
            !=========================================
            !if( x<=1.0d0 .and. y<=1.0d0 .and. z<=1.0d0 )then
            !    p%of(id)%loc%phi%now(i,j,k) = 1.0_8
            !else
            !    p%of(id)%loc%phi%now(i,j,k) = -1.0_8
            !end if
            
            ! bubble burst
            !=========================================
            ! if( z<-2.0d0 )then
                ! if( x**2.0d0 + y**2.0d0 + (z+3.2d0)**2.0d0 < 1.0d0 )then
                    ! p%of(id)%loc%phi%now(i,j,k) = dsqrt( x**2.0d0 + y**2.0d0 + (z+3.2d0)**2.0d0 ) - 1.0d0
                ! else
                    ! p%of(id)%loc%phi%now(i,j,k) = min( dsqrt( x**2.0d0 + y**2.0d0 + (z+3.2d0)**2.0d0 ) - 1.0d0, -2.0-z )
                ! endif
            ! else
                ! p%of(id)%loc%phi%now(i,j,k) = -2.0d0-z
            ! endif
            
            ! Milk Crown
            !=========================================
            !if( z <= 0.1876d0 )then    
            !   p%of(id)%loc%phi%now(i,j,k) = -z + 0.1876d0
            !else if( dsqrt(x**2.0d0+y**2.0d0+(z-0.8d0)**2.0d0)<=0.5d0 )then
            !   p%of(id)%loc%phi%now(i,j,k) = -dsqrt(x**2.0d0+y**2.0d0+(z-0.8d0)**2.0d0)+0.5d0
            !else 
            !   p%of(id)%loc%phi%now(i,j,k)= MAX(-z+0.1876d0,-dsqrt(x**2.0d0+y**2.0d0+(z-0.8d0)**2.0d0)+0.5d0 )
            !end if
            
            ! Gravity Current
            !========================================
            p%of(id)%loc%c%now(i,j,k) = (-x+2.0d0) / p%glb%ls_wid 
            
            if( p%of(id)%loc%c%now(i,j,k) > 1.0d0 - 1.0d-12 )then
                p%of(id)%loc%c%now(i,j,k) = 1.0d0
            else if( p%of(id)%loc%c%now(i,j,k) < -1.0d0 + 1.0d-12 )then
                p%of(id)%loc%c%now(i,j,k)  = 0.0d0
            else
                p%of(id)%loc%c%now(i,j,k)  = 0.5d0 * (1.0d0 + p%of(id)%loc%c%now(i,j,k) + dsin(pi*p%of(id)%loc%c%now(i,j,k)) / pi )
            endif
            
            !x = 0.5d0*( p%glb%x(i)+p%glb%x(i+1) )
            !y = 0.5d0*( p%glb%y(j)+p%glb%y(j+1) )
            !x = p%glb%x(i)
            !y = p%glb%y(j)
            !z = 0.5d0*( p%glb%z(k)+p%glb%z(k+1) )
            
            !if( dsqrt(x**2.0d0+y**2.0d0+(z-0.8d0)**2.0d0)<=0.5d0 )then
            !   p%of(id)%loc%vel%z%now(i,j,k) = -1.0d0
            !endif 
            
            p%of(id)%loc%vel%x%now(i,j,k) = 0.0_8
            p%of(id)%loc%vel%y%now(i,j,k) = 0.0_8
            p%of(id)%loc%vel%z%now(i,j,k) = 0.0_8
            
        end do
        end do
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%phi%now)
        call p%of(id)%bc(0,p%of(id)%loc%vof%now)
        call p%of(id)%bc(0,p%of(id)%loc%c%now)
        call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now,p%of(id)%loc%vel%z%now)

    enddo
    !$omp end parallel do
    
    write(*,*)"Init data finish"
    
    call pt%phi%sync;write(*,*)"phi sync"
    call pt%vel%sync;write(*,*)"vel sync"
    call pt%vof%sync;write(*,*)"vof sync"
    call pt%c%sync;write(*,*)"c sync"
    write(*,*)"data sync finish"
    
    !write(*,*)"start redistancing"
    !call level_set_rk3_redis(0)
    
    write(*,*)"find cell center velocity"
    call p%node_vel
    call pt%nvel%sync
    call ns_init
    
    !---------------------------------------------------

    !$omp parallel do 
    do id = 0, p%glb%threads-1
        call p%of(id)%loc%init  
    enddo
    !$omp end parallel do

    write(*,*)"ns ab setup"
    call ns_ab_setup

    call p%ls_mv(.true.)
    call p%total_c(.true.)
    call p%total_e(.true.)

    write(*,*) "plotting"
    call plot
    call p%plot

    p%glb%ls_adv = 0.0d0
    p%glb%ls_red = 0.0d0
    p%glb%ppe    = 0.0d0
    p%glb%ns     = 0.0d0
    p%glb%syn    = 0.0d0
    
    write(*,*)"rho difference:",(( p%glb%re * p%glb%mu_1 ) / ( p%glb%rho_1 * p%glb%L ))**2.0d0 / (p%glb%g*p%glb%L)
 
end subroutine
