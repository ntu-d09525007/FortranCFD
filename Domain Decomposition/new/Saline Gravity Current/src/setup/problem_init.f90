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
    !call plt%init(p)
    !write(*,*)"plotter init finish"
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
        
            p%of(id)%loc%solid%now(i,j,k) = 0.0d0
        
            do ii = 1, ug
            do jj = 1, ug
            do kk = 1, ug
                
                x = 0.5d0*( p%glb%x(i)+p%glb%x(i-1) ) + real(ii,8)*p%glb%dx/real(ug,8)
                y = 0.5d0*( p%glb%y(j)+p%glb%y(j-1) ) + real(jj,8)*p%glb%dy/real(ug,8)
                z = 0.5d0*( p%glb%z(k)+p%glb%z(k-1) ) + real(kk,8)*p%glb%dz/real(ug,8)         
                
                if( z < p%glb%zend - x*p%glb%slope )then
                    p%of(id)%loc%solid%now(i,j,k) = p%of(id)%loc%solid%now(i,j,k) + 1.0d0/real(ug,8)**3.0d0
                endif
                
            end do
            end do
            end do
            
            p%of(id)%loc%phi%now(i,j,k) = 2.0d0*p%of(id)%loc%solid%now(i,j,k)-1.0d0
            
            p%of(id)%loc%vel%x%now(i,j,k) = 0.0_8
            p%of(id)%loc%vel%y%now(i,j,k) = 0.0_8
            p%of(id)%loc%vel%z%now(i,j,k) = 0.0_8
            
        end do
        end do
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%phi%now)
        call p%of(id)%bc(0,p%of(id)%loc%solid%now)
    enddo
    !$omp end parallel do
    
    call pt%phi%sync
    call pt%vel%sync
    call pt%solid%sync
    
    !====================== Smooth Solid ==============================
    
    call level_set_rk3_redis(0)
    call p%ls_funs
    
    !$omp parallel do private(i,j,k)
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks-p%glb%ghc, p%of(id)%loc%ke+p%glb%ghc
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
            p%of(id)%loc%solid%now(i,j,k) = p%of(id)%loc%heavy%now(i,j,k) 
        enddo
        enddo
        enddo
    
        call p%of(id)%bc(0,p%of(id)%loc%solid%now)
        
    enddo
    
    call pt%solid%sync
    
    !====================== Smooth C ==============================
    
    !$omp parallel do private(i,j,k,x,y,z,heavy)
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            x = p%glb%x(i)
            y = p%glb%y(j)
            z = p%glb%z(k)
            
            ! if( p%of(id)%loc%solid%now(i,j,k) < 0.5d0 .and. x>p%glb%h0/p%glb%H/p%glb%slope )then
                ! p%of(id)%loc%s%now(i,j,k) = 1.0d0-z
            ! else
                ! p%of(id)%loc%s%now(i,j,k) = 0.0d0
            ! endif
            
            ! if( p%of(id)%loc%solid%now(i,j,k) < 0.5d0 .and. x<p%glb%h0/p%glb%H/p%glb%slope)then
                ! p%of(id)%loc%c%now(i,j,k) = 1.0d0
            ! else
                ! p%of(id)%loc%c%now(i,j,k) = -1.0d0
            ! endif
            
            !p%of(id)%loc%phi%now(i,j,k) = p%of(id)%loc%c%now(i,j,k) 
            
            p%of(id)%loc%s%now(i,j,k) = 1.0d0-z
            p%of(id)%loc%c%now(i,j,k) = (p%glb%h0/p%glb%H/p%glb%slope - x) / p%glb%ls_wid
            
            if( p%of(id)%loc%c%now(i,j,k)>1.0d0 ) then
                p%of(id)%loc%c%now(i,j,k) = 1.0d0
            else if ( p%of(id)%loc%c%now(i,j,k) < -1.0d0 )then
                p%of(id)%loc%c%now(i,j,k) = 0.0d0
            else
                p%of(id)%loc%c%now(i,j,k) = 0.5d0 * (1.0_8 + p%of(id)%loc%c%now(i,j,k) + dsin(pi*p%of(id)%loc%c%now(i,j,k)) / pi )
            endif
            
        end do
        end do
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%phi%now)
        call p%of(id)%bc(0,p%of(id)%loc%c%now)
        call p%of(id)%bc(0,p%of(id)%loc%s%now)

    enddo
    !$omp end parallel do
    
    ! write(*,*)"start redistancing"
    ! call level_set_rk3_redis(0)
    
    ! call p%ls_funs
    
    ! !$omp parallel do private(i,j,jj,k,x,y,z,heavy)
    ! do id = 0, p%glb%threads-1
        
        ! do k = p%of(id)%loc%ks-p%glb%ghc, p%of(id)%loc%ke+p%glb%ghc
        ! do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        ! do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
            ! p%of(id)%loc%c%now(i,j,k) = p%of(id)%loc%heavy%now(i,j,k) 
        ! enddo
        ! enddo
        ! enddo
        
        ! call p%of(id)%bc(0,p%of(id)%loc%c%now)
        
        ! j=p%of(id)%loc%js
        
        ! do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        ! do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            ! x = p%glb%x(i)
            ! z = p%glb%z(k)
            
            ! heavy = 1.0d0 - z
            ! y = (x - heavy / p%glb%slope)/p%glb%ls_wid
            
            ! if( z<p%glb%zend-p%glb%h0/p%glb%H )then
                ! if( abs(y) < 1.0d0 )then
                    ! p%of(id)%loc%s%now(i,j,k) = 0.5d0 * heavy * (1.0_8 + y + dsin(pi*y) / pi )
                ! endif
            ! else
                ! y = ( x - p%glb%h0/p%glb%H/p%glb%slope)/p%glb%ls_wid
                ! if( abs(y) < 1.0d0 )then
                    ! p%of(id)%loc%s%now(i,j,k) = 0.5d0 * heavy * (1.0_8 + y + dsin(pi*y) / pi )
                ! endif
            ! endif
            
            ! do jj = p%of(id)%loc%js+1, p%of(id)%loc%je
                ! p%of(id)%loc%s%now(i,jj,k) = p%of(id)%loc%s%now(i,j,k) 
            ! enddo
            
        ! enddo
        ! enddo
        
        ! call p%of(id)%bc(0,p%of(id)%loc%s%now)
        
    ! enddo
    ! !$omp end parallel do
    
    call pt%s%sync
    call pt%c%sync 

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

    call p%ls_mv
    call p%total_c
    p%glb%ivol = p%glb%vol
    p%glb%imass = p%glb%mass
    p%glb%ivolv = p%glb%volv
    p%glb%imassv = p%glb%massv
    p%glb%icsum = p%glb%csum
    p%glb%issum = p%glb%ssum
    call p%sync

    write(*,*) "plotting"
    call plot
	call p%plot

    p%glb%ls_adv = 0.0d0
    p%glb%ls_red = 0.0d0
    p%glb%ppe    = 0.0d0
    p%glb%ns     = 0.0d0
    p%glb%syn    = 0.0d0

        
end subroutine
