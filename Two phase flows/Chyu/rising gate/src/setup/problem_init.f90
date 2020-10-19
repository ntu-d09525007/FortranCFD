subroutine problem_init()
use all
!$ use omp_lib
implicit none
integer :: id, i, j, k, ug, ii,jj,kk
real(8) :: x, y, z, err
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
        
    call p%show
    !---------------------------------------------------
    
    id=0

    ug=30
    !$omp parallel do collapse(3), private(ii,jj,kk,x,y,z)
    do k = p%loc%ks, p%loc%ke
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
    
        p%loc%vof%now(i,j,k) = 0.0d0
        p%loc%ibm%solid%now(i,j,k) = 0.0d0
        p%loc%vel%x%now(i,j,k) = 0.0d0
        p%loc%vel%y%now(i,j,k) = 0.0d0
        p%loc%vel%z%now(i,j,k) = 0.0d0
    
        do ii = 1, ug
        do jj = 1, ug
        do kk = 1, ug
            
            x = 0.5d0*( p%glb%x(i,j,k)+p%glb%x(i-1,j,k) ) + real(ii,8)*p%glb%dx/real(ug,8)
            y = 0.5d0*( p%glb%y(i,j,k)+p%glb%y(i,j-1,k) ) + real(jj,8)*p%glb%dy/real(ug,8)
            z = 0.5d0*( p%glb%z(i,j,k)+p%glb%z(i,j,k-1) ) + real(kk,8)*p%glb%dz/real(ug,8)
            
            ! dambreak -- partial failure
            !=========================================
            ! if( x>1.0d0 .and. x<1.0d0+2.0d0*p%glb%dx .and. abs(y-1.0d0)>0.2d0 )then
            !     p%loc%ibm%solid%now(i,j,k) = p%loc%ibm%solid%now(i,j,k) + 1.0d0/real(ug,8)**3.0d0
            ! end if
            
            ! dambreak -- rising gate
            !========================================
            if( x>5.0d0/3.0d0 .and. x<5.0d0/3.0d0+2.0d0*p%glb%dx .and. z<4.0d0/3.0d0 )then
                p%loc%ibm%solid%now(i,j,k) = p%loc%ibm%solid%now(i,j,k) + 1.0d0/real(ug,8)**3.0d0
            endif

        end do
        end do
        end do
    
        x = p%glb%x(i,j,k)
        y = p%glb%y(i,j,k)
        z = p%glb%z(i,j,k)
        
        ! vortex
        !=========================================
        !p%loc%phi%now(i,j,k) = - dsqrt( (x-0.35d0)**2.0d0 + (y-0.35d0)**2.0d0 + (z-0.35d0)**2.0d0 ) + 0.15d0
        
        ! dambreak -- drybed
        !=========================================
        if( x<=5.0d0/3.0d0 .and. z<=1.0d0 )then
            p%loc%phi%now(i,j,k) = 1.0_8
        else
            p%loc%phi%now(i,j,k) = -1.0_8
        end if

        ! dambreak -- wetbed
        !=========================================
        ! if( x<=2.53d0 .and. z<=1.0 )then
        !     p%loc%phi%now(i,j,k) = 1.0_8
        ! else if ( x>2.53d0 .and. z<=0.12d0 )then
        !     p%loc%phi%now(i,j,k) = 1.0_8
        ! else 
        !     p%loc%phi%now(i,j,k) = -1.0_8
        ! end if

        ! dambreak -- partial failure
        !=========================================
        ! if( x<1.0d0 .and. z<0.6d0 )then
        !     p%loc%phi%now(i,j,k) = 1.0_8
        ! else 
        !     p%loc%phi%now(i,j,k) = -1.0_8
        ! endif
        
    end do
    end do
    end do
    !$omp end parallel do

    p%loc%ibm%z = 0.0d0
    p%loc%ibm%w = 0.408d0 / p%glb%U
    
    call bc(p%loc%phi%now)
    call bc(p%loc%vof%now)
    call bc(p%loc%ibm%solid%now)
    call velbc(p%loc%vel%x%now,p%loc%vel%y%now,p%loc%vel%z%now)
    
    write(*,*)"Init data finish"

    write(*,*)"start redistancing"
    call level_set_rk3_redis(0)
    
    write(*,*)"find cell center velocity"
    call node_vel
    call ns_init
    call p%loc%init
    !--------------------------------------------------- 

    write(*,*)"ns ab setup"
    call ns_ab_setup

    call ls_mv
    p%glb%ivol = p%glb%vol
    p%glb%imass = p%glb%mass
    p%glb%ivolv = p%glb%volv
    p%glb%imassv = p%glb%massv


    write(*,*) "plotting"
    call plot
    call p%plot

    p%glb%ls_adv = 0.0d0
    p%glb%ls_red = 0.0d0
    p%glb%ppe    = 0.0d0
    p%glb%ns     = 0.0d0
    p%glb%syn    = 0.0d0
        
end subroutine
