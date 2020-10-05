subroutine problem_init()
use all
!$ use omp_lib
implicit none
integer :: i, j, ug, ii,jj
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
    write(*,*)"data init finish"
        
    call p%show
    !---------------------------------------------------

    ug=30
    !$omp parallel do collapse(2), private(ii,jj,x,y)
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
    
        p%loc%vof%now(i,j) = 0.0d0
    
        ! do ii = 1, ug
        ! do jj = 1, ug
            
        !     x = 0.5d0*( p%glb%x(i,j)+p%glb%x(i-1,j) ) + real(ii,8)*p%glb%dx/real(ug,8)
        !     y = 0.5d0*( p%glb%y(i,j)+p%glb%y(i,j-1) ) + real(jj,8)*p%glb%dy/real(ug,8)
            
        ! end do
        ! end do
    
        x = p%glb%x(i,j)
        y = p%glb%y(i,j)

        ! dambreak -- wetbed
        !=========================================
        if( x<=1.0d0 .and. y<=1.0 )then
            p%loc%phi%now(i,j) = 1.0_8
        else 
            p%loc%phi%now(i,j) = -1.0_8
        end if
        
        p%loc%vel%x%now(i,j) = 0.0_8
        p%loc%vel%y%now(i,j) = 0.0_8
        
    end do
    end do
    !$omp end parallel do
    
    call bc(p%loc%phi%now)
    call bc(p%loc%vof%now)
    call velbc(p%loc%vel%x%now,p%loc%vel%y%now)
    
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
