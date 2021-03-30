subroutine problem_init()
use all
!$ use omp_lib
implicit none
integer :: id, i, j, k, ug, ii,jj,kk
real(8) :: x, y, z, err
CHARACTER(100) :: NAME_OF_FILE
real(8) :: kx, kz, kh
    
    NAME_OF_FILE="default.txt"
    
    ! WRITE(*,*)"============================================"
    ! WRITE(*,*)'Returning the files in directory "/input" '
    ! WRITE(*,*)
    ! CALL SYSTEM("ls ./input")
    ! WRITE(*,*)
    ! WRITE(*,*)"============================================"
    ! WRITE(*,*)'Selet an input file (<d> : default.txt), with full name'
    ! WRITE(*,*)
    ! READ(*,*)NAME_OF_FILE
    !if( name_of_file == 'd' )name_of_file = "default.txt"
    
    call p%init( "./input/"//trim(name_of_file) )
    call pt%init(p)
        
    call p%show
    
    kh = p%wa%k * abs(p%glb%zstart)

    ug=30
    !$omp parallel do private(i,j,k,ii,jj,kk,x,y,z,kx,kz)
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            x = p%glb%x(i,j,k)
            y = p%glb%y(i,j,k)
            z = p%glb%z(i,j,k)

            kx = p%wa%k * x
            kz = p%wa%k * z

            if( z <= p%wa%L * Stokes_wave_interface(kx,p%wa%steepness,kh)  )then
                p%of(id)%loc%phi%now(i,j,k) = 1.0d0
                p%of(id)%loc%vel%x%now(i,j,k) = p%wa%U * Stokes_wave_u(kx,p%wa%steepness,kh,kz)
                p%of(id)%loc%vel%z%now(i,j,k) = p%wa%U * Stokes_wave_v(kx,p%wa%steepness,kh,kz)
            else
                p%of(id)%loc%phi%now(i,j,k) = -1.0d0
                p%of(id)%loc%vel%x%now(i,j,k) = 0.0d0
                p%of(id)%loc%vel%z%now(i,j,k) = 0.0d0
            endif
            
            p%of(id)%loc%vel%y%now(i,j,k) = 0.0_8
            
        end do
        end do
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%phi%now)
        call p%of(id)%bc(0,p%of(id)%loc%vof%now)

        call p%of(id)%bc(0,p%of(id)%loc%vel%x%now)
        call p%of(id)%bc(0,p%of(id)%loc%vel%y%now)
        call p%of(id)%bc(0,p%of(id)%loc%vel%z%now)
    enddo
    !$omp end parallel do

    call pt%vel%sync
    call pt%phi%sync
    call pt%vof%sync

    call level_set_rk3_redis(0)

    call p%node_vel
    call pt%nvel%sync
    call ns_init
    !---------------------------------------------------

    !$omp parallel do 
    do id = 0, p%glb%threads-1
        call p%of(id)%loc%init  
    enddo
    !$omp end parallel do

    call ns_ab_setup

    call p%ls_mv
    p%glb%ivol = p%glb%vol
    p%glb%imass = p%glb%mass
    p%glb%ivolv = p%glb%volv
    p%glb%imassv = p%glb%massv
    call p%sync

    write(*,*) "plotting"
    call plot
    call p%plot
    
    p%glb%ls_adv = 0.0d0
    p%glb%ls_red = 0.0d0
    p%glb%ppe    = 0.0d0
    p%glb%ns     = 0.0d0
    p%glb%syn    = 0.0d0

    p%glb%iter = 0
    p%glb%time = 0.0d0

end subroutine
