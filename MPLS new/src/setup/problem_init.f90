subroutine problem_init()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,ug,ii,jj
real(8) :: x,y,err
real(8) :: kx,ky
CHARACTER(100) :: NAME_OF_FILE
    
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
    
    ug=30
    !$omp parallel do private(i,j,ii,jj,x,y)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            x = p%glb%x(i,j)
            y = p%glb%y(i,j)

            p%of(id)%loc%phi%now(i,j) = - dsqrt((x-0.5d0)**2.0d0+(y-0.75d0)**2.0d0)+0.15d0

            p%of(id)%loc%vel%x%now(i,j) = 0.0_8
            p%of(id)%loc%vel%y%now(i,j) = 0.0_8

        end do
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%phi%now)
        call p%of(id)%bc(0,p%of(id)%loc%vof%now)

        call p%of(id)%bc(0,p%of(id)%loc%vel%x%now)
        call p%of(id)%bc(0,p%of(id)%loc%vel%y%now)

    enddo
    !$omp end parallel do

    call pt%vel%sync
    call pt%phi%sync
    call pt%vof%sync

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
