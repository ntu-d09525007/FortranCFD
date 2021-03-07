module tree
!$ use omp_lib
use branches
implicit none

type wave
real(8) :: phase_speed, wavelength, amplitude
real(8) :: k, w
real(8) :: L, U
end type wave

type filemanager
integer :: ls_mv, damdata
end type filemanager

type multigrid
integer :: n, nx, ny
integer,dimension(:,:),allocatable :: node
integer,dimension(:),allocatable :: ipiv, i, j
real(8),dimension(:,:),allocatable :: A, AA
real(8),dimension(:),allocatable :: B, work, error, sol
end type multigrid

type manager
type(wave) :: wa
type(global) :: glb
type(filemanager) :: fil
type(job),allocatable :: of(:)
type(multigrid) :: mg
contains
procedure show => manager_show
procedure read => manager_read
procedure init => manager_init
procedure mg_setup => manager_mg_setup
procedure mg_solve_exact => manager_mg_solve_exact
procedure mg_find_error => manager_mg_find_error
procedure mg_solve_correct => manager_mg_solve_correct
procedure sync => manager_sync
procedure ls_mv => manager_ls_mv
procedure ls_funs => manager_ls_funs
procedure rho_mu => manager_rho_mu
procedure surface_norms => manager_surface_norms
procedure surface_norms2 => manager_surface_norms_sec
procedure node_vel => manager_node_vel
procedure curv => manager_curv
procedure switch => manager_switch
procedure plot => manager_plot
end type manager

contains 

include './tree_multigrid.f90'
include './tree_twophaseflows.f90'
include './tree_plot.f90'

subroutine manager_show(p)
implicit none
class(manager) :: p
integer :: id, level

    write(*,*)" --- Problem information --- "
    write(*,'(A30,I5)')"Number of computing threads:",p%glb%threads
    write(*,'(A20)')p%glb%name
    write(*,'(A20,F15.8)')"dx:",p%glb%dx
    write(*,'(A20,F15.8)')"dy:",p%glb%dy

    write(*,'(A20,F15.8)')"dt:",p%glb%dt
    write(*,'(A20,F15.4)')"Re:",p%glb%re
    write(*,'(A20,F15.4)')"We:",p%glb%we
    write(*,'(A20,F15.4)')"Fr:",p%glb%fr
    write(*,'(A20,F10.4)')"Density ratio:",p%glb%rho_12
    write(*,'(A20,F10.4)')"Viscosity ratio:",p%glb%mu_12
    write(*,'(A20,I5,A3,I5)')"Grids:",p%glb%node_x,"x",p%glb%node_y
    write(*,'(A20,I5,A3,I5)')"Threads Grid:",p%glb%grid_x,"x",p%glb%grid_y
    write(*,'(A20,I5)')"Multigrid level:",p%glb%level
    write(*,'(A20,I8)')"Overlap layer",p%of(id)%glb%ghc
    write(*,'(A20,L5)')"X Periodic:",p%glb%xper
    write(*,'(A20,L5)')"Y Periodic:",p%glb%yper

    ! write(*,*)" --- SubDomain Information  --- "
    ! do id = 0, p%glb%threads-1
    !    write(*,'("ID ",I2,": (",I2,",",I2,",",I2,")")')ID,p%of(id)%loc%idx,p%of(id)%loc%idy,p%of(id)%loc%idz
    !    write(*,'(A20,I4,A3,I4)')"X index:",p%of(id)%loc%is,"~",p%of(id)%loc%ie
    !    write(*,'(A20,I4,A3,I4)')"Y index:",p%of(id)%loc%js,"~",p%of(id)%loc%je
    !    write(*,'(A20,I4,A3,I4)')"Z index:",p%of(id)%loc%ks,"~",p%of(id)%loc%ke
    !    write(*,*)""
    ! end do

    ! write(*,'(A)')"Multigrid information"
    ! do id  = 0, p%glb%threads-1
    !     write(*,'("ID ",I2," :",I3,"x",I3,"x",I3)')id,p%of(id)%loc%ie-p%of(id)%loc%is+1,p%of(id)%loc%je-p%of(id)%loc%js+1,p%of(id)%loc%ke-p%of(id)%loc%ks+1
    !     write(*,*)"-----------------------"
    !     do level = 1, p%glb%level
    !         write(*,'("Level: ",I2," ",I5,"x",I5,"x",I5)')level,p%of(id)%loc%mg(level)%nx,p%of(id)%loc%mg(level)%ny,p%of(id)%loc%mg(level)%nz
    !         write(*,'("dx=",ES11.4,",dy=",ES11.4,",dz=",ES11.4)')p%of(id)%loc%mg(level)%dx,p%of(id)%loc%mg(level)%dy,p%of(id)%loc%mg(level)%dz
    !     end do
    !     write(*,*)"======================="
    ! end do
    
end subroutine

subroutine manager_read(p,path)
implicit none
class(manager) :: p
character(*) :: path
integer :: x,y

 open(unit=526,file=trim(path),status='old')
 
 read(526,*)
 read(526,*)p%glb%method
 read(526,*)
 read(526,*)p%glb%name
 read(526,*)
 read(526,*)p%glb%grid_x, p%glb%grid_y
 read(526,*)
 read(526,*)p%glb%level
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
 read(526,*)p%glb%t_tol, p%glb%t_w
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
 read(526,*)p%glb%gx, p%glb%gy
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
 read(526,*)
 read(526,*)x,y
 read(526,*)
 read(526,*)p%wa%phase_speed, p%wa%wavelength, p%wa%amplitude
 close(unit=526)

 if( x==1 )then
    p%glb%xper = .true.
 else
    p%glb%xper = .false.
 endif

 if( y==1 )then
    p%glb%yper = .true.
 else
    p%glb%yper = .false.
 endif
 
end subroutine

subroutine manager_init(p,path)
implicit none
class(manager) :: p
character(*) :: path
integer :: max_threads
integer :: i, j, id
real(8) :: mag

    call p%read(path)

    write(*,*)"finish read data from file"
    
    p%fil%ls_mv = 15
    open(unit=p%fil%ls_mv,file="./out/"//trim(p%glb%name)//"_MVloss.plt")
    write(p%fil%ls_mv,*)'variables = "T" "Loss of mass" "Loss of Volume" '

    p%fil%damdata = 16
    open(unit=p%fil%damdata,file="./out/"//trim(p%glb%name)//"_DamData.plt")
    write(p%fil%ls_mv,*)'variables = "T" "Damfront" "Wall" '
    
    p%glb%threads = p%glb%grid_x * p%glb%grid_y 
    
    allocate( p%of(0:p%glb%threads-1))

    call omp_set_dynamic(.false.)
    call omp_set_num_threads(min(omp_get_max_threads(),p%glb%threads))

    write(*,*)"finish allocating nmumber of jobs"
    
    p%glb%node_x = p%glb%ug * ( p%glb%xend - p%glb%xstart )
    p%glb%node_y = p%glb%ug * ( p%glb%yend - p%glb%ystart )

    allocate( p%glb%x(0:p%glb%node_x+1,0:p%glb%node_y+1),&
             &p%glb%y(0:p%glb%node_x+1,0:p%glb%node_y+1))

    write(*,*)"finish allocating public grids"
    
    !$omp parallel do
    do id = 0, p%glb%threads-1
        allocate( p%of(id)%glb%x(0:p%glb%node_x+1,0:p%glb%node_y+1),&
                 &p%of(id)%glb%y(0:p%glb%node_x+1,0:p%glb%node_y+1))
    enddo
    !$omp end parallel do
    
    write(*,*)"finish allocating private grids"
    
    p%glb%dx = ( p%glb%xend - p%glb%xstart ) / p%glb%node_x
    p%glb%dy = ( p%glb%yend - p%glb%ystart ) / p%glb%node_y
  
    !$omp parallel do
    do i = 1, p%glb%node_x
        p%glb%x(i,:) = p%glb%xstart + (i-0.5)*p%glb%dx
    enddo
    !$omp end parallel do
    
    !$omp parallel do
    do j = 1, p%glb%node_y
        p%glb%y(:,j) = p%glb%ystart + (j-0.5)*p%glb%dy
    enddo
    !$omp end parallel do
    
    p%glb%x(0,:)=p%glb%xstart; p%glb%x(p%glb%node_x+1,:)=p%glb%xend
    p%glb%y(:,0)=p%glb%ystart; p%glb%y(:,p%glb%node_y+1)=p%glb%yend

    write(*,*)"finish assigning grid data"
    
    p%glb%dt = p%glb%dt * p%glb%dx
    p%glb%rdt = p%glb%rdt * p%glb%dx
    p%glb%ls_wid = p%glb%ls_wid * p%glb%dx
    
    p%glb%p_tol = 0.1_8 ** p%glb%p_tol
    p%glb%t_tol = 0.1_8 ** p%glb%t_tol

    p%glb%p_b = 0.1_8 ** p%glb%p_b

    write(*,*)"finsih numeric setting"

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
        p%glb%we = p%glb%rho_1 * p%glb%u**2.0d0 * p%glb%L / p%glb%sigma
        p%glb%fr = p%glb%u**2.0d0 / ( p%glb%g * p%glb%L ) 
    endif

    write(*,*)"finish calculating problem parameters"
    
    call p%sync

    write(*,*)"finish assinging problem parameters to job"
    
    !$omp parallel do
    do id = 0, p%glb%threads-1
        call p%of(id)%init(id,p%glb%grid_x,p%glb%grid_y)
    enddo
    !$omp end parallel do

    write(*,*)"finish initializing job data"

    allocate( p%glb%id(0:p%glb%grid_x-1,0:p%glb%grid_y-1) )

    do id = 0, p%glb%threads-1
        p%glb%id(p%of(id)%loc%idx,p%of(id)%loc%idy) = id
    enddo

    p%glb%time = 0.0_8
    p%glb%iter = 0
    p%glb%pid = 0
    
    p%glb%ls_adv = 0.0d0
    p%glb%ls_red = 0.0d0
    p%glb%ppe    = 0.0d0
    p%glb%ns     = 0.0d0
    p%glb%syn    = 0.0d0
    
    mag = dsqrt(p%glb%gx**2.0d0+p%glb%gy**2.0d0)
    
    p%glb%gx = p%glb%gx / mag
    p%glb%gy = p%glb%gy / mag
    
    call system_clock( count_rate=p%glb%cpurate )
    
    if( p%glb%level > 0)then
        call p%mg_setup
        write(*,*)"finish initializing multigrid exact solver"
    endif

    ! ----- Wave papmeters setup --------
    p%wa%k = 2.0d0 * dacos(-1.0d0) * p%glb%L / p%wa%wavelength
    p%wa%w = p%wa%phase_speed / p%glb%U * p%wa%k
    p%wa%L = p%wa%amplitude / p%glb%L 
    p%wa%U = p%wa%L * p%wa%w
    
end subroutine


subroutine manager_sync(p)
implicit none
class(manager) :: p
integer :: i, id

!$omp parallel do
do id = 0, p%glb%threads-1
    p%of(id)%glb = p%glb
enddo
!$omp end parallel do
 
end subroutine

subroutine manager_node_vel(p)
implicit none
class(manager) :: p
integer :: id,i,j

!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        p%of(id)%loc%nvel%x%now(i,j) = 0.5d0 * ( p%of(id)%loc%vel%x%now(i-1,j) + p%of(id)%loc%vel%x%now(i,j) )
        p%of(id)%loc%nvel%y%now(i,j) = 0.5d0 * ( p%of(id)%loc%vel%y%now(i,j-1) + p%of(id)%loc%vel%y%now(i,j) )
    end do
    end do

    call p%of(id)%nvelbc(p%of(id)%loc%nvel%x%now,p%of(id)%loc%nvel%y%now)
 
enddo       
!$omp end parallel do

end subroutine

subroutine manager_switch(p)
implicit none
class(manager) :: p
integer :: id

!$omp parallel do
do id = 0, p%glb%threads-1
    
    call p%of(id)%loc%phi%switch
    call p%of(id)%loc%vof%switch
    
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
    
enddo      
!$omp end parallel do

end subroutine

end module tree
