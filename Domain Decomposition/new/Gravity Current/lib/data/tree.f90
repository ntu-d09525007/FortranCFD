module tree
!$ use omp_lib
use branches
implicit none

type filemanager
integer :: ls_mv
integer :: c, energy
end type filemanager

type multigrid
integer :: n, nx, ny, nz
integer,dimension(:,:,:),allocatable :: node
integer,dimension(:),allocatable :: ipiv, i, j, k
real(8),dimension(:,:),allocatable :: A, AA
real(8),dimension(:),allocatable :: B, work, error, sol
end type multigrid

type manager
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
procedure total_c => manager_total_c
procedure total_e => manager_total_e
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
    write(*,'(A20,F15.8)')"dz:",p%glb%dz
    write(*,'(A20,F15.8)')"dt:",p%glb%dt
    write(*,'(A20,F15.4)')"Re:",p%glb%re
    write(*,'(A20,F15.4)')"We:",p%glb%we
    write(*,'(A20,F15.4)')"Fr:",p%glb%fr
    write(*,'(A20,F10.4)')"Density ratio:",p%glb%rho_12
    write(*,'(A20,F10.4)')"Viscosity ratio:",p%glb%mu_12
    write(*,'(A20,I5,A3,I5,A3,I5)')"Grids:",p%glb%node_x,"x",p%glb%node_y,"x",p%glb%node_z
    write(*,'(A20,I5,A3,I5,A3,I5)')"Threads Grid:",p%glb%grid_x,"x",p%glb%grid_y,"x",p%glb%grid_z
    write(*,'(A20,I5)')"Multigrid level:",p%glb%level
    
end subroutine

subroutine manager_read(p,path)
implicit none
class(manager) :: p
character(*) :: path

 open(unit=526,file=trim(path),status='old')
 
 read(526,*)
 read(526,*)p%glb%method
 read(526,*)
 read(526,*)p%glb%name
 read(526,*)
 read(526,*)p%glb%grid_x, p%glb%grid_y, p%glb%grid_z
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
 read(526,*)p%glb%zstart, p%glb%zend
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
 read(526,*)p%glb%gx,p%glb%gy,p%glb%gz
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
 read(526,*)p%glb%wbc(1), p%glb%wbc(2) 
 read(526,*)
 read(526,*)p%glb%us_c
 
 close(unit=526)
 
end subroutine

subroutine manager_init(p,path)
implicit none
class(manager) :: p
character(*) :: path
integer :: max_threads
integer :: i, j, k, id
real(8) :: mag

    call p%read(path)

    write(*,*)"finish read data from file"
    
    !p%fil%ls_mv = 15
    !open(unit=p%fil%ls_mv,file="./out/"//trim(p%glb%name)//"_MVloss.plt")
    !write(p%fil%ls_mv,*)'variables = "T" "Loss of mass" "Loss of Volume" '

    p%fil%c = 16
    open(unit=p%fil%c,file="./out/"//trim(p%glb%name)//"_Closs.plt")
    write(p%fil%c,*)'variables = "T" "Loss of concentration(%)" '

    p%fil%energy = 17
    open(unit=p%fil%energy,file="./out/"//trim(p%glb%name)//"_Energy.plt")
    write(p%fil%energy,*)'variables = "T" "Kine" "Poten" "Diss" "Sed" '
    
    p%glb%threads = p%glb%grid_x * p%glb%grid_y * p%glb%grid_z
    
    allocate( p%of(0:p%glb%threads-1))

    call omp_set_dynamic(.false.)
    call omp_set_num_threads(min(omp_get_max_threads(),p%glb%threads))

    write(*,*)"finish allocating nmumber of jobs"
    
    p%glb%node_x = p%glb%ug * ( p%glb%xend - p%glb%xstart )
    p%glb%node_y = p%glb%ug * ( p%glb%yend - p%glb%ystart )
    p%glb%node_z = p%glb%ug * ( p%glb%zend - p%glb%zstart )
    
    allocate( p%glb%x(0:p%glb%node_x+1,0:p%glb%node_y+1,0:p%glb%node_z+1), &
              p%glb%y(0:p%glb%node_x+1,0:p%glb%node_y+1,0:p%glb%node_z+1), &
              p%glb%z(0:p%glb%node_x+1,0:p%glb%node_y+1,0:p%glb%node_z+1) )

    write(*,*)"finish allocating public grids"
    
    !$omp parallel do
    do id = 0, p%glb%threads-1
        allocate( p%of(id)%glb%x(0:p%glb%node_x+1,0:p%glb%node_y+1,0:p%glb%node_z+1), &
                  p%of(id)%glb%y(0:p%glb%node_x+1,0:p%glb%node_y+1,0:p%glb%node_z+1), &
                  p%of(id)%glb%z(0:p%glb%node_x+1,0:p%glb%node_y+1,0:p%glb%node_z+1)  )
    enddo
    !$omp end parallel do
    
    write(*,*)"finish allocating private grids"
    
    p%glb%dx = ( p%glb%xend - p%glb%xstart ) / p%glb%node_x
    p%glb%dy = ( p%glb%yend - p%glb%ystart ) / p%glb%node_y
    p%glb%dz = ( p%glb%zend - p%glb%zstart ) / p%glb%node_z
        
    !$omp parallel do
    do i = 1, p%glb%node_x
        p%glb%x(i,:,:) = p%glb%xstart + (i-0.5)*p%glb%dx
    enddo
    !$omp end parallel do
    
    !$omp parallel do
    do j = 1, p%glb%node_y
        p%glb%y(:,j,:) = p%glb%ystart + (j-0.5)*p%glb%dy
    enddo
    !$omp end parallel do
    
    !$omp parallel do
    do k = 1, p%glb%node_z
        p%glb%z(:,:,k) = p%glb%zstart + (k-0.5)*p%glb%dz
    enddo
    !$omp end parallel do
    
    p%glb%x(0,:,:)=p%glb%xstart; p%glb%x(p%glb%node_x+1,:,:)=p%glb%xend
    p%glb%y(:,0,:)=p%glb%ystart; p%glb%y(:,p%glb%node_y+1,:)=p%glb%yend
    p%glb%z(:,:,0)=p%glb%zstart; p%glb%z(:,:,p%glb%node_z+1)=p%glb%zend

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
        call p%of(id)%init(id,p%glb%grid_x,p%glb%grid_y,p%glb%grid_z)
    enddo
    !$omp end parallel do

    write(*,*)"finish initializing job data"
    
    p%glb%time = 0.0_8
    p%glb%iter = 0
    p%glb%pid = 0
    
    p%glb%ls_adv = 0.0d0
    p%glb%ls_red = 0.0d0
    p%glb%ppe    = 0.0d0
    p%glb%ns     = 0.0d0
    p%glb%syn    = 0.0d0
    
    p%glb%e_kine = 0.0d0
    p%glb%e_diss = 0.0d0
    p%glb%e_poten = 0.0d0
    p%glb%e_sed  = 0.0d0

    mag = dsqrt(p%glb%gx**2.0d0+p%glb%gy**2.0d0+p%glb%gz**2.0d0)
    
    p%glb%gx = p%glb%gx / mag
    p%glb%gy = p%glb%gy / mag
    p%glb%gz = p%glb%gz / mag
    
    call system_clock( count_rate=p%glb%cpurate )
    
	if( p%glb%level > 0)then
        call p%mg_setup
        write(*,*)"finish initializing multigrid exact solver"
	endif

end subroutine

subroutine manager_total_e(p,init)
implicit none
class(manager) :: p
integer :: id,i,j,k
real(8) :: kine, poten, diss, sed, dv
real(8) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
logical :: init

dv = p%glb%dx * p%glb%dy * p%glb%dz

kine=0.0d0;poten=0.0d0;diss=0.0d0;sed=0.0d0
!$omp parallel do reduction(+:kine,poten,diss,sed), private(i,j,k,ux,uy,uz,vx,vy,vz,wx,wy,wz)
do id = 0, p%glb%threads-1

    do k = p%of(id)%loc%ks-p%glb%ghc, p%of(id)%loc%ke+p%glb%ghc
    do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc

        call p%of(id)%loc%ccdsolvers%x%solve("uccd",p%of(id)%loc%nvel%x%old(:,j,k),&
            p%of(id)%loc%vel_ten%xx(:,j,k),p%of(id)%loc%vel_ten%xxx(:,j,k),p%of(id)%loc%nvel%x%tmp(:,j,k))
        call p%of(id)%loc%ccdsolvers%x%solve("uccd",p%of(id)%loc%nvel%y%old(:,j,k),&
            p%of(id)%loc%vel_ten%yx(:,j,k),p%of(id)%loc%vel_ten%yxx(:,j,k),p%of(id)%loc%tdata%x%s1(:,j,k))
        call p%of(id)%loc%ccdsolvers%x%solve("uccd",p%of(id)%loc%nvel%z%old(:,j,k),&
            p%of(id)%loc%vel_ten%zx(:,j,k),p%of(id)%loc%vel_ten%zxx(:,j,k),p%of(id)%loc%tdata%x%s2(:,j,k))

    end do
    end do

    do k = p%of(id)%loc%ks-p%glb%ghc, p%of(id)%loc%ke+p%glb%ghc
    do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc

        call p%of(id)%loc%ccdsolvers%y%solve("uccd",p%of(id)%loc%nvel%x%old(i,:,k),&
            p%of(id)%loc%vel_ten%xy(i,:,k),p%of(id)%loc%vel_ten%xyy(i,:,k),p%of(id)%loc%tdata%y%s1(i,:,k))
        call p%of(id)%loc%ccdsolvers%y%solve("uccd",p%of(id)%loc%nvel%y%old(i,:,k),&
            p%of(id)%loc%vel_ten%yy(i,:,k),p%of(id)%loc%vel_ten%yyy(i,:,k), p%of(id)%loc%nvel%y%tmp(i,:,k))
        call p%of(id)%loc%ccdsolvers%y%solve("uccd",p%of(id)%loc%nvel%z%old(i,:,k),&
            p%of(id)%loc%vel_ten%zy(i,:,k),p%of(id)%loc%vel_ten%zyy(i,:,k),p%of(id)%loc%tdata%y%s2(i,:,k))

    end do
    end do

    do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
    do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc

        call p%of(id)%loc%ccdsolvers%z%solve("uccd",p%of(id)%loc%nvel%x%old(i,j,:),&
            p%of(id)%loc%vel_ten%xz(i,j,:),p%of(id)%loc%vel_ten%xzz(i,j,:),p%of(id)%loc%tdata%z%s1(i,j,:))
        call p%of(id)%loc%ccdsolvers%z%solve("uccd",p%of(id)%loc%nvel%y%old(i,j,:),&
            p%of(id)%loc%vel_ten%yz(i,j,:),p%of(id)%loc%vel_ten%yzz(i,j,:),p%of(id)%loc%tdata%z%s2(i,j,:))
        call p%of(id)%loc%ccdsolvers%z%solve("uccd",p%of(id)%loc%nvel%z%old(i,j,:),&
            p%of(id)%loc%vel_ten%zz(i,j,:),p%of(id)%loc%vel_ten%zzz(i,j,:), p%of(id)%loc%nvel%z%tmp(i,j,:))
                        
    end do
    end do

    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie

        ux=p%of(id)%loc%vel_ten%xx(i,j,k);uy=p%of(id)%loc%vel_ten%xy(i,j,k);uz=p%of(id)%loc%vel_ten%xz(i,j,k)
        vx=p%of(id)%loc%vel_ten%yx(i,j,k);vy=p%of(id)%loc%vel_ten%yy(i,j,k);vz=p%of(id)%loc%vel_ten%yz(i,j,k)
        wx=p%of(id)%loc%vel_ten%zx(i,j,k);wy=p%of(id)%loc%vel_ten%zy(i,j,k);wz=p%of(id)%loc%vel_ten%zz(i,j,k)

        kine=kine+(p%of(id)%loc%nvel%x%now(i,j,k)**2.0d0+p%of(id)%loc%nvel%y%now(i,j,k)**2.0d0+p%of(id)%loc%nvel%z%now(i,j,k)**2.0d0)/2.0d0
        poten=poten+p%glb%z(i,j,k)*p%of(id)%loc%c%now(i,j,k)
        diss=diss+(ux**2.0d0+vy**2.0d0+wz**2.0d0+((vx+uy)**2.0d0+(wy+vz)**2.0d0+(uz+wx)**2.0d0)/2.0d0)
        sed=sed+p%glb%us_c*p%of(id)%loc%c%now(i,j,k)
        
    enddo
    enddo
    enddo

enddo
!$omp end parallel do

kine=kine*dv
poten=poten*dv
diss=2.0d0*diss*dv/p%glb%re
sed=sed*dv

p%glb%e_kine=kine
p%glb%e_poten=poten
p%glb%e_diss=p%glb%e_diss+diss*p%glb%dt
p%glb%e_sed=p%glb%e_sed+sed*p%glb%dt

p%glb%e_t = p%glb%e_kine+p%glb%e_poten+p%glb%e_diss+p%glb%e_sed

if(init)p%glb%e_t0=p%glb%e_t

end subroutine

subroutine manager_total_c(p,flag)
implicit none
class(manager) :: p
integer :: id, i, j, k
real(8) :: dv, csum
logical :: flag

    dv = p%glb%dx * p%glb%dy * p%glb%dz

    csum=0.0d0
    !$omp parallel do private(i,j,k), reduction(+:csum)    
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            
            csum = csum + p%of(id)%loc%c%now(i,j,k)*dv
        
        enddo
        enddo
        enddo
    
    enddo
    !$omp end parallel do
    
    p%glb%csum = csum
    if(flag)p%glb%icsum=csum

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
integer :: id,i,j,k

    !$omp parallel do private(i,j,k)
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            p%of(id)%loc%nvel%x%now(i,j,k) = 0.5d0 * ( p%of(id)%loc%vel%x%now(i-1,j,k) + p%of(id)%loc%vel%x%now(i,j,k) )
            p%of(id)%loc%nvel%y%now(i,j,k) = 0.5d0 * ( p%of(id)%loc%vel%y%now(i,j-1,k) + p%of(id)%loc%vel%y%now(i,j,k) )
            p%of(id)%loc%nvel%z%now(i,j,k) = 0.5d0 * ( p%of(id)%loc%vel%z%now(i,j,k-1) + p%of(id)%loc%vel%z%now(i,j,k) )
        end do
        end do
        end do
    
        call p%of(id)%nvelbc(p%of(id)%loc%nvel%x%now,p%of(id)%loc%nvel%y%now,p%of(id)%loc%nvel%z%now)
     
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
        call p%of(id)%loc%c%switch
        
    enddo      
    !$omp end parallel do

end subroutine



end module tree
