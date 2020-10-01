module tree
!$ use omp_lib
use branches
implicit none

type filemanager
integer :: ls_mv
end type filemanager

type multigrid
integer :: n, nx, ny, nz
integer,dimension(:,:,:),allocatable :: node
integer,dimension(:),allocatable :: ipiv, i, j, k
real(8),dimension(:,:),allocatable :: A, AA
real(8),dimension(:),allocatable :: B, work
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
procedure sync => manager_sync
procedure ls_mv => manager_ls_mv
procedure ls_funs => manager_ls_funs
procedure rho_mu => manager_rho_mu
procedure surface_norms => manager_surface_norms
procedure surface_norms2 => manager_surface_norms_sec
procedure node_vel => manager_node_vel
procedure curv => manager_curv
procedure switch => manager_switch
procedure vortex => manager_vortex_identification
procedure plot => manager_plot
procedure plots => manager_plots
end type manager

contains 

subroutine manager_show(p)
implicit none
class(manager) :: p
integer :: id, level

    write(*,*)" --- Problem information --- "
    write(*,'(A30,I5)')"Number of computing threads:",p%glb%threads
    write(*,'(A30,I5)')"Number of ploting threads:",p%glb%pthreads
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
    
    !write(*,'(A20,I8)')"Overlap layer",p%of(id)%glb%ghc
    !write(*,*)" --- SubDomain Information  --- "
    !do id = 0, p%glb%threads-1
    !   write(*,'("ID ",I2,": (",I1,",",I1,",",I1,")")')ID,p%of(id)%loc%idx,p%of(id)%loc%idy,p%of(id)%loc%idz
    !   write(*,'(A20,I4,A3,I4)')"X index:",p%of(id)%loc%is,"~",p%of(id)%loc%ie
    !   write(*,'(A20,I4,A3,I4)')"Y index:",p%of(id)%loc%js,"~",p%of(id)%loc%je
    !   write(*,'(A20,I4,A3,I4)')"Z index:",p%of(id)%loc%ks,"~",p%of(id)%loc%ke
    !   write(*,*)""
    !end do
    
    ! write(*,'(A)')"Multigrid information"
    ! do id  = 0, p%glb%threads-1
        ! write(*,'("ID ",I2," :",I3,"x",I3,"x",I3)')id,p%of(id)%loc%ie-p%of(id)%loc%is+1,p%of(id)%loc%je-p%of(id)%loc%js+1,p%of(id)%loc%ke-p%of(id)%loc%ks+1
        ! write(*,*)"-----------------------"
        ! do level = 1, p%glb%level
            ! write(*,'("Level: ",I2," ",I5,"x",I5,"x",I5)')level,p%of(id)%loc%mg(level)%nx,p%of(id)%loc%mg(level)%ny,p%of(id)%loc%mg(level)%nz
            ! write(*,'("dx=",ES11.4,",dy=",ES11.4,",dz=",ES11.4)')p%of(id)%loc%mg(level)%dx,p%of(id)%loc%mg(level)%dy,p%of(id)%loc%mg(level)%dz
        ! end do
        ! write(*,*)"======================="
    ! end do
    
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
 read(526,*)p%glb%pthreads
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
    
    p%fil%ls_mv = 15
    open(unit=p%fil%ls_mv,file="./out/"//trim(p%glb%name)//"_MVloss.plt")
    write(p%fil%ls_mv,*)'variables = "T" "Loss of mass" "Loss of Volume" '
    
    p%glb%threads = p%glb%grid_x * p%glb%grid_y * p%glb%grid_z
    
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

    write(*,*)"finish allocating nmumber of jobs"
    
    p%glb%node_x = p%glb%ug * ( p%glb%xend - p%glb%xstart )
    p%glb%node_y = p%glb%ug * ( p%glb%yend - p%glb%ystart )
    p%glb%node_z = p%glb%ug * ( p%glb%zend - p%glb%zstart )
    
    allocate( p%glb%x(0:p%glb%node_x+1), p%glb%y(0:p%glb%node_y+1), p%glb%z(0:p%glb%node_z+1) )

    write(*,*)"finish allocating public grids"
    
    !$omp parallel private(id), num_threads(p%glb%threads)
        id = 0
        !$ id = omp_get_thread_num()
        allocate( p%of(id)%glb%x(0:p%glb%node_x+1), p%of(id)%glb%y(0:p%glb%node_y+1), p%of(id)%glb%z(0:p%glb%node_z+1) )
    !$omp end parallel
    
    write(*,*)"finish allocating private grids"
    
    p%glb%dx = ( p%glb%xend - p%glb%xstart ) / p%glb%node_x
    p%glb%dy = ( p%glb%yend - p%glb%ystart ) / p%glb%node_y
    p%glb%dz = ( p%glb%zend - p%glb%zstart ) / p%glb%node_z
        
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
    
    !$omp parallel do
    do k = 1, p%glb%node_z
        p%glb%z(k) = p%glb%zstart + (k-0.5)*p%glb%dz
    enddo
    !$omp end parallel do
    
    p%glb%x(0)=p%glb%xstart; p%glb%x(p%glb%node_x+1)=p%glb%xend
    p%glb%y(0)=p%glb%ystart; p%glb%y(p%glb%node_y+1)=p%glb%yend
    p%glb%z(0)=p%glb%zstart; p%glb%z(p%glb%node_z+1)=p%glb%zend

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
    
    !$omp parallel private(id), num_threads(p%glb%threads)
        id=0
        !$ id = omp_get_thread_num()
        call p%of(id)%init(id,p%glb%grid_x,p%glb%grid_y,p%glb%grid_z)
    !$omp end parallel

    write(*,*)"finish initializing job data"
    
    p%glb%time = 0.0_8
    p%glb%iter = 0
    p%glb%pid = 0
    
    p%glb%ls_adv = 0.0d0
    p%glb%ls_red = 0.0d0
    p%glb%ppe    = 0.0d0
    p%glb%ns     = 0.0d0
    p%glb%syn    = 0.0d0
    
    mag = dsqrt(p%glb%gx**2.0d0+p%glb%gy**2.0d0+p%glb%gz**2.0d0)
    
    p%glb%gx = p%glb%gx / mag
    p%glb%gy = p%glb%gy / mag
    p%glb%gz = p%glb%gz / mag
    
    call system_clock( count_rate=p%glb%cpurate )
    
    call p%mg_setup
    write(*,*)"finish initializing multigrid exact solver"

end subroutine

subroutine manager_mg_setup(p)
implicit none
class(manager) :: p
integer :: i,j,k,id,level,ind,ii,jj,kk
integer :: nx,ny,nz,n,info
integer :: imax,imin,jmax,jmin,kmax,kmin
real(8) :: dx,dy,dz

    !write(*,*)"start multigrid setup"

    level = p%glb%level
    
    p%mg%n = p%glb%grid_x*p%glb%grid_y*p%glb%grid_z*p%of(0)%loc%mg(level)%n
    p%mg%nx = p%glb%grid_x*p%of(0)%loc%mg(level)%nx
    p%mg%ny = p%glb%grid_y*p%of(0)%loc%mg(level)%ny
    p%mg%nz = p%glb%grid_z*p%of(0)%loc%mg(level)%nz
    
    write(*,'(A,4I6)')"Local A size",p%of(0)%loc%mg(level)%n,p%of(0)%loc%mg(level)%nx,p%of(0)%loc%mg(level)%ny,p%of(0)%loc%mg(level)%nz
    write(*,'(A,4I6)')"Global A size",p%mg%n,p%mg%nx,p%mg%ny,p%mg%nz
    
    allocate(p%mg%A(p%mg%n,p%mg%n),p%mg%AA(p%mg%n,p%mg%n),p%mg%B(p%mg%n),p%mg%ipiv(p%mg%n))
    allocate(p%mg%i(p%mg%n),p%mg%j(p%mg%n),p%mg%k(p%mg%n),p%mg%node(p%mg%nx,p%mg%ny,p%mg%nz))
    !allocate(p%mg%iA(p%mg%n,p%mg%n),p%mg%work(p%mg%n))
    
    dx = p%of(0)%loc%mg(level)%dx
    dy = p%of(0)%loc%mg(level)%dy
    dz = p%of(0)%loc%mg(level)%dz
    
    do i = 1, p%mg%n
    do j = 1, p%mg%n
        p%mg%A(i,j)=0.0d0
    enddo
    enddo
    
    imax=0; imin=10000
    jmax=0; jmin=10000
    kmax=0; kmin=10000
    do i = 1, p%mg%n    
        id = (i-1) / p%of(0)%loc%mg(level)%n
        ! global index of local A
        ind = i - id*p%of(0)%loc%mg(level)%n     
        
        ! global (i,j,k) of global A
        ii = p%of(id)%loc%mg(level)%i(ind) + p%of(0)%loc%mg(level)%nx*p%of(id)%loc%idx 
        jj = p%of(id)%loc%mg(level)%j(ind) + p%of(0)%loc%mg(level)%ny*p%of(id)%loc%idy
        kk = p%of(id)%loc%mg(level)%k(ind) + p%of(0)%loc%mg(level)%nz*p%of(id)%loc%idz
        
        p%mg%i(i) = ii
        p%mg%j(i) = jj
        p%mg%k(i) = kk
        p%mg%node(ii,jj,kk) = i
        
        p%mg%A(i,i) = -2.0d0/dx**2.0d0-2.0d0/dy**2.0d0-2.0d0/dz**2.0d0
        
        imax=max(imax,ii); imin=min(imin,ii)
        jmax=max(jmax,jj); jmin=min(jmin,jj)
        kmax=max(kmax,kk); kmin=min(kmin,kk)
    end do
    
    write(*,'(A15,2I10)')"I matrix",imin,imax
    write(*,'(A15,2I10)')"J matrix",jmin,jmax
    write(*,'(A15,2I10)')"K matrix",kmin,kmax
    
    do i = 1, p%mg%n
    
        if(p%mg%i(i)>1)then
            p%mg%A(i, p%mg%node(p%mg%i(i)-1,p%mg%j(i),p%mg%k(i)) ) = 1.0d0/dx**2.0d0
        else
            p%mg%A(i,i) = p%mg%A(i,i) + 1.0d0/dx**2.0d0
        endif
        
        if(p%mg%i(i)<p%mg%nx)then
            p%mg%A(i, p%mg%node(p%mg%i(i)+1,p%mg%j(i),p%mg%k(i)) ) = 1.0d0/dx**2.0d0
        else
            p%mg%A(i,i) = p%mg%A(i,i) + 1.0d0/dx**2.0d0
        endif
        
        if(p%mg%j(i)>1)then
            p%mg%A(i, p%mg%node(p%mg%i(i),p%mg%j(i)-1,p%mg%k(i)) ) = 1.0d0/dy**2.0d0
        else
            p%mg%A(i,i) = p%mg%A(i,i) + 1.0d0/dy**2.0d0
        endif
        
        if(p%mg%j(i)<p%mg%ny)then
            p%mg%A(i, p%mg%node(p%mg%i(i),p%mg%j(i)+1,p%mg%k(i)) ) = 1.0d0/dy**2.0d0
        else
            p%mg%A(i,i) = p%mg%A(i,i) + 1.0d0/dy**2.0d0
        endif   

        if(p%mg%k(i)>1)then
            p%mg%A(i, p%mg%node(p%mg%i(i),p%mg%j(i),p%mg%k(i)-1) ) = 1.0d0/dz**2.0d0
        else
            p%mg%A(i,i) = p%mg%A(i,i) + 1.0d0/dz**2.0d0
        endif
        
        if(p%mg%k(i)<p%mg%nz)then
            p%mg%A(i, p%mg%node(p%mg%i(i),p%mg%j(i),p%mg%k(i)+1) ) = 1.0d0/dz**2.0d0
        else
            p%mg%A(i,i) = p%mg%A(i,i) + 1.0d0/dz**2.0d0
        endif   
        
    end do
    
end subroutine

subroutine manager_mg_solve_exact(p,output)
implicit none
class(manager) :: p
integer :: i,j,k,id,level
integer :: ii,jj,kk,ind,info
real(8) :: error, tmp
logical :: output

level = p%glb%level

!$omp parallel do private(id,ind,ii,jj,kk), num_threads(p%glb%threads)
do i = 1, p%mg%n
    id =  (i-1) / p%of(0)%loc%mg(level)%n
    ind = i - id*p%of(0)%loc%mg(level)%n
    ii = p%of(id)%loc%mg(level)%i(ind)
    jj = p%of(id)%loc%mg(level)%j(ind)
    kk = p%of(id)%loc%mg(level)%k(ind)
    p%mg%B(i) = p%of(id)%loc%mg(level)%src(ii,jj,kk)
    !call random_number(p%of(id)%loc%mg(level)%sol(ii,jj,kk))
enddo
!$omp end parallel do 

!$omp parallel do num_threads(p%glb%threads)
do i = 1, p%mg%n
do j = 1, p%mg%n
    p%mg%AA(i,j) = p%mg%A(i,j)
end do
end do
!$omp end parallel do

! do i = 1, p%mg%n
! p%mg%B(i)=0.0d0
! do j = 1, p%mg%n
    ! id =  (j-1) / p%of(0)%loc%mg(level)%n
    ! ind = j - id*p%of(0)%loc%mg(level)%n
    ! ii = p%of(id)%loc%mg(level)%i(ind)
    ! jj = p%of(id)%loc%mg(level)%j(ind)
    ! kk = p%of(id)%loc%mg(level)%k(ind)
    ! p%mg%B(i)=p%mg%B(i)+p%mg%A(i,j)*p%of(id)%loc%mg(level)%sol(ii,jj,kk)
! enddo
! id =  (i-1) / p%of(0)%loc%mg(level)%n
! ind = i - id*p%of(0)%loc%mg(level)%n
! ii = p%of(id)%loc%mg(level)%i(ind)
! jj = p%of(id)%loc%mg(level)%j(ind)
! kk = p%of(id)%loc%mg(level)%k(ind)
! p%of(id)%loc%mg(level)%src(ii,jj,kk) = p%mg%B(i)
! enddo

call dgesv(p%mg%n,1,p%mg%AA,p%mg%n,p%mg%ipiv,p%mg%B,p%mg%n,info)
if(info.ne.0)stop "Something Wrong with Multigrid Solver!!"

error=0.0d0
!$omp parallel do private(id,ind,j,ii,jj,kk,tmp), num_threads(p%glb%threads), reduction(max:error)
do i = 1, p%mg%n
    id = (i-1) / p%of(0)%loc%mg(level)%n
    ind = i - id*p%of(0)%loc%mg(level)%n
    ii = p%of(id)%loc%mg(level)%i(ind)
    jj = p%of(id)%loc%mg(level)%j(ind)
    kk = p%of(id)%loc%mg(level)%k(ind)
    !error = max(error, abs(p%of(id)%loc%mg(level)%sol(ii,jj,kk)-p%mg%B(i)))
    p%of(id)%loc%mg(level)%sol(ii,jj,kk) = p%mg%B(i) 
    tmp=-p%of(id)%loc%mg(level)%src(ii,jj,kk)
    do j = 1, p%mg%n
        id = (i-1) / p%of(0)%loc%mg(level)%n
        ind = i - id*p%of(0)%loc%mg(level)%n
        ii = p%of(id)%loc%mg(level)%i(ind)
        jj = p%of(id)%loc%mg(level)%j(ind)
        kk = p%of(id)%loc%mg(level)%k(ind)
        tmp = tmp + p%mg%A(i,j)*p%mg%B(j)
    enddo
    error = max(error,abs(tmp))
enddo
!$omp end parallel do 

if(output)write(*,*)"Exact solver error:",error

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
integer :: id, i, j, k
real(8) :: x, heavy, hp, pi, eps


    eps = 1.0d-12
    pi = dacos(-1.0_8)

    !$omp parallel private(id,i,j,k,x,heavy,hp), num_threads(p%glb%threads), shared(eps,pi)
    
        id=0
        !$ id = omp_get_thread_num()
        
        do k = p%of(id)%loc%ks-p%glb%ghc, p%of(id)%loc%ke+p%glb%ghc
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        
            x = p%of(id)%loc%phi%now(i,j,k) / p%glb%ls_wid
            
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
            
            heavy = max(min(heavy,1.0d0),0.0d0) 
            
            p%of(id)%loc%heavy%now(i,j,k) = heavy
            p%of(id)%loc%delta%now(i,j,k) = hp
            p%of(id)%loc%sign%now(i,j,k) = 2.0_8*heavy-1.0_8
            
        end do
        end do
        end do
    
    !$omp end parallel

end subroutine

subroutine manager_rho_mu(p)
implicit none
class(manager) :: p
integer :: id,i,j,k
real(8) :: heavy

    call p%ls_funs
    
    !$omp parallel private(id,i,j,k,heavy), num_threads(p%glb%threads)
    
        id=0
        !$ id = omp_get_thread_num()

        do k = p%of(id)%loc%ks-p%glb%ghc, p%of(id)%loc%ke+p%glb%ghc        
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
            
            heavy = p%of(id)%loc%vof%now(i,j,k)
            if( p%glb%method .ne. 3)heavy = p%of(id)%loc%heavy%now(i,j,k)
            
            p%of(id)%loc%rho%now(i,j,k) = heavy + p%glb%rho_12 * (1.0_8 - heavy )
            p%of(id)%loc%mu%now(i,j,k)  = heavy + p%glb%mu_12  * (1.0_8 - heavy )
            
        end do
        end do
        end do
        
    !$omp end parallel
    
    
end subroutine

subroutine manager_ls_mv(p)
implicit none
class(manager) :: p
integer :: id, i, j, k
real(8) :: mass, vol, rho
real(8) :: dv

    dv = p%glb%dx * p%glb%dy * p%glb%dz

    !call p%rho_mu

    !===========================  LS 
    
    mass = 0.0_8; vol=0.0_8
    
    !$omp parallel private(id,i,j,k,rho), num_threads(p%glb%threads), reduction(+:mass,vol)
        
        id=0
        !$ id = omp_get_thread_num()
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            
            !rho = p%of(id)%loc%heavy%now(i,j,k) + p%glb%rho_12 * (1.0d0 - p%of(id)%loc%heavy%now(i,j,k))
            !mass = mass + rho*p%of(id)%loc%heavy%now(i,j,k)*dv
            !vol = vol + p%of(id)%loc%heavy%now(i,j,k)*dv
            
            vol = vol + p%of(id)%loc%phi%now(i,j,k)*dv
        
        enddo
        enddo
        enddo
    
    !$omp end parallel
    
    p%glb%mass = mass
    p%glb%vol = vol
    
    !===========================  VOF 

    mass = 0.0_8; vol=0.0_8
    
    !$omp parallel private(id,i,j,k,rho), num_threads(p%glb%threads), reduction(+:mass,vol)
        
        id=0
        !$ id = omp_get_thread_num()
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            
            rho = p%of(id)%loc%vof%now(i,j,k) + p%glb%rho_12 * (1.0d0-p%of(id)%loc%vof%now(i,j,k))
            mass = mass + rho*p%of(id)%loc%vof%now(i,j,k)*dv
            vol = vol + p%of(id)%loc%vof%now(i,j,k)*dv
        
        enddo
        enddo
        enddo
    
    !$omp end parallel
    
    p%glb%massv = mass
    p%glb%volv = vol
    
    call p%sync


end subroutine

subroutine manager_surface_norms(p)
implicit none
class(manager) :: p
integer :: id,I,J,k  
    
    !$omp parallel private(id,i,j,k), num_threads(p%glb%threads)
    
        id=0
        !$ id = omp_get_thread_num()
        
        do k = p%of(id)%loc%ks-p%glb%ghc, p%of(id)%loc%ke+p%glb%ghc
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        
            call p%of(id)%loc%ccd%x%solve_fixed_central(15.0_8/16.0_8, p%of(id)%loc%phi%now(:,j,k),&
                                                                      &p%of(id)%loc%normals%x%now(:,j,k),&
                                                                      &p%of(id)%loc%normals%xx%now(:,j,k) )
                                                                    
        enddo
        enddo
        
        do k = p%of(id)%loc%ks-p%glb%ghc, p%of(id)%loc%ke+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        
            call p%of(id)%loc%ccd%y%solve_fixed_central(15.0_8/16.0_8, p%of(id)%loc%phi%now(i,:,k),&
                                                                      &p%of(id)%loc%normals%y%now(i,:,k),&
                                                                      &p%of(id)%loc%normals%yy%now(i,:,k) )
                                                                            
        enddo
        enddo
 
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        
            call p%of(id)%loc%ccd%z%solve_fixed_central(15.0_8/16.0_8, p%of(id)%loc%phi%now(i,j,:),&
                                                                      &p%of(id)%loc%normals%z%now(i,j,:),&
                                                                      &p%of(id)%loc%normals%zz%now(i,j,:) )
                                                                            
        enddo
        enddo       
                
        !===========================================
        
        do k = p%of(id)%loc%ks-p%glb%ghc, p%of(id)%loc%ke+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        
            call p%of(id)%loc%ccd%y%solve_fixed_central(15.0_8/16.0_8, p%of(id)%loc%normals%x%now(i,:,k),&
                                                                      &p%of(id)%loc%normals%xy%now(i,:,k),&
                                                                      &p%of(id)%loc%normals%curv%now(i,:,k) )
                                                                            
        enddo
        enddo

        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        
            call p%of(id)%loc%ccd%z%solve_fixed_central(15.0_8/16.0_8, p%of(id)%loc%normals%x%now(i,j,:),&
                                                                      &p%of(id)%loc%normals%xz%now(i,j,:),&
                                                                      &p%of(id)%loc%normals%curv%now(i,j,:) )
                                                                            
        enddo
        enddo 

        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        
            call p%of(id)%loc%ccd%z%solve_fixed_central(15.0_8/16.0_8, p%of(id)%loc%normals%y%now(i,j,:),&
                                                                      &p%of(id)%loc%normals%yz%now(i,j,:),&
                                                                      &p%of(id)%loc%normals%curv%now(i,j,:) )
                                                                            
        enddo
        enddo

        !===========================================
        
        call p%of(id)%bc(0,p%of(id)%loc%normals%x%now)
        call p%of(id)%bc(0,p%of(id)%loc%normals%xx%now)
        
        call p%of(id)%bc(0,p%of(id)%loc%normals%y%now)
        call p%of(id)%bc(0,p%of(id)%loc%normals%yy%now)
       
        call p%of(id)%bc(0,p%of(id)%loc%normals%z%now)
        call p%of(id)%bc(0,p%of(id)%loc%normals%zz%now)
         
        call p%of(id)%bc(0,p%of(id)%loc%normals%xy%now)
        call p%of(id)%bc(0,p%of(id)%loc%normals%xz%now)
        call p%of(id)%bc(0,p%of(id)%loc%normals%yz%now)
                                        
    !$omp end parallel
    
end subroutine

subroutine manager_surface_norms_sec(p)
implicit none
class(manager) :: p
integer :: id,I,J,k  
    
    !$omp parallel private(id,i,j,k), num_threads(p%glb%threads)
    
        id=0
        !$ id = omp_get_thread_num()
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            p%of(id)%loc%normals%x%now(i,j,k) = 0.5d0*( p%of(id)%loc%phi%now(i+1,j,k) - p%of(id)%loc%phi%now(i-1,j,k) )/p%glb%dx
            p%of(id)%loc%normals%y%now(i,j,k) = 0.5d0*( p%of(id)%loc%phi%now(i,j+1,k) - p%of(id)%loc%phi%now(i,j-1,k) )/p%glb%dy
            p%of(id)%loc%normals%z%now(i,j,k) = 0.5d0*( p%of(id)%loc%phi%now(i,j,k+1) - p%of(id)%loc%phi%now(i,j,k-1) )/p%glb%dz
        
            p%of(id)%loc%normals%xx%now(i,j,k) = ( p%of(id)%loc%phi%now(i+1,j,k) - 2.0d0*p%of(id)%loc%phi%now(i,j,k) - p%of(id)%loc%phi%now(i-1,j,k) )/p%glb%dx**2.0d0
            p%of(id)%loc%normals%yy%now(i,j,k) = ( p%of(id)%loc%phi%now(i,j+1,k) - 2.0d0*p%of(id)%loc%phi%now(i,j,k) - p%of(id)%loc%phi%now(i,j-1,k) )/p%glb%dy**2.0d0
            p%of(id)%loc%normals%zz%now(i,j,k) = ( p%of(id)%loc%phi%now(i,j,k+1) - 2.0d0*p%of(id)%loc%phi%now(i,j,k) - p%of(id)%loc%phi%now(i,j,k-1) )/p%glb%dz**2.0d0
            
            p%of(id)%loc%normals%xy%now(i,j,k) = ( p%of(id)%loc%phi%now(i+1,j+1,k)+p%of(id)%loc%phi%now(i-1,j-1,k) &
                                               & - p%of(id)%loc%phi%now(i-1,j+1,k)-p%of(id)%loc%phi%now(i+1,j-1,k) )/(4.0d0*p%glb%dx*p%glb%dy)
            p%of(id)%loc%normals%xz%now(i,j,k) = ( p%of(id)%loc%phi%now(i+1,j,k+1)+p%of(id)%loc%phi%now(i-1,j,k-1) &
                                               & - p%of(id)%loc%phi%now(i-1,j,k+1)-p%of(id)%loc%phi%now(i+1,j,k-1) )/(4.0d0*p%glb%dx*p%glb%dz)
            p%of(id)%loc%normals%yz%now(i,j,k) = ( p%of(id)%loc%phi%now(i,j+1,k+1)+p%of(id)%loc%phi%now(i,j-1,k-1) &
                                               & - p%of(id)%loc%phi%now(i,j+1,k-1)-p%of(id)%loc%phi%now(i,j-1,k+1) )/(4.0d0*p%glb%dz*p%glb%dy)
        end do
        end do
        end do

        !===========================================
        
        call p%of(id)%bc(0,p%of(id)%loc%normals%x%now)
        call p%of(id)%bc(0,p%of(id)%loc%normals%xx%now)
        
        call p%of(id)%bc(0,p%of(id)%loc%normals%y%now)
        call p%of(id)%bc(0,p%of(id)%loc%normals%yy%now)
       
        call p%of(id)%bc(0,p%of(id)%loc%normals%z%now)
        call p%of(id)%bc(0,p%of(id)%loc%normals%zz%now)
         
        call p%of(id)%bc(0,p%of(id)%loc%normals%xy%now)
        call p%of(id)%bc(0,p%of(id)%loc%normals%xz%now)
        call p%of(id)%bc(0,p%of(id)%loc%normals%yz%now)
                                    
    !$omp end parallel
    
end subroutine

subroutine manager_curv(p)
implicit none
class(manager) :: p
integer :: id,i,j,k
real(8) :: fx,fxx,fy,fyy,fz,fzz,fxy,fxz,fyz

    !$omp parallel private(id,i,j,k,fx,fxx,fy,fyy,fz,fzz,fxy,fxz,fyz), num_threads(p%glb%threads)
    
        id=0
        !$ id = omp_get_thread_num()
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            fx  = p%of(id)%loc%normals%x%now(i,j,k)
            fy  = p%of(id)%loc%normals%y%now(i,j,k)
            fz  = p%of(id)%loc%normals%z%now(i,j,k)
            
            fxx = p%of(id)%loc%normals%xx%now(i,j,k)
            fyy = p%of(id)%loc%normals%yy%now(i,j,k)
            fzz = p%of(id)%loc%normals%zz%now(i,j,k)
            
            fxy = p%of(id)%loc%normals%xy%now(i,j,k)
            fxz = p%of(id)%loc%normals%xz%now(i,j,k)
            fyz = p%of(id)%loc%normals%yz%now(i,j,k)
            
            p%of(id)%loc%normals%curv%now(i,j,k) = (fxx*(fy**2+fz**2)+fyy*(fx**2+fz**2)+fzz*(fx**2+fy**2) &
                                            & - 2.0d0*(fxy*fx*fy+fxz*fx*fz+fyz*fy*fz)) / (fx**2+fy**2+fz**2+1.0d-12)**1.5d0         
        end do
        end do
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%normals%curv%now)
    
    !$omp end parallel 

end subroutine

subroutine manager_node_vel(p)
implicit none
class(manager) :: p
integer :: id,i,j,k

    !$omp parallel private(id,i,j,k), num_threads(p%glb%threads)
    
        id=0
        !$ id = omp_get_thread_num()
        
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
            
    !$omp end parallel

end subroutine

subroutine manager_vortex_identification(p)
implicit none
class(manager) :: p
integer :: id,i,j,k
real(8) :: ux,uy,uz,vx,vy,vz,wx,wy,wz,l1,l2,l3
real(8) :: R

    !$omp parallel private(id,i,j,k,ux,uy,uz,vx,vy,vz,wx,wy,wz,l1,l2,l3,R), num_threads(p%glb%threads)
    
        id=0
        !$ id = omp_get_thread_num()
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            ux = 0.5d0*( p%of(id)%loc%nvel%x%now(i+1,j,k)-p%of(id)%loc%nvel%x%now(i-1,j,k) )/p%glb%dx
            uy = 0.5d0*( p%of(id)%loc%nvel%x%now(i,j+1,k)-p%of(id)%loc%nvel%x%now(i,j-1,k) )/p%glb%dy
            uz = 0.5d0*( p%of(id)%loc%nvel%x%now(i,j,k+1)-p%of(id)%loc%nvel%x%now(i,j,k-1) )/p%glb%dz
            
            vx = 0.5d0*( p%of(id)%loc%nvel%y%now(i+1,j,k)-p%of(id)%loc%nvel%y%now(i-1,j,k) )/p%glb%dx
            vy = 0.5d0*( p%of(id)%loc%nvel%y%now(i,j+1,k)-p%of(id)%loc%nvel%y%now(i,j-1,k) )/p%glb%dy
            vz = 0.5d0*( p%of(id)%loc%nvel%y%now(i,j,k+1)-p%of(id)%loc%nvel%y%now(i,j,k-1) )/p%glb%dz
            
            wx = 0.5d0*( p%of(id)%loc%nvel%z%now(i+1,j,k)-p%of(id)%loc%nvel%z%now(i-1,j,k) )/p%glb%dx
            wy = 0.5d0*( p%of(id)%loc%nvel%z%now(i,j+1,k)-p%of(id)%loc%nvel%z%now(i,j-1,k) )/p%glb%dy
            wz = 0.5d0*( p%of(id)%loc%nvel%z%now(i,j,k+1)-p%of(id)%loc%nvel%z%now(i,j,k-1) )/p%glb%dz
                    
            call eigenvalue_of_tensor(Ux,Uy,Uz,Vx,Vy,Vz,Wx,Wy,Wz,L1,L2,L3)
            
            R = Ux*Vy*Wz - Ux*Vz*Wy - Uy*Vx*Wz + Uy*Vz*Wx + Uz*Vx*Wy - Uz*Vy*Wx
            
            p%of(id)%loc%q_cri%now(i,j,k) = -0.5d0*(L1+L2+L3)
            p%of(id)%loc%q_cri%now(i,j,k) = (p%of(id)%loc%q_cri%now(i,j,k)/3.0d0)**3.0d0 + (R/2.0d0)**2.0d0
            p%of(id)%loc%lam2_cri%now(i,j,k) = L1+L2+L3 - MAX(MAX(L1,L2),L3) - MIN(MIN(L1,L2),L3)
            
            p%of(id)%loc%vort%x%now(i,j,k) = Wy-Vz
            p%of(id)%loc%vort%y%now(i,j,k) = Uz-Wx
            p%of(id)%loc%vort%z%now(i,j,k) = Vx-Uy
            
            p%of(id)%loc%lamb%x%now(i,j,k) = p%of(id)%loc%nvel%y%now(i,j,k)*p%of(id)%loc%vort%z%now(i,j,k) &
                                        &  - p%of(id)%loc%nvel%z%now(i,j,k)*p%of(id)%loc%vort%y%now(i,j,k)
                                        
            p%of(id)%loc%lamb%y%now(i,j,k) = p%of(id)%loc%nvel%z%now(i,j,k)*p%of(id)%loc%vort%x%now(i,j,k) &
                                        &  - p%of(id)%loc%nvel%x%now(i,j,k)*p%of(id)%loc%vort%z%now(i,j,k)
                                        
            p%of(id)%loc%lamb%z%now(i,j,k) = p%of(id)%loc%nvel%x%now(i,j,k)*p%of(id)%loc%vort%y%now(i,j,k) &
                                        &  - p%of(id)%loc%nvel%y%now(i,j,k)*p%of(id)%loc%vort%x%now(i,j,k)
                                        
            p%of(id)%loc%helicity%now(i,j,k) = p%of(id)%loc%nvel%x%now(i,j,k)*p%of(id)%loc%vort%x%now(i,j,k) &
                                            &+ p%of(id)%loc%nvel%y%now(i,j,k)*p%of(id)%loc%vort%y%now(i,j,k) &
                                            &+ p%of(id)%loc%nvel%z%now(i,j,k)*p%of(id)%loc%vort%z%now(i,j,k) 
        end do
        end do
        end do
        
        
    !$omp end parallel

end subroutine

SUBROUTINE eigenvalue_of_tensor(Ux,Uy,Uz,Vx,Vy,Vz,Wx,Wy,Wz,L1,L2,L3)
implicit none
real(8) :: Ux,Uy,Uz,Vx,Vy,Vz,Wx,Wy,Wz
real(8) :: a, b, c, d, e, f, p, q
real(8) :: cfB, cfC, cfD
real(8) :: theta, L1, L2, L3, eps

 eps = 1.0d-12

 a = Ux**2.0d0 + Uy*Vx + Uz*Wx
 b = (Ux*Uy)/2.0d0 + (Ux*Vx)/2.0d0 + (Uy*Vy)/2.0d0 + (Vx*Vy)/2.0d0 + (Uz*Wy)/2.0d0 + (Vz*Wx)/2.0d0
 c = (Ux*Uz)/2.0d0 + (Uy*Vz)/2.0d0 + (Ux*Wx)/2.0d0 + (Uz*Wz)/2.0d0 + (Vx*Wy)/2.0d0 + (Wx*Wz)/2.0d0
 d = Vy**2.0d0 + Uy*Vx + Vz*Wy
 e = (Uz*Vx)/2.0d0 + (Uy*Wx)/2.0d0 + (Vy*Vz)/2.0d0 + (Vy*Wy)/2.0d0 + (Vz*Wz)/2.0d0 + (Wy*Wz)/2.0d0
 f = Wz**2.0d0 + Uz*Wx + Vz*Wy
 
 cfB = - (a + d + f)
 cfC = - (b**2.0d0 + c**2.0d0 + e**2.0d0 - a*d - a*f - d*f)
 cfD = f*b**2.0d0 - 2.0d0*b*c*e + d*c**2.0d0 + a*e**2.0d0 - a*d*f
 
 p = cfC - cfB**2.0d0/3.0
 q = cfD - cfB*cfC/3.0d0 + 2.0d0*cfB**3.0d0/27.0d0
 
 theta = dacos(-0.5*q / ( dsqrt((-p/3.0d0)**3.0d0) + eps ))
 
 L1 = -cfB/3.0 + dsqrt(-p/3.0)* 2.0*dcos(theta/3.0)
 L2 = -cfB/3.0 - dsqrt(-p/3.0)*( dcos(theta/3.0) + dsqrt(3.0_8)*dsin(theta/3.0) )
 L3 = -cfB/3.0 - dsqrt(-p/3.0)*( dcos(theta/3.0) - dsqrt(3.0_8)*dsin(theta/3.0) )
 
end subroutine

subroutine manager_plots(p)
implicit none
class(manager) :: p
integer :: ix,iy,iz,i,j,k,id
character(6) :: name
real(8) :: x,y,z

 if ( abs(p%glb%time - p%glb%pid * p%glb%t2p) > p%glb%dt ) return
 
 call p%vortex
 
 write(name,'(I2.2,"_",I3.3)')0,p%glb%pid
 open(unit=777,file="./out/"//trim(p%glb%name)//'_'//name//".vtk")  
 WRITE(777,'(A)')"# vtk DataFile Version 3.0"
 write(777,'(A)')"vtk TEST"
 WRITE(777,'(A)')"ASCII"
 WRITE(777,'(A)')"DATASET STRUCTURED_POINTS"
 WRITE(777,'(A,3I6)')"DIMENSIONS ",p%glb%node_x,p%glb%node_y,p%glb%node_z
 WRITE(777,'(A,3ES15.4)')"SPACING ",p%glb%dx, p%glb%dy, p%glb%dz  
 WRITE(777,'(A,3ES15.4)')"ORIGIN ",p%glb%xstart-p%glb%dx/2.0,p%glb%ystart-p%glb%dy/2.0,p%glb%zstart-p%glb%dz/2.0
 WRITE(777,'(A,I12)')"POINT_DATA ",(p%glb%node_x)*(p%glb%node_y)*(p%glb%node_z)
 
 100 format(F20.10)
 101 format(3F20.10)
 
 write(777,'(A)')"SCALARS C FLOAT"
 write(777,'(A)')"LOOKUP_TABLE DEFAULT"
 do iz=0, p%glb%grid_z-1
 id = iz*p%glb%grid_x*p%glb%grid_y 
 do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do iy=0, p%glb%grid_y-1
    do j = p%of(id+iy*p%glb%grid_x)%loc%js, p%of(id+iy*p%glb%grid_x)%loc%je
        do ix=0, p%glb%grid_x-1
        do i = p%of(id+iy*p%glb%grid_x+ix)%loc%is, p%of(id+iy*p%glb%grid_x+ix)%loc%ie
            write(777,100)p%of(id+iy*p%glb%grid_x+ix)%loc%phi%now(i,j,k)
        enddo
        enddo       
    enddo
    enddo   
 enddo
 enddo

 write(777,'(A)')"SCALARS Q FLOAT"
 write(777,'(A)')"LOOKUP_TABLE DEFAULT"
 do iz=0, p%glb%grid_z-1
 id = iz*p%glb%grid_x*p%glb%grid_y 
 do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do iy=0, p%glb%grid_y-1
    do j = p%of(id+iy*p%glb%grid_x)%loc%js, p%of(id+iy*p%glb%grid_x)%loc%je
        do ix=0, p%glb%grid_x-1
        do i = p%of(id+iy*p%glb%grid_x+ix)%loc%is, p%of(id+iy*p%glb%grid_x+ix)%loc%ie
            write(777,100)p%of(id+iy*p%glb%grid_x+ix)%loc%q_cri%now(i,j,k)
        enddo
        enddo       
    enddo
    enddo   
 enddo
 enddo

 write(777,'(A)')"SCALARS Helicity FLOAT"
 write(777,'(A)')"LOOKUP_TABLE DEFAULT"
 do iz=0, p%glb%grid_z-1
 id = iz*p%glb%grid_x*p%glb%grid_y 
 do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do iy=0, p%glb%grid_y-1
    do j = p%of(id+iy*p%glb%grid_x)%loc%js, p%of(id+iy*p%glb%grid_x)%loc%je
        do ix=0, p%glb%grid_x-1
        do i = p%of(id+iy*p%glb%grid_x+ix)%loc%is, p%of(id+iy*p%glb%grid_x+ix)%loc%ie
            write(777,100)p%of(id+iy*p%glb%grid_x+ix)%loc%helicity%now(i,j,k)
        enddo
        enddo       
    enddo
    enddo   
 enddo
 enddo

 write(777,'(A)')"VECTORS Velocity FLOAT"
 do iz=0, p%glb%grid_z-1
 id = iz*p%glb%grid_x*p%glb%grid_y 
 do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do iy=0, p%glb%grid_y-1
    do j = p%of(id+iy*p%glb%grid_x)%loc%js, p%of(id+iy*p%glb%grid_x)%loc%je
        do ix=0, p%glb%grid_x-1
        do i = p%of(id+iy*p%glb%grid_x+ix)%loc%is, p%of(id+iy*p%glb%grid_x+ix)%loc%ie
            write(777,101)p%of(id+iy*p%glb%grid_x+ix)%loc%nvel%x%now(i,j,k),&
                      &   p%of(id+iy*p%glb%grid_x+ix)%loc%nvel%y%now(i,j,k),&
                      &   p%of(id+iy*p%glb%grid_x+ix)%loc%nvel%z%now(i,j,k)
        enddo
        enddo       
    enddo
    enddo   
 enddo
 enddo

 write(777,'(A)')"VECTORS Vorticity FLOAT"
 do iz=0, p%glb%grid_z-1
 id = iz*p%glb%grid_x*p%glb%grid_y 
 do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do iy=0, p%glb%grid_y-1
    do j = p%of(id+iy*p%glb%grid_x)%loc%js, p%of(id+iy*p%glb%grid_x)%loc%je
        do ix=0, p%glb%grid_x-1
        do i = p%of(id+iy*p%glb%grid_x+ix)%loc%is, p%of(id+iy*p%glb%grid_x+ix)%loc%ie
            write(777,101)p%of(id+iy*p%glb%grid_x+ix)%loc%vort%x%now(i,j,k),&
                      &   p%of(id+iy*p%glb%grid_x+ix)%loc%vort%y%now(i,j,k),&
                      &   p%of(id+iy*p%glb%grid_x+ix)%loc%vort%z%now(i,j,k)
        enddo
        enddo       
    enddo
    enddo   
 enddo
 enddo

 write(777,'(A)')"VECTORS Lamb FLOAT"
 do iz=0, p%glb%grid_z-1
 id = iz*p%glb%grid_x*p%glb%grid_y 
 do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do iy=0, p%glb%grid_y-1
    do j = p%of(id+iy*p%glb%grid_x)%loc%js, p%of(id+iy*p%glb%grid_x)%loc%je
        do ix=0, p%glb%grid_x-1
        do i = p%of(id+iy*p%glb%grid_x+ix)%loc%is, p%of(id+iy*p%glb%grid_x+ix)%loc%ie
            write(777,101)p%of(id+iy*p%glb%grid_x+ix)%loc%lamb%x%now(i,j,k),&
                      &   p%of(id+iy*p%glb%grid_x+ix)%loc%lamb%y%now(i,j,k),&
                      &   p%of(id+iy*p%glb%grid_x+ix)%loc%lamb%z%now(i,j,k)
        enddo
        enddo       
    enddo
    enddo   
 enddo
 enddo 
 close(777)
 
 p%glb%pid = p%glb%pid + 1
 call p%sync
 
end subroutine

subroutine manager_plot(p)
implicit none
class(manager) :: p
integer :: ix,iy,iz,i,j,k,id
character(6) :: name
real(8) :: x,y,z

 if ( abs(p%glb%time - p%glb%pid * p%glb%t2p) > p%glb%dt ) return

 !$omp parallel private(id,i,j,k,name,x,y,z), num_threads(p%glb%threads)
 
    id=0
    !$ id = omp_get_thread_num()
    
    write(name,'(I2.2,"_",I3.3)')id,p%glb%pid
    open(unit=777+id,file="./out/"//trim(p%glb%name)//'_'//name//".vtk")
    
    WRITE(id+777,'(A)')"# vtk DataFile Version 3.0"
    write(id+777,'(A)')"vtk TEST"
    WRITE(id+777,'(A)')"ASCII"
    WRITE(id+777,'(A)')"DATASET STRUCTURED_POINTS"
    
    WRITE(id+777,'(A,3I6)')"DIMENSIONS ",(p%of(id)%loc%ie-p%of(id)%loc%is+2),(p%of(id)%loc%je-p%of(id)%loc%js+2),(p%of(id)%loc%ke-p%of(id)%loc%ks+2)
    
    WRITE(id+777,'(A,3ES15.4)')"SPACING ",p%glb%dx, p%glb%dy, p%glb%dz

    x = p%glb%x(p%of(id)%loc%is-1)
    y = p%glb%y(p%of(id)%loc%js-1)
    z = p%glb%z(p%of(id)%loc%ks-1)
    
    if( p%of(id)%loc%idx==0 ) x = x - 0.5d0*p%glb%dx
    if( p%of(id)%loc%idy==0 ) y = y - 0.5d0*p%glb%dy
    if( p%of(id)%loc%idz==0 ) z = z - 0.5d0*p%glb%dz
    
    WRITE(id+777,'(A,3ES15.4)')"ORIGIN ",x,y,z
        
    WRITE(id+777,'(A,I12)')"POINT_DATA ",(p%of(id)%loc%ie-p%of(id)%loc%is+2)*(p%of(id)%loc%je-p%of(id)%loc%js+2)*(p%of(id)%loc%ke-p%of(id)%loc%ks+2)
    
    write(id+777,'(A)')"SCALARS phi FLOAT"
    write(id+777,'(A)')"LOOKUP_TABLE DEFAULT"
    
    do k = p%of(id)%loc%ks-1, p%of(id)%loc%ke
    do j = p%of(id)%loc%js-1, p%of(id)%loc%je
    do i = p%of(id)%loc%is-1, p%of(id)%loc%ie
        write(id+777,*)p%of(id)%loc%phi%now(i,j,k)
    enddo
    enddo
    enddo
    
    close(id+777)
    
 !$omp end parallel 
 
 p%glb%pid = p%glb%pid + 1
call p%sync

end subroutine

end module tree
