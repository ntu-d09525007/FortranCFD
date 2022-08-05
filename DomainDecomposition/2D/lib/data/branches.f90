module branches
use roots
use ccd_solvers
use time_solver
use mutligrid_root
implicit none

type global
character(20) :: name
integer :: level
integer :: method, mpls
integer :: how_to_paras
integer :: grid_x, grid_y, threads ! numbers of CPUs 
integer :: node_x, node_y
integer :: ug ! grid per unit length
integer :: ghc ! ghost cell
integer :: ubc(2), vbc(2)
integer :: pid, iter, piter
integer(8) :: cpurate
real(8) :: xstart, xend, ystart, yend
real(8) :: time, t2s, t2p ! time to stop/plot
real(8) :: dx, dy, dt, rdt
real(8) :: mass, vol, imass, ivol
real(8) :: massv, volv, imassv, ivolv
real(8) :: re, we, fr
real(8) :: t_w, t_tol
real(8) :: p_w1, p_w2, p_tol, p_b
REAL(8) :: ls_wid
real(8) :: L, T, U ! characteristic length, time, velocity
real(8) :: G, sigma ! gravity, surface tension
real(8) :: gx,gy
real(8) :: mu_1, mu_2, rho_1, rho_2
real(8) :: mu_12, rho_12
real(8) :: vel_div, vel_sdiv, ns_linf, ns_l2f, ppe_linf
real(8) :: btn_sf, btn_g
real(8) :: ls_adv, ls_red, ppe, ns, syn
real(8),dimension(:,:),allocatable :: x, y, zeros
integer,allocatable :: id(:,:)
logical :: xper, yper
end type global

type local
integer :: id, idx, idy
integer :: is, ie, js, je
type(dum_matrices)  :: coe
type(time_recorded) :: heavy, delta, grad, sign
type(time_recorded) :: phi, p, rho, mu, vof
type(time_recorded_derivatives) :: normals
type(time_recorded_vec) :: vel, nvel, velsrc
type(tsolver_data) :: tdata
type(ccd_manager) :: ccdsolvers
!-------------------------------------------
type(multigrid_root),dimension(:),allocatable :: mg
contains
procedure init => job_loc_init
end type local

type job
type(global) :: glb
type(local)  :: loc
contains
procedure init => job_init
procedure bc => job_bc
procedure velbc => job_vel_bc
procedure nvelbc => job_nvel_bc
procedure find_stag_vel => job_find_stag_vel
end type job

contains

include './branches_velbc.f90'

subroutine job_init(p,id,gx,gy)
implicit none
class(job) :: p
integer :: gx,gy,gz
integer :: id, IS, IE, JS, JE
integer :: nx, ny, nz, level

P%LOC%ID = ID

p%loc%idy = id / gx
p%loc%idx = (id - p%loc%idy*gx)

P%LOC%IS =  P%LOC%IDx * P%GLB%NODE_X / gX + 1
P%LOC%IE = (P%LOC%IDx+1)*P%GLB%NODE_X / gx

P%LOC%JS =  P%LOC%IDy * P%GLB%NODE_Y / gY + 1
P%LOC%JE = (P%LOC%IDy+1)*P%GLB%NODE_Y / gy

IF(P%LOC%IE-P%LOC%IS+1<P%GLB%GHC .OR. P%LOC%JE-P%LOC%JS+1<P%GLB%GHC )THEN
    WRITE(*,*)"Data Synchorziation will fail for this number of threads and grids settings."
    pause
ENDIF

IS = P%LOC%IS - P%GLB%GHC; IE = P%LOC%IE + P%GLB%GHC
JS = P%LOC%JS - P%GLB%GHC; JE = P%LOC%JE + P%GLB%GHC

! coefficients matrices
CALL P%LOC%COE%ALLOC(IS,IE,JS,JE)

! level set method
CALL P%LOC%PHI%ALLOC(IS,IE,JS,JE)
CALL P%LOC%VOF%ALLOC(IS,IE,JS,JE)
CALL P%LOC%heavy%ALLOC(IS,IE,JS,JE)
CALL P%LOC%sign%ALLOC(IS,IE,JS,JE)
CALL P%LOC%delta%ALLOC(IS,IE,JS,JE)
CALL P%LOC%grad%ALLOC(IS,IE,JS,JE)

CALL P%LOC%NORMALS%ALLOC(IS,IE,JS,JE)

! fluid varaible
CALL P%LOC%P%ALLOC(IS,IE,JS,JE)
CALL P%LOC%RHO%ALLOC(IS,IE,JS,JE)
CALL P%LOC%MU%ALLOC(IS,IE,JS,JE)

CALL P%LOC%VEL%ALLOC(IS,IE,JS,JE)
CALL P%LOC%NVEL%ALLOC(IS,IE,JS,JE)
CALL P%LOC%VELSRC%ALLOC(IS,IE,JS,JE)

! solvers
CALL P%LOC%tdata%ALLOC(IS,IE,JS,JE,P%GLB%DT,P%GLB%t_w,P%GLB%GHC)
CALL P%LOC%ccdsolvers%init(IS,IE,JS,JE,P%GLB%DX,P%GLB%DY,P%GLB%DT)

if( p%glb%level>0 )then

! multigrid 
allocate(p%loc%mg(p%glb%level))
nx = (p%loc%ie-p%loc%is+1)
ny = (p%loc%je-p%loc%js+1)

do level = 1, p%glb%level

    p%loc%mg(level)%idx = p%loc%idx
    p%loc%mg(level)%idy = p%loc%idy

    p%loc%mg(level)%gx = p%glb%grid_x
    p%loc%mg(level)%gy = p%glb%grid_y

    if(mod(nx,2**(level-1)).eq.0 .and. nx/2**(level-1)>0)then
        p%loc%mg(level)%nx = nx/2**(level-1)
        p%loc%mg(level)%dx = p%glb%dx*2.0**(level-1) 
    else
        write(*,'(A,I2,I3)')"Multigrid level is too high for x, current level:",level,nx,2**(level-1)
        stop
    endif

    if(mod(ny,2**(level-1)).eq.0 .and. ny/2**(level-1)>0)then
        p%loc%mg(level)%ny = ny/2**(level-1)
        p%loc%mg(level)%dy = p%glb%dy*2.0**(level-1) 
    else
        write(*,'(A,I2,I3)')"Multigrid level is too high for y, current level:",level,ny,2**(level-1)
        stop
    endif
    
    call p%loc%mg(level)%init
    
enddo

endif

end subroutine

subroutine job_loc_init(pp)
implicit none
class(local) :: pp

    call pp%phi%init
    call pp%vof%init
    
    call pp%heavy%init
    call pp%delta%init
    call pp%grad%init
    call pp%sign%init
    
    call pp%p%init
    call pp%rho%init
    call pp%mu%init
    call pp%normals%init
    
    call pp%vel%init
    call pp%nvel%init
    call pp%velsrc%init

end subroutine

subroutine job_bc(p,order,dat)
implicit none
class(job) :: p
real(8), dimension(p%loc%is-p%glb%ghc:p%loc%ie+p%glb%ghc,&
                  &p%loc%js-p%glb%ghc:p%loc%je+p%glb%ghc) :: dat
integer :: order, i, j

 if( order == 0 ) then
    
    do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
    do i = 1, p%glb%ghc
        dat(p%loc%is-i,j) = dat(p%loc%is,j)
        dat(p%loc%ie+i,j) = dat(p%loc%ie,j)
    enddo
    enddo
    
    do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
    do j = 1, p%glb%ghc
        dat(i,p%loc%js-j) = dat(i,p%loc%js)
        dat(i,p%loc%je+j) = dat(i,p%loc%je)
    enddo
    enddo
        
 else 
 
    write(*,*)" Wrong order for boundary condition, stop the program."
    stop
    
 endif
    

end subroutine

subroutine job_find_stag_vel(p,u,v,us,vs)
implicit none
class(job) :: p
real(8), dimension(p%loc%is-p%glb%ghc:p%loc%ie+p%glb%ghc,&
                  &p%loc%js-p%glb%ghc:p%loc%je+p%glb%ghc) :: u,v,us,vs
integer :: i, j

!===============================================================
!  Us_t + Us*U_x +  V*U_y +  W*U_z = 0
!  Vs_t +  U*V_x + Vs*V_y + WW*V_z = 0
!  Ws_t + UU*W_x + VV*W_y + Ws*W_z = 0
!===============================================================

do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie

    V(i,j) = 0.25d0*( Vs(i,j)+Vs(i,j-1)+Vs(i+1,j)+Vs(i+1,j-1) )
    U(i,j) = 0.25d0*( Us(i,j)+Us(i-1,j)+Us(i,j+1)+Us(i-1,j+1) )

end do
end do

end subroutine

end module branches


