module branches
use roots
use ccd_solvers
use srk6_solver
implicit none

type global
character(20) :: name
integer :: how_to_paras
integer :: threads ! numbers of CPUs 
integer :: pthreads ! number of pieces of plot
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
real(8) :: re, we, fr
real(8) :: srk6_w, srk6_tol
real(8) :: p_w1, p_w2, p_tol, p_b
REAL(8) :: ls_wid
real(8) :: L, T, U ! characteristic length, time, velocity
real(8) :: G, sigma ! gravity, surface tension
real(8) :: mu_1, mu_2, rho_1, rho_2
real(8) :: mu_12, rho_12
real(8) :: vel_div, ns_linf, ns_l2f, ppe_linf
real(8) :: btn_sf, btn_g
real(8) :: ls_adv, ls_red, ppe, ns 
real(8),dimension(:),allocatable :: x, y
end type global

type local
integer :: id
integer :: is, ie, js, je
type(dum_matrices)  :: coe
type(time_recorded) :: heavy, delta, grad, sign
type(time_recorded) :: phi, p, rho, mu
type(time_recorded_derivatives) :: normals
type(time_recorded_vec) :: vel, nvel, velsrc
type(srk6_data) :: srk6
type(CCD_ROOTS) :: ccd
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

subroutine job_init(p,id)
implicit none
class(job) :: p
integer :: id, IS, IE, JS, JE

P%LOC%ID = ID
P%LOC%IS = 1; P%LOC%IE = P%GLB%NODE_X
P%LOC%JS = 1; P%LOC%JE = P%GLB%NODE_Y
	
P%LOC%IS =  ID*P%GLB%NODE_X / P%GLB%THREADS + 1
P%LOC%IE = (ID+1)*P%GLB%NODE_X / P%GLB%THREADS

IS = P%LOC%IS - P%GLB%GHC; IE = P%LOC%IE + P%GLB%GHC
JS = P%LOC%JS - P%GLB%GHC; JE = P%LOC%JE + P%GLB%GHC

! coefficients matrices
CALL P%LOC%COE%ALLOC(IS,IE,JS,JE)

! level set method
CALL P%LOC%PHI%ALLOC(IS,IE,JS,JE)
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

!========================================================
CALL P%LOC%SRK6%ALLOC(IS,IE,JS,JE,P%GLB%DT,P%GLB%srk6_w,P%GLB%GHC)

CALL P%LOC%CCD%ALLOC(IS,IE,JS,JE,P%GLB%DX,P%GLB%DY,P%GLB%DT)

end subroutine

subroutine job_loc_init(pp)
implicit none
class(local) :: pp

	call pp%phi%init
	
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
real(8), dimension(p%loc%is-p%glb%ghc:p%loc%ie+p%glb%ghc,p%loc%js-p%glb%ghc:p%loc%je+p%glb%ghc) :: dat
integer :: order, i, j

 if( order == 0 ) then
	
	do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
	do i = 1, p%glb%ghc
		dat(p%loc%is-i,j) = dat(p%loc%is,j)
		dat(p%loc%ie+i,j) = dat(p%loc%ie,j)
	enddo
	enddo
	
	do j = 1, p%glb%ghc
	do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
		dat(i,p%loc%js-j) = dat(i,p%loc%js)
		dat(i,p%loc%je+j) = dat(i,p%loc%je)
	enddo
	enddo
	
 else if( order == 1)then

	do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
	do i = 1, p%glb%ghc
		dat(p%loc%is-i,j) = 2.0_8*dat(p%loc%is-i+1,j)-dat(p%loc%is-i+2,j)
		dat(p%loc%ie+i,j) = 2.0_8*dat(p%loc%ie+i-1,j)-dat(p%loc%ie+i-2,j)
	enddo
	enddo
	
	do j = 1, p%glb%ghc
	do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
		dat(i,p%loc%js-j) = 2.0_8*dat(i,p%loc%js-j+1) - dat(i,p%loc%js-j+2)
		dat(i,p%loc%je+j) = 2.0_8*dat(i,p%loc%je+j-1) - dat(i,p%loc%je+j-2)
	enddo
	enddo

 else if( order == 2)then

	do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
	do i = 1, p%glb%ghc
		dat(p%loc%is-i,j) = 3.0_8*dat(p%loc%is-i+1,j)-3.0_8*dat(p%loc%is-i+2,j)+dat(p%loc%is-i+3,j)
		dat(p%loc%ie+i,j) = 3.0_8*dat(p%loc%ie+i-1,j)-3.0_8*dat(p%loc%ie+i-2,j)+dat(p%loc%ie+i-3,j)
	enddo
	enddo
	
	do j = 1, p%glb%ghc
	do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
		dat(i,p%loc%js-j) = 3.0_8*dat(i,p%loc%js-j+1) - 3.0_8*dat(i,p%loc%js-j+2) + dat(i,p%loc%js-j+3)
		dat(i,p%loc%je+j) = 3.0_8*dat(i,p%loc%je+j-1) - 3.0_8*dat(i,p%loc%je+j-2) + dat(i,p%loc%je+j-3)
	enddo
	enddo
	
 else 
 
	write(*,*)" Wrong order for boundary condition, stop the program."
	stop
	
 endif
	

end subroutine

subroutine job_vel_bc(p,u,v)
! doi.org/10.1063/1.1761178
implicit none
class(job) :: p
real(8), dimension(p%loc%is-p%glb%ghc:p%loc%ie+p%glb%ghc,p%loc%js-p%glb%ghc:p%loc%je+p%glb%ghc) :: u, v
integer :: i, j


	!==========================================
	!  X-direction velocity boundary condition
	!==========================================

	if( p%loc%id == 0 )then
	
		if( p%glb%ubc(1) == 1 )then
		
			do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc		
			do i = 1, p%glb%ghc
				u(p%loc%is-i,j) = - u(p%loc%is-2+i,j)
				v(p%loc%is-i,j) = - v(p%loc%is-1+i,j)
			end do	
			u(p%loc%is-1,j) = 0.0_8
			end do
		
		else if ( p%glb%ubc(1) == 2)then

			do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc		
			do i = 1, p%glb%ghc
				u(p%loc%is-i,j) = u(p%loc%is-2+i,j)
				v(p%loc%is-i,j) = v(p%loc%is-1+i,j)
			end do	
			u(p%loc%is-1,j) = 0.0_8
			end do
					
		endif
	
	else if ( p%loc%id == p%glb%threads-1 )then
	
		if( p%glb%ubc(2) == 1 )then
		
			do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc		
			do i = 1, p%glb%ghc
				u(p%loc%ie+i,j) = - u(p%loc%ie-i,j)
				v(p%loc%ie+i,j) = - v(p%loc%ie+1-i,j)
			end do	
			u(p%loc%ie,j) = 0.0_8
			end do
		
		else if ( p%glb%ubc(2) == 2)then

			do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc		
			do i = 1, p%glb%ghc
				u(p%loc%ie+i,j) = u(p%loc%ie-i,j)
				v(p%loc%ie+i,j) = v(p%loc%ie+1-i,j)
			end do	
			u(p%loc%ie,j) = 0.0_8
			end do
					
		endif
			
	endif

	!==========================================
	!  Y-direction velocity boundary condition
	!==========================================
	
	if( p%glb%vbc(1) == 1 )then
		
		do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
		do j = 1, p%glb%ghc
			u(i,p%loc%js-j) = - u(i,p%loc%js-1+j)
			v(i,p%loc%js-j) = - v(i,p%loc%js-2+j)
		enddo
		v(i,p%loc%js-1) = 0.0_8
		enddo
				
	else if ( p%glb%vbc(1) == 2 )then

		do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
		do j = 1, p%glb%ghc
			u(i,p%loc%js-j) = u(i,p%loc%js-1+j)
			v(i,p%loc%js-j) = v(i,p%loc%js-2+j)
		enddo
		v(i,p%loc%js-1) = 0.0_8
		enddo
					
	endif

	if( p%glb%vbc(2) == 1 )then
		
		do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
		do j = 1, p%glb%ghc
			u(i,p%loc%je+j) = - u(i,p%loc%je+1-j)
			v(i,p%loc%je+j) = - v(i,p%loc%je-j)
		enddo
		v(i,p%loc%je) = 0.0_8
		enddo
				
	else if ( p%glb%vbc(2) == 2 )then

		do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
		do j = 1, p%glb%ghc
			u(i,p%loc%je+j) = u(i,p%loc%je+1-j)
			v(i,p%loc%je+j) = v(i,p%loc%je-j)
		enddo
		v(i,p%loc%je) = 0.0_8
		enddo
					
	endif
	
end subroutine

subroutine job_nvel_bc(p,u,v)
! doi.org/10.1063/1.1761178
implicit none
class(job) :: p
real(8), dimension(p%loc%is-p%glb%ghc:p%loc%ie+p%glb%ghc,p%loc%js-p%glb%ghc:p%loc%je+p%glb%ghc) :: u, v
integer :: i, j

	!==========================================
	!  X-direction velocity boundary condition
	!==========================================

	if( p%loc%id == 0 )then
	
		if( p%glb%ubc(1) == 1 )then
		
			do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc		
			do i = 1, p%glb%ghc
				u(p%loc%is-i,j) = - u(p%loc%is-1+i,j)
				v(p%loc%is-i,j) = - v(p%loc%is-1+i,j)
			end do	
			end do
		
		else if ( p%glb%ubc(1) == 2)then

			do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc		
			do i = 1, p%glb%ghc
				u(p%loc%is-i,j) = u(p%loc%is-1+i,j)
				v(p%loc%is-i,j) = v(p%loc%is-1+i,j)
			end do	
			end do
					
		endif
	
	else if ( p%loc%id == p%glb%threads-1 )then
	
		if( p%glb%ubc(2) == 1 )then
		
			do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc		
			do i = 1, p%glb%ghc
				u(p%loc%ie+i,j) = - u(p%loc%ie+1-i,j)
				v(p%loc%ie+i,j) = - v(p%loc%ie+1-i,j)
			end do	
			end do
		
		else if ( p%glb%ubc(2) == 2)then

			do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc		
			do i = 1, p%glb%ghc
				u(p%loc%ie+i,j) = u(p%loc%ie+1-i,j)
				v(p%loc%ie+i,j) = v(p%loc%ie+1-i,j)
			end do	
			end do
					
		endif
			
	endif

	!==========================================
	!  Y-direction velocity boundary condition
	!==========================================
	
	if( p%glb%vbc(1) == 1 )then
		
		do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
		do j = 1, p%glb%ghc
			u(i,p%loc%js-j) = - u(i,p%loc%js-1+j)
			v(i,p%loc%js-j) = - v(i,p%loc%js-1+j)
		enddo
		enddo
				
	else if ( p%glb%vbc(1) == 2 )then

		do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
		do j = 1, p%glb%ghc
			u(i,p%loc%js-j) = u(i,p%loc%js-1+j)
			v(i,p%loc%js-j) = v(i,p%loc%js-1+j)
		enddo
		enddo
					
	endif

	if( p%glb%vbc(2) == 1 )then
		
		do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
		do j = 1, p%glb%ghc
			u(i,p%loc%je+j) = - u(i,p%loc%je+1-j)
			v(i,p%loc%je+j) = - v(i,p%loc%je+1-j)
		enddo
		enddo
				
	else if ( p%glb%vbc(2) == 2 )then

		do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
		do j = 1, p%glb%ghc
			u(i,p%loc%je+j) = u(i,p%loc%je+1-j)
			v(i,p%loc%je+j) = v(i,p%loc%je+1-j)
		enddo
		enddo
					
	endif
	
end subroutine

subroutine job_find_stag_vel(p,u,v,us,vs)
implicit none
class(job) :: p
real(8), dimension(p%loc%is-p%glb%ghc:p%loc%ie+p%glb%ghc,p%loc%js-p%glb%ghc:p%loc%je+p%glb%ghc) :: u, v, us, vs
integer :: i, j

	!===============================================================
	!  Us_t + Us*U_x +  V*U_y = 0
	!  Vs_t +  U*V_x + Vs*V_y = 0
	!===============================================================
 
	do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
	do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
		v(i,j) = 0.25_8*( vs(i,j)+vs(i,j-1)+vs(i+1,j)+vs(i+1,j-1) )
		u(i,j) = 0.25_8*( us(i,j)+us(i,j+1)+us(i-1,j)+us(i-1,j+1) )
	end do
	end do
 
	!==========================================
	!  X-direction velocity boundary condition
	!==========================================

	if( p%loc%id == 0 )then
	
		if( p%glb%ubc(1) == 1 )then
		
			do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc		
			do i = 1, p%glb%ghc
				u(p%loc%is-i,j) = - u(p%loc%is-1+i,j)
				v(p%loc%is-i,j) = - v(p%loc%is-2+i,j)
			end do	
			v(p%loc%is-1,j) = 0.0_8
			end do
		
		else if ( p%glb%ubc(1) == 2)then

			do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc		
			do i = 1, p%glb%ghc
				u(p%loc%is-i,j) = u(p%loc%is-1+i,j)
				v(p%loc%is-i,j) = v(p%loc%is-2+i,j)
			end do	
			!v(p%loc%is-1,j) = 0.0_8
			end do
					
		endif
	
	else if ( p%loc%id == p%glb%threads-1 )then
	
		if( p%glb%ubc(2) == 1 )then
		
			do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc		
			do i = 1, p%glb%ghc
				u(p%loc%ie+i,j) = - u(p%loc%ie+1-i,j)
				v(p%loc%ie+i,j) = - v(p%loc%ie-i,j)
			end do	
			v(p%loc%ie,j) = 0.0_8
			end do
		
		else if ( p%glb%ubc(2) == 2)then

			do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc		
			do i = 1, p%glb%ghc
				u(p%loc%ie+i,j) = u(p%loc%ie+1-i,j)
				v(p%loc%ie+i,j) = v(p%loc%ie-i,j)
			end do	
			end do
					
		endif
			
	endif

	!==========================================
	!  Y-direction velocity boundary condition
	!==========================================
	
	if( p%glb%vbc(1) == 1 )then
		
		do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
		do j = 1, p%glb%ghc
			u(i,p%loc%js-j) = - u(i,p%loc%js-2+j)
		end do
		u(i,p%loc%js-1) = 0.0_8
		end do
				
	else if ( p%glb%vbc(1) == 2 )then

		do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
		do j = 1, p%glb%ghc
			u(i,p%loc%js-j) = u(i,p%loc%js-2+j)
		end do
		end do
					
	endif

	if( p%glb%vbc(2) == 1 )then
		
		do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
		do j = 1, p%glb%ghc
			u(i,p%loc%je+j) = - u(i,p%loc%je-j)
			v(i,p%loc%je+j) = - v(i,p%loc%je+1-j)
		end do
		u(i,p%loc%je) = 0.0_8
		end do
				
	else if ( p%glb%vbc(2) == 2 )then

		do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
		do j = 1, p%glb%ghc
			u(i,p%loc%je+j) = u(i,p%loc%je-j)
			v(i,p%loc%je+j) = v(i,p%loc%je+1-j)
		end do
		end do
					
	endif
	
end subroutine

end module branches


