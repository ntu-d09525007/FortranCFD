subroutine ppe_mg_solver(input)
use all
!$ use omp_lib
implicit none
integer :: iter, level, id, input
integer :: i,j,k
integer(8) :: cpustart, cpuend


call system_clock(cpustart)

call ppe_mg_solver_src(input)
     
do iter = 1, 200

    p%glb%piter=p%glb%piter+1
    
    call multigrid_residual(1,.true.)
    
    if(iter==1)then
        call multigrid_full_V_cycle(.false.,5)
    else
        call multigrid_v_cycle(.false.)
    endif
    
    !call multigrid_W_cycle(.false.)
    
    call multigrid_residual(1,.false.)
    
    if( p%of(0)%loc%mg(1)%l2norm .lt. p%glb%p_tol )exit
    
    if(mod(iter,50).eq.0)write(*,*)"Final:",iter,p%of(0)%loc%mg(1)%l2norm
    
enddo

p%glb%ppe_linf = p%of(0)%loc%mg(1)%l2norm

!$omp parallel do private(i,j,k)
do id = 0, p%glb%threads-1
    
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie 
        p%of(id)%loc%p%now(i,j,k) = p%of(id)%loc%mg(1)%sol(i-p%of(id)%loc%is+1,&
                                                        &  j-p%of(id)%loc%js+1,&
                                                        &  k-p%of(id)%loc%ks+1)
                                                            
    enddo
    enddo
    enddo
    
    call p%of(id)%bc(0,p%of(id)%loc%p%now)

enddo    
!$omp end parallel do

call pt%p%sync
call ppe_mg_correction

call system_clock(cpuend)
p%glb%ppe = p%glb%ppe + real(cpuend-cpustart,kind=8)/real(p%glb%cpurate,kind=8)
    
end subroutine

subroutine ppe_mg_solver_init
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k

    !$omp parallel do private(i,j,k)
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            p%of(id)%loc%coe%r(i,j,k) = 1.0d0 - 2.0d0*p%of(id)%glb%rho_12 &
            & / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i+1,j,k))         
            p%of(id)%loc%coe%l(i,j,k) = 1.0d0 - 2.0d0*p%of(id)%glb%rho_12 &
            & / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i-1,j,k))
            p%of(id)%loc%coe%f(i,j,k) = 1.0d0 - 2.0d0*p%of(id)%glb%rho_12 &
            & / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j+1,k))
            p%of(id)%loc%coe%b(i,j,k) = 1.0d0 - 2.0d0*p%of(id)%glb%rho_12 &
            & / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j-1,k))
            p%of(id)%loc%coe%u(i,j,k) = 1.0d0 - 2.0d0*p%of(id)%glb%rho_12 &
            & / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j,k+1))
            p%of(id)%loc%coe%d(i,j,k) = 1.0d0 - 2.0d0*p%of(id)%glb%rho_12 &
            & / (p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j,k-1))
                        
            p%of(id)%loc%coe%r(i,j,k) = p%of(id)%loc%coe%r(i,j,k) / p%glb%dx**2.0d0
            p%of(id)%loc%coe%l(i,j,k) = p%of(id)%loc%coe%l(i,j,k) / p%glb%dx**2.0d0
            p%of(id)%loc%coe%f(i,j,k) = p%of(id)%loc%coe%f(i,j,k) / p%glb%dy**2.0d0
            p%of(id)%loc%coe%b(i,j,k) = p%of(id)%loc%coe%b(i,j,k) / p%glb%dy**2.0d0
            p%of(id)%loc%coe%u(i,j,k) = p%of(id)%loc%coe%u(i,j,k) / p%glb%dz**2.0d0
            p%of(id)%loc%coe%d(i,j,k) = p%of(id)%loc%coe%d(i,j,k) / p%glb%dz**2.0d0
            
            p%of(id)%loc%coe%c(i,j,k) = - ( p%of(id)%loc%coe%r(i,j,k) + p%of(id)%loc%coe%l(i,j,k) + &
                                    &       p%of(id)%loc%coe%f(i,j,k) + p%of(id)%loc%coe%b(i,j,k) + &
                                    &       p%of(id)%loc%coe%u(i,j,k) + p%of(id)%loc%coe%d(i,j,k) )
                            
            if( i==1 )then
                p%of(id)%loc%coe%c(i,j,k)=p%of(id)%loc%coe%c(i,j,k)+p%of(id)%loc%coe%l(i,j,k)
                p%of(id)%loc%coe%l(i,j,k)=0.0d0
            endif
            
            if( i==p%glb%node_x )then
                p%of(id)%loc%coe%c(i,j,k)=p%of(id)%loc%coe%c(i,j,k)+p%of(id)%loc%coe%r(i,j,k)
                p%of(id)%loc%coe%r(i,j,k)=0.0d0
            endif
            
            if( j==1 )then
                p%of(id)%loc%coe%c(i,j,k)=p%of(id)%loc%coe%c(i,j,k)+p%of(id)%loc%coe%b(i,j,k)
                p%of(id)%loc%coe%b(i,j,k)=0.0d0
            endif
            
            if( j==p%glb%node_y )then
                p%of(id)%loc%coe%c(i,j,k)=p%of(id)%loc%coe%c(i,j,k)+p%of(id)%loc%coe%f(i,j,k)
                p%of(id)%loc%coe%f(i,j,k)=0.0d0
            endif
            
            if( k==1 )then
                p%of(id)%loc%coe%c(i,j,k)=p%of(id)%loc%coe%c(i,j,k)+p%of(id)%loc%coe%d(i,j,k)
                p%of(id)%loc%coe%d(i,j,k)=0.0d0
            endif
            
            if( k==p%glb%node_z )then
                p%of(id)%loc%coe%c(i,j,k)=p%of(id)%loc%coe%c(i,j,k)+p%of(id)%loc%coe%u(i,j,k)
                p%of(id)%loc%coe%u(i,j,k)=0.0d0
            endif
            
        end do
        end do
        end do
        
    enddo
    !$omp end parallel do
    
end subroutine

subroutine ppe_mg_solver_src(iter)
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k,iter
real(8) :: pi,x,y

!$omp parallel do private(i,j,k)
do id = 0, p%glb%threads-1
    
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
        p%of(id)%loc%coe%src(i,j,k) = ( ( p%of(id)%loc%vel%x%now(i,j,k) - p%of(id)%loc%vel%x%now(i-1,j,k) ) / p%glb%dx + &
                                    &   ( p%of(id)%loc%vel%y%now(i,j,k) - p%of(id)%loc%vel%y%now(i,j-1,k) ) / p%glb%dy + & 
                                    &   ( p%of(id)%loc%vel%z%now(i,j,k) - p%of(id)%loc%vel%z%now(i,j,k-1) ) / p%glb%dz ) / p%glb%dt   
    enddo
    enddo
    enddo
    
    do k = 1, p%of(id)%loc%mg(1)%nz
    do j = 1, p%of(id)%loc%mg(1)%ny
    do i = 1, p%of(id)%loc%mg(1)%nx
        
        p%of(id)%loc%mg(1)%src(i,j,k) = p%of(id)%loc%coe%src(i+p%of(id)%loc%is-1,j+p%of(id)%loc%js-1,k+p%of(id)%loc%ks-1)
        p%of(id)%loc%mg(1)%sol(i,j,k) = p%of(id)%loc%p%tmp(i+p%of(id)%loc%is-1,j+p%of(id)%loc%js-1,k+p%of(id)%loc%ks-1)
        
    enddo
    enddo
    enddo

enddo    
!$omp end parallel do

    call multigrid_sync(1)

end subroutine

subroutine ppe_mg_correction
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: src

    !$omp parallel do private(i,j,k,src)
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            
            src = ( p%of(id)%loc%p%now(i+1,j,k) - p%of(id)%loc%p%now(i,j,k) )/p%glb%dx
            p%of(id)%loc%vel%x%now(i,j,k) = p%of(id)%loc%vel%x%now(i,j,k) - src * p%glb%dt
            
            src = ( p%of(id)%loc%p%now(i,j+1,k) - p%of(id)%loc%p%now(i,j,k) )/p%glb%dy
            p%of(id)%loc%vel%y%now(i,j,k) = p%of(id)%loc%vel%y%now(i,j,k) - src * p%glb%dt
            
            src = ( p%of(id)%loc%p%now(i,j,k+1) - p%of(id)%loc%p%now(i,j,k) )/p%glb%dz
            p%of(id)%loc%vel%z%now(i,j,k) = p%of(id)%loc%vel%z%now(i,j,k) - src * p%glb%dt
            
        enddo
        enddo
        enddo
        
        call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now,p%of(id)%loc%vel%z%now)
    
    enddo
    !$omp end parallel do
    
    call pt%vel%sync
    
end subroutine
