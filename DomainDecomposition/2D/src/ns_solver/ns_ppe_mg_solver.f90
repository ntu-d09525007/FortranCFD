subroutine ppe_mg_solver(input)
use all
!$ use omp_lib
implicit none
integer :: iter, id, initer
integer :: i,j,k,input
integer(8) :: cpustart, cpuend
real(8) :: grow

call system_clock(cpustart)

call ppe_mg_solver_init
call ppe_mg_solver_src(input)

initer=0

110 do iter = 1, 100

    p%glb%piter=p%glb%piter+1

    call multigrid_residual(1,.true.)
    
    !call multigrid_full_V_cycle(.false.,5)
    !call multigrid_w_cycle(.false.)
    call multigrid_v_cycle(.false.)   
    
    call multigrid_residual(1,.false.)

    grow = p%of(0)%loc%mg(1)%l2norm / p%of(0)%loc%mg(1)%l2norm0
    
    if( grow > 1.0d0 .or. abs(grow-1.0d0)<0.01 )then
        if( initer < 5 )then
            initer=initer+1
            call ppe_mg_get_pressure
            call ppe_mg_solver_src(0,1.5d0)
            goto 110
        else
            if( p%of(0)%loc%mg(1)%l2norm < p%glb%p_tol )then
                goto 115
            endif
        endif
    else if ( p%of(0)%loc%mg(1)%l2norm < p%glb%p_tol*0.001d0 ) then
        goto 115
    endif
   
    !if(mod(iter,50).eq.0)then
    !    write(*,*)"Final:",initer,iter,p%of(0)%loc%mg(1)%l2norm
    !endif
    
enddo

115 p%glb%ppe_linf = p%of(0)%loc%mg(1)%l2norm

call ppe_mg_get_pressure
call ppe_mg_correction

call system_clock(cpuend)
p%glb%ppe = p%glb%ppe + real(cpuend-cpustart,kind=8)/real(p%glb%cpurate,kind=8)

end subroutine

subroutine ppe_mg_solver_init
use all
!$ use omp_lib
implicit none
integer :: id,i,j

    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            p%of(id)%loc%coe%r(i,j) = 1.0d0 - 2.0d0*p%of(id)%glb%rho_12 &
            & / (p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i+1,j))         
            p%of(id)%loc%coe%l(i,j) = 1.0d0 - 2.0d0*p%of(id)%glb%rho_12 &
            & / (p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i-1,j))
            p%of(id)%loc%coe%f(i,j) = 1.0d0 - 2.0d0*p%of(id)%glb%rho_12 &
            & / (p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i,j+1))
            p%of(id)%loc%coe%b(i,j) = 1.0d0 - 2.0d0*p%of(id)%glb%rho_12 &
            & / (p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i,j-1))
                        
            p%of(id)%loc%coe%r(i,j) = p%of(id)%loc%coe%r(i,j) / p%glb%dx**2.0d0
            p%of(id)%loc%coe%l(i,j) = p%of(id)%loc%coe%l(i,j) / p%glb%dx**2.0d0
            p%of(id)%loc%coe%f(i,j) = p%of(id)%loc%coe%f(i,j) / p%glb%dy**2.0d0
            p%of(id)%loc%coe%b(i,j) = p%of(id)%loc%coe%b(i,j) / p%glb%dy**2.0d0

            p%of(id)%loc%coe%c(i,j) = - ( p%of(id)%loc%coe%r(i,j) + p%of(id)%loc%coe%l(i,j) + &
                                    &       p%of(id)%loc%coe%f(i,j) + p%of(id)%loc%coe%b(i,j) )
                            
            ! if( i==1 )then
            !     p%of(id)%loc%coe%c(i,j)=p%of(id)%loc%coe%c(i,j)+p%of(id)%loc%coe%l(i,j)
            !     p%of(id)%loc%coe%l(i,j)=0.0d0
            ! endif
            
            ! if( i==p%glb%node_x )then
            !     p%of(id)%loc%coe%c(i,j)=p%of(id)%loc%coe%c(i,j)+p%of(id)%loc%coe%r(i,j)
            !     p%of(id)%loc%coe%r(i,j)=0.0d0
            ! endif
            
            ! if( j==1 )then
            !     p%of(id)%loc%coe%c(i,j)=p%of(id)%loc%coe%c(i,j)+p%of(id)%loc%coe%b(i,j)
            !     p%of(id)%loc%coe%b(i,j)=0.0d0
            ! endif
            
            ! if( j==p%glb%node_y )then
            !     p%of(id)%loc%coe%c(i,j)=p%of(id)%loc%coe%c(i,j)+p%of(id)%loc%coe%f(i,j)
            !     p%of(id)%loc%coe%f(i,j)=0.0d0
            ! endif

        end do
        end do
        
    enddo
    !$omp end parallel do
    
end subroutine

subroutine ppe_mg_solver_src(iter,w)
use all
!$ use omp_lib
implicit none
integer :: id,i,j,iter
real(8) :: pi,x,y
real(8), optional :: w

if( present(w) )then

    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js-p%of(id)%glb%ghc, p%of(id)%loc%je+p%of(id)%glb%ghc
        do i = p%of(id)%loc%is-p%of(id)%glb%ghc, p%of(id)%loc%ie+p%of(id)%glb%ghc   
            p%of(id)%loc%p%tmp(i,j) = w*p%of(id)%loc%p%now(i,j) + (1.0d0-w)*p%of(id)%loc%p%tmp(i,j)
        enddo
        enddo
        
    enddo
    !$omp end parallel do

else

    if(iter>3)then
        !$omp parallel do private(i,j)
        do id = 0, p%glb%threads-1
            
            do j = p%of(id)%loc%js-p%of(id)%glb%ghc, p%of(id)%loc%je+p%of(id)%glb%ghc
            do i = p%of(id)%loc%is-p%of(id)%glb%ghc, p%of(id)%loc%ie+p%of(id)%glb%ghc   
                p%of(id)%loc%p%tmp(i,j) = 2.0d0*p%of(id)%loc%p%old(i,j) - p%of(id)%loc%p%old2(i,j)        
            enddo
            enddo
            
        enddo
        !$omp end parallel do
    else
        !$omp parallel do private(i,j)
        do id = 0, p%glb%threads-1
            
            do j = p%of(id)%loc%js-p%of(id)%glb%ghc, p%of(id)%loc%je+p%of(id)%glb%ghc
            do i = p%of(id)%loc%is-p%of(id)%glb%ghc, p%of(id)%loc%ie+p%of(id)%glb%ghc   
                p%of(id)%loc%p%tmp(i,j) = p%of(id)%loc%p%now(i,j) 
            enddo
            enddo
            
        enddo
        !$omp end parallel do
    endif

endif

!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
        p%of(id)%loc%coe%src(i,j) = ( ( p%of(id)%loc%vel%x%now(i,j) - p%of(id)%loc%vel%x%now(i-1,j) ) / p%glb%dx + &
                                    &   ( p%of(id)%loc%vel%y%now(i,j) - p%of(id)%loc%vel%y%now(i,j-1) ) / p%glb%dy ) / p%glb%dt
                                        
        p%of(id)%loc%coe%src(i,j) = p%of(id)%loc%coe%src(i,j) * p%of(id)%glb%rho_12 + &
                                    & p%of(id)%loc%coe%c(i,j) * p%of(id)%loc%p%tmp(i,j)   + &
                                    & p%of(id)%loc%coe%r(i,j) * p%of(id)%loc%p%tmp(i+1,j) + p%of(id)%loc%coe%l(i,j) * p%of(id)%loc%p%tmp(i-1,j) + &
                                    & p%of(id)%loc%coe%f(i,j) * p%of(id)%loc%p%tmp(i,j+1) + p%of(id)%loc%coe%b(i,j) * p%of(id)%loc%p%tmp(i,j-1) 
        
    enddo
    enddo
    
    do j = 1, p%of(id)%loc%mg(1)%ny
    do i = 1, p%of(id)%loc%mg(1)%nx
        
        p%of(id)%loc%mg(1)%src(i,j) = p%of(id)%loc%coe%src(i+p%of(id)%loc%is-1,j+p%of(id)%loc%js-1)
        p%of(id)%loc%mg(1)%sol(i,j) = p%of(id)%loc%p%tmp(i+p%of(id)%loc%is-1,j+p%of(id)%loc%js-1)
        
    enddo
    enddo

enddo    
!$omp end parallel do

call pt%mg%sync(1)

end subroutine

subroutine ppe_mg_correction
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: src

!$omp parallel do private(i,j,src)
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
        src = 1.0d0 / p%glb%rho_12 * ( p%of(id)%loc%p%now(i+1,j) - p%of(id)%loc%p%now(i,j) )/p%glb%dx
        src = src + ( 2.0d0/(p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i+1,j)) - 1.0d0/p%glb%rho_12 ) * &
        & ( p%of(id)%loc%p%tmp(i+1,j) - p%of(id)%loc%p%tmp(i,j) )/p%glb%dx
        
        p%of(id)%loc%vel%x%now(i,j) = p%of(id)%loc%vel%x%now(i,j) - src * p%glb%dt
        
        src = 1.0d0 / p%glb%rho_12 * ( p%of(id)%loc%p%now(i,j+1) - p%of(id)%loc%p%now(i,j) )/p%glb%dy
        src = src + ( 2.0d0/(p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i,j+1)) - 1.0d0/p%glb%rho_12 ) * &
        & ( p%of(id)%loc%p%tmp(i,j+1) - p%of(id)%loc%p%tmp(i,j) )/p%glb%dy
        
        p%of(id)%loc%vel%y%now(i,j) = p%of(id)%loc%vel%y%now(i,j) - src * p%glb%dt
        
    enddo
    enddo
    
enddo
!$omp end parallel do

call ns_velbc
    
end subroutine

subroutine ppe_mg_get_pressure()
use all 
!$use omp_lib
implicit none
integer :: id,i,j

!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie 
        p%of(id)%loc%p%now(i,j) = p%of(id)%loc%mg(1)%sol(i-p%of(id)%loc%is+1,j-p%of(id)%loc%js+1)
                                                            
    enddo
    enddo
    
    call p%of(id)%bc(0,p%of(id)%loc%p%now)

enddo    
!$omp end parallel do

call pt%p%sync

end subroutine
