subroutine ppe_sor_init
use all
!$ use omp_lib
implicit none
integer :: i,j

    !$omp parallel do collapse(2)  
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
    
        p%loc%coe%r(i,j) = 2.0d0 / (p%loc%rho%now(i,j)+p%loc%rho%now(i+1,j))
        p%loc%coe%l(i,j) = 2.0d0 / (p%loc%rho%now(i,j)+p%loc%rho%now(i-1,j))
        p%loc%coe%f(i,j) = 2.0d0 / (p%loc%rho%now(i,j)+p%loc%rho%now(i,j+1))
        p%loc%coe%b(i,j) = 2.0d0 / (p%loc%rho%now(i,j)+p%loc%rho%now(i,j-1))
        
        p%loc%coe%r(i,j) = p%loc%coe%r(i,j) / p%glb%dx**2.0d0
        p%loc%coe%l(i,j) = p%loc%coe%l(i,j) / p%glb%dx**2.0d0
        p%loc%coe%f(i,j) = p%loc%coe%f(i,j) / p%glb%dy**2.0d0
        p%loc%coe%b(i,j) = p%loc%coe%b(i,j) / p%glb%dy**2.0d0
        
        p%loc%coe%c(i,j) = - ( p%loc%coe%r(i,j) + p%loc%coe%l(i,j) + &
                                &p%loc%coe%f(i,j) + p%loc%coe%b(i,j) )
                        
        if( i==1 )then
            p%loc%coe%c(i,j)=p%loc%coe%c(i,j)+p%loc%coe%l(i,j)
            p%loc%coe%l(i,j)=0.0d0
        endif
        
        if( i==p%glb%node_x )then
            p%loc%coe%c(i,j)=p%loc%coe%c(i,j)+p%loc%coe%r(i,j)
            p%loc%coe%r(i,j)=0.0d0
        endif
        
        if( j==1 )then
            p%loc%coe%c(i,j)=p%loc%coe%c(i,j)+p%loc%coe%b(i,j)
            p%loc%coe%b(i,j)=0.0d0
        endif
        
        if( j==p%glb%node_y )then
            p%loc%coe%c(i,j)=p%loc%coe%c(i,j)+p%loc%coe%f(i,j)
            p%loc%coe%f(i,j)=0.0d0
        endif
        
    end do
    end do
    !$omp end parallel do
    
end subroutine

subroutine ppe_sor_solver(tol)
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k,iter
integer(8) :: cpustart, cpuend
real(8) :: sump, err, w, pcal, tol

    call system_clock(cpustart)

    !$omp parallel do collapse(2)    
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
    
        p%loc%coe%src(i,j) = ( ( p%loc%vel%x%now(i,j) - p%loc%vel%x%now(i-1,j) ) / p%glb%dx + &
                                 ( p%loc%vel%y%now(i,j) - p%loc%vel%y%now(i,j-1) ) / p%glb%dy ) / p%glb%dt
                                                
    end do
    end do 
    !$omp end parallel do
    
    w = p%glb%p_w1
    
do
    
    p%glb%piter=p%glb%piter+1

    !$omp parallel do collapse(2)
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
    
        p%loc%p%tmp(i,j) = p%loc%p%now(i,j)
        
        p%loc%p%now(i,j) = p%loc%coe%src(i,j) - p%loc%coe%r(i,j)*p%loc%p%now(i+1,j) &
                                          &   - p%loc%coe%l(i,j)*p%loc%p%now(i-1,j) &
                                          &   - p%loc%coe%f(i,j)*p%loc%p%now(i,j+1) &
                                          &   - p%loc%coe%b(i,j)*p%loc%p%now(i,j-1)
                                                        
        p%loc%p%now(i,j) = p%loc%p%now(i,j) / p%loc%coe%c(i,j)  

    end do
    end do
    !$omp end parallel do
    
    sump=0.0d0
    !$omp parallel do collapse(2), reduction(+:sump)    
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
           
        sump = sump + p%loc%p%now(i,j)

    end do
    end do
    !$omp end parallel do
    
    sump = sump / ( p%glb%node_x * p%glb%node_y )

    err=0.0d0    
    !$omp parallel do collapse(2), reduction(max:err)  
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
        
        p%loc%p%now(i,j) = p%loc%p%now(i,j) - sump

        p%loc%p%now(i,j) = w * p%loc%p%now(i,j) + (1.0d0-w)*p%loc%p%tmp(i,j)

        err = max( err,abs(p%loc%p%now(i,j)-p%loc%p%tmp(i,j)) )

    end do
    end do
    !$omp end parallel do

    call bc(p%loc%p%now)

    if( err < tol ) exit
    !if( err < p%glb%p_b ) w = p%glb%p_w2
    if( err > 10 .and. p%glb%piter > 100000 )then
        write(*,*)"The solution can not converge in PPE :",err
        stop
    end if
    
    if( mod(p%glb%piter,5000) .eq. 0 )then
        write(*,'("PPE iter:",I8,",error:",ES15.7)')p%glb%piter,err
    endif
    
end do

    p%glb%ppe_linf = err
    
    call system_clock(cpuend)
    p%glb%ppe = p%glb%ppe + real(cpuend-cpustart,kind=8)/real(p%glb%cpurate,kind=8)
    
    call ns_momentum_correction
    
end subroutine

subroutine ns_momentum_correction
use all
!$ use omp_lib
implicit none
integer :: i,j
real(8) :: px,py,rho,ux,vy

    !$omp parallel do collapse(2), private(px,py,rho)   
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
        
        rho = (p%loc%rho%now(i,j)+p%loc%rho%now(i+1,j))/2.0d0
        px  = (p%loc%p%now(i+1,j)-p%loc%p%now(i,j)) / p%glb%dx  
        p%loc%vel%x%now(i,j) = p%loc%vel%x%now(i,j) - p%glb%dt*px/rho
        
        rho = (p%loc%rho%now(i,j)+p%loc%rho%now(i,j+1))/2.0d0
        py  = (p%loc%p%now(i,j+1)-p%loc%p%now(i,j)) / p%glb%dy        
        p%loc%vel%y%now(i,j) = p%loc%vel%y%now(i,j) - p%glb%dt*py/rho
        
    end do
    end do
    !$omp end parallel do

    call velbc(p%loc%vel%x%now,p%loc%vel%y%now)
  
end subroutine
