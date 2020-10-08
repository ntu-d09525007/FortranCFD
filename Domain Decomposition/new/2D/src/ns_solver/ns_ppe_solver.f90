subroutine ppe_sor_init
use all
!$ use omp_lib
implicit none
integer :: id,i,j

    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            p%of(id)%loc%coe%r(i,j) = 2.0d0 / (p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i+1,j))
            p%of(id)%loc%coe%l(i,j) = 2.0d0 / (p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i-1,j))
            p%of(id)%loc%coe%f(i,j) = 2.0d0 / (p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i,j+1))
            p%of(id)%loc%coe%b(i,j) = 2.0d0 / (p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i,j-1))
  
            p%of(id)%loc%coe%r(i,j) = p%of(id)%loc%coe%r(i,j) / p%glb%dx**2.0d0
            p%of(id)%loc%coe%l(i,j) = p%of(id)%loc%coe%l(i,j) / p%glb%dx**2.0d0
            p%of(id)%loc%coe%f(i,j) = p%of(id)%loc%coe%f(i,j) / p%glb%dy**2.0d0
            p%of(id)%loc%coe%b(i,j) = p%of(id)%loc%coe%b(i,j) / p%glb%dy**2.0d0
  
            p%of(id)%loc%coe%c(i,j) = - ( p%of(id)%loc%coe%r(i,j) + p%of(id)%loc%coe%l(i,j) + &
                                    &     p%of(id)%loc%coe%f(i,j) + p%of(id)%loc%coe%b(i,j) )
                            
            if( i==1 )then
                p%of(id)%loc%coe%c(i,j)=p%of(id)%loc%coe%c(i,j)+p%of(id)%loc%coe%l(i,j)
                p%of(id)%loc%coe%l(i,j)=0.0d0
            endif
            
            if( i==p%glb%node_x )then
                p%of(id)%loc%coe%c(i,j)=p%of(id)%loc%coe%c(i,j)+p%of(id)%loc%coe%r(i,j)
                p%of(id)%loc%coe%r(i,j)=0.0d0
            endif
            
            if( j==1 )then
                p%of(id)%loc%coe%c(i,j)=p%of(id)%loc%coe%c(i,j)+p%of(id)%loc%coe%b(i,j)
                p%of(id)%loc%coe%b(i,j)=0.0d0
            endif
            
            if( j==p%glb%node_y )then
                p%of(id)%loc%coe%c(i,j)=p%of(id)%loc%coe%c(i,j)+p%of(id)%loc%coe%f(i,j)
                p%of(id)%loc%coe%f(i,j)=0.0d0
            endif
            
        end do
        end do
        
    enddo
    !$omp end parallel do
    
end subroutine

subroutine ppe_sor_solver(tol)
use all
!$ use omp_lib
implicit none
integer :: id,i,j,iter
integer(8) :: cpustart, cpuend
real(8) :: sump, err, w, pcal, tol

    call system_clock(cpustart)

    call ppe_sor_init

    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            p%of(id)%loc%coe%src(i,j) = ( ( p%of(id)%loc%vel%x%now(i,j) - p%of(id)%loc%vel%x%now(i-1,j) ) / p%glb%dx + &
                                          ( p%of(id)%loc%vel%y%now(i,j) - p%of(id)%loc%vel%y%now(i,j-1) ) / p%glb%dy ) / p%glb%dt
                                                    
        end do
        end do
       
    enddo   
    !$omp end parallel do
    
    w = p%glb%p_w1
    
do
    
    p%glb%piter=p%glb%piter+1

    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1

        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            p%of(id)%loc%p%tmp(i,j) = p%of(id)%loc%p%now(i,j)
            
            p%of(id)%loc%p%now(i,j) = p%of(id)%loc%coe%src(i,j) - p%of(id)%loc%coe%r(i,j)*p%of(id)%loc%p%now(i+1,j) &
                                                            &   - p%of(id)%loc%coe%l(i,j)*p%of(id)%loc%p%now(i-1,j) &
                                                            &   - p%of(id)%loc%coe%f(i,j)*p%of(id)%loc%p%now(i,j+1) &
                                                            &   - p%of(id)%loc%coe%b(i,j)*p%of(id)%loc%p%now(i,j-1)
                                                            
            p%of(id)%loc%p%now(i,j) = p%of(id)%loc%p%now(i,j) / p%of(id)%loc%coe%c(i,j)  
  
        end do
        end do
    
    enddo
    !$omp end parallel do
    
    sump=0.0d0
    !$omp parallel do private(i,j), reduction(+:sump)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
               
            sump = sump + p%of(id)%loc%p%now(i,j)

        end do
        end do
    
    enddo
    !$omp end parallel do
    
    sump = sump / ( p%glb%node_x * p%glb%node_y )

    err=0.0d0    
    !$omp parallel do private(i,j), reduction(max:err)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            
            p%of(id)%loc%p%now(i,j) = p%of(id)%loc%p%now(i,j) - sump

            p%of(id)%loc%p%now(i,j) = w * p%of(id)%loc%p%now(i,j) + (1.0d0-w)*p%of(id)%loc%p%tmp(i,j)

            err = max( err,abs(p%of(id)%loc%p%now(i,j)-p%of(id)%loc%p%tmp(i,j)) )

        end do
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%p%now)
     
    enddo       
    !$omp end parallel do

    call pt%p%sync
 
    if( err < tol ) exit
    !if( err < p%glb%p_b ) w = p%glb%p_w2
    if( err > 10 .and. p%glb%piter > 10000 )then
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
integer :: id,i,j
real(8) :: px,py,rho,ux,vy

    !$omp parallel do private(i,j,px,py,rho)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            
            rho = (p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i+1,j))/2.0d0
            px  = (p%of(id)%loc%p%now(i+1,j)-p%of(id)%loc%p%now(i,j)) / p%glb%dx  
            p%of(id)%loc%vel%x%now(i,j) = p%of(id)%loc%vel%x%now(i,j) - p%glb%dt*px/rho
            
            rho = (p%of(id)%loc%rho%now(i,j)+p%of(id)%loc%rho%now(i,j+1))/2.0d0
            py  = (p%of(id)%loc%p%now(i,j+1)-p%of(id)%loc%p%now(i,j)) / p%glb%dy        
            p%of(id)%loc%vel%y%now(i,j) = p%of(id)%loc%vel%y%now(i,j) - p%glb%dt*py/rho

        end do
        end do
    
        call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now)
      
    enddo       
    !$omp end parallel do
    
    call pt%vel%sync
    
end subroutine
