subroutine vof_wlic_solver()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,ii,ib
real(8) :: eps,isgn,alpha,beta,w,a1,a3,a4,a5,xc
integer(8) :: cpustart, cpuend

    call system_clock(cpustart)

    eps = 1.0d-12

    call surface_norms2

    ! x direction 
    !$omp parallel do collapse(2), private(ii,ib,isgn,alpha,beta,w,a1,a3,a4,a5,xc)
    do j = p%loc%js-1, p%loc%je
    do i = p%loc%is-1, p%loc%ie
        
        if( p%loc%vel%x%old(i,j) > 0.0d0 )then
            ii = i 
            isgn = 1.0d0
        else
            ii = i+1
            isgn = 0.0d0
        endif
        
        beta = 2.3d0
        
        if( p%loc%vof%now(ii,j) > 1.0d0-eps .or. p%loc%vof%now(ii,j) < eps )then
           
           p%loc%tdata%x%l1(i,j) = p%loc%vof%now(ii,j) * p%loc%vel%x%old(i,j) * p%glb%dt
          
        else
        
            ib = max(1,ii-1)
            
            if( p%loc%vof%now(ib,j) < p%loc%vof%now(ii+1,j) )then
                alpha = 1.0d0
            else
                alpha = -1.0d0
            endif
            
            a1 = dexp( beta*( 2.0d0 * p%loc%vof%now(ii,j) - 1.0d0 ) / alpha )
            a3 = dexp( beta )
            xc = 0.5d0 / beta * dlog( (a3*a3-a1*a3)/(a1*a3-1.0d0) )
            a4 = dcosh( beta * ( isgn - p%loc%vel%x%old(i,j)*p%glb%dt/p%glb%dx - xc ) )
            a5 = dcosh( beta * ( isgn - xc ) )
            
            w  = abs(p%loc%normals%x%now(ii,j)) + abs(p%loc%normals%y%now(ii,j)) 
            w  = abs(p%loc%normals%x%now(ii,j)) / w
        
            p%loc%tdata%x%l1(i,j) = 0.5d0*( p%loc%vel%x%old(i,j)*p%glb%dt - alpha*p%glb%dx/beta*dlog(a4/a5) )
        
            p%loc%tdata%x%l1(i,j) = p%loc%tdata%x%l1(i,j)*w + ( 1.0d0 - w ) * p%loc%vof%now(ii,j) * p%loc%vel%x%old(i,j) * p%glb%dt
                                            
        endif
        
    end do
    end do
    !$omp end parallel do
    
    !$omp parallel do collapse(2)
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
        
        p%loc%vof%now(i,j) = p%loc%vof%now(i,j) - ( p%loc%tdata%x%l1(i,j) - p%loc%tdata%x%l1(i-1,j) ) / p%glb%dx + &
                           & p%loc%vof%old(i,j) * p%glb%dt * ( p%loc%vel%x%old(i,j) - p%loc%vel%x%old(i-1,j) ) / p%glb%dx
        
    end do
    end do 
    !$omp end parallel do
    
    call bc(p%loc%vof%now)

    ! y direction 
    !$omp parallel do collapse(2), private(ib,isgn,alpha,beta,w,a1,a3,a4,a5,xc)
    do j = p%loc%js-1, p%loc%je
    do i = p%loc%is-1, p%loc%ie
        
        if( p%loc%vel%y%old(i,j) > 0.0d0 )then
            ii = j 
            isgn = 1.0d0
        else
            ii = j+1
            isgn = 0.0d0
        endif
        
        beta = 2.3d0
        
        if( p%loc%vof%now(i,ii) > 1.0d0-eps .or. p%loc%vof%now(i,ii) < eps )then
           
           p%loc%tdata%y%l1(i,j) = p%loc%vof%now(i,ii) * p%loc%vel%y%old(i,j) * p%glb%dt
          
        else
        
            ib = max(1,ii-1)
            
            if( p%loc%vof%now(i,ib) < p%loc%vof%now(i,ii+1) )then
                alpha = 1.0d0
            else
                alpha = -1.0d0
            endif
            
            a1 = dexp( beta*( 2.0d0 * p%loc%vof%now(i,ii) - 1.0d0 ) / alpha )
            a3 = dexp( beta )
            xc = 0.5d0 / beta * dlog( ( a3*a3-a1*a3)/(a1*a3-1.0d0) )
            a4 = dcosh( beta * ( isgn - p%loc%vel%y%old(i,j)*p%glb%dt/p%glb%dy - xc ) )
            a5 = dcosh( beta * ( isgn - xc ) )
            
            w  = abs(p%loc%normals%x%now(i,ii)) + abs(p%loc%normals%y%now(i,ii)) 
            w  = abs(p%loc%normals%y%now(i,ii)) / w
        
            p%loc%tdata%y%l1(i,j) = 0.5d0*( p%loc%vel%y%old(i,j)*p%glb%dt - alpha*p%glb%dy/beta*dlog(a4/a5) )
        
            p%loc%tdata%y%l1(i,j) = p%loc%tdata%y%l1(i,j)*w + ( 1.0d0 - w ) * p%loc%vof%now(i,ii) * p%loc%vel%y%old(i,j) * p%glb%dt
                                            
        endif
        
    end do
    end do
    !$omp end parallel do
    
    !$omp parallel do collapse(2)
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
        
        p%loc%vof%now(i,j) = p%loc%vof%now(i,j) - ( p%loc%tdata%y%l1(i,j) - p%loc%tdata%y%l1(i,j-1) ) / p%glb%dy + &
                           & p%loc%vof%old(i,j) * p%glb%dt * ( p%loc%vel%y%old(i,j) - p%loc%vel%y%old(i,j-1) ) / p%glb%dy
        
        if(p%loc%vof%now(i,j)<eps)p%loc%vof%now(i,j)=0.0d0
        if(p%loc%vof%now(i,j)>1.0-eps)p%loc%vof%now(i,j)=1.0d0

    end do
    end do 
    !$omp end parallel do
    
    call bc(p%loc%vof%now)

    call system_clock(cpuend)
    p%glb%ls_adv = p%glb%ls_adv + real(cpuend-cpustart,kind=8) / real( p%glb%cpurate, kind=8 )
    
end subroutine

subroutine clsvof_recon()
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: time
integer(8) :: cpustart, cpuend

    call system_clock(cpustart)
    
    !$omp parallel do collapse(2)
    do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
    do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
        p%loc%phi%tmp(i,j) = p%loc%phi%now(i,j)   
        p%loc%phi%now(i,j) = 2.0d0*p%loc%vof%now(i,j) - 1.0d0
    end do
    end do
    !$omp end parallel do
    
    call level_set_redis_init(0)
    
    time = 0.0d0
    
do

    time = time + p%glb%rdt
    
    call level_set_rk3_redis_solver(0)
    
    if( time>3.0d0*p%glb%dx ) exit
    
end do 

    !$omp parallel do collapse(2)
    do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
    do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
        if( abs(p%loc%phi%tmp(i,j)) < p%glb%ls_wid )p%loc%phi%tmp(i,j)=p%loc%phi%now(i,j)
        p%loc%phi%now(i,j) = p%loc%phi%tmp(i,j)
    end do
    end do
    !$omp end parallel do

    call system_clock(cpuend)
    p%glb%ls_red = p%glb%ls_red + real(cpuend-cpustart,kind=8)/real(p%glb%cpurate,kind=8)   
    
    call level_set_rk3_redis(1)

end subroutine





