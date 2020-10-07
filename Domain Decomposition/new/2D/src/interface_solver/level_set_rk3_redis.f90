subroutine level_set_rk3_redis(btn)
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k,btn,iter
real(8) :: time, error, timestop
integer(8) :: cpustart, cpuend

    call system_clock(cpustart)

    call level_set_redis_init(btn)
    
    iter = 0
    time = 0.0_8

    if( btn .eq. 0  )then
        timestop = 1.5d0*max(p%glb%xend-p%glb%xstart,p%glb%yend-p%glb%ystart,p%glb%zend-p%glb%zstart)
    else
        timestop = 3.0d0 * p%glb%dx
    end if
    
do

    iter = iter + 1
    time = time + p%glb%rdt

    !$omp parallel do private(i,j,k)
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            p%of(id)%loc%phi%tmp(i,j,k) = p%of(id)%loc%phi%now(i,j,k)
        end do
        end do 
        end do
        
    enddo   
    !$omp end parallel do
    
    call level_set_rk3_redis_solver(btn)
    
    error=0.0_8
    
    !$omp parallel do private(i,j,k), reduction(max:error)
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            error = max( error, abs(p%of(id)%loc%phi%tmp(i,j,k)-p%of(id)%loc%phi%now(i,j,k)) )
        end do
        end do 
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%phi%now)
        
    enddo
    !$omp end parallel do
    
    call pt%phi%sync
    
    if( time>timestop .or.  error.le.1.0d-8) exit
    
    if( mod(iter,100) .eq. 0 )then
        write(*,'("LS Init:",I8,F8.5,ES15.4)')iter,time,error
    end if
    
end do 

    call system_clock(cpuend)
    p%glb%ls_red = p%glb%ls_red + real(btn,kind=8)*real(cpuend-cpustart,kind=8)/real(p%glb%cpurate,kind=8)

end subroutine

subroutine level_set_redis_init(btn)
use all
!$ use omp_lib
implicit none
integer :: btn,id,i,j,k
real(8) :: grad 

    !$omp parallel do private(i,j,k)
    do id = 0, p%glb%threads-1
        call p%of(id)%bc(0,p%of(id)%loc%phi%now)
    enddo
    !$omp end parallel do
    
    call pt%phi%sync
    call p%ls_funs

    !$omp parallel do private(i,j,k)
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            p%of(id)%loc%sign%tmp(i,j,k) = p%of(id)%loc%sign%now(i,j,k)
        end do
        end do 
        end do
        
    enddo
    !$omp end parallel do
        
    if( btn.eq.0) call level_set_redis_stable()
        

end subroutine

subroutine level_set_redis_stable()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k,btn
real(8) :: grad

!$omp parallel do private(i,j,k)
do id = 0, p%glb%threads-1
    
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        p%of(id)%loc%sign%tmp(i,j,k) = p%of(id)%loc%phi%now(i,j,k)
    end do
    end do 
    end do

enddo
!$omp end parallel do

call level_set_redis_gradient()

grad = 0.0_8

!$omp parallel do private(i,j,k), reduction(max:grad)
do id = 0, p%glb%threads-1
                
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        grad = max( grad, p%of(id)%loc%grad%now(i,j,k) )
    end do
    end do 
    end do

enddo
!$omp end parallel do

!$omp parallel do private(i,j,k)
do id = 0, p%glb%threads-1
    
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie         
        p%of(id)%loc%phi%now(i,j,k) = p%of(id)%loc%phi%now(i,j,k) / grad            
    end do
    end do
    end do
    
    call p%of(id)%bc(0,p%of(id)%loc%phi%now)
    
enddo
!$omp end parallel do

call pt%phi%sync


end subroutine

subroutine level_set_rk3_redis_solver(btn)
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k,btn
real(8) :: src

call level_set_redis_gradient
call level_set_redis_lambda(btn)

!$omp parallel do private(i,j,k,src)
do id = 0, p%glb%threads-1
    
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
    
        p%of(id)%loc%tdata%x%l1(i,j,k) = (p%of(id)%loc%sign%tmp(i,j,k)-p%of(id)%loc%tdata%x%s1(i,j,k))*p%of(id)%loc%grad%now(i,j,k)-p%of(id)%loc%sign%tmp(i,j,k)
        
        src = p%of(id)%loc%tdata%x%l1(i,j,k)

        p%of(id)%loc%phi%now(i,j,k) = p%of(id)%loc%phi%now(i,j,k) - p%glb%rdt * src
        
    end do
    end do 
    end do
    
    call p%of(id)%bc(0,p%of(id)%loc%phi%now)

enddo
!$omp end parallel do

call pt%phi%sync

call level_set_redis_gradient
call level_set_redis_lambda(btn)

!$omp parallel do private(i,j,k,src)
do id = 0, p%glb%threads-1
    
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
            
        p%of(id)%loc%tdata%x%l2(i,j,k) = (p%of(id)%loc%sign%tmp(i,j,k)-p%of(id)%loc%tdata%x%s1(i,j,k))*p%of(id)%loc%grad%now(i,j,k)-p%of(id)%loc%sign%tmp(i,j,k)
        
        src = ( -3.0_8*p%of(id)%loc%tdata%x%l1(i,j,k)+p%of(id)%loc%tdata%x%l2(i,j,k) ) / 4.0_8
        
        p%of(id)%loc%phi%now(i,j,k) = p%of(id)%loc%phi%now(i,j,k) - p%glb%rdt * src
        
    end do
    end do 
    end do
    
    call p%of(id)%bc(0,p%of(id)%loc%phi%now)

enddo
!$omp end parallel do

call pt%phi%sync    

call level_set_redis_gradient
call level_set_redis_lambda(btn)

!$omp parallel do private(i,j,k,src)
do id = 0, p%glb%threads-1
    
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
    
        p%of(id)%loc%tdata%x%l3(i,j,k) = (p%of(id)%loc%sign%tmp(i,j,k)-p%of(id)%loc%tdata%x%s1(i,j,k))*p%of(id)%loc%grad%now(i,j,k)-p%of(id)%loc%sign%tmp(i,j,k)
        
        src = ( -p%of(id)%loc%tdata%x%l1(i,j,k)-p%of(id)%loc%tdata%x%l2(i,j,k)+8.0_8*p%of(id)%loc%tdata%x%l3(i,j,k) ) / 12.0_8
        
        p%of(id)%loc%phi%now(i,j,k) = p%of(id)%loc%phi%now(i,j,k) - p%glb%rdt * src
        
    end do
    end do 
    end do
    
    call p%of(id)%bc(0,p%of(id)%loc%phi%now)
    
enddo
!$omp end parallel do

call pt%phi%sync    
    
end subroutine

subroutine level_set_redis_gradient()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: upp,upm,ump,umm,vpp,vpm,vmp,vmm,wpp,wpm,wmp,wmm

!$omp parallel do private(i,j,k)
do id = 0, p%glb%threads-1
    
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie     
        p%of(id)%loc%normals%x%now(i,j,k) = (p%of(id)%loc%phi%now(i,j,k)-p%of(id)%loc%phi%now(i-1,j,k))/p%glb%dx
        p%of(id)%loc%normals%y%now(i,j,k) = (p%of(id)%loc%phi%now(i,j,k)-p%of(id)%loc%phi%now(i,j-1,k))/p%glb%dy
        p%of(id)%loc%normals%z%now(i,j,k) = (p%of(id)%loc%phi%now(i,j,k)-p%of(id)%loc%phi%now(i,j,k-1))/p%glb%dz        
    end do
    end do 
    end do
    
    call p%of(id)%bc(0,p%of(id)%loc%normals%x%now)
    call p%of(id)%bc(0,p%of(id)%loc%normals%y%now)
    call p%of(id)%bc(0,p%of(id)%loc%normals%z%now)

enddo
!$omp end parallel do

call pt%normals%sync

!$omp parallel do private(i,j,k,upp,upm,ump,umm,vpp,vpm,vmp,vmm,wpp,wpm,wmp,wmm)
do id = 0, p%glb%threads-1
    
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je     
        call wenojs_flux(p%of(id)%loc%normals%x%now(:,j,k),p%of(id)%loc%tdata%x%s1(:,j,k),p%of(id)%loc%tdata%x%s2(:,j,k),&
                        &p%of(id)%loc%is,p%of(id)%loc%ie,p%glb%ghc)                         
    end do 
    end do
    
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do i = p%of(id)%loc%is, p%of(id)%loc%ie 
        call wenojs_flux(p%of(id)%loc%normals%y%now(i,:,k),p%of(id)%loc%tdata%y%s1(i,:,k),p%of(id)%loc%tdata%y%s2(i,:,k),&
                        &p%of(id)%loc%js,p%of(id)%loc%je,p%glb%ghc)                         
    end do
    end do

    do j = p%of(id)%loc%js, p%of(id)%loc%je 
    do i = p%of(id)%loc%is, p%of(id)%loc%ie 
        call wenojs_flux(p%of(id)%loc%normals%z%now(i,j,:),p%of(id)%loc%tdata%z%s1(i,j,:),p%of(id)%loc%tdata%z%s2(i,j,:),&
                        &p%of(id)%loc%ks,p%of(id)%loc%ke,p%glb%ghc)                         
    end do
    end do
    
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke 
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
        upm=-MIN(p%of(id)%loc%tdata%x%s1(i,j,k),0.0_8)
        upp= MAX(p%of(id)%loc%tdata%x%s1(i,j,k),0.0_8)
        umm=-MIN(p%of(id)%loc%tdata%x%s2(i,j,k),0.0_8)
        ump= MAX(p%of(id)%loc%tdata%x%s2(i,j,k),0.0_8)
        
        vpm=-MIN(p%of(id)%loc%tdata%y%s1(i,j,k),0.0_8)
        vpp= MAX(p%of(id)%loc%tdata%y%s1(i,j,k),0.0_8)
        vmm=-MIN(p%of(id)%loc%tdata%y%s2(i,j,k),0.0_8)
        vmp= MAX(p%of(id)%loc%tdata%y%s2(i,j,k),0.0_8)
        
        wpm=-MIN(p%of(id)%loc%tdata%z%s1(i,j,k),0.0_8)
        wpp= MAX(p%of(id)%loc%tdata%z%s1(i,j,k),0.0_8)
        wmm=-MIN(p%of(id)%loc%tdata%z%s2(i,j,k),0.0_8)
        wmp= MAX(p%of(id)%loc%tdata%z%s2(i,j,k),0.0_8)
        
        if( p%of(id)%loc%sign%tmp(i,j,k) >= 0.0_8 )then
            p%of(id)%loc%grad%now(i,j,k) = dsqrt( MAX(upm,ump)**2.0d0 + MAX(vpm,vmp)**2.0d0 + MAX(wpm,wmp)**2.0d0  )
        else 
            p%of(id)%loc%grad%now(i,j,k) = dsqrt( MAX(upp,umm)**2.0d0 + MAX(vpp,vmm)**2.0d0 + MAX(wpp,wmm)**2.0d0 )
        end if
    
    end do
    end do
    end do 

enddo   
!$omp end parallel do

end subroutine

subroutine wenojs_flux(f,fp,fm,is,ie,ghc)
implicit none
integer :: i, is, ie, ghc
real(8),dimension(is-ghc:ie+ghc) :: f, fp, fm
real(8) :: a1,a2,a3,b1,b2,b3,w1,w2,w3,eps

EPS = 1.0D-10

do i = is-1, ie
    
    b1 = 13.0_8*(f(i-2)-2.0_8*f(i-1)+f(i))**2 + 3.0_8*(f(i-2)-4.0_8*f(i-1)+3.0_8*f(i))**2
    b2 = 13.0_8*(f(i-1)-2.0_8*f(i)+f(i+1))**2 + 3.0_8*(f(i-1)-f(i+1))**2
    b3 = 13.0_8*(f(i)-2.0_8*f(i+1)+f(i+2))**2 + 3.0_8*(3.0_8*f(i)-4.0_8*f(i+1)+f(i+2))**2
    
    a1 = 1.0_8/(EPS+b1)**2
    a2 = 6.0_8/(EPS+b2)**2
    a3 = 3.0_8/(EPS+b3)**2
    
    w1 = a1/(a1+a2+a3)
    w2 = a2/(a1+a2+a3)
    w3 = a3/(a1+a2+a3)
    
    fm(i) = w1/3.0_8*f(i-2) - (7.0_8*w1+w2)/6.0_8*f(i-1) + (11.0_8*w1+5.0_8*w2+2.0_8*w3)/6.0_8*f(i) &
            + (2.0_8*w2+5.0_8*w3)/6.0_8*f(i+1) - w3/6.0_8*f(i+2)
        
    b3 = 13.0_8*(f(i-1)-2.0_8*f(i)  +f(i+1))**2 + 3.0_8*(f(i-1)-4.0_8*f(i)+3.0_8*f(i+1))**2
    b2 = 13.0_8*(f(i)  -2.0_8*f(i+1)+f(i+2))**2 + 3.0_8*(f(i)-f(i+2))**2
    b1 = 13.0_8*(f(i+1)-2.0_8*f(i+2)+f(i+3))**2 + 3.0_8*(3.0_8*f(i+1)-4.0_8*f(i+2)+f(i+3))**2
    
    a1 = 1.0_8/(EPS+b1)**2
    a2 = 6.0_8/(EPS+b2)**2
    a3 = 3.0_8/(EPS+b3)**2  
    
    w1 = a1 / (a1+a2+a3)
    w2 = a2 / (a1+a2+a3)
    w3 = a3 / (a1+a2+a3)
    
    fp(i) =  w3*(-f(i-1)+5.0_8*f(i)+2.0_8*f(i+1))/6.0_8 &
            +w2*(2.0_8*f(i)+5.0_8*f(i+1)-f(i+2))/6.0_8 &
            +w1*(11.0_8*f(i+1)-7.0_8*f(i+2)+2.0_8*f(i+3))/6.0_8 
    
    
end do

end subroutine

subroutine level_set_redis_lambda(btn)
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k,btn,ii,jj,kk
real(8) :: a,b,lam

if( btn==0 )then

    !$omp parallel do private(i,j,k)
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie 
            p%of(id)%loc%tdata%x%s1(i,j,k) = 0.0d0      
        end do
        end do 
        end do
        
    enddo
    !$omp end parallel do   
    
    return
    
endif

call p%ls_funs

!$omp parallel do private(i,j,k,lam)
do id = 0, p%glb%threads-1
    
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
        lam = p%of(id)%loc%delta%now(i,j,k)
        !lam = (2.0d0*(1.0d0-p%glb%rho_12)*p%of(id)%loc%heavy%now(i,j,k)+p%glb%rho_12)*p%of(id)%loc%delta%now(i,j,k)
        
        p%of(id)%loc%tdata%x%s2(i,j,k) = lam*p%of(id)%loc%sign%tmp(i,j,k)*( p%of(id)%loc%grad%now(i,j,k) - 1.0d0 ) 
        p%of(id)%loc%tdata%x%s3(i,j,k) = p%of(id)%loc%grad%now(i,j,k)*p%of(id)%loc%delta%now(i,j,k)*lam
            
    end do 
    end do
    end do
        
    call p%of(id)%bc(0,p%of(id)%loc%tdata%x%s2)
    call p%of(id)%bc(0,p%of(id)%loc%tdata%x%s3)

enddo       
!$omp end parallel do
    
call pt%tdatax%sync

!$omp parallel do private(i,j,k,ii,jj,kk,a,b)
do id = 0, p%glb%threads-1
        
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
            
        a = 51.0d0*p%of(id)%loc%tdata%x%s2(i,j,k)
        b = 51.0d0*p%of(id)%loc%tdata%x%s3(i,j,k)
            
        do kk = -1, 1
        do jj = -1, 1
        do ii = -1, 1
            a = a + p%of(id)%loc%tdata%x%s2(i+ii,j+jj,k+kk)
            b = b + p%of(id)%loc%tdata%x%s3(i+ii,j+jj,k+kk)
        end do
        end do 
        end do

        p%of(id)%loc%tdata%x%s1(i,j,k) = 0.0d0
        if( abs(b)>1.0d-12 )p%of(id)%loc%tdata%x%s1(i,j,k) = a/b*p%of(id)%loc%delta%now(i,j,k)
        
    end do
    end do 
    end do

enddo
!$omp end parallel do

end subroutine
