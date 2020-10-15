subroutine level_set_rk3_redis(btn)
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k,btn,iter
real(8) :: time, error, timestop
integer(8) :: cpustart, cpuend

    id=0

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

    !$omp parallel do collapse(3)   
    do k = p%loc%ks, p%loc%ke
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
        p%loc%phi%tmp(i,j,k) = p%loc%phi%now(i,j,k)
    end do
    end do 
    end do
    !$omp end parallel do 
           
    call level_set_rk3_redis_solver(btn)
    
    error=0.0_8
    
    !$omp parallel do collapse(3) 
    do k = p%loc%ks, p%loc%ke
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
        error = max( error, abs(p%loc%phi%tmp(i,j,k)-p%loc%phi%now(i,j,k)) )
    end do
    end do 
    end do
    !$omp end parallel do 
    
    if( time>timestop .or.  error.le.1.0d-5) exit
    
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

id=0

call ls_funs

!$omp parallel do collapse(3) 
do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie
    p%loc%sign%tmp(i,j,k) = p%loc%sign%now(i,j,k)
end do
end do 
end do
    
if( btn.eq.0) call level_set_redis_stable()
        

end subroutine

subroutine level_set_redis_stable()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k,btn
real(8) :: grad

id=0

!$omp parallel do collapse(3) 
do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie
    p%loc%sign%tmp(i,j,k) = p%loc%phi%now(i,j,k)
end do
end do 
end do
!$omp end parallel do 

call level_set_redis_gradient()

grad = 0.0_8
!$omp parallel do collapse(3), reduction(max:grad)
do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie
    grad = max( grad, p%loc%grad%now(i,j,k) )
end do
end do 
end do
!$omp end parallel do

!$omp parallel do collapse(3) 
do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie         
    p%loc%phi%now(i,j,k) = p%loc%phi%now(i,j,k) / grad            
end do
end do
end do
!$omp end parallel do 

call bc(p%loc%phi%now)

end subroutine

subroutine level_set_rk3_redis_solver(btn)
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k,btn
real(8) :: src

id=0

call level_set_redis_gradient
call level_set_redis_lambda(btn)

!$omp parallel do collapse(3), private(src)     
do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie

    p%loc%tdata%x%l1(i,j,k) = (p%loc%sign%tmp(i,j,k)-p%loc%tdata%x%s1(i,j,k))*p%loc%grad%now(i,j,k)-p%loc%sign%tmp(i,j,k)
    
    src = p%loc%tdata%x%l1(i,j,k)

    p%loc%phi%now(i,j,k) = p%loc%phi%now(i,j,k) - p%glb%rdt * src
    
end do
end do 
end do
!$omp end parallel do 

call bc(p%loc%phi%now)

call level_set_redis_gradient
call level_set_redis_lambda(btn)

!$omp parallel do collapse(3), private(src)  
do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie
        
    p%loc%tdata%x%l2(i,j,k) = (p%loc%sign%tmp(i,j,k)-p%loc%tdata%x%s1(i,j,k))*p%loc%grad%now(i,j,k)-p%loc%sign%tmp(i,j,k)
    
    src = ( -3.0_8*p%loc%tdata%x%l1(i,j,k)+p%loc%tdata%x%l2(i,j,k) ) / 4.0_8
    
    p%loc%phi%now(i,j,k) = p%loc%phi%now(i,j,k) - p%glb%rdt * src
    
end do
end do 
end do
!$omp end parallel do 

call bc(p%loc%phi%now)
 
call level_set_redis_gradient
call level_set_redis_lambda(btn)

!$omp parallel do collapse(3), private(src)    
do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie

    p%loc%tdata%x%l3(i,j,k) = (p%loc%sign%tmp(i,j,k)-p%loc%tdata%x%s1(i,j,k))*p%loc%grad%now(i,j,k)-p%loc%sign%tmp(i,j,k)
    
    src = ( -p%loc%tdata%x%l1(i,j,k)-p%loc%tdata%x%l2(i,j,k)+8.0_8*p%loc%tdata%x%l3(i,j,k) ) / 12.0_8
    
    p%loc%phi%now(i,j,k) = p%loc%phi%now(i,j,k) - p%glb%rdt * src
    
end do
end do 
end do
!$omp end parallel do 

call bc(p%loc%phi%now)
    
end subroutine

subroutine level_set_redis_gradient()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: upp,upm,ump,umm,vpp,vpm,vmp,vmm,wpp,wpm,wmp,wmm

!$omp parallel do collapse(3) 
do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie     
    p%loc%normals%x%now(i,j,k) = (p%loc%phi%now(i,j,k)-p%loc%phi%now(i-1,j,k))/p%glb%dx
    p%loc%normals%y%now(i,j,k) = (p%loc%phi%now(i,j,k)-p%loc%phi%now(i,j-1,k))/p%glb%dy
    p%loc%normals%z%now(i,j,k) = (p%loc%phi%now(i,j,k)-p%loc%phi%now(i,j,k-1))/p%glb%dz        
end do
end do 
end do
!$omp end parallel do

call bc(p%loc%normals%x%now)
call bc(p%loc%normals%y%now)
call bc(p%loc%normals%z%now)

!$omp parallel do collapse(2)   
do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je     
    call wenojs_flux(p%loc%normals%x%now(:,j,k),p%loc%tdata%x%s1(:,j,k),p%loc%tdata%x%s2(:,j,k),&
                    &p%loc%is,p%loc%ie,p%glb%ghc)                         
end do 
end do
!$omp end parallel do

!$omp parallel do collapse(2)
do k = p%loc%ks, p%loc%ke
do i = p%loc%is, p%loc%ie 
    call wenojs_flux(p%loc%normals%y%now(i,:,k),p%loc%tdata%y%s1(i,:,k),p%loc%tdata%y%s2(i,:,k),&
                    &p%loc%js,p%loc%je,p%glb%ghc)                         
end do
end do
!$omp end parallel do

!$omp parallel do collapse(2)
do j = p%loc%js, p%loc%je 
do i = p%loc%is, p%loc%ie 
    call wenojs_flux(p%loc%normals%z%now(i,j,:),p%loc%tdata%z%s1(i,j,:),p%loc%tdata%z%s2(i,j,:),&
                    &p%loc%ks,p%loc%ke,p%glb%ghc)                         
end do
end do
!$omp end parallel do 

!$omp parallel do collapse(3), private(upp,upm,ump,umm,vpp,vpm,vmp,vmm,wpp,wpm,wmp,wmm)
do k = p%loc%ks, p%loc%ke 
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie
    
    upm=-MIN(p%loc%tdata%x%s1(i,j,k),0.0_8)
    upp= MAX(p%loc%tdata%x%s1(i,j,k),0.0_8)
    umm=-MIN(p%loc%tdata%x%s2(i,j,k),0.0_8)
    ump= MAX(p%loc%tdata%x%s2(i,j,k),0.0_8)
    
    vpm=-MIN(p%loc%tdata%y%s1(i,j,k),0.0_8)
    vpp= MAX(p%loc%tdata%y%s1(i,j,k),0.0_8)
    vmm=-MIN(p%loc%tdata%y%s2(i,j,k),0.0_8)
    vmp= MAX(p%loc%tdata%y%s2(i,j,k),0.0_8)
    
    wpm=-MIN(p%loc%tdata%z%s1(i,j,k),0.0_8)
    wpp= MAX(p%loc%tdata%z%s1(i,j,k),0.0_8)
    wmm=-MIN(p%loc%tdata%z%s2(i,j,k),0.0_8)
    wmp= MAX(p%loc%tdata%z%s2(i,j,k),0.0_8)
    
    if( p%loc%sign%tmp(i,j,k) >= 0.0_8 )then
        p%loc%grad%now(i,j,k) = dsqrt( MAX(upm,ump)**2.0d0 + MAX(vpm,vmp)**2.0d0 + MAX(wpm,wmp)**2.0d0  )
    else 
        p%loc%grad%now(i,j,k) = dsqrt( MAX(upp,umm)**2.0d0 + MAX(vpp,vmm)**2.0d0 + MAX(wpp,wmm)**2.0d0 )
    end if

end do
end do
end do 
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

    !$omp parallel do collapse(3)
    do k = p%loc%ks, p%loc%ke
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie 
        p%loc%tdata%x%s1(i,j,k) = 0.0d0      
    end do
    end do 
    end do
    !$omp end parallel do
        
    return
    
endif

call ls_funs

!$omp parallel do collapse(3), private(lam)  
do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie
    
    !lam = p%loc%delta%now(i,j,k)
    lam = (2.0d0*(1.0d0-p%glb%rho_12)*p%loc%heavy%now(i,j,k)+p%glb%rho_12)*p%loc%delta%now(i,j,k)
    
    p%loc%tdata%x%s2(i,j,k) = lam*p%loc%sign%tmp(i,j,k)*( p%loc%grad%now(i,j,k) - 1.0d0 ) 
    p%loc%tdata%x%s3(i,j,k) = p%loc%grad%now(i,j,k)*p%loc%delta%now(i,j,k)*lam
        
end do 
end do
end do
!$omp end parallel do

call bc(p%loc%tdata%x%s2)
call bc(p%loc%tdata%x%s3)

!$omp parallel do collapse(3), private(ii,jj,kk,a,b)     
do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie
        
    a = 51.0d0*p%loc%tdata%x%s2(i,j,k)
    b = 51.0d0*p%loc%tdata%x%s3(i,j,k)
        
    do kk = -1, 1
    do jj = -1, 1
    do ii = -1, 1
        a = a + p%loc%tdata%x%s2(i+ii,j+jj,k+kk)
        b = b + p%loc%tdata%x%s3(i+ii,j+jj,k+kk)
    end do
    end do 
    end do

    p%loc%tdata%x%s1(i,j,k) = 0.0d0
    if( abs(b)>1.0d-12 )p%loc%tdata%x%s1(i,j,k) = a/b*p%loc%delta%now(i,j,k)
    
end do
end do 
end do
!$omp end parallel do

end subroutine
