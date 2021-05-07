subroutine level_set_rk3_redis(btn,T)
use all
!$ use omp_lib
implicit none
integer :: id,i,j,btn,iter
real(8) :: time, error, timestop
real(8), optional :: T
integer(8) :: cpustart, cpuend

    call system_clock(cpustart)

    call level_set_redis_init(btn)
    
    iter = 0
    time = 0.0_8

    if( btn .eq. 0  )then
        timestop = 1.5d0*max(p%glb%xend-p%glb%xstart,p%glb%yend-p%glb%ystart)
    else 
        timestop = 2.5d0 * p%glb%dx
    end if

    if( present(T) ) timestop = T
    
do

    iter = iter + 1
    time = time + p%glb%rdt

    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            p%of(id)%loc%phi%tmp(i,j) = p%of(id)%loc%phi%now(i,j)
        end do 
        end do
        
    enddo   
    !$omp end parallel do
    
    call level_set_rk3_redis_solver(btn)
    
    error=0.0_8
    
    !$omp parallel do private(i,j), reduction(max:error)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            error = max( error, abs(p%of(id)%loc%phi%tmp(i,j)-p%of(id)%loc%phi%now(i,j)) )
        end do 
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%phi%now)
        
    enddo
    !$omp end parallel do
    
    call pt%phi%sync
    
    if( time>timestop .or.  error.le.1.0d-8) exit
    
    if( mod(iter,500) .eq. 0 )then
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
integer :: btn,id,i,j
real(8) :: grad 

    !$omp parallel do 
    do id = 0, p%glb%threads-1
        call p%of(id)%bc(0,p%of(id)%loc%phi%now)
    enddo
    !$omp end parallel do
    
    call pt%phi%sync
    call p%ls_funs

    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            p%of(id)%loc%sign%tmp(i,j) = p%of(id)%loc%sign%now(i,j)
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
integer :: id,i,j,btn
real(8) :: grad

!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        p%of(id)%loc%sign%tmp(i,j) = p%of(id)%loc%phi%now(i,j)
    end do 
    end do

enddo
!$omp end parallel do

call level_set_redis_gradient()

grad = 0.0_8

!$omp parallel do private(i,j), reduction(max:grad)
do id = 0, p%glb%threads-1
                
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        grad = max( grad, p%of(id)%loc%grad%now(i,j) )
    end do 
    end do

enddo
!$omp end parallel do

!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie         
        p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%phi%now(i,j) / grad            
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
integer :: id,i,j,btn
real(8) :: src

call level_set_redis_gradient
call level_set_redis_lambda(btn)

!$omp parallel do private(i,j,src)
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
    
        p%of(id)%loc%tdata%x%l1(i,j) = (p%of(id)%loc%sign%tmp(i,j)-p%of(id)%loc%tdata%x%s1(i,j))*p%of(id)%loc%grad%now(i,j)-p%of(id)%loc%sign%tmp(i,j)
        
        src = p%of(id)%loc%tdata%x%l1(i,j)

        p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%phi%now(i,j) - p%glb%rdt * src
        
    end do 
    end do
    
    call p%of(id)%bc(0,p%of(id)%loc%phi%now)

enddo
!$omp end parallel do

call pt%phi%sync

call level_set_redis_gradient
call level_set_redis_lambda(btn)

!$omp parallel do private(i,j,src)
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
            
        p%of(id)%loc%tdata%x%l2(i,j) = (p%of(id)%loc%sign%tmp(i,j)-p%of(id)%loc%tdata%x%s1(i,j))*p%of(id)%loc%grad%now(i,j)-p%of(id)%loc%sign%tmp(i,j)
        
        src = ( -3.0_8*p%of(id)%loc%tdata%x%l1(i,j)+p%of(id)%loc%tdata%x%l2(i,j) ) / 4.0_8
        
        p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%phi%now(i,j) - p%glb%rdt * src
        
    end do 
    end do
    
    call p%of(id)%bc(0,p%of(id)%loc%phi%now)

enddo
!$omp end parallel do

call pt%phi%sync    

call level_set_redis_gradient
call level_set_redis_lambda(btn)

!$omp parallel do private(i,j,src)
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
    
        p%of(id)%loc%tdata%x%l3(i,j) = (p%of(id)%loc%sign%tmp(i,j)-p%of(id)%loc%tdata%x%s1(i,j))*p%of(id)%loc%grad%now(i,j)-p%of(id)%loc%sign%tmp(i,j)
        
        src = ( -p%of(id)%loc%tdata%x%l1(i,j)-p%of(id)%loc%tdata%x%l2(i,j)+8.0_8*p%of(id)%loc%tdata%x%l3(i,j) ) / 12.0_8
        
        p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%phi%now(i,j) - p%glb%rdt * src
        
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
integer :: id,i,j
real(8) :: upp,upm,ump,umm,vpp,vpm,vmp,vmm
real(8) :: a,b

!$omp parallel do 
do id = 0, p%glb%threads-1
    call p%of(id)%bc(0,p%of(id)%loc%phi%now)
enddo
!$omp end parallel do

call pt%phi%sync

!$omp parallel do private(i,j,a,b,upp,upm,ump,umm,vpp,vpm,vmp,vmm)
do id = 0, p%glb%threads-1

    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie

        a=1.0/(12.0*p%glb%dx)*(-(p%of(id)%loc%phi%now(i-1,j)-p%of(id)%loc%phi%now(i-2,j)) &
                           +7.0*(p%of(id)%loc%phi%now(i,j)  -p%of(id)%loc%phi%now(i-1,j)) &
                           +7.0*(p%of(id)%loc%phi%now(i+1,j)-p%of(id)%loc%phi%now(i,j)) &
                               -(p%of(id)%loc%phi%now(i+2,j)-p%of(id)%loc%phi%now(i+1,j)))

        b=1.0/(12.0*p%glb%dy)*(-(p%of(id)%loc%phi%now(i,j-1)-p%of(id)%loc%phi%now(i,j-2)) &
                           +7.0*(p%of(id)%loc%phi%now(i,j)-p%of(id)%loc%phi%now(i,j-1)) &
                           +7.0*(p%of(id)%loc%phi%now(i,j+1)-p%of(id)%loc%phi%now(i,j)) &
                               -(p%of(id)%loc%phi%now(i,j+2)-p%of(id)%loc%phi%now(i,j+1)))

        p%of(id)%loc%tdata%x%s1(i,j)=a &
               +1.0/p%glb%dx*phyn((p%of(id)%loc%phi%now(i+3,j)-2.0*p%of(id)%loc%phi%now(i+2,j)+p%of(id)%loc%phi%now(i+1,j)), &
                                  (p%of(id)%loc%phi%now(i+2,j)-2.0*p%of(id)%loc%phi%now(i+1,j)+p%of(id)%loc%phi%now(i  ,j)), &
                                  (p%of(id)%loc%phi%now(i+1,j)-2.0*p%of(id)%loc%phi%now(i  ,j)+p%of(id)%loc%phi%now(i-1,j)), &
                                  (p%of(id)%loc%phi%now(i  ,j)-2.0*p%of(id)%loc%phi%now(i-1,j)+p%of(id)%loc%phi%now(i-2,j)))

        p%of(id)%loc%tdata%x%s2(i,j)=a &
               -1.0/p%glb%dx*phyn((p%of(id)%loc%phi%now(i-3,j)-2.0*p%of(id)%loc%phi%now(i-2,j)+p%of(id)%loc%phi%now(i-1,j)), &
                                  (p%of(id)%loc%phi%now(i-2,j)-2.0*p%of(id)%loc%phi%now(i-1,j)+p%of(id)%loc%phi%now(i  ,j)), &
                                  (p%of(id)%loc%phi%now(i-1,j)-2.0*p%of(id)%loc%phi%now(i  ,j)+p%of(id)%loc%phi%now(i+1,j)), &
                                  (p%of(id)%loc%phi%now(i  ,j)-2.0*p%of(id)%loc%phi%now(i+1,j)+p%of(id)%loc%phi%now(i+2,j)))

        p%of(id)%loc%tdata%y%s1(i,j)=b &
               +1.0/p%glb%dy*phyn((p%of(id)%loc%phi%now(i,j+3)-2.0*p%of(id)%loc%phi%now(i,j+2)+p%of(id)%loc%phi%now(i,j+1)), &
                                  (p%of(id)%loc%phi%now(i,j+2)-2.0*p%of(id)%loc%phi%now(i,j+1)+p%of(id)%loc%phi%now(i,j  )), &
                                  (p%of(id)%loc%phi%now(i,j+1)-2.0*p%of(id)%loc%phi%now(i,j  )+p%of(id)%loc%phi%now(i,j-1)), &
                                  (p%of(id)%loc%phi%now(i,j  )-2.0*p%of(id)%loc%phi%now(i,j-1)+p%of(id)%loc%phi%now(i,j-2)))

        p%of(id)%loc%tdata%y%s2(i,j)=b &
               -1.0/p%glb%dy*phyn((p%of(id)%loc%phi%now(i,j-3)-2.0*p%of(id)%loc%phi%now(i,j-2)+p%of(id)%loc%phi%now(i,j-1)), &
                                  (p%of(id)%loc%phi%now(i,j-2)-2.0*p%of(id)%loc%phi%now(i,j-1)+p%of(id)%loc%phi%now(i,j  )), &
                                  (p%of(id)%loc%phi%now(i,j-1)-2.0*p%of(id)%loc%phi%now(i,j  )+p%of(id)%loc%phi%now(i,j+1)), &
                                  (p%of(id)%loc%phi%now(i,j  )-2.0*p%of(id)%loc%phi%now(i,j+1)+p%of(id)%loc%phi%now(i,j+2)))

        upm=-MIN(p%of(id)%loc%tdata%x%s1(i,j),0.0_8)
        upp= MAX(p%of(id)%loc%tdata%x%s1(i,j),0.0_8)
        umm=-MIN(p%of(id)%loc%tdata%x%s2(i,j),0.0_8)
        ump= MAX(p%of(id)%loc%tdata%x%s2(i,j),0.0_8)
        
        vpm=-MIN(p%of(id)%loc%tdata%y%s1(i,j),0.0_8)
        vpp= MAX(p%of(id)%loc%tdata%y%s1(i,j),0.0_8)
        vmm=-MIN(p%of(id)%loc%tdata%y%s2(i,j),0.0_8)
        vmp= MAX(p%of(id)%loc%tdata%y%s2(i,j),0.0_8)
        
        if( p%of(id)%loc%sign%tmp(i,j) >= 0.0_8 )then
            p%of(id)%loc%grad%now(i,j) = dsqrt( MAX(upm,ump)**2.0d0 + MAX(vpm,vmp)**2.0d0 )
        else 
            p%of(id)%loc%grad%now(i,j) = dsqrt( MAX(upp,umm)**2.0d0 + MAX(vpp,vmm)**2.0d0 )
        end if
    
    end do
    end do 

enddo   
!$omp end parallel do

end subroutine

subroutine level_set_redis_lambda(btn)
use all
!$ use omp_lib
implicit none
integer :: id,i,j,btn,ii,jj
real(8) :: a,b,lam

if( btn==0 )then

    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie 
            p%of(id)%loc%tdata%x%s1(i,j) = 0.0d0      
        end do 
        end do
        
    enddo
    !$omp end parallel do   
    
    return
    
endif

call p%ls_funs

!$omp parallel do private(i,j,lam)
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
        lam = p%of(id)%loc%delta%now(i,j)
        !lam = (2.0d0*(1.0d0-p%glb%rho_12)*p%of(id)%loc%heavy%now(i,j)+p%glb%rho_12)*p%of(id)%loc%delta%now(i,j)
        
        p%of(id)%loc%tdata%x%s2(i,j) = lam*p%of(id)%loc%sign%tmp(i,j)*( p%of(id)%loc%grad%now(i,j) - 1.0d0 ) 
        p%of(id)%loc%tdata%x%s3(i,j) = p%of(id)%loc%grad%now(i,j)*p%of(id)%loc%delta%now(i,j)*lam
            
    end do
    end do
        
    call p%of(id)%bc(0,p%of(id)%loc%tdata%x%s2)
    call p%of(id)%bc(0,p%of(id)%loc%tdata%x%s3)

enddo       
!$omp end parallel do
    
call pt%tdatax%sync

!$omp parallel do private(i,j,ii,jj,a,b)
do id = 0, p%glb%threads-1
        
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
            
        a = 16.0d0*p%of(id)%loc%tdata%x%s2(i,j)
        b = 16.0d0*p%of(id)%loc%tdata%x%s3(i,j)
            
        do jj = -1, 1
        do ii = -1, 1
            a = a + p%of(id)%loc%tdata%x%s2(i+ii,j+jj)
            b = b + p%of(id)%loc%tdata%x%s3(i+ii,j+jj)
        end do
        end do

        p%of(id)%loc%tdata%x%s1(i,j) = 0.0d0
        if( abs(b)>1.0d-12 )p%of(id)%loc%tdata%x%s1(i,j) = a/b*p%of(id)%loc%delta%now(i,j)
        
    end do 
    end do

enddo
!$omp end parallel do

end subroutine

function phyn(a,b,c,d)
implicit none
real(8) :: a,b,c,d,phyn
real(8) :: is0,is1,is2,alp0,alp1,alp2,w0,w2
real(8) :: eps
eps=1.0d-6
is0=13.0d0*(a-b)**2.0d0+3.0d0*(a-3.0d0*b)**2.0d0
is1=13.0d0*(b-c)**2.0d0+3.0d0*(b+c)**2.0d0
is2=13.0d0*(c-d)**2.0d0+3.0d0*(3.0d0*c-d)**2.0d0
alp0=1.0d0/(eps+is0)**2.0d0
alp1=6.0d0/(eps+is1)**2.0d0
alp2=3.0d0/(eps+is2)**2.0d0
w0=alp0/(alp0+alp1+alp2)
w2=alp2/(alp0+alp1+alp2)
phyn=w0/3.0d0*(a-2.0d0*b+c)+(w2-0.5d0)/6.0d0*(b-2.0d0*c+d)
end function
