subroutine ns_ab_setup
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k

    call ns_ab_adv_source_sec

end subroutine

subroutine ns_ab_solver_mg
use all
implicit none
integer :: iter,initer,relax_iter
integer :: id,i,j,k

relax_iter = 0

! if( p%glb%iter<relax_iter)then
    ! call ppe_sor_init
! else
    ! call ppe_mg_solver_init
! endif

call ns_ab_adv_source
    
iter=0
p%glb%piter=0
    
do 
    
    iter=iter+1 
    
    call ns_linearize
    call ns_ab_diff_source
    call ns_ab_predictor   
 
    !if( p%glb%iter<relax_iter )then
    !    call ppe_sor_solver(1.0d-5)
    !else
        call ppe_mg_solver(iter)
    !endif


    if(iter>5)exit
        
end do

    
end subroutine

subroutine ns_ab_solver_SOR
use all
implicit none
integer :: iter,initer
integer :: id,i,j,k
real(8) :: tol

call ppe_sor_init
call ns_ab_adv_source
    
iter=0
p%glb%piter=0
    
do 
    
    iter=iter+1 
    
    call ns_linearize
    call ns_ab_diff_source
    call ns_ab_predictor
    call ppe_sor_solver(p%glb%p_tol)

    if(iter>5)exit
        
end do

end subroutine

subroutine ns_ab_predictor
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: src

!$omp parallel private(id,i,j,k,src), num_threads(p%glb%threads)
    
    id=0
    !$ id = omp_get_thread_num()
        
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
            
        src = 1.50*p%of(id)%loc%velsrc%x%now(i,j,k) - 0.5d0*p%of(id)%loc%velsrc%x%old(i,j,k) 
        src = src + p%of(id)%loc%velsrc%x%tmp(i,j,k)            
        p%of(id)%loc%vel%x%now(i,j,k) = p%of(id)%loc%vel%x%old(i,j,k) + p%glb%dt * src
        
        src = 1.50*p%of(id)%loc%velsrc%y%now(i,j,k) - 0.5d0*p%of(id)%loc%velsrc%y%old(i,j,k) 
        src = src + p%of(id)%loc%velsrc%y%tmp(i,j,k)    
        p%of(id)%loc%vel%y%now(i,j,k) = p%of(id)%loc%vel%y%old(i,j,k) + p%glb%dt * src
            
        src = 1.50*p%of(id)%loc%velsrc%z%now(i,j,k) - 0.5d0*p%of(id)%loc%velsrc%z%old(i,j,k) 
        src = src + p%of(id)%loc%velsrc%z%tmp(i,j,k)    
        p%of(id)%loc%vel%z%now(i,j,k) = p%of(id)%loc%vel%z%old(i,j,k) + p%glb%dt * src
            
    end do
    end do 
    end do
        
    call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now,p%of(id)%loc%vel%z%now)
        
!$omp end parallel
    
call pt%vel%sync
    
end subroutine

subroutine ns_ab_adv_source
implicit none

call ns_ab_adv_source_sec
!call ns_ab_adv_source_quick
    
end subroutine

subroutine ns_ab_adv_source_sec
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: u,v,w,xp,xm,yp,ym,zp,zm

!$omp parallel private(id,i,j,k), num_threads(p%glb%threads)
    
    id=0
    !$ id = omp_get_thread_num()
        
    call p%of(id)%find_stag_vel( p%of(id)%loc%tdata%x%s1, p%of(id)%loc%tdata%y%s1, p%of(id)%loc%tdata%z%s1, &
                                &p%of(id)%loc%tdata%x%s2, p%of(id)%loc%tdata%y%s2, p%of(id)%loc%tdata%z%s2, &
                                &p%of(id)%loc%vel%x%old, p%of(id)%loc%vel%y%old, p%of(id)%loc%vel%z%old )
        
!$omp end parallel
    
call pt%tdatax%sync
call pt%tdatay%sync
call pt%tdataz%sync
    
!$omp parallel private(id,i,j,k,u,v,w,xp,xm,yp,ym,zp,zm), num_threads(p%glb%threads)
        
    id=0
    !$ id = omp_get_thread_num()
        
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
        u = p%of(id)%loc%vel%x%old(i,j,k)
        v = p%of(id)%loc%tdata%y%s1(i,j,k)
        w = p%of(id)%loc%tdata%z%s1(i,j,k)
            
        xp = 0.5d0*(-p%of(id)%loc%vel%x%old(i+2,j,k)+4.0d0*p%of(id)%loc%vel%x%old(i+1,j,k)-3.0d0*p%of(id)%loc%vel%x%old(i,j,k))/p%glb%dx 
        xm = 0.5d0*( p%of(id)%loc%vel%x%old(i-2,j,k)-4.0d0*p%of(id)%loc%vel%x%old(i-1,j,k)+3.0d0*p%of(id)%loc%vel%x%old(i,j,k))/p%glb%dx
        
        yp = 0.5d0*(-p%of(id)%loc%vel%x%old(i,j+2,k)+4.0d0*p%of(id)%loc%vel%x%old(i,j+1,k)-3.0d0*p%of(id)%loc%vel%x%old(i,j,k))/p%glb%dy
        ym = 0.5d0*( p%of(id)%loc%vel%x%old(i,j-2,k)-4.0d0*p%of(id)%loc%vel%x%old(i,j-1,k)+3.0d0*p%of(id)%loc%vel%x%old(i,j,k))/p%glb%dy
            
        zp = 0.5d0*(-p%of(id)%loc%vel%x%old(i,j,k+2)+4.0d0*p%of(id)%loc%vel%x%old(i,j,k+1)-3.0d0*p%of(id)%loc%vel%x%old(i,j,k))/p%glb%dz
        zm = 0.5d0*( p%of(id)%loc%vel%x%old(i,j,k-2)-4.0d0*p%of(id)%loc%vel%x%old(i,j,k-1)+3.0d0*p%of(id)%loc%vel%x%old(i,j,k))/p%glb%dz
        
        p%of(id)%loc%velsrc%x%now(i,j,k) = - ((u+abs(u))*xm+(u-abs(u))*xp)/2.0d0 &
                                        &  - ((v+abs(v))*ym+(v-abs(v))*yp)/2.0d0 &
                                        &  - ((w+abs(w))*zm+(w-abs(w))*zp)/2.0d0  
                
        !-----------------------------------------------------------
            
        u = p%of(id)%loc%tdata%x%s1(i,j,k)
        v = p%of(id)%loc%vel%y%old(i,j,k)
        w = p%of(id)%loc%tdata%z%s2(i,j,k)
            
        xp = 0.5d0*(-p%of(id)%loc%vel%y%old(i+2,j,k)+4.0d0*p%of(id)%loc%vel%y%old(i+1,j,k)-3.0d0*p%of(id)%loc%vel%y%old(i,j,k))/p%glb%dx 
        xm = 0.5d0*( p%of(id)%loc%vel%y%old(i-2,j,k)-4.0d0*p%of(id)%loc%vel%y%old(i-1,j,k)+3.0d0*p%of(id)%loc%vel%y%old(i,j,k))/p%glb%dx
            
        yp = 0.5d0*(-p%of(id)%loc%vel%y%old(i,j+2,k)+4.0d0*p%of(id)%loc%vel%y%old(i,j+1,k)-3.0d0*p%of(id)%loc%vel%y%old(i,j,k))/p%glb%dy
        ym = 0.5d0*( p%of(id)%loc%vel%y%old(i,j-2,k)-4.0d0*p%of(id)%loc%vel%y%old(i,j-1,k)+3.0d0*p%of(id)%loc%vel%y%old(i,j,k))/p%glb%dy
            
        zp = 0.5d0*(-p%of(id)%loc%vel%y%old(i,j,k+2)+4.0d0*p%of(id)%loc%vel%y%old(i,j,k+1)-3.0d0*p%of(id)%loc%vel%y%old(i,j,k))/p%glb%dz
        zm = 0.5d0*( p%of(id)%loc%vel%y%old(i,j,k-2)-4.0d0*p%of(id)%loc%vel%y%old(i,j,k-1)+3.0d0*p%of(id)%loc%vel%y%old(i,j,k))/p%glb%dz
            
        p%of(id)%loc%velsrc%y%now(i,j,k) = - ((u+abs(u))*xm+(u-abs(u))*xp)/2.0d0 &
                                        &  - ((v+abs(v))*ym+(v-abs(v))*yp)/2.0d0 &
                                        &  - ((w+abs(w))*zm+(w-abs(w))*zp)/2.0d0
            
        !-----------------------------------------------------------
            
        u = p%of(id)%loc%tdata%x%s2(i,j,k)
        v = p%of(id)%loc%tdata%y%s2(i,j,k)
        w = p%of(id)%loc%vel%z%old(i,j,k)
            
        xp = 0.5d0*(-p%of(id)%loc%vel%z%old(i+2,j,k)+4.0d0*p%of(id)%loc%vel%z%old(i+1,j,k)-3.0d0*p%of(id)%loc%vel%z%old(i,j,k))/p%glb%dx 
        xm = 0.5d0*( p%of(id)%loc%vel%z%old(i-2,j,k)-4.0d0*p%of(id)%loc%vel%z%old(i-1,j,k)+3.0d0*p%of(id)%loc%vel%z%old(i,j,k))/p%glb%dx
            
        yp = 0.5d0*(-p%of(id)%loc%vel%z%old(i,j+2,k)+4.0d0*p%of(id)%loc%vel%z%old(i,j+1,k)-3.0d0*p%of(id)%loc%vel%z%old(i,j,k))/p%glb%dy
        ym = 0.5d0*( p%of(id)%loc%vel%z%old(i,j-2,k)-4.0d0*p%of(id)%loc%vel%z%old(i,j-1,k)+3.0d0*p%of(id)%loc%vel%z%old(i,j,k))/p%glb%dy
            
        zp = 0.5d0*(-p%of(id)%loc%vel%z%old(i,j,k+2)+4.0d0*p%of(id)%loc%vel%z%old(i,j,k+1)-3.0d0*p%of(id)%loc%vel%z%old(i,j,k))/p%glb%dz
        zm = 0.5d0*( p%of(id)%loc%vel%z%old(i,j,k-2)-4.0d0*p%of(id)%loc%vel%z%old(i,j,k-1)+3.0d0*p%of(id)%loc%vel%z%old(i,j,k))/p%glb%dz
            
        p%of(id)%loc%velsrc%z%now(i,j,k) = - ((u+abs(u))*xm+(u-abs(u))*xp)/2.0d0 &
                                        &  - ((v+abs(v))*ym+(v-abs(v))*yp)/2.0d0 &
                                        &  - ((w+abs(w))*zm+(w-abs(w))*zp)/2.0d0
                                          
    end do
    end do 
    end do
        
!$omp end parallel
    
end subroutine

subroutine ns_ab_adv_source_quick
!--------------------------------------
! u*Px = [ u_{i}+u_{i+1} ] * [ P_{i+1/2}-P_{i-1/2} ] / dx
!--------------------------------------
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: uh, vh, wh

    !$omp parallel private(id,i,j,k,uh,vh,wh), num_threads(p%glb%threads)
    
        id=0
        !$ id = omp_get_thread_num()
        
        do k = p%of(id)%loc%ks-1, p%of(id)%loc%ke
        do j = p%of(id)%loc%js-1, p%of(id)%loc%je
        do i = p%of(id)%loc%is-1, p%of(id)%loc%ie
        
            uh = 0.5d0*( p%of(id)%loc%vel%x%old(i,j,k) + p%of(id)%loc%vel%x%old(i+1,j,k) )
            vh = 0.5d0*( p%of(id)%loc%vel%y%old(i,j,k) + p%of(id)%loc%vel%y%old(i+1,j,k) )
            wh = 0.5d0*( p%of(id)%loc%vel%z%old(i,j,k) + p%of(id)%loc%vel%z%old(i+1,j,k) )
            
            if( uh>=0.0d0 )then
                p%of(id)%loc%tdata%x%s1(i,j,k) = ( -p%of(id)%loc%vel%x%old(i-1,j,k)+6.0d0*p%of(id)%loc%vel%x%old(i,j,k)+3.0d0*p%of(id)%loc%vel%x%old(i+1,j,k) )/8.0d0
            else
                p%of(id)%loc%tdata%x%s1(i,j,k) = ( -p%of(id)%loc%vel%x%old(i+2,j,k)+6.0d0*p%of(id)%loc%vel%x%old(i+1,j,k)+3.0d0*p%of(id)%loc%vel%x%old(i,j,k) )/8.0d0
            endif
            
            if( vh>=0.0d0 )then
                p%of(id)%loc%tdata%x%s2(i,j,k) = ( -p%of(id)%loc%vel%x%old(i,j-1,k)+6.0d0*p%of(id)%loc%vel%x%old(i,j,k)+3.0d0*p%of(id)%loc%vel%x%old(i,j+1,k) )/8.0d0
            else
                p%of(id)%loc%tdata%x%s2(i,j,k) = ( -p%of(id)%loc%vel%x%old(i,j+2,k)+6.0d0*p%of(id)%loc%vel%x%old(i,j+1,k)+3.0d0*p%of(id)%loc%vel%x%old(i,j,k) )/8.0d0
            endif
            
            if( wh>=0.0d0 )then
                p%of(id)%loc%tdata%x%s3(i,j,k) = ( -p%of(id)%loc%vel%x%old(i,j,k-1)+6.0d0*p%of(id)%loc%vel%x%old(i,j,k)+3.0d0*p%of(id)%loc%vel%x%old(i,j,k+1) )/8.0d0
            else
                p%of(id)%loc%tdata%x%s3(i,j,k) = ( -p%of(id)%loc%vel%x%old(i,j,k+2)+6.0d0*p%of(id)%loc%vel%x%old(i,j,k+1)+3.0d0*p%of(id)%loc%vel%x%old(i,j,k) )/8.0d0
            endif
            
            !-----------------------------------------------------------
            
            uh = 0.5d0*( p%of(id)%loc%vel%x%old(i,j,k) + p%of(id)%loc%vel%x%old(i,j+1,k) )
            vh = 0.5d0*( p%of(id)%loc%vel%y%old(i,j,k) + p%of(id)%loc%vel%y%old(i,j+1,k) )
            wh = 0.5d0*( p%of(id)%loc%vel%z%old(i,j,k) + p%of(id)%loc%vel%z%old(i,j+1,k) )
            
            if( uh>=0.0d0 )then
                p%of(id)%loc%tdata%y%s1(i,j,k) = ( -p%of(id)%loc%vel%y%old(i-1,j,k)+6.0d0*p%of(id)%loc%vel%y%old(i,j,k)+3.0d0*p%of(id)%loc%vel%y%old(i+1,j,k) )/8.0d0
            else
                p%of(id)%loc%tdata%y%s1(i,j,k) = ( -p%of(id)%loc%vel%y%old(i+2,j,k)+6.0d0*p%of(id)%loc%vel%y%old(i+1,j,k)+3.0d0*p%of(id)%loc%vel%y%old(i,j,k) )/8.0d0
            endif
            
            if( vh>=0.0d0 )then
                p%of(id)%loc%tdata%y%s2(i,j,k) = ( -p%of(id)%loc%vel%y%old(i,j-1,k)+6.0d0*p%of(id)%loc%vel%y%old(i,j,k)+3.0d0*p%of(id)%loc%vel%y%old(i,j+1,k) )/8.0d0
            else
                p%of(id)%loc%tdata%y%s2(i,j,k) = ( -p%of(id)%loc%vel%y%old(i,j+2,k)+6.0d0*p%of(id)%loc%vel%y%old(i,j+1,k)+3.0d0*p%of(id)%loc%vel%y%old(i,j,k) )/8.0d0
            endif
            
            if( wh>=0.0d0 )then
                p%of(id)%loc%tdata%y%s3(i,j,k) = ( -p%of(id)%loc%vel%y%old(i,j,k-1)+6.0d0*p%of(id)%loc%vel%y%old(i,j,k)+3.0d0*p%of(id)%loc%vel%y%old(i,j,k+1) )/8.0d0
            else
                p%of(id)%loc%tdata%y%s3(i,j,k) = ( -p%of(id)%loc%vel%y%old(i,j,k+2)+6.0d0*p%of(id)%loc%vel%y%old(i,j,k+1)+3.0d0*p%of(id)%loc%vel%y%old(i,j,k) )/8.0d0
            endif
            
            !-----------------------------------------------------------

            uh = 0.5d0*( p%of(id)%loc%vel%x%old(i,j,k) + p%of(id)%loc%vel%x%old(i,j,k+1) )
            vh = 0.5d0*( p%of(id)%loc%vel%y%old(i,j,k) + p%of(id)%loc%vel%y%old(i,j,k+1) )
            wh = 0.5d0*( p%of(id)%loc%vel%z%old(i,j,k) + p%of(id)%loc%vel%z%old(i,j,k+1) )
            
            if( uh>=0.0d0 )then
                p%of(id)%loc%tdata%z%s1(i,j,k) = ( -p%of(id)%loc%vel%z%old(i-1,j,k)+6.0d0*p%of(id)%loc%vel%z%old(i,j,k)+3.0d0*p%of(id)%loc%vel%z%old(i+1,j,k) )/8.0d0
            else
                p%of(id)%loc%tdata%z%s1(i,j,k) = ( -p%of(id)%loc%vel%z%old(i+2,j,k)+6.0d0*p%of(id)%loc%vel%z%old(i+1,j,k)+3.0d0*p%of(id)%loc%vel%z%old(i,j,k) )/8.0d0
            endif
            
            if( vh>=0.0d0 )then
                p%of(id)%loc%tdata%z%s2(i,j,k) = ( -p%of(id)%loc%vel%z%old(i,j-1,k)+6.0d0*p%of(id)%loc%vel%z%old(i,j,k)+3.0d0*p%of(id)%loc%vel%z%old(i,j+1,k) )/8.0d0
            else
                p%of(id)%loc%tdata%z%s2(i,j,k) = ( -p%of(id)%loc%vel%z%old(i,j+2,k)+6.0d0*p%of(id)%loc%vel%z%old(i,j+1,k)+3.0d0*p%of(id)%loc%vel%z%old(i,j,k) )/8.0d0
            endif
            
            if( wh>=0.0d0 )then
                p%of(id)%loc%tdata%z%s3(i,j,k) = ( -p%of(id)%loc%vel%z%old(i,j,k-1)+6.0d0*p%of(id)%loc%vel%z%old(i,j,k)+3.0d0*p%of(id)%loc%vel%z%old(i,j,k+1) )/8.0d0
            else
                p%of(id)%loc%tdata%z%s3(i,j,k) = ( -p%of(id)%loc%vel%z%old(i,j,k+2)+6.0d0*p%of(id)%loc%vel%z%old(i,j,k+1)+3.0d0*p%of(id)%loc%vel%z%old(i,j,k) )/8.0d0
            endif
            
        end do
        end do
        end do
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            uh = 0.5d0*( p%of(id)%loc%vel%x%old(i,j,k) + p%of(id)%loc%vel%x%old(i+1,j,k) )
            vh = 0.5d0*( p%of(id)%loc%vel%y%old(i,j,k) + p%of(id)%loc%vel%y%old(i+1,j,k) )
            wh = 0.5d0*( p%of(id)%loc%vel%z%old(i,j,k) + p%of(id)%loc%vel%z%old(i+1,j,k) )
        
            p%of(id)%loc%velsrc%x%now(i,j,k) = - uh*( p%of(id)%loc%tdata%x%s1(i,j,k)-p%of(id)%loc%tdata%x%s1(i-1,j,k) )/p%glb%dx &
                                              &- vh*( p%of(id)%loc%tdata%x%s2(i,j,k)-p%of(id)%loc%tdata%x%s2(i,j-1,k) )/p%glb%dy &
                                              &- wh*( p%of(id)%loc%tdata%x%s3(i,j,k)-p%of(id)%loc%tdata%x%s3(i,j,k-1) )/p%glb%dz
            
            uh = 0.5d0*( p%of(id)%loc%vel%x%old(i,j,k) + p%of(id)%loc%vel%x%old(i,j+1,k) )
            vh = 0.5d0*( p%of(id)%loc%vel%y%old(i,j,k) + p%of(id)%loc%vel%y%old(i,j+1,k) )
            wh = 0.5d0*( p%of(id)%loc%vel%z%old(i,j,k) + p%of(id)%loc%vel%z%old(i,j+1,k) )
            
            p%of(id)%loc%velsrc%y%now(i,j,k) = - uh*( p%of(id)%loc%tdata%y%s1(i,j,k)-p%of(id)%loc%tdata%y%s1(i-1,j,k) )/p%glb%dx &
                                              &- vh*( p%of(id)%loc%tdata%y%s2(i,j,k)-p%of(id)%loc%tdata%y%s2(i,j-1,k) )/p%glb%dy &
                                              &- wh*( p%of(id)%loc%tdata%y%s3(i,j,k)-p%of(id)%loc%tdata%y%s3(i,j,k-1) )/p%glb%dz
            
            uh = 0.5d0*( p%of(id)%loc%vel%x%old(i,j,k) + p%of(id)%loc%vel%x%old(i,j,k+1) )
            vh = 0.5d0*( p%of(id)%loc%vel%y%old(i,j,k) + p%of(id)%loc%vel%y%old(i,j,k+1) )
            wh = 0.5d0*( p%of(id)%loc%vel%z%old(i,j,k) + p%of(id)%loc%vel%z%old(i,j,k+1) )
            
            p%of(id)%loc%velsrc%z%now(i,j,k) = - uh*( p%of(id)%loc%tdata%z%s1(i,j,k)-p%of(id)%loc%tdata%z%s1(i-1,j,k) )/p%glb%dx &
                                              &- vh*( p%of(id)%loc%tdata%z%s2(i,j,k)-p%of(id)%loc%tdata%z%s2(i,j-1,k) )/p%glb%dy &
                                              &- wh*( p%of(id)%loc%tdata%z%s3(i,j,k)-p%of(id)%loc%tdata%z%s3(i,j,k-1) )/p%glb%dz
                                              
        end do
        end do
        end do
    
    !$omp end parallel
    
end subroutine

subroutine ns_ab_diff_source
implicit none

    call u_source
    call v_source
    call w_source 
    
end subroutine

subroutine u_source()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: rho,mu,delta,xx,yy,zz
real(8) :: ux,uy,uz,vx,wx,phix,phiy,phiz,dif_x,dif_y,dif_z,curv

    !$omp parallel private(id,i,j,k,rho,mu,delta,xx,yy,zz), &
    !$omp& private(ux,uy,uz,vx,wx,phix,phiy,phiz,dif_x,dif_y,dif_z,curv), num_threads(p%glb%threads)
    
        id=0
        !$ id = omp_get_thread_num()
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            !rho = 0.5d0*(p%of(id)%loc%rho%old(i,j,k)+p%of(id)%loc%rho%old(i+1,j,k))
            !mu = 0.5d0*(p%of(id)%loc%mu%old(i,j,k)+p%of(id)%loc%mu%old(i+1,j,k))
            !delta = 0.5d0*(p%of(id)%loc%delta%old(i,j,k)+p%of(id)%loc%delta%old(i+1,j,k))
            !curv = (p%of(id)%loc%normals%curv%old(i,j,k)+p%of(id)%loc%normals%curv%old(i+1,j,k))/2.0d0
            
            xx = (p%of(id)%loc%vel%x%old(I+1,J,K)-2.0d0*p%of(id)%loc%vel%x%old(I,J,K)+p%of(id)%loc%vel%x%old(I-1,J,K))/p%glb%dx**2.0d0
            yy = (p%of(id)%loc%vel%x%old(I,J+1,K)-2.0d0*p%of(id)%loc%vel%x%old(I,J,K)+p%of(id)%loc%vel%x%old(I,J-1,K))/p%glb%dy**2.0d0
            zz = (p%of(id)%loc%vel%x%old(I,J,K+1)-2.0d0*p%of(id)%loc%vel%x%old(I,J,K)+p%of(id)%loc%vel%x%old(I,J,K-1))/p%glb%dz**2.0d0
            
            !ux = 0.5d0*( p%of(id)%loc%vel%x%old(i+1,j,k)-p%of(id)%loc%vel%x%old(i-1,j,k) )/p%glb%dx
            !uy = 0.5d0*( p%of(id)%loc%vel%x%old(i,j+1,k)-p%of(id)%loc%vel%x%old(i,j-1,k) )/p%glb%dy
            !uz = 0.5d0*( p%of(id)%loc%vel%x%old(i,j,k+1)-p%of(id)%loc%vel%x%old(i,j,k-1) )/p%glb%dz
    
            !vx = 0.5d0*( p%of(id)%loc%vel%y%old(i+1,j,k)-p%of(id)%loc%vel%y%old(i,j,k)+p%of(id)%loc%vel%y%old(i+1,j-1,k)-p%of(id)%loc%vel%y%old(i,j-1,k) )/p%glb%dx
            !wx = 0.5d0*( p%of(id)%loc%vel%z%old(i+1,j,k)-p%of(id)%loc%vel%z%old(i,j,k)+p%of(id)%loc%vel%z%old(i+1,j,k-1)-p%of(id)%loc%vel%z%old(i,j,k-1) )/p%glb%dx
    
            !phix = (p%of(id)%loc%phi%old(i+1,j,k)-p%of(id)%loc%phi%old(i,j,k))/p%glb%dx
            !phiy = 0.25d0*(p%of(id)%loc%phi%old(i+1,j+1,k)-p%of(id)%loc%phi%old(i+1,j-1,k)+p%of(id)%loc%phi%old(i,j+1,k)-p%of(id)%loc%phi%old(i,j-1,k))/p%glb%dy
            !phiz = 0.25d0*(p%of(id)%loc%phi%old(i+1,j,k+1)-p%of(id)%loc%phi%old(i+1,j,k-1)+p%of(id)%loc%phi%old(i,j,k+1)-p%of(id)%loc%phi%old(i,j,k-1))/p%glb%dz
            
            !dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*2.0d0*ux/(rho*p%glb%re)
            !dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*(uy+vx)/(rho*p%glb%re)
            !dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*(uz+wx)/(rho*p%glb%re)
            
			dif_x = xx/p%glb%re 
            dif_y = yy/p%glb%re 
            dif_z = zz/p%glb%re
			
            p%of(id)%loc%velsrc%x%tmp(i,j,k) = dif_x + dif_y + dif_z !+ p%glb%gx*p%glb%btn_g / p%glb%fr - p%glb%btn_sf*curv*delta*phix / (p%glb%we*rho)
            
            ! =================================================
            
            !rho = 0.5d0*(p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i+1,j,k))
            !mu = 0.5d0*(p%of(id)%loc%mu%now(i,j,k)+p%of(id)%loc%mu%now(i+1,j,k))
            !delta = 0.5d0*(p%of(id)%loc%delta%now(i,j,k)+p%of(id)%loc%delta%now(i+1,j,k))
            !curv = (p%of(id)%loc%normals%curv%now(i,j,k)+p%of(id)%loc%normals%curv%now(i+1,j,k))/2.0d0
            
            xx = (p%of(id)%loc%vel%x%tmp(I+1,J,K)-2.0d0*p%of(id)%loc%vel%x%tmp(I,J,K)+p%of(id)%loc%vel%x%tmp(I-1,J,K))/p%glb%dx**2.0d0
            yy = (p%of(id)%loc%vel%x%tmp(I,J+1,K)-2.0d0*p%of(id)%loc%vel%x%tmp(I,J,K)+p%of(id)%loc%vel%x%tmp(I,J-1,K))/p%glb%dy**2.0d0
            zz = (p%of(id)%loc%vel%x%tmp(I,J,K+1)-2.0d0*p%of(id)%loc%vel%x%tmp(I,J,K)+p%of(id)%loc%vel%x%tmp(I,J,K-1))/p%glb%dz**2.0d0
            
            !ux = 0.5d0*( p%of(id)%loc%vel%x%tmp(i+1,j,k)-p%of(id)%loc%vel%x%tmp(i-1,j,k) )/p%glb%dx
            !uy = 0.5d0*( p%of(id)%loc%vel%x%tmp(i,j+1,k)-p%of(id)%loc%vel%x%tmp(i,j-1,k) )/p%glb%dy
            !uz = 0.5d0*( p%of(id)%loc%vel%x%tmp(i,j,k+1)-p%of(id)%loc%vel%x%tmp(i,j,k-1) )/p%glb%dz
    
            !vx = 0.5d0*( p%of(id)%loc%vel%y%tmp(i+1,j,k)-p%of(id)%loc%vel%y%tmp(i,j,k)+p%of(id)%loc%vel%y%tmp(i+1,j-1,k)-p%of(id)%loc%vel%y%tmp(i,j-1,k) )/p%glb%dx
            !wx = 0.5d0*( p%of(id)%loc%vel%z%tmp(i+1,j,k)-p%of(id)%loc%vel%z%tmp(i,j,k)+p%of(id)%loc%vel%z%tmp(i+1,j,k-1)-p%of(id)%loc%vel%z%tmp(i,j,k-1) )/p%glb%dx
    
            !phix = (p%of(id)%loc%phi%now(i+1,j,k)-p%of(id)%loc%phi%now(i,j,k))/p%glb%dx
            !phiy = 0.25d0*(p%of(id)%loc%phi%now(i+1,j+1,k)-p%of(id)%loc%phi%now(i+1,j-1,k)+p%of(id)%loc%phi%now(i,j+1,k)-p%of(id)%loc%phi%now(i,j-1,k))/p%glb%dy
            !phiz = 0.25d0*(p%of(id)%loc%phi%now(i+1,j,k+1)-p%of(id)%loc%phi%now(i+1,j,k-1)+p%of(id)%loc%phi%now(i,j,k+1)-p%of(id)%loc%phi%now(i,j,k-1))/p%glb%dz
            
            !dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*2.0d0*ux/(rho*p%glb%re)
            !dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*(uy+vx)/(rho*p%glb%re)
            !dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*(uz+wx)/(rho*p%glb%re)
			
            dif_x = xx/p%glb%re 
            dif_y = yy/p%glb%re 
            dif_z = zz/p%glb%re
			
            !p%of(id)%loc%velsrc%x%tmp(i,j,k) = ( p%of(id)%loc%velsrc%x%tmp(i,j,k) + &
            !& dif_x + dif_y + dif_z + p%glb%gx*p%glb%btn_g / p%glb%fr - p%glb%btn_sf*curv*delta*phix / (p%glb%we*rho) )/2.0d0
			
			p%of(id)%loc%velsrc%x%tmp(i,j,k) = ( p%of(id)%loc%velsrc%x%tmp(i,j,k) + dif_x + dif_y + dif_z )/2.0d0
            
        end do
        end do
        end do 
        
    !$omp end parallel 

end subroutine

subroutine v_source()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: rho,mu,delta,xx,yy,zz
real(8) :: vx,vy,vz,wy,uy,phix,phiy,phiz,dif_x,dif_y,dif_z,curv

    !$omp parallel private(id,i,j,k,rho,mu,delta,xx,yy,zz), &
    !$omp& private(vx,vy,vz,wy,uy,phix,phiy,phiz,dif_x,dif_y,dif_z,curv), num_threads(p%glb%threads)
    
        id=0
        !$ id = omp_get_thread_num()
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            !rho = 0.5d0*(p%of(id)%loc%rho%old(i,j,k)+p%of(id)%loc%rho%old(i,j+1,k))
            !mu = 0.5d0*(p%of(id)%loc%mu%old(i,j,k)+p%of(id)%loc%mu%old(i,j+1,k))
            !delta = 0.5d0*(p%of(id)%loc%delta%old(i,j,k)+p%of(id)%loc%delta%old(i,j+1,k))
            !curv = (p%of(id)%loc%normals%curv%old(i,j,k)+p%of(id)%loc%normals%curv%old(i,j+1,k))/2.0d0
            
            xx = (p%of(id)%loc%vel%y%old(I+1,J,K)-2.0d0*p%of(id)%loc%vel%y%old(I,J,K)+p%of(id)%loc%vel%y%old(I-1,J,K))/p%glb%dx**2.0d0
            yy = (p%of(id)%loc%vel%y%old(I,J+1,K)-2.0d0*p%of(id)%loc%vel%y%old(I,J,K)+p%of(id)%loc%vel%y%old(I,J-1,K))/p%glb%dy**2.0d0
            zz = (p%of(id)%loc%vel%y%old(I,J,K+1)-2.0d0*p%of(id)%loc%vel%y%old(I,J,K)+p%of(id)%loc%vel%y%old(I,J,K-1))/p%glb%dz**2.0d0
            
            !vx = 0.5d0*( p%of(id)%loc%vel%y%old(i+1,j,k)-p%of(id)%loc%vel%y%old(i-1,j,k) )/p%glb%dx
            !vy = 0.5d0*( p%of(id)%loc%vel%y%old(i,j+1,k)-p%of(id)%loc%vel%y%old(i,j-1,k) )/p%glb%dy
            !vz = 0.5d0*( p%of(id)%loc%vel%y%old(i,j,k+1)-p%of(id)%loc%vel%y%old(i,j,k-1) )/p%glb%dz
    
            !uy = 0.5d0*( p%of(id)%loc%vel%x%old(i,j+1,k)-p%of(id)%loc%vel%x%old(i,j,k)+p%of(id)%loc%vel%x%old(i-1,j+1,k)-p%of(id)%loc%vel%x%old(i-1,j,k) )/p%glb%dy
            !wy = 0.5d0*( p%of(id)%loc%vel%z%old(i,j+1,k)-p%of(id)%loc%vel%z%old(i,j,k)+p%of(id)%loc%vel%z%old(i,j+1,k-1)-p%of(id)%loc%vel%z%old(i,j,k-1) )/p%glb%dy
    
            !phix = 0.25d0*(p%of(id)%loc%phi%old(i+1,j,k)-p%of(id)%loc%phi%old(i-1,j,k)+p%of(id)%loc%phi%old(i+1,j+1,k)-p%of(id)%loc%phi%old(i-1,j+1,k))/p%glb%dx
            !phiy = ( p%of(id)%loc%phi%old(i,j+1,k)-p%of(id)%loc%phi%old(i,j,k) )/p%glb%dy
            !phiz = 0.25d0*(p%of(id)%loc%phi%old(i,j+1,k+1)-p%of(id)%loc%phi%old(i,j+1,k-1)+p%of(id)%loc%phi%old(i,j,k+1)-p%of(id)%loc%phi%old(i,j,k-1))/p%glb%dz
            
            !dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*(uy+vx)/(rho*p%glb%re)
            !dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*2.0d0*vy/(rho*p%glb%re)
            !dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*(wy+vz)/(rho*p%glb%re)
            
			dif_x = xx/p%glb%re 
            dif_y = yy/p%glb%re 
            dif_z = zz/p%glb%re
			
            p%of(id)%loc%velsrc%y%tmp(i,j,k) = dif_x + dif_y + dif_z !+ p%glb%gy*p%glb%btn_g / p%glb%fr - p%glb%btn_sf*curv*delta*phiy / (p%glb%we*rho)
            
            ! =========================================================================
            
            !rho = 0.5d0*(p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j+1,k))
            !mu = 0.5d0*(p%of(id)%loc%mu%now(i,j,k)+p%of(id)%loc%mu%now(i,j+1,k))
            !delta = 0.5d0*(p%of(id)%loc%delta%now(i,j,k)+p%of(id)%loc%delta%now(i,j+1,k))
            !curv = (p%of(id)%loc%normals%curv%now(i,j,k)+p%of(id)%loc%normals%curv%now(i,j+1,k))/2.0d0
            
            xx = (p%of(id)%loc%vel%y%tmp(I+1,J,K)-2.0d0*p%of(id)%loc%vel%y%tmp(I,J,K)+p%of(id)%loc%vel%y%tmp(I-1,J,K))/p%glb%dx**2.0d0
            yy = (p%of(id)%loc%vel%y%tmp(I,J+1,K)-2.0d0*p%of(id)%loc%vel%y%tmp(I,J,K)+p%of(id)%loc%vel%y%tmp(I,J-1,K))/p%glb%dy**2.0d0
            zz = (p%of(id)%loc%vel%y%tmp(I,J,K+1)-2.0d0*p%of(id)%loc%vel%y%tmp(I,J,K)+p%of(id)%loc%vel%y%tmp(I,J,K-1))/p%glb%dz**2.0d0
            
            !vx = 0.5d0*( p%of(id)%loc%vel%y%tmp(i+1,j,k)-p%of(id)%loc%vel%y%tmp(i-1,j,k) )/p%glb%dx
            !vy = 0.5d0*( p%of(id)%loc%vel%y%tmp(i,j+1,k)-p%of(id)%loc%vel%y%tmp(i,j-1,k) )/p%glb%dy
            !vz = 0.5d0*( p%of(id)%loc%vel%y%tmp(i,j,k+1)-p%of(id)%loc%vel%y%tmp(i,j,k-1) )/p%glb%dz
    
            !uy = 0.5d0*( p%of(id)%loc%vel%x%tmp(i,j+1,k)-p%of(id)%loc%vel%x%tmp(i,j,k)+p%of(id)%loc%vel%x%tmp(i-1,j+1,k)-p%of(id)%loc%vel%x%tmp(i-1,j,k) )/p%glb%dy
            !wy = 0.5d0*( p%of(id)%loc%vel%z%tmp(i,j+1,k)-p%of(id)%loc%vel%z%tmp(i,j,k)+p%of(id)%loc%vel%z%tmp(i,j+1,k-1)-p%of(id)%loc%vel%z%tmp(i,j,k-1) )/p%glb%dy
    
            !phix = 0.25d0*(p%of(id)%loc%phi%now(i+1,j,k)-p%of(id)%loc%phi%now(i-1,j,k)+p%of(id)%loc%phi%now(i+1,j+1,k)-p%of(id)%loc%phi%now(i-1,j+1,k))/p%glb%dx
            !phiy = ( p%of(id)%loc%phi%now(i,j+1,k)-p%of(id)%loc%phi%now(i,j,k) )/p%glb%dy
            !phiz = 0.25d0*(p%of(id)%loc%phi%now(i,j+1,k+1)-p%of(id)%loc%phi%now(i,j+1,k-1)+p%of(id)%loc%phi%now(i,j,k+1)-p%of(id)%loc%phi%now(i,j,k-1))/p%glb%dz
            
            !dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*(uy+vx)/(rho*p%glb%re)
            !dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*2.0d0*vy/(rho*p%glb%re)
            !dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*(wy+vz)/(rho*p%glb%re)
            
			dif_x = xx/p%glb%re 
            dif_y = yy/p%glb%re 
            dif_z = zz/p%glb%re
			
            !p%of(id)%loc%velsrc%y%tmp(i,j,k) = ( p%of(id)%loc%velsrc%y%tmp(i,j,k) + &
            !& dif_x + dif_y + dif_z + p%glb%gy*p%glb%btn_g / p%glb%fr - p%glb%btn_sf*curv*delta*phiy / (p%glb%we*rho) )/2.0d0
			
			p%of(id)%loc%velsrc%y%tmp(i,j,k) = ( p%of(id)%loc%velsrc%y%tmp(i,j,k) +  dif_x + dif_y + dif_z )/2.0d0
            
        end do
        end do
        end do 
        
    !$omp end parallel 

end subroutine

subroutine w_source()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: rho,mu,delta,xx,yy,zz,cnew,cold
real(8) :: wx,wy,wz,uz,vz,phix,phiy,phiz,dif_x,dif_y,dif_z,curv

    !$omp parallel private(id,i,j,k,rho,mu,delta,xx,yy,zz), &
    !$omp& private(wx,wy,wz,uz,vz,phix,phiy,phiz,dif_x,dif_y,dif_z,curv,cnew,cold), num_threads(p%glb%threads)
    
        id=0
        !$ id = omp_get_thread_num()
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            !rho = 0.5d0*(p%of(id)%loc%rho%old(i,j,k)+p%of(id)%loc%rho%old(i,j,k+1))
            !mu = 0.5d0*(p%of(id)%loc%mu%old(i,j,k)+p%of(id)%loc%mu%old(i,j,k+1))
            !delta = 0.5d0*(p%of(id)%loc%delta%old(i,j,k)+p%of(id)%loc%delta%old(i,j,k+1))
            !curv = (p%of(id)%loc%normals%curv%old(i,j,k)+p%of(id)%loc%normals%curv%old(i,j+1,k))/2.0d0
            
            xx = (p%of(id)%loc%vel%z%old(I+1,J,K)-2.0d0*p%of(id)%loc%vel%z%old(I,J,K)+p%of(id)%loc%vel%z%old(I-1,J,K))/p%glb%dx**2.0d0
            yy = (p%of(id)%loc%vel%z%old(I,J+1,K)-2.0d0*p%of(id)%loc%vel%z%old(I,J,K)+p%of(id)%loc%vel%z%old(I,J-1,K))/p%glb%dy**2.0d0
            zz = (p%of(id)%loc%vel%z%old(I,J,K+1)-2.0d0*p%of(id)%loc%vel%z%old(I,J,K)+p%of(id)%loc%vel%z%old(I,J,K-1))/p%glb%dz**2.0d0
            
            !wx = 0.5d0*( p%of(id)%loc%vel%z%old(i+1,j,k)-p%of(id)%loc%vel%z%old(i-1,j,k) )/p%glb%dx
            !wy = 0.5d0*( p%of(id)%loc%vel%z%old(i,j+1,k)-p%of(id)%loc%vel%z%old(i,j-1,k) )/p%glb%dy
            !wz = 0.5d0*( p%of(id)%loc%vel%z%old(i,j,k+1)-p%of(id)%loc%vel%z%old(i,j,k-1) )/p%glb%dz
    
            !uz = 0.5d0*( p%of(id)%loc%vel%x%old(i,j,k+1)-p%of(id)%loc%vel%x%old(i,j,k)+p%of(id)%loc%vel%x%old(i-1,j,k+1)-p%of(id)%loc%vel%x%old(i-1,j,k) )/p%glb%dz
            !vz = 0.5d0*( p%of(id)%loc%vel%y%old(i,j,k+1)-p%of(id)%loc%vel%y%old(i,j,k)+p%of(id)%loc%vel%y%old(i,j-1,k+1)-p%of(id)%loc%vel%y%old(i,j-1,k) )/p%glb%dz
    
            !phix = 0.25d0*( p%of(id)%loc%phi%old(i+1,j,k+1)-p%of(id)%loc%phi%old(i-1,j,k+1) + p%of(id)%loc%phi%old(i+1,j,k)-p%of(id)%loc%phi%old(i-1,j,k) )/p%glb%dx
            !phiy = 0.25d0*( p%of(id)%loc%phi%old(i,j+1,k+1)-p%of(id)%loc%phi%old(i,j-1,k+1) + p%of(id)%loc%phi%old(i,j+1,k)-p%of(id)%loc%phi%old(i,j-1,k) )/p%glb%dy
            !phiz = ( p%of(id)%loc%phi%old(i,j,k+1)-p%of(id)%loc%phi%old(i,j,k) )/p%glb%dz
            
            !dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*(uz+wx)/(rho*p%glb%re)
            !dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*(vz+wy)/(rho*p%glb%re)
            !dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*2.0d0*wz/(rho*p%glb%re)
            
			dif_x = xx/p%glb%re 
            dif_y = yy/p%glb%re 
            dif_z = zz/p%glb%re
			
            p%of(id)%loc%velsrc%z%tmp(i,j,k) = dif_x + dif_y + dif_z !+ p%glb%gz*p%glb%btn_g / p%glb%fr - p%glb%btn_sf*curv*delta*phiz / (p%glb%we*rho)
            
            ! ==========================================================================
            
            !rho = 0.5d0*(p%of(id)%loc%rho%now(i,j,k)+p%of(id)%loc%rho%now(i,j,k+1))
            !mu = 0.5d0*(p%of(id)%loc%mu%now(i,j,k)+p%of(id)%loc%mu%now(i,j,k+1))
            !delta = 0.5d0*(p%of(id)%loc%delta%now(i,j,k)+p%of(id)%loc%delta%now(i,j,k+1))
            !curv = (p%of(id)%loc%normals%curv%now(i,j,k)+p%of(id)%loc%normals%curv%now(i,j+1,k))/2.0d0
            
            xx = (p%of(id)%loc%vel%z%tmp(I+1,J,K)-2.0d0*p%of(id)%loc%vel%z%tmp(I,J,K)+p%of(id)%loc%vel%z%tmp(I-1,J,K))/p%glb%dx**2.0d0
            yy = (p%of(id)%loc%vel%z%tmp(I,J+1,K)-2.0d0*p%of(id)%loc%vel%z%tmp(I,J,K)+p%of(id)%loc%vel%z%tmp(I,J-1,K))/p%glb%dy**2.0d0
            zz = (p%of(id)%loc%vel%z%tmp(I,J,K+1)-2.0d0*p%of(id)%loc%vel%z%tmp(I,J,K)+p%of(id)%loc%vel%z%tmp(I,J,K-1))/p%glb%dz**2.0d0
            
            !wx = 0.5d0*( p%of(id)%loc%vel%z%tmp(i+1,j,k)-p%of(id)%loc%vel%z%tmp(i-1,j,k) )/p%glb%dx
            !wy = 0.5d0*( p%of(id)%loc%vel%z%tmp(i,j+1,k)-p%of(id)%loc%vel%z%tmp(i,j-1,k) )/p%glb%dy
            !wz = 0.5d0*( p%of(id)%loc%vel%z%tmp(i,j,k+1)-p%of(id)%loc%vel%z%tmp(i,j,k-1) )/p%glb%dz
    
            !uz = 0.5d0*( p%of(id)%loc%vel%x%tmp(i,j,k+1)-p%of(id)%loc%vel%x%tmp(i,j,k)+p%of(id)%loc%vel%x%tmp(i-1,j,k+1)-p%of(id)%loc%vel%x%tmp(i-1,j,k) )/p%glb%dz
            !vz = 0.5d0*( p%of(id)%loc%vel%y%tmp(i,j,k+1)-p%of(id)%loc%vel%y%tmp(i,j,k)+p%of(id)%loc%vel%y%tmp(i,j-1,k+1)-p%of(id)%loc%vel%y%tmp(i,j-1,k) )/p%glb%dz
    
            !phix = 0.25d0*( p%of(id)%loc%phi%now(i+1,j,k+1)-p%of(id)%loc%phi%now(i-1,j,k+1) + p%of(id)%loc%phi%now(i+1,j,k)-p%of(id)%loc%phi%now(i-1,j,k) )/p%glb%dx
            !phiy = 0.25d0*( p%of(id)%loc%phi%now(i,j+1,k+1)-p%of(id)%loc%phi%now(i,j-1,k+1) + p%of(id)%loc%phi%now(i,j+1,k)-p%of(id)%loc%phi%now(i,j-1,k) )/p%glb%dy
            !phiz = ( p%of(id)%loc%phi%now(i,j,k+1)-p%of(id)%loc%phi%now(i,j,k) )/p%glb%dz
            
            !dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*(uz+wx)/(rho*p%glb%re)
            !dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*(vz+wy)/(rho*p%glb%re)
            !dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*2.0d0*wz/(rho*p%glb%re)
            
			dif_x = xx/p%glb%re 
            dif_y = yy/p%glb%re 
            dif_z = zz/p%glb%re
			
            !p%of(id)%loc%velsrc%z%tmp(i,j,k) = ( p%of(id)%loc%velsrc%z%tmp(i,j,k) + &
            !& dif_x + dif_y + dif_z + p%glb%gz*p%glb%btn_g / p%glb%fr - p%glb%btn_sf*curv*delta*phiz / (p%glb%we*rho) ) /2.0d0
			
			p%of(id)%loc%velsrc%z%tmp(i,j,k) = ( p%of(id)%loc%velsrc%z%tmp(i,j,k) + dif_x + dif_y + dif_z  ) /2.0d0
            
			! ==========================================================================
			
			cnew = 0.5d0*( p%of(id)%loc%phi%now(i,j,k)+p%of(id)%loc%phi%now(i,j,k+1) )
			cold = 0.5d0*( p%of(id)%loc%phi%old(i,j,k)+p%of(id)%loc%phi%old(i,j,k+1) )
			
			p%of(id)%loc%velsrc%z%tmp(i,j,k) = p%of(id)%loc%velsrc%z%tmp(i,j,k) - 0.5d0*(cnew+cold)
			
        end do
        end do
        end do 
        
    !$omp end parallel 

end subroutine
