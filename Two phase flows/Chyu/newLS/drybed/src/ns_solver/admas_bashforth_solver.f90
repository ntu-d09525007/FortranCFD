subroutine ns_ab_setup
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k

    call ns_ab_adv_source

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

    if(iter>5)exit
        
end do

call ns_check_convergence_vel
call ppe_sor_solver(p%glb%p_tol)

end subroutine

subroutine ns_ab_predictor
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: src

!$omp parallel do collapse(3), private(src)      
do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie
        
    src = 1.50*p%loc%velsrc%x%now(i,j,k) - 0.5d0*p%loc%velsrc%x%old(i,j,k) 
    src = src + p%loc%velsrc%x%tmp(i,j,k)            
    p%loc%vel%x%now(i,j,k) = p%loc%vel%x%old(i,j,k) + p%glb%dt * src
    
    src = 1.50*p%loc%velsrc%y%now(i,j,k) - 0.5d0*p%loc%velsrc%y%old(i,j,k) 
    src = src + p%loc%velsrc%y%tmp(i,j,k)    
    p%loc%vel%y%now(i,j,k) = p%loc%vel%y%old(i,j,k) + p%glb%dt * src
        
    src = 1.50*p%loc%velsrc%z%now(i,j,k) - 0.5d0*p%loc%velsrc%z%old(i,j,k) 
    src = src + p%loc%velsrc%z%tmp(i,j,k)    
    p%loc%vel%z%now(i,j,k) = p%loc%vel%z%old(i,j,k) + p%glb%dt * src
        
end do
end do 
end do
!$omp end parallel do
    
call velbc(p%loc%vel%x%now,p%loc%vel%y%now,p%loc%vel%z%now)
    
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
   
call find_stag_vel( p%loc%tdata%x%s1, p%loc%tdata%y%s1, p%loc%tdata%z%s1, &
                &p%loc%tdata%x%s2, p%loc%tdata%y%s2, p%loc%tdata%z%s2, &
                &p%loc%vel%x%old, p%loc%vel%y%old, p%loc%vel%z%old )
    
!$omp parallel do collapse(3), private(u,v,w,xp,xm,yp,ym,zp,zm)
do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie
    
    u = p%loc%vel%x%old(i,j,k)
    v = p%loc%tdata%y%s1(i,j,k)
    w = p%loc%tdata%z%s1(i,j,k)
        
    xp = 0.5d0*(-p%loc%vel%x%old(i+2,j,k)+4.0d0*p%loc%vel%x%old(i+1,j,k)-3.0d0*p%loc%vel%x%old(i,j,k))/p%glb%dx 
    xm = 0.5d0*( p%loc%vel%x%old(i-2,j,k)-4.0d0*p%loc%vel%x%old(i-1,j,k)+3.0d0*p%loc%vel%x%old(i,j,k))/p%glb%dx
    
    yp = 0.5d0*(-p%loc%vel%x%old(i,j+2,k)+4.0d0*p%loc%vel%x%old(i,j+1,k)-3.0d0*p%loc%vel%x%old(i,j,k))/p%glb%dy
    ym = 0.5d0*( p%loc%vel%x%old(i,j-2,k)-4.0d0*p%loc%vel%x%old(i,j-1,k)+3.0d0*p%loc%vel%x%old(i,j,k))/p%glb%dy
        
    zp = 0.5d0*(-p%loc%vel%x%old(i,j,k+2)+4.0d0*p%loc%vel%x%old(i,j,k+1)-3.0d0*p%loc%vel%x%old(i,j,k))/p%glb%dz
    zm = 0.5d0*( p%loc%vel%x%old(i,j,k-2)-4.0d0*p%loc%vel%x%old(i,j,k-1)+3.0d0*p%loc%vel%x%old(i,j,k))/p%glb%dz
    
    p%loc%velsrc%x%now(i,j,k) = - ((u+abs(u))*xm+(u-abs(u))*xp)/2.0d0 &
                                    &  - ((v+abs(v))*ym+(v-abs(v))*yp)/2.0d0 &
                                    &  - ((w+abs(w))*zm+(w-abs(w))*zp)/2.0d0  
            
    !-----------------------------------------------------------
        
    u = p%loc%tdata%x%s1(i,j,k)
    v = p%loc%vel%y%old(i,j,k)
    w = p%loc%tdata%z%s2(i,j,k)
        
    xp = 0.5d0*(-p%loc%vel%y%old(i+2,j,k)+4.0d0*p%loc%vel%y%old(i+1,j,k)-3.0d0*p%loc%vel%y%old(i,j,k))/p%glb%dx 
    xm = 0.5d0*( p%loc%vel%y%old(i-2,j,k)-4.0d0*p%loc%vel%y%old(i-1,j,k)+3.0d0*p%loc%vel%y%old(i,j,k))/p%glb%dx
        
    yp = 0.5d0*(-p%loc%vel%y%old(i,j+2,k)+4.0d0*p%loc%vel%y%old(i,j+1,k)-3.0d0*p%loc%vel%y%old(i,j,k))/p%glb%dy
    ym = 0.5d0*( p%loc%vel%y%old(i,j-2,k)-4.0d0*p%loc%vel%y%old(i,j-1,k)+3.0d0*p%loc%vel%y%old(i,j,k))/p%glb%dy
        
    zp = 0.5d0*(-p%loc%vel%y%old(i,j,k+2)+4.0d0*p%loc%vel%y%old(i,j,k+1)-3.0d0*p%loc%vel%y%old(i,j,k))/p%glb%dz
    zm = 0.5d0*( p%loc%vel%y%old(i,j,k-2)-4.0d0*p%loc%vel%y%old(i,j,k-1)+3.0d0*p%loc%vel%y%old(i,j,k))/p%glb%dz
        
    p%loc%velsrc%y%now(i,j,k) = - ((u+abs(u))*xm+(u-abs(u))*xp)/2.0d0 &
                                    &  - ((v+abs(v))*ym+(v-abs(v))*yp)/2.0d0 &
                                    &  - ((w+abs(w))*zm+(w-abs(w))*zp)/2.0d0
        
    !-----------------------------------------------------------
        
    u = p%loc%tdata%x%s2(i,j,k)
    v = p%loc%tdata%y%s2(i,j,k)
    w = p%loc%vel%z%old(i,j,k)
        
    xp = 0.5d0*(-p%loc%vel%z%old(i+2,j,k)+4.0d0*p%loc%vel%z%old(i+1,j,k)-3.0d0*p%loc%vel%z%old(i,j,k))/p%glb%dx 
    xm = 0.5d0*( p%loc%vel%z%old(i-2,j,k)-4.0d0*p%loc%vel%z%old(i-1,j,k)+3.0d0*p%loc%vel%z%old(i,j,k))/p%glb%dx
        
    yp = 0.5d0*(-p%loc%vel%z%old(i,j+2,k)+4.0d0*p%loc%vel%z%old(i,j+1,k)-3.0d0*p%loc%vel%z%old(i,j,k))/p%glb%dy
    ym = 0.5d0*( p%loc%vel%z%old(i,j-2,k)-4.0d0*p%loc%vel%z%old(i,j-1,k)+3.0d0*p%loc%vel%z%old(i,j,k))/p%glb%dy
        
    zp = 0.5d0*(-p%loc%vel%z%old(i,j,k+2)+4.0d0*p%loc%vel%z%old(i,j,k+1)-3.0d0*p%loc%vel%z%old(i,j,k))/p%glb%dz
    zm = 0.5d0*( p%loc%vel%z%old(i,j,k-2)-4.0d0*p%loc%vel%z%old(i,j,k-1)+3.0d0*p%loc%vel%z%old(i,j,k))/p%glb%dz
        
    p%loc%velsrc%z%now(i,j,k) = - ((u+abs(u))*xm+(u-abs(u))*xp)/2.0d0 &
                                    &  - ((v+abs(v))*ym+(v-abs(v))*yp)/2.0d0 &
                                    &  - ((w+abs(w))*zm+(w-abs(w))*zp)/2.0d0
                                      
end do
end do 
end do
!$omp end parallel do
    
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

!$omp parallel do collapse(3), private(uh,vh,wh)        
do k = p%loc%ks-1, p%loc%ke
do j = p%loc%js-1, p%loc%je
do i = p%loc%is-1, p%loc%ie

    uh = 0.5d0*( p%loc%vel%x%old(i,j,k) + p%loc%vel%x%old(i+1,j,k) )
    vh = 0.5d0*( p%loc%vel%y%old(i,j,k) + p%loc%vel%y%old(i+1,j,k) )
    wh = 0.5d0*( p%loc%vel%z%old(i,j,k) + p%loc%vel%z%old(i+1,j,k) )
    
    if( uh>=0.0d0 )then
        p%loc%tdata%x%s1(i,j,k) = ( -p%loc%vel%x%old(i-1,j,k)+6.0d0*p%loc%vel%x%old(i,j,k)+3.0d0*p%loc%vel%x%old(i+1,j,k) )/8.0d0
    else
        p%loc%tdata%x%s1(i,j,k) = ( -p%loc%vel%x%old(i+2,j,k)+6.0d0*p%loc%vel%x%old(i+1,j,k)+3.0d0*p%loc%vel%x%old(i,j,k) )/8.0d0
    endif
    
    if( vh>=0.0d0 )then
        p%loc%tdata%x%s2(i,j,k) = ( -p%loc%vel%x%old(i,j-1,k)+6.0d0*p%loc%vel%x%old(i,j,k)+3.0d0*p%loc%vel%x%old(i,j+1,k) )/8.0d0
    else
        p%loc%tdata%x%s2(i,j,k) = ( -p%loc%vel%x%old(i,j+2,k)+6.0d0*p%loc%vel%x%old(i,j+1,k)+3.0d0*p%loc%vel%x%old(i,j,k) )/8.0d0
    endif
    
    if( wh>=0.0d0 )then
        p%loc%tdata%x%s3(i,j,k) = ( -p%loc%vel%x%old(i,j,k-1)+6.0d0*p%loc%vel%x%old(i,j,k)+3.0d0*p%loc%vel%x%old(i,j,k+1) )/8.0d0
    else
        p%loc%tdata%x%s3(i,j,k) = ( -p%loc%vel%x%old(i,j,k+2)+6.0d0*p%loc%vel%x%old(i,j,k+1)+3.0d0*p%loc%vel%x%old(i,j,k) )/8.0d0
    endif
    
    !-----------------------------------------------------------
    
    uh = 0.5d0*( p%loc%vel%x%old(i,j,k) + p%loc%vel%x%old(i,j+1,k) )
    vh = 0.5d0*( p%loc%vel%y%old(i,j,k) + p%loc%vel%y%old(i,j+1,k) )
    wh = 0.5d0*( p%loc%vel%z%old(i,j,k) + p%loc%vel%z%old(i,j+1,k) )
    
    if( uh>=0.0d0 )then
        p%loc%tdata%y%s1(i,j,k) = ( -p%loc%vel%y%old(i-1,j,k)+6.0d0*p%loc%vel%y%old(i,j,k)+3.0d0*p%loc%vel%y%old(i+1,j,k) )/8.0d0
    else
        p%loc%tdata%y%s1(i,j,k) = ( -p%loc%vel%y%old(i+2,j,k)+6.0d0*p%loc%vel%y%old(i+1,j,k)+3.0d0*p%loc%vel%y%old(i,j,k) )/8.0d0
    endif
    
    if( vh>=0.0d0 )then
        p%loc%tdata%y%s2(i,j,k) = ( -p%loc%vel%y%old(i,j-1,k)+6.0d0*p%loc%vel%y%old(i,j,k)+3.0d0*p%loc%vel%y%old(i,j+1,k) )/8.0d0
    else
        p%loc%tdata%y%s2(i,j,k) = ( -p%loc%vel%y%old(i,j+2,k)+6.0d0*p%loc%vel%y%old(i,j+1,k)+3.0d0*p%loc%vel%y%old(i,j,k) )/8.0d0
    endif
    
    if( wh>=0.0d0 )then
        p%loc%tdata%y%s3(i,j,k) = ( -p%loc%vel%y%old(i,j,k-1)+6.0d0*p%loc%vel%y%old(i,j,k)+3.0d0*p%loc%vel%y%old(i,j,k+1) )/8.0d0
    else
        p%loc%tdata%y%s3(i,j,k) = ( -p%loc%vel%y%old(i,j,k+2)+6.0d0*p%loc%vel%y%old(i,j,k+1)+3.0d0*p%loc%vel%y%old(i,j,k) )/8.0d0
    endif
    
    !-----------------------------------------------------------

    uh = 0.5d0*( p%loc%vel%x%old(i,j,k) + p%loc%vel%x%old(i,j,k+1) )
    vh = 0.5d0*( p%loc%vel%y%old(i,j,k) + p%loc%vel%y%old(i,j,k+1) )
    wh = 0.5d0*( p%loc%vel%z%old(i,j,k) + p%loc%vel%z%old(i,j,k+1) )
    
    if( uh>=0.0d0 )then
        p%loc%tdata%z%s1(i,j,k) = ( -p%loc%vel%z%old(i-1,j,k)+6.0d0*p%loc%vel%z%old(i,j,k)+3.0d0*p%loc%vel%z%old(i+1,j,k) )/8.0d0
    else
        p%loc%tdata%z%s1(i,j,k) = ( -p%loc%vel%z%old(i+2,j,k)+6.0d0*p%loc%vel%z%old(i+1,j,k)+3.0d0*p%loc%vel%z%old(i,j,k) )/8.0d0
    endif
    
    if( vh>=0.0d0 )then
        p%loc%tdata%z%s2(i,j,k) = ( -p%loc%vel%z%old(i,j-1,k)+6.0d0*p%loc%vel%z%old(i,j,k)+3.0d0*p%loc%vel%z%old(i,j+1,k) )/8.0d0
    else
        p%loc%tdata%z%s2(i,j,k) = ( -p%loc%vel%z%old(i,j+2,k)+6.0d0*p%loc%vel%z%old(i,j+1,k)+3.0d0*p%loc%vel%z%old(i,j,k) )/8.0d0
    endif
    
    if( wh>=0.0d0 )then
        p%loc%tdata%z%s3(i,j,k) = ( -p%loc%vel%z%old(i,j,k-1)+6.0d0*p%loc%vel%z%old(i,j,k)+3.0d0*p%loc%vel%z%old(i,j,k+1) )/8.0d0
    else
        p%loc%tdata%z%s3(i,j,k) = ( -p%loc%vel%z%old(i,j,k+2)+6.0d0*p%loc%vel%z%old(i,j,k+1)+3.0d0*p%loc%vel%z%old(i,j,k) )/8.0d0
    endif
    
end do
end do
end do
!$omp end parallel do

!$omp parallel do collapse(3)
do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie

    uh = 0.5d0*( p%loc%vel%x%old(i,j,k) + p%loc%vel%x%old(i+1,j,k) )
    vh = 0.5d0*( p%loc%vel%y%old(i,j,k) + p%loc%vel%y%old(i+1,j,k) )
    wh = 0.5d0*( p%loc%vel%z%old(i,j,k) + p%loc%vel%z%old(i+1,j,k) )

    p%loc%velsrc%x%now(i,j,k) = - uh*( p%loc%tdata%x%s1(i,j,k)-p%loc%tdata%x%s1(i-1,j,k) )/p%glb%dx &
                                      &- vh*( p%loc%tdata%x%s2(i,j,k)-p%loc%tdata%x%s2(i,j-1,k) )/p%glb%dy &
                                      &- wh*( p%loc%tdata%x%s3(i,j,k)-p%loc%tdata%x%s3(i,j,k-1) )/p%glb%dz
    
    uh = 0.5d0*( p%loc%vel%x%old(i,j,k) + p%loc%vel%x%old(i,j+1,k) )
    vh = 0.5d0*( p%loc%vel%y%old(i,j,k) + p%loc%vel%y%old(i,j+1,k) )
    wh = 0.5d0*( p%loc%vel%z%old(i,j,k) + p%loc%vel%z%old(i,j+1,k) )
    
    p%loc%velsrc%y%now(i,j,k) = - uh*( p%loc%tdata%y%s1(i,j,k)-p%loc%tdata%y%s1(i-1,j,k) )/p%glb%dx &
                                      &- vh*( p%loc%tdata%y%s2(i,j,k)-p%loc%tdata%y%s2(i,j-1,k) )/p%glb%dy &
                                      &- wh*( p%loc%tdata%y%s3(i,j,k)-p%loc%tdata%y%s3(i,j,k-1) )/p%glb%dz
    
    uh = 0.5d0*( p%loc%vel%x%old(i,j,k) + p%loc%vel%x%old(i,j,k+1) )
    vh = 0.5d0*( p%loc%vel%y%old(i,j,k) + p%loc%vel%y%old(i,j,k+1) )
    wh = 0.5d0*( p%loc%vel%z%old(i,j,k) + p%loc%vel%z%old(i,j,k+1) )
    
    p%loc%velsrc%z%now(i,j,k) = - uh*( p%loc%tdata%z%s1(i,j,k)-p%loc%tdata%z%s1(i-1,j,k) )/p%glb%dx &
                                      &- vh*( p%loc%tdata%z%s2(i,j,k)-p%loc%tdata%z%s2(i,j-1,k) )/p%glb%dy &
                                      &- wh*( p%loc%tdata%z%s3(i,j,k)-p%loc%tdata%z%s3(i,j,k-1) )/p%glb%dz
                                      
end do
end do
end do
!$omp end parallel do
    
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

    !$omp parallel do collapse(3), private(i,j,k,rho,mu,delta,xx,yy,zz), &
    !$omp& private(ux,uy,uz,vx,wx,phix,phiy,phiz,dif_x,dif_y,dif_z,curv)
    do k = p%loc%ks, p%loc%ke
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
    
        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i+1,j,k))
        mu = 0.5d0*(p%loc%mu%old(i,j,k)+p%loc%mu%old(i+1,j,k))
        delta = 0.5d0*(p%loc%delta%old(i,j,k)+p%loc%delta%old(i+1,j,k))
        curv = (p%loc%normals%curv%old(i,j,k)+p%loc%normals%curv%old(i+1,j,k))/2.0d0
        
        xx = (p%loc%vel%x%old(I+1,J,K)-2.0d0*p%loc%vel%x%old(I,J,K)+p%loc%vel%x%old(I-1,J,K))/p%glb%dx**2.0d0
        yy = (p%loc%vel%x%old(I,J+1,K)-2.0d0*p%loc%vel%x%old(I,J,K)+p%loc%vel%x%old(I,J-1,K))/p%glb%dy**2.0d0
        zz = (p%loc%vel%x%old(I,J,K+1)-2.0d0*p%loc%vel%x%old(I,J,K)+p%loc%vel%x%old(I,J,K-1))/p%glb%dz**2.0d0
        
        ux = 0.5d0*( p%loc%vel%x%old(i+1,j,k)-p%loc%vel%x%old(i-1,j,k) )/p%glb%dx
        uy = 0.5d0*( p%loc%vel%x%old(i,j+1,k)-p%loc%vel%x%old(i,j-1,k) )/p%glb%dy
        uz = 0.5d0*( p%loc%vel%x%old(i,j,k+1)-p%loc%vel%x%old(i,j,k-1) )/p%glb%dz

        vx = 0.5d0*( p%loc%vel%y%old(i+1,j,k)-p%loc%vel%y%old(i,j,k)+p%loc%vel%y%old(i+1,j-1,k)-p%loc%vel%y%old(i,j-1,k) )/p%glb%dx
        wx = 0.5d0*( p%loc%vel%z%old(i+1,j,k)-p%loc%vel%z%old(i,j,k)+p%loc%vel%z%old(i+1,j,k-1)-p%loc%vel%z%old(i,j,k-1) )/p%glb%dx

        phix = (p%loc%phi%old(i+1,j,k)-p%loc%phi%old(i,j,k))/p%glb%dx
        phiy = 0.25d0*(p%loc%phi%old(i+1,j+1,k)-p%loc%phi%old(i+1,j-1,k)+p%loc%phi%old(i,j+1,k)-p%loc%phi%old(i,j-1,k))/p%glb%dy
        phiz = 0.25d0*(p%loc%phi%old(i+1,j,k+1)-p%loc%phi%old(i+1,j,k-1)+p%loc%phi%old(i,j,k+1)-p%loc%phi%old(i,j,k-1))/p%glb%dz
        
        dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*2.0d0*ux/(rho*p%glb%re)
        dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*(uy+vx)/(rho*p%glb%re)
        dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*(uz+wx)/(rho*p%glb%re)
        
        p%loc%velsrc%x%tmp(i,j,k) = dif_x + dif_y + dif_z &
        & + p%glb%gx*p%glb%btn_g / p%glb%fr &
        & - p%glb%btn_sf*curv*delta*phix / (p%glb%we*rho) 
        
        ! =================================================
        
        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i+1,j,k))
        mu = 0.5d0*(p%loc%mu%old(i,j,k)+p%loc%mu%old(i+1,j,k))
        delta = 0.5d0*(p%loc%delta%now(i,j,k)+p%loc%delta%now(i+1,j,k))
        curv = (p%loc%normals%curv%now(i,j,k)+p%loc%normals%curv%now(i+1,j,k))/2.0d0
        
        xx = (p%loc%vel%x%tmp(I+1,J,K)-2.0d0*p%loc%vel%x%tmp(I,J,K)+p%loc%vel%x%tmp(I-1,J,K))/p%glb%dx**2.0d0
        yy = (p%loc%vel%x%tmp(I,J+1,K)-2.0d0*p%loc%vel%x%tmp(I,J,K)+p%loc%vel%x%tmp(I,J-1,K))/p%glb%dy**2.0d0
        zz = (p%loc%vel%x%tmp(I,J,K+1)-2.0d0*p%loc%vel%x%tmp(I,J,K)+p%loc%vel%x%tmp(I,J,K-1))/p%glb%dz**2.0d0
        
        ux = 0.5d0*( p%loc%vel%x%tmp(i+1,j,k)-p%loc%vel%x%tmp(i-1,j,k) )/p%glb%dx
        uy = 0.5d0*( p%loc%vel%x%tmp(i,j+1,k)-p%loc%vel%x%tmp(i,j-1,k) )/p%glb%dy
        uz = 0.5d0*( p%loc%vel%x%tmp(i,j,k+1)-p%loc%vel%x%tmp(i,j,k-1) )/p%glb%dz

        vx = 0.5d0*( p%loc%vel%y%tmp(i+1,j,k)-p%loc%vel%y%tmp(i,j,k)+p%loc%vel%y%tmp(i+1,j-1,k)-p%loc%vel%y%tmp(i,j-1,k) )/p%glb%dx
        wx = 0.5d0*( p%loc%vel%z%tmp(i+1,j,k)-p%loc%vel%z%tmp(i,j,k)+p%loc%vel%z%tmp(i+1,j,k-1)-p%loc%vel%z%tmp(i,j,k-1) )/p%glb%dx

        phix = (p%loc%phi%now(i+1,j,k)-p%loc%phi%now(i,j,k))/p%glb%dx
        phiy = 0.25d0*(p%loc%phi%now(i+1,j+1,k)-p%loc%phi%now(i+1,j-1,k)+p%loc%phi%now(i,j+1,k)-p%loc%phi%now(i,j-1,k))/p%glb%dy
        phiz = 0.25d0*(p%loc%phi%now(i+1,j,k+1)-p%loc%phi%now(i+1,j,k-1)+p%loc%phi%now(i,j,k+1)-p%loc%phi%now(i,j,k-1))/p%glb%dz
        
        dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*2.0d0*ux/(rho*p%glb%re)
        dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*(uy+vx)/(rho*p%glb%re)
        dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*(uz+wx)/(rho*p%glb%re)
        
        p%loc%velsrc%x%tmp(i,j,k) = ( p%loc%velsrc%x%tmp(i,j,k) + &
        & dif_x + dif_y + dif_z &
        & + p%glb%gx*p%glb%btn_g / p%glb%fr &
        & - p%glb%btn_sf*curv*delta*phix / (p%glb%we*rho) )/2.0d0
        
    end do
    end do
    end do 
    !$omp end parallel do

end subroutine

subroutine v_source()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: rho,mu,delta,xx,yy,zz
real(8) :: vx,vy,vz,wy,uy,phix,phiy,phiz,dif_x,dif_y,dif_z,curv

    !$omp parallel do collapse(3), private(i,j,k,rho,mu,delta,xx,yy,zz), &
    !$omp& private(vx,vy,vz,wy,uy,phix,phiy,phiz,dif_x,dif_y,dif_z,curv)
    do k = p%loc%ks, p%loc%ke
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
    
        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i,j+1,k))
        mu = 0.5d0*(p%loc%mu%old(i,j,k)+p%loc%mu%old(i,j+1,k))
        delta = 0.5d0*(p%loc%delta%old(i,j,k)+p%loc%delta%old(i,j+1,k))
        curv = (p%loc%normals%curv%old(i,j,k)+p%loc%normals%curv%old(i,j+1,k))/2.0d0
        
        xx = (p%loc%vel%y%old(I+1,J,K)-2.0d0*p%loc%vel%y%old(I,J,K)+p%loc%vel%y%old(I-1,J,K))/p%glb%dx**2.0d0
        yy = (p%loc%vel%y%old(I,J+1,K)-2.0d0*p%loc%vel%y%old(I,J,K)+p%loc%vel%y%old(I,J-1,K))/p%glb%dy**2.0d0
        zz = (p%loc%vel%y%old(I,J,K+1)-2.0d0*p%loc%vel%y%old(I,J,K)+p%loc%vel%y%old(I,J,K-1))/p%glb%dz**2.0d0
        
        vx = 0.5d0*( p%loc%vel%y%old(i+1,j,k)-p%loc%vel%y%old(i-1,j,k) )/p%glb%dx
        vy = 0.5d0*( p%loc%vel%y%old(i,j+1,k)-p%loc%vel%y%old(i,j-1,k) )/p%glb%dy
        vz = 0.5d0*( p%loc%vel%y%old(i,j,k+1)-p%loc%vel%y%old(i,j,k-1) )/p%glb%dz

        uy = 0.5d0*( p%loc%vel%x%old(i,j+1,k)-p%loc%vel%x%old(i,j,k)+p%loc%vel%x%old(i-1,j+1,k)-p%loc%vel%x%old(i-1,j,k) )/p%glb%dy
        wy = 0.5d0*( p%loc%vel%z%old(i,j+1,k)-p%loc%vel%z%old(i,j,k)+p%loc%vel%z%old(i,j+1,k-1)-p%loc%vel%z%old(i,j,k-1) )/p%glb%dy

        phix = 0.25d0*(p%loc%phi%old(i+1,j,k)-p%loc%phi%old(i-1,j,k)+p%loc%phi%old(i+1,j+1,k)-p%loc%phi%old(i-1,j+1,k))/p%glb%dx
        phiy = ( p%loc%phi%old(i,j+1,k)-p%loc%phi%old(i,j,k) )/p%glb%dy
        phiz = 0.25d0*(p%loc%phi%old(i,j+1,k+1)-p%loc%phi%old(i,j+1,k-1)+p%loc%phi%old(i,j,k+1)-p%loc%phi%old(i,j,k-1))/p%glb%dz
        
        dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*(uy+vx)/(rho*p%glb%re)
        dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*2.0d0*vy/(rho*p%glb%re)
        dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*(wy+vz)/(rho*p%glb%re)
        
        p%loc%velsrc%y%tmp(i,j,k) = dif_x + dif_y + dif_z &
        & + p%glb%gy*p%glb%btn_g / p%glb%fr &
        & - p%glb%btn_sf*curv*delta*phiy / (p%glb%we*rho)
        
        ! =========================================================================
        
        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i,j+1,k))
        mu = 0.5d0*(p%loc%mu%old(i,j,k)+p%loc%mu%old(i,j+1,k))
        delta = 0.5d0*(p%loc%delta%now(i,j,k)+p%loc%delta%now(i,j+1,k))
        curv = (p%loc%normals%curv%now(i,j,k)+p%loc%normals%curv%now(i,j+1,k))/2.0d0
        
        xx = (p%loc%vel%y%tmp(I+1,J,K)-2.0d0*p%loc%vel%y%tmp(I,J,K)+p%loc%vel%y%tmp(I-1,J,K))/p%glb%dx**2.0d0
        yy = (p%loc%vel%y%tmp(I,J+1,K)-2.0d0*p%loc%vel%y%tmp(I,J,K)+p%loc%vel%y%tmp(I,J-1,K))/p%glb%dy**2.0d0
        zz = (p%loc%vel%y%tmp(I,J,K+1)-2.0d0*p%loc%vel%y%tmp(I,J,K)+p%loc%vel%y%tmp(I,J,K-1))/p%glb%dz**2.0d0
        
        vx = 0.5d0*( p%loc%vel%y%tmp(i+1,j,k)-p%loc%vel%y%tmp(i-1,j,k) )/p%glb%dx
        vy = 0.5d0*( p%loc%vel%y%tmp(i,j+1,k)-p%loc%vel%y%tmp(i,j-1,k) )/p%glb%dy
        vz = 0.5d0*( p%loc%vel%y%tmp(i,j,k+1)-p%loc%vel%y%tmp(i,j,k-1) )/p%glb%dz

        uy = 0.5d0*( p%loc%vel%x%tmp(i,j+1,k)-p%loc%vel%x%tmp(i,j,k)+p%loc%vel%x%tmp(i-1,j+1,k)-p%loc%vel%x%tmp(i-1,j,k) )/p%glb%dy
        wy = 0.5d0*( p%loc%vel%z%tmp(i,j+1,k)-p%loc%vel%z%tmp(i,j,k)+p%loc%vel%z%tmp(i,j+1,k-1)-p%loc%vel%z%tmp(i,j,k-1) )/p%glb%dy

        phix = 0.25d0*(p%loc%phi%now(i+1,j,k)-p%loc%phi%now(i-1,j,k)+p%loc%phi%now(i+1,j+1,k)-p%loc%phi%now(i-1,j+1,k))/p%glb%dx
        phiy = ( p%loc%phi%now(i,j+1,k)-p%loc%phi%now(i,j,k) )/p%glb%dy
        phiz = 0.25d0*(p%loc%phi%now(i,j+1,k+1)-p%loc%phi%now(i,j+1,k-1)+p%loc%phi%now(i,j,k+1)-p%loc%phi%now(i,j,k-1))/p%glb%dz
        
        dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*(uy+vx)/(rho*p%glb%re)
        dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*2.0d0*vy/(rho*p%glb%re)
        dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*(wy+vz)/(rho*p%glb%re)
        
        p%loc%velsrc%y%tmp(i,j,k) = ( p%loc%velsrc%y%tmp(i,j,k) + &
        & dif_x + dif_y + dif_z + &
        & p%glb%gy*p%glb%btn_g / p%glb%fr &
        & - p%glb%btn_sf*curv*delta*phiy / (p%glb%we*rho) )/2.0d0
        
    end do
    end do
    end do 
    !$omp end parallel do

end subroutine

subroutine w_source()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: rho,mu,delta,xx,yy,zz
real(8) :: wx,wy,wz,uz,vz,phix,phiy,phiz,dif_x,dif_y,dif_z,curv

    !$omp parallel do collapse(3), private(i,j,k,rho,mu,delta,xx,yy,zz), &
    !$omp& private(wx,wy,wz,uz,vz,phix,phiy,phiz,dif_x,dif_y,dif_z,curv)
    do k = p%loc%ks, p%loc%ke
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
    
        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i,j,k+1))
        mu = 0.5d0*(p%loc%mu%old(i,j,k)+p%loc%mu%old(i,j,k+1))
        delta = 0.5d0*(p%loc%delta%old(i,j,k)+p%loc%delta%old(i,j,k+1))
        curv = (p%loc%normals%curv%old(i,j,k)+p%loc%normals%curv%old(i,j+1,k))/2.0d0
        
        xx = (p%loc%vel%z%old(I+1,J,K)-2.0d0*p%loc%vel%z%old(I,J,K)+p%loc%vel%z%old(I-1,J,K))/p%glb%dx**2.0d0
        yy = (p%loc%vel%z%old(I,J+1,K)-2.0d0*p%loc%vel%z%old(I,J,K)+p%loc%vel%z%old(I,J-1,K))/p%glb%dy**2.0d0
        zz = (p%loc%vel%z%old(I,J,K+1)-2.0d0*p%loc%vel%z%old(I,J,K)+p%loc%vel%z%old(I,J,K-1))/p%glb%dz**2.0d0
        
        wx = 0.5d0*( p%loc%vel%z%old(i+1,j,k)-p%loc%vel%z%old(i-1,j,k) )/p%glb%dx
        wy = 0.5d0*( p%loc%vel%z%old(i,j+1,k)-p%loc%vel%z%old(i,j-1,k) )/p%glb%dy
        wz = 0.5d0*( p%loc%vel%z%old(i,j,k+1)-p%loc%vel%z%old(i,j,k-1) )/p%glb%dz

        uz = 0.5d0*( p%loc%vel%x%old(i,j,k+1)-p%loc%vel%x%old(i,j,k)+p%loc%vel%x%old(i-1,j,k+1)-p%loc%vel%x%old(i-1,j,k) )/p%glb%dz
        vz = 0.5d0*( p%loc%vel%y%old(i,j,k+1)-p%loc%vel%y%old(i,j,k)+p%loc%vel%y%old(i,j-1,k+1)-p%loc%vel%y%old(i,j-1,k) )/p%glb%dz

        phix = 0.25d0*( p%loc%phi%old(i+1,j,k+1)-p%loc%phi%old(i-1,j,k+1) + p%loc%phi%old(i+1,j,k)-p%loc%phi%old(i-1,j,k) )/p%glb%dx
        phiy = 0.25d0*( p%loc%phi%old(i,j+1,k+1)-p%loc%phi%old(i,j-1,k+1) + p%loc%phi%old(i,j+1,k)-p%loc%phi%old(i,j-1,k) )/p%glb%dy
        phiz = ( p%loc%phi%old(i,j,k+1)-p%loc%phi%old(i,j,k) )/p%glb%dz
        
        dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*(uz+wx)/(rho*p%glb%re)
        dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*(vz+wy)/(rho*p%glb%re)
        dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*2.0d0*wz/(rho*p%glb%re)
        
        p%loc%velsrc%z%tmp(i,j,k) = dif_x + dif_y + dif_z &
        & + p%glb%gz*p%glb%btn_g / p%glb%fr &
        & - p%glb%btn_sf*curv*delta*phiz / (p%glb%we*rho) 
        
        ! ==========================================================================
        
        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i,j,k+1))
        mu = 0.5d0*(p%loc%mu%old(i,j,k)+p%loc%mu%old(i,j,k+1))
        delta = 0.5d0*(p%loc%delta%now(i,j,k)+p%loc%delta%now(i,j,k+1))
        curv = (p%loc%normals%curv%now(i,j,k)+p%loc%normals%curv%now(i,j+1,k))/2.0d0
        
        xx = (p%loc%vel%z%tmp(I+1,J,K)-2.0d0*p%loc%vel%z%tmp(I,J,K)+p%loc%vel%z%tmp(I-1,J,K))/p%glb%dx**2.0d0
        yy = (p%loc%vel%z%tmp(I,J+1,K)-2.0d0*p%loc%vel%z%tmp(I,J,K)+p%loc%vel%z%tmp(I,J-1,K))/p%glb%dy**2.0d0
        zz = (p%loc%vel%z%tmp(I,J,K+1)-2.0d0*p%loc%vel%z%tmp(I,J,K)+p%loc%vel%z%tmp(I,J,K-1))/p%glb%dz**2.0d0
        
        wx = 0.5d0*( p%loc%vel%z%tmp(i+1,j,k)-p%loc%vel%z%tmp(i-1,j,k) )/p%glb%dx
        wy = 0.5d0*( p%loc%vel%z%tmp(i,j+1,k)-p%loc%vel%z%tmp(i,j-1,k) )/p%glb%dy
        wz = 0.5d0*( p%loc%vel%z%tmp(i,j,k+1)-p%loc%vel%z%tmp(i,j,k-1) )/p%glb%dz

        uz = 0.5d0*( p%loc%vel%x%tmp(i,j,k+1)-p%loc%vel%x%tmp(i,j,k)+p%loc%vel%x%tmp(i-1,j,k+1)-p%loc%vel%x%tmp(i-1,j,k) )/p%glb%dz
        vz = 0.5d0*( p%loc%vel%y%tmp(i,j,k+1)-p%loc%vel%y%tmp(i,j,k)+p%loc%vel%y%tmp(i,j-1,k+1)-p%loc%vel%y%tmp(i,j-1,k) )/p%glb%dz

        phix = 0.25d0*( p%loc%phi%now(i+1,j,k+1)-p%loc%phi%now(i-1,j,k+1) + p%loc%phi%now(i+1,j,k)-p%loc%phi%now(i-1,j,k) )/p%glb%dx
        phiy = 0.25d0*( p%loc%phi%now(i,j+1,k+1)-p%loc%phi%now(i,j-1,k+1) + p%loc%phi%now(i,j+1,k)-p%loc%phi%now(i,j-1,k) )/p%glb%dy
        phiz = ( p%loc%phi%now(i,j,k+1)-p%loc%phi%now(i,j,k) )/p%glb%dz
        
        dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*(uz+wx)/(rho*p%glb%re)
        dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*(vz+wy)/(rho*p%glb%re)
        dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*2.0d0*wz/(rho*p%glb%re)
        
        p%loc%velsrc%z%tmp(i,j,k) = ( p%loc%velsrc%z%tmp(i,j,k) + &
        & dif_x + dif_y + dif_z + &
        & p%glb%gz*p%glb%btn_g / p%glb%fr &
        & - p%glb%btn_sf*curv*delta*phiz / (p%glb%we*rho) ) /2.0d0
        
    end do
    end do
    end do 
    !$omp end parallel do

end subroutine
