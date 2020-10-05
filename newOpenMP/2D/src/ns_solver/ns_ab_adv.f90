subroutine ns_ab_adv_source_sec
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: u,v,xp,xm,yp,ym
    
!$omp parallel do collapse(2), private(u,v,xp,xm,yp,ym)
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie
    
    u = p%loc%vel%x%tmp(i,j)
    v = p%loc%tdata%y%s1(i,j)
        
    xp = 0.5d0*(-p%loc%vel%x%old(i+2,j)+4.0d0*p%loc%vel%x%old(i+1,j)-3.0d0*p%loc%vel%x%old(i,j))/p%glb%dx 
    xm = 0.5d0*( p%loc%vel%x%old(i-2,j)-4.0d0*p%loc%vel%x%old(i-1,j)+3.0d0*p%loc%vel%x%old(i,j))/p%glb%dx
    
    yp = 0.5d0*(-p%loc%vel%x%old(i,j+2)+4.0d0*p%loc%vel%x%old(i,j+1)-3.0d0*p%loc%vel%x%old(i,j))/p%glb%dy
    ym = 0.5d0*( p%loc%vel%x%old(i,j-2)-4.0d0*p%loc%vel%x%old(i,j-1)+3.0d0*p%loc%vel%x%old(i,j))/p%glb%dy
            
    p%loc%velsrc%x%now(i,j) = - ((u+abs(u))*xm+(u-abs(u))*xp)/2.0d0 &
                           &  - ((v+abs(v))*ym+(v-abs(v))*yp)/2.0d0 
            
    !-----------------------------------------------------------
        
    u = p%loc%tdata%x%s1(i,j)
    v = p%loc%vel%y%tmp(i,j)
        
    xp = 0.5d0*(-p%loc%vel%y%old(i+2,j)+4.0d0*p%loc%vel%y%old(i+1,j)-3.0d0*p%loc%vel%y%old(i,j))/p%glb%dx 
    xm = 0.5d0*( p%loc%vel%y%old(i-2,j)-4.0d0*p%loc%vel%y%old(i-1,j)+3.0d0*p%loc%vel%y%old(i,j))/p%glb%dx
        
    yp = 0.5d0*(-p%loc%vel%y%old(i,j+2)+4.0d0*p%loc%vel%y%old(i,j+1)-3.0d0*p%loc%vel%y%old(i,j))/p%glb%dy
    ym = 0.5d0*( p%loc%vel%y%old(i,j-2)-4.0d0*p%loc%vel%y%old(i,j-1)+3.0d0*p%loc%vel%y%old(i,j))/p%glb%dy
     
    p%loc%velsrc%y%now(i,j) = - ((u+abs(u))*xm+(u-abs(u))*xp)/2.0d0 &
                           &  - ((v+abs(v))*ym+(v-abs(v))*yp)/2.0d0
                                      
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
integer :: id,i,j
real(8) :: uh, vh

!$omp parallel do collapse(2), private(uh,vh)
do j = p%loc%js-1, p%loc%je
do i = p%loc%is-1, p%loc%ie

    uh = 0.5d0*( p%loc%vel%x%tmp(i,j) + p%loc%vel%x%tmp(i+1,j) )
    vh = 0.5d0*( p%loc%vel%y%tmp(i,j) + p%loc%vel%y%tmp(i+1,j) )
    
    if( uh>=0.0d0 )then
        p%loc%tdata%x%s1(i,j) = ( -p%loc%vel%x%old(i-1,j)+6.0d0*p%loc%vel%x%old(i,j)+3.0d0*p%loc%vel%x%old(i+1,j) )/8.0d0
    else
        p%loc%tdata%x%s1(i,j) = ( -p%loc%vel%x%old(i+2,j)+6.0d0*p%loc%vel%x%old(i+1,j)+3.0d0*p%loc%vel%x%old(i,j) )/8.0d0
    endif
    
    if( vh>=0.0d0 )then
        p%loc%tdata%x%s2(i,j) = ( -p%loc%vel%x%old(i,j-1)+6.0d0*p%loc%vel%x%old(i,j)+3.0d0*p%loc%vel%x%old(i,j+1) )/8.0d0
    else
        p%loc%tdata%x%s2(i,j) = ( -p%loc%vel%x%old(i,j+2)+6.0d0*p%loc%vel%x%old(i,j+1)+3.0d0*p%loc%vel%x%old(i,j) )/8.0d0
    endif
    
    !-----------------------------------------------------------
    
    uh = 0.5d0*( p%loc%vel%x%tmp(i,j) + p%loc%vel%x%tmp(i,j+1) )
    vh = 0.5d0*( p%loc%vel%y%tmp(i,j) + p%loc%vel%y%tmp(i,j+1) )
    
    if( uh>=0.0d0 )then
        p%loc%tdata%y%s1(i,j) = ( -p%loc%vel%y%old(i-1,j)+6.0d0*p%loc%vel%y%old(i,j)+3.0d0*p%loc%vel%y%old(i+1,j) )/8.0d0
    else
        p%loc%tdata%y%s1(i,j) = ( -p%loc%vel%y%old(i+2,j)+6.0d0*p%loc%vel%y%old(i+1,j)+3.0d0*p%loc%vel%y%old(i,j) )/8.0d0
    endif
    
    if( vh>=0.0d0 )then
        p%loc%tdata%y%s2(i,j) = ( -p%loc%vel%y%old(i,j-1)+6.0d0*p%loc%vel%y%old(i,j)+3.0d0*p%loc%vel%y%old(i,j+1) )/8.0d0
    else
        p%loc%tdata%y%s2(i,j) = ( -p%loc%vel%y%old(i,j+2)+6.0d0*p%loc%vel%y%old(i,j+1)+3.0d0*p%loc%vel%y%old(i,j) )/8.0d0
    endif
    
end do
end do
!$omp end parallel do

!$omp parallel do collapse(2)
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie

    uh = 0.5d0*( p%loc%vel%x%tmp(i,j) + p%loc%vel%x%tmp(i+1,j) )
    vh = 0.5d0*( p%loc%vel%y%tmp(i,j) + p%loc%vel%y%tmp(i+1,j) )

    p%loc%velsrc%x%now(i,j) = - uh*( p%loc%tdata%x%s1(i,j)-p%loc%tdata%x%s1(i-1,j) )/p%glb%dx &
                               &- vh*( p%loc%tdata%x%s2(i,j)-p%loc%tdata%x%s2(i,j-1) )/p%glb%dy
    
    uh = 0.5d0*( p%loc%vel%x%tmp(i,j) + p%loc%vel%x%tmp(i,j+1) )
    vh = 0.5d0*( p%loc%vel%y%tmp(i,j) + p%loc%vel%y%tmp(i,j+1) )
    
    p%loc%velsrc%y%now(i,j) = - uh*( p%loc%tdata%y%s1(i,j)-p%loc%tdata%y%s1(i-1,j) )/p%glb%dx &
                               &- vh*( p%loc%tdata%y%s2(i,j)-p%loc%tdata%y%s2(i,j-1) )/p%glb%dy
                               
end do
end do
!$omp end parallel do
    
end subroutine