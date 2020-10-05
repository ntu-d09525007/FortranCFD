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

!$omp parallel do private(i,j,k)
do id = 0, p%glb%threads-1
        
    call p%of(id)%find_stag_vel( p%of(id)%loc%tdata%x%s1, p%of(id)%loc%tdata%y%s1, p%of(id)%loc%tdata%z%s1, &
                                &p%of(id)%loc%tdata%x%s2, p%of(id)%loc%tdata%y%s2, p%of(id)%loc%tdata%z%s2, &
                                &p%of(id)%loc%vel%x%old, p%of(id)%loc%vel%y%old, p%of(id)%loc%vel%z%old )
enddo       
!$omp end parallel do
    
call pt%tdatax%sync
call pt%tdatay%sync
call pt%tdataz%sync
    
!$omp parallel do private(i,j,k,u,v,w,xp,xm,yp,ym,zp,zm)
do id = 0, p%glb%threads-1
        
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

enddo   
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

    !$omp parallel do private(i,j,k,uh,vh,wh)
    do id = 0, p%glb%threads-1
        
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
    
    enddo
    !$omp end parallel do
    
end subroutine

subroutine ns_ab_uccd()
use all
implicit none
integer :: id,i,j,k

!$omp parallel do
do id = 0, p%glb%threads-1

enddo
!$omp end parallel do

end subroutine