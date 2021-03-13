subroutine ns_ab_adv_source
implicit none

call ns_ab_adv_source_sec
!call ns_ab_adv_source_quick
!call ns_ab_adv_source_uccd
    
end subroutine

subroutine ns_ab_adv_source_sec
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: u,v,xp,xm,yp,ym

!$omp parallel do 
do id = 0, p%glb%threads-1
        
    call p%of(id)%find_stag_vel( p%of(id)%loc%tdata%x%s1, p%of(id)%loc%tdata%y%s1, &
                                &p%of(id)%loc%vel%x%old , p%of(id)%loc%vel%y%old )
enddo       
!$omp end parallel do
    
call pt%tdatax%sync
call pt%tdatay%sync
    
!$omp parallel do private(i,j,u,v,xp,xm,yp,ym)
do id = 0, p%glb%threads-1
        
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
        u = p%of(id)%loc%vel%x%old(i,j)
        v = p%of(id)%loc%tdata%y%s1(i,j)

        xp = 0.5d0*(-p%of(id)%loc%vel%x%old(i+2,j)+4.0d0*p%of(id)%loc%vel%x%old(i+1,j)-3.0d0*p%of(id)%loc%vel%x%old(i,j))/p%glb%dx 
        xm = 0.5d0*( p%of(id)%loc%vel%x%old(i-2,j)-4.0d0*p%of(id)%loc%vel%x%old(i-1,j)+3.0d0*p%of(id)%loc%vel%x%old(i,j))/p%glb%dx
        
        yp = 0.5d0*(-p%of(id)%loc%vel%x%old(i,j+2)+4.0d0*p%of(id)%loc%vel%x%old(i,j+1)-3.0d0*p%of(id)%loc%vel%x%old(i,j))/p%glb%dy
        ym = 0.5d0*( p%of(id)%loc%vel%x%old(i,j-2)-4.0d0*p%of(id)%loc%vel%x%old(i,j-1)+3.0d0*p%of(id)%loc%vel%x%old(i,j))/p%glb%dy
            
        p%of(id)%loc%velsrc%x%now(i,j) = - ((u+abs(u))*xm+(u-abs(u))*xp)/2.0d0 &
                                        &  - ((v+abs(v))*ym+(v-abs(v))*yp)/2.0d0 
                
        !-----------------------------------------------------------
            
        u = p%of(id)%loc%tdata%x%s1(i,j)
        v = p%of(id)%loc%vel%y%old(i,j)
            
        xp = 0.5d0*(-p%of(id)%loc%vel%y%old(i+2,j)+4.0d0*p%of(id)%loc%vel%y%old(i+1,j)-3.0d0*p%of(id)%loc%vel%y%old(i,j))/p%glb%dx 
        xm = 0.5d0*( p%of(id)%loc%vel%y%old(i-2,j)-4.0d0*p%of(id)%loc%vel%y%old(i-1,j)+3.0d0*p%of(id)%loc%vel%y%old(i,j))/p%glb%dx
            
        yp = 0.5d0*(-p%of(id)%loc%vel%y%old(i,j+2)+4.0d0*p%of(id)%loc%vel%y%old(i,j+1)-3.0d0*p%of(id)%loc%vel%y%old(i,j))/p%glb%dy
        ym = 0.5d0*( p%of(id)%loc%vel%y%old(i,j-2)-4.0d0*p%of(id)%loc%vel%y%old(i,j-1)+3.0d0*p%of(id)%loc%vel%y%old(i,j))/p%glb%dy
            
        p%of(id)%loc%velsrc%y%now(i,j) = - ((u+abs(u))*xm+(u-abs(u))*xp)/2.0d0 &
                                        &  - ((v+abs(v))*ym+(v-abs(v))*yp)/2.0d0
            
    end do 
    end do

enddo   
!$omp end parallel do
    
end subroutine