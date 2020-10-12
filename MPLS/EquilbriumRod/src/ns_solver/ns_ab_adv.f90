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
real(8) :: u,v,xp,xm,yp,ym,zp,zm

!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1
        
    call p%of(id)%find_stag_vel( p%of(id)%loc%tdata%x%s1, p%of(id)%loc%tdata%y%s1,&
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

subroutine ns_ab_adv_source_quick
!--------------------------------------
! u*Px = [ u_{i}+u_{i+1} ] * [ P_{i+1/2}-P_{i-1/2} ] / dx
!--------------------------------------
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: uh,vh

    !$omp parallel do private(i,j,uh,vh)
    do id = 0, p%glb%threads-1

        do j = p%of(id)%loc%js-1, p%of(id)%loc%je
        do i = p%of(id)%loc%is-1, p%of(id)%loc%ie
        
            uh = 0.5d0*( p%of(id)%loc%vel%x%old(i,j) + p%of(id)%loc%vel%x%old(i+1,j) )
            vh = 0.5d0*( p%of(id)%loc%vel%y%old(i,j) + p%of(id)%loc%vel%y%old(i+1,j) )

            if( uh>=0.0d0 )then
                p%of(id)%loc%tdata%x%s1(i,j) = ( -p%of(id)%loc%vel%x%old(i-1,j)+6.0d0*p%of(id)%loc%vel%x%old(i,j)+3.0d0*p%of(id)%loc%vel%x%old(i+1,j) )/8.0d0
            else
                p%of(id)%loc%tdata%x%s1(i,j) = ( -p%of(id)%loc%vel%x%old(i+2,j)+6.0d0*p%of(id)%loc%vel%x%old(i+1,j)+3.0d0*p%of(id)%loc%vel%x%old(i,j) )/8.0d0
            endif
            
            if( vh>=0.0d0 )then
                p%of(id)%loc%tdata%x%s2(i,j) = ( -p%of(id)%loc%vel%x%old(i,j-1)+6.0d0*p%of(id)%loc%vel%x%old(i,j)+3.0d0*p%of(id)%loc%vel%x%old(i,j+1) )/8.0d0
            else
                p%of(id)%loc%tdata%x%s2(i,j) = ( -p%of(id)%loc%vel%x%old(i,j+2)+6.0d0*p%of(id)%loc%vel%x%old(i,j+1)+3.0d0*p%of(id)%loc%vel%x%old(i,j) )/8.0d0
            endif
            
            !-----------------------------------------------------------
            
            uh = 0.5d0*( p%of(id)%loc%vel%x%old(i,j) + p%of(id)%loc%vel%x%old(i,j+1) )
            vh = 0.5d0*( p%of(id)%loc%vel%y%old(i,j) + p%of(id)%loc%vel%y%old(i,j+1) )

            if( uh>=0.0d0 )then
                p%of(id)%loc%tdata%y%s1(i,j) = ( -p%of(id)%loc%vel%y%old(i-1,j)+6.0d0*p%of(id)%loc%vel%y%old(i,j)+3.0d0*p%of(id)%loc%vel%y%old(i+1,j) )/8.0d0
            else
                p%of(id)%loc%tdata%y%s1(i,j) = ( -p%of(id)%loc%vel%y%old(i+2,j)+6.0d0*p%of(id)%loc%vel%y%old(i+1,j)+3.0d0*p%of(id)%loc%vel%y%old(i,j) )/8.0d0
            endif
            
            if( vh>=0.0d0 )then
                p%of(id)%loc%tdata%y%s2(i,j) = ( -p%of(id)%loc%vel%y%old(i,j-1)+6.0d0*p%of(id)%loc%vel%y%old(i,j)+3.0d0*p%of(id)%loc%vel%y%old(i,j+1) )/8.0d0
            else
                p%of(id)%loc%tdata%y%s2(i,j) = ( -p%of(id)%loc%vel%y%old(i,j+2)+6.0d0*p%of(id)%loc%vel%y%old(i,j+1)+3.0d0*p%of(id)%loc%vel%y%old(i,j) )/8.0d0
            endif
            
        end do
        end do
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            uh = 0.5d0*( p%of(id)%loc%vel%x%old(i,j) + p%of(id)%loc%vel%x%old(i+1,j) )
            vh = 0.5d0*( p%of(id)%loc%vel%y%old(i,j) + p%of(id)%loc%vel%y%old(i+1,j) )

            p%of(id)%loc%velsrc%x%now(i,j) = - uh*( p%of(id)%loc%tdata%x%s1(i,j)-p%of(id)%loc%tdata%x%s1(i-1,j) )/p%glb%dx &
                                            &- vh*( p%of(id)%loc%tdata%x%s2(i,j)-p%of(id)%loc%tdata%x%s2(i,j-1) )/p%glb%dy 
                      
            uh = 0.5d0*( p%of(id)%loc%vel%x%old(i,j) + p%of(id)%loc%vel%x%old(i,j+1) )
            vh = 0.5d0*( p%of(id)%loc%vel%y%old(i,j) + p%of(id)%loc%vel%y%old(i,j+1) )
 
            p%of(id)%loc%velsrc%y%now(i,j) = - uh*( p%of(id)%loc%tdata%y%s1(i,j)-p%of(id)%loc%tdata%y%s1(i-1,j) )/p%glb%dx &
                                            &- vh*( p%of(id)%loc%tdata%y%s2(i,j)-p%of(id)%loc%tdata%y%s2(i,j-1) )/p%glb%dy 
                                             
        end do
        end do
    
    enddo
    !$omp end parallel do
    
end subroutine

subroutine ns_ab_adv_source_uccd()
use all
implicit none
integer :: id,i,j
real(8) :: u,v

!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1
        
    call p%of(id)%find_stag_vel( p%of(id)%loc%tdata%x%s1, p%of(id)%loc%tdata%y%s1, &
                                &p%of(id)%loc%vel%x%old, p%of(id)%loc%vel%y%old )
enddo       
!$omp end parallel do

call pt%tdatax%sync
call pt%tdatay%sync

!$omp parallel do private(u,v,i,j)
do id = 0, p%glb%threads-1

    do j = p%of(id)%loc%js, p%of(id)%loc%je

        call p%of(id)%loc%ccdsolvers%x%solve("uccd",p%of(id)%loc%vel%x%old(:,j),&
                                            &p%of(id)%loc%vel_ten%xx(:,j),&
                                            &p%of(id)%loc%vel_ten%xxx(:,j),&
                                            &p%of(id)%loc%vel%x%old(:,j)) 

        call p%of(id)%loc%ccdsolvers%x%solve("uccd",p%of(id)%loc%vel%y%old(:,j),&
                                            &p%of(id)%loc%vel_ten%yx(:,j),&
                                            &p%of(id)%loc%vel_ten%yxx(:,j),&
                                            &p%of(id)%loc%tdata%x%s1(:,j)) 

    enddo

enddo
!$omp end parallel do

!$omp parallel do private(u,v,i,j)
do id = 0, p%glb%threads-1

    do j = p%of(id)%loc%js, p%of(id)%loc%je

        call p%of(id)%loc%ccdsolvers%y%solve("uccd",p%of(id)%loc%vel%x%old(i,:),&
                                            &p%of(id)%loc%vel_ten%xy(i,:),&
                                            &p%of(id)%loc%vel_ten%xyy(i,:),&
                                            &p%of(id)%loc%tdata%y%s1(i,:))

        call p%of(id)%loc%ccdsolvers%y%solve("uccd",p%of(id)%loc%vel%y%old(i,:),&
                                            &p%of(id)%loc%vel_ten%yy(i,:),&
                                            &p%of(id)%loc%vel_ten%yyy(i,:),&
                                            &p%of(id)%loc%vel%y%old(i,:))

    enddo

enddo
!$omp end parallel do


!$omp parallel do private(u,v,i,j)
do id = 0, p%glb%threads-1

    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie

        u = p%of(id)%loc%vel%x%old(i,j)
        v = p%of(id)%loc%tdata%y%s1(i,j)

        p%of(id)%loc%velsrc%x%now(i,j) = - u * p%of(id)%loc%vel_ten%xx(i,j) &
                                      &  - v * p%of(id)%loc%vel_ten%xy(i,j)

        u = p%of(id)%loc%tdata%x%s1(i,j)
        v = p%of(id)%loc%vel%y%old(i,j)

        p%of(id)%loc%velsrc%y%now(i,j) = - u * p%of(id)%loc%vel_ten%yx(i,j) &
                                      &  - v * p%of(id)%loc%vel_ten%yy(i,j)

    enddo
    enddo

enddo
!$omp end parallel do

end subroutine