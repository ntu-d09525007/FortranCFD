subroutine ns_split_solver
use all
!$ use omp_lib
implicit none

p%glb%piter=0

call ppe_sor_init
call ns_split_adv
call ns_split_diff
call ns_check_convergence_vel
call ppe_sor_solver(p%glb%p_tol)
    
end subroutine

subroutine ns_split_adv
use all
!$ use omp_lib
implicit none
integer :: i,j,id,iter
integer(8) :: cpustart, cpuend
real(8) :: err


    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1
        
        call p%of(id)%loc%tdata%init(.true.,p%of(id)%loc%vel%x%old,p%of(id)%loc%vel%y%old)

        call p%of(id)%find_stag_vel( p%of(id)%loc%velsrc%x%now, p%of(id)%loc%velsrc%y%now,  &
                                    &p%of(id)%loc%vel%x%old   , p%of(id)%loc%vel%y%old  )

    enddo           
    !$omp end parallel do
    
    call pt%velsrc%sync
    
    iter=0
    
    do 
    
        iter = iter + 1
        
        call ns_split_adv_source()
        
        err=0.0_8  
        !$omp parallel do reduction(max:err)
        do id = 0, p%glb%threads-1
            call p%of(id)%loc%tdata%solve_srk6(err)
        enddo
        !$omp end parallel do
    
        if( err < p%glb%t_tol )exit
        
        if( mod(iter,100) == 0)write(*,*)"NS adv solver,",iter,err
        
    end do
    
    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1
        
        call p%of(id)%loc%tdata%final_srk6()

        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        
            p%of(id)%loc%vel%x%now(i,j) = p%of(id)%loc%tdata%x%target(i,j)
            p%of(id)%loc%vel%y%now(i,j) = p%of(id)%loc%tdata%y%target(i,j)

        end do
        end do

        call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now)

    enddo
    !$omp end parallel do

    call pt%vel%sync
        
end subroutine

subroutine ns_split_adv_source
use all
!$ use omp_lib
implicit none
integer :: id,i,j

    !$omp parallel do 
    do id = 0, p%glb%threads-1

        call p%of(id)%bc(0,p%of(id)%loc%tdata%x%s1);call p%of(id)%bc(0,p%of(id)%loc%tdata%x%s2);call p%of(id)%bc(0,p%of(id)%loc%tdata%x%s3)
        call p%of(id)%bc(0,p%of(id)%loc%tdata%y%s1);call p%of(id)%bc(0,p%of(id)%loc%tdata%y%s2);call p%of(id)%bc(0,p%of(id)%loc%tdata%y%s3)

    enddo           
    !$omp end parallel do
   
    call pt%tdatax%sync
    call pt%tdatay%sync
    
    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1

        do j = p%of(id)%loc%js, p%of(id)%loc%je

            ! ux ------------------------------------------------------------------------
            call p%of(id)%loc%ccdsolvers%x%solve("srkccd",p%of(id)%loc%tdata%x%s1(:,j),&
                                                &p%of(id)%loc%tdata%x%ss1(:,j),&
                                                &p%of(id)%loc%vel%x%tmp(:,j),&
                                                &p%of(id)%loc%vel%x%old(:,j))   

            call p%of(id)%loc%ccdsolvers%x%solve("srkccd",p%of(id)%loc%tdata%x%s2(:,j),&
                                                &p%of(id)%loc%tdata%x%ss2(:,j),&
                                                &p%of(id)%loc%vel%x%tmp(:,j),&
                                                &p%of(id)%loc%vel%x%old(:,j)) 
                                                
            call p%of(id)%loc%ccdsolvers%x%solve("srkccd",p%of(id)%loc%tdata%x%s3(:,j),&
                                                &p%of(id)%loc%tdata%x%ss3(:,j),&
                                                &p%of(id)%loc%vel%x%tmp(:,j),&
                                                &p%of(id)%loc%vel%x%old(:,j)) 
            
            ! vx ------------------------------------------------------------------------
            call p%of(id)%loc%ccdsolvers%x%solve("srkccd",p%of(id)%loc%tdata%y%s1(:,j),&
                                                &p%of(id)%loc%tdata%y%ss1(:,j),&
                                                &p%of(id)%loc%vel%x%tmp(:,j),&
                                                &p%of(id)%loc%velsrc%x%now(:,j))

            call p%of(id)%loc%ccdsolvers%x%solve("srkccd",p%of(id)%loc%tdata%y%s2(:,j),&
                                                &p%of(id)%loc%tdata%y%ss2(:,j),&
                                                &p%of(id)%loc%vel%x%tmp(:,j),&
                                                &p%of(id)%loc%velsrc%x%now(:,j))

            call p%of(id)%loc%ccdsolvers%x%solve("srkccd",p%of(id)%loc%tdata%y%s3(:,j),&
                                                &p%of(id)%loc%tdata%y%ss3(:,j),&
                                                &p%of(id)%loc%vel%x%tmp(:,j),&
                                                &p%of(id)%loc%velsrc%x%now(:,j))
        end do
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            p%of(id)%loc%tdata%x%l1(i,j) = - p%of(id)%loc%vel%x%old(i,j)*p%of(id)%loc%tdata%x%ss1(i,j)
            p%of(id)%loc%tdata%x%l2(i,j) = - p%of(id)%loc%vel%x%old(i,j)*p%of(id)%loc%tdata%x%ss2(i,j)
            p%of(id)%loc%tdata%x%l3(i,j) = - p%of(id)%loc%vel%x%old(i,j)*p%of(id)%loc%tdata%x%ss3(i,j)

            p%of(id)%loc%tdata%y%l1(i,j) = - p%of(id)%loc%velsrc%x%now(i,j)*p%of(id)%loc%tdata%y%ss1(i,j)
            p%of(id)%loc%tdata%y%l2(i,j) = - p%of(id)%loc%velsrc%x%now(i,j)*p%of(id)%loc%tdata%y%ss2(i,j)
            p%of(id)%loc%tdata%y%l3(i,j) = - p%of(id)%loc%velsrc%x%now(i,j)*p%of(id)%loc%tdata%y%ss3(i,j)

        end do
        end do
            
        do i = p%of(id)%loc%is, p%of(id)%loc%ie

            ! uy ------------------------------------------------------------------------
            call p%of(id)%loc%ccdsolvers%y%solve("srkccd",p%of(id)%loc%tdata%x%s1(i,:),&
                                                &p%of(id)%loc%tdata%x%ss1(i,:),&
                                                &p%of(id)%loc%vel%y%tmp(i,:),&
                                                &p%of(id)%loc%velsrc%y%now(i,:))

            call p%of(id)%loc%ccdsolvers%y%solve("srkccd",p%of(id)%loc%tdata%x%s2(i,:),&
                                                &p%of(id)%loc%tdata%x%ss2(i,:),&
                                                &p%of(id)%loc%vel%y%tmp(i,:),&
                                                &p%of(id)%loc%velsrc%y%now(i,:))

            call p%of(id)%loc%ccdsolvers%y%solve("srkccd",p%of(id)%loc%tdata%x%s3(i,:),&
                                                &p%of(id)%loc%tdata%x%ss3(i,:),&
                                                &p%of(id)%loc%vel%y%tmp(i,:),&
                                                &p%of(id)%loc%velsrc%y%now(i,:))

            ! vy ------------------------------------------------------------------------
            call p%of(id)%loc%ccdsolvers%y%solve("srkccd",p%of(id)%loc%tdata%y%s1(i,:),&
                                                &p%of(id)%loc%tdata%y%ss1(i,:),&
                                                &p%of(id)%loc%vel%y%tmp(i,:),&
                                                &p%of(id)%loc%vel%y%old(i,:))

            call p%of(id)%loc%ccdsolvers%y%solve("srkccd",p%of(id)%loc%tdata%y%s2(i,:),&
                                                &p%of(id)%loc%tdata%y%ss2(i,:),&
                                                &p%of(id)%loc%vel%y%tmp(i,:),&
                                                &p%of(id)%loc%vel%y%old(i,:))

            call p%of(id)%loc%ccdsolvers%y%solve("srkccd",p%of(id)%loc%tdata%y%s3(i,:),&
                                                &p%of(id)%loc%tdata%y%ss3(i,:),&
                                                &p%of(id)%loc%vel%y%tmp(i,:),&
                                                &p%of(id)%loc%vel%y%old(i,:))

        end do

        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            p%of(id)%loc%tdata%x%l1(i,j) = p%of(id)%loc%tdata%x%l1(i,j) - p%of(id)%loc%velsrc%y%now(i,j)*p%of(id)%loc%tdata%x%ss1(i,j)
            p%of(id)%loc%tdata%x%l2(i,j) = p%of(id)%loc%tdata%x%l2(i,j) - p%of(id)%loc%velsrc%y%now(i,j)*p%of(id)%loc%tdata%x%ss2(i,j)
            p%of(id)%loc%tdata%x%l3(i,j) = p%of(id)%loc%tdata%x%l3(i,j) - p%of(id)%loc%velsrc%y%now(i,j)*p%of(id)%loc%tdata%x%ss3(i,j)

            p%of(id)%loc%tdata%y%l1(i,j) = p%of(id)%loc%tdata%y%l1(i,j) - p%of(id)%loc%vel%y%old(i,j)*p%of(id)%loc%tdata%y%ss1(i,j)
            p%of(id)%loc%tdata%y%l2(i,j) = p%of(id)%loc%tdata%y%l2(i,j) - p%of(id)%loc%vel%y%old(i,j)*p%of(id)%loc%tdata%y%ss2(i,j)
            p%of(id)%loc%tdata%y%l3(i,j) = p%of(id)%loc%tdata%y%l3(i,j) - p%of(id)%loc%vel%y%old(i,j)*p%of(id)%loc%tdata%y%ss3(i,j)

        end do
        end do
     
    enddo           
    !$omp end parallel do
    
end subroutine

subroutine ns_split_diff
use all
!$ use omp_lib
implicit none
integer :: id,i,j,iter
real(8) :: err

!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1
   
   do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
   do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        p%of(id)%loc%vel%x%old2(i,j) = p%of(id)%loc%vel%x%now(i,j)
        p%of(id)%loc%vel%y%old2(i,j) = p%of(id)%loc%vel%y%now(i,j)
   enddo
   enddo
        
enddo    
!$omp end parallel do

!==============================================================================

do iter = 1, 5

call ns_linearize

!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1

   do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
   do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        p%of(id)%loc%vel%x%now(i,j) = 0.5d0*( p%of(id)%loc%vel%x%now(i,j) + p%of(id)%loc%vel%x%old2(i,j) ) 
        p%of(id)%loc%vel%y%now(i,j) = 0.5d0*( p%of(id)%loc%vel%y%now(i,j) + p%of(id)%loc%vel%y%old2(i,j) )
   enddo
   enddo

    call p%of(id)%find_stag_vel( p%of(id)%loc%velsrc%x%now, p%of(id)%loc%velsrc%y%now,  &
                                &p%of(id)%loc%vel%x%now   , p%of(id)%loc%vel%y%now  )

enddo
!$omp end parallel do

call pt%velsrc%sync

!$omp parallel do
do id = 0, p%glb%threads-1
    
    call ns_split_diff_source( p%of(id), &
                              &p%of(id)%loc%velsrc%x%now, p%of(id)%loc%velsrc%y%now, &
                              &p%of(id)%loc%vel%x%now   , p%of(id)%loc%vel%y%now   , &
                              &p%of(id)%loc%velsrc%x%tmp, p%of(id)%loc%velsrc%y%tmp, &
                              &p%of(id)%loc%rho%old,p%of(id)%loc%mu%old,p%of(id)%loc%delta%old,p%of(id)%loc%normals%curv%old)
    
enddo    
!$omp end parallel do


!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
        p%of(id)%loc%vel%x%now(i,j) = p%of(id)%loc%vel%x%old2(i,j) + p%of(id)%loc%velsrc%x%tmp(i,j) * p%glb%dt                      
        p%of(id)%loc%vel%y%now(i,j) = p%of(id)%loc%vel%y%old2(i,j) + p%of(id)%loc%velsrc%y%tmp(i,j) * p%glb%dt                      

    end do
    end do
        
    call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now)
    
enddo
!$omp end parallel do
    
call pt%vel%sync


end do

end subroutine

subroutine ns_split_diff_source(q,u,v,us,vs,sx,sy,irho,imu,idelta,icurv)
use all
implicit none
type(job) :: q
real(8), dimension(q%loc%is-q%glb%ghc:q%loc%ie+q%glb%ghc,&
                  &q%loc%js-q%glb%ghc:q%loc%je+q%glb%ghc) :: u,v,us,vs,sx,sy,irho,imu,idelta,icurv
real(8) :: rho,mu,delta,curv,xx,yy
real(8) :: ux,uy,vx,vy,phix,phiy
real(8) :: dif_x,dif_y             
integer :: i,j

! ---------------------------------------------- Calculate derivatives

    do j = q%loc%js, q%loc%je
        call q%loc%ccdsolvers%x%solve("srkccd",us(:,j),us(:,j),q%loc%tdata%x%s1(:,j),q%loc%tdata%x%ss1(:,j))
        call q%loc%ccdsolvers%x%solve("srkccd",vs(:,j),vs(:,j),q%loc%tdata%y%s1(:,j),q%loc%tdata%y%ss1(:,j))

        call q%loc%ccdsolvers%x%solve("srkccd",v(:,j),v(:,j),q%loc%tdata%y%l1(:,j))
    end do

    do i = q%loc%is, q%loc%ie
        call q%loc%ccdsolvers%y%solve("srkccd",us(i,:),us(i,:),q%loc%tdata%x%s2(i,:),q%loc%tdata%x%ss2(i,:))
        call q%loc%ccdsolvers%y%solve("srkccd",vs(i,:),vs(i,:),q%loc%tdata%y%s2(i,:),q%loc%tdata%y%ss2(i,:))

        call q%loc%ccdsolvers%y%solve("srkccd",u(i,:),u(i,:),q%loc%tdata%x%l1(i,:))
    end do
! ----------------------------------------------
    
    do j = q%loc%js, q%loc%je
    do i = q%loc%is, q%loc%ie
        
        rho = (irho(i,j)+irho(i+1,j))/2.0d0
        mu  = (imu(i,j)+imu(i+1,j))/2.0d0
        delta = (idelta(i,j)+idelta(i+1,j))/2.0d0
        curv = (icurv(i,j)+icurv(i+1,j))/2.0d0
        
        xx = q%loc%tdata%x%ss1(i,j); yy = q%loc%tdata%x%ss2(i,j); 
        ux = q%loc%tdata%x%s1(i,j) ; uy = q%loc%tdata%x%s2(i,j) ; 
        vx = q%loc%tdata%y%l1(i,j)
        
        phix = 0.5d0*( q%loc%normals%x%old(i,j)+q%loc%normals%x%old(i+1,j) )
        phiy = 0.5d0*( q%loc%normals%y%old(i,j)+q%loc%normals%y%old(i+1,j) )
 
        dif_x = mu/rho*xx/q%glb%re + (1.0d0-q%glb%mu_12)*delta*phix*2.0d0*ux/(rho*q%glb%re)
        dif_y = mu/rho*yy/q%glb%re + (1.0d0-q%glb%mu_12)*delta*phiy*(uy+vx)/(rho*q%glb%re)
     
        sx(i,j) = dif_x + dif_y + q%glb%gx*q%glb%btn_g/q%glb%fr - q%glb%btn_sf*curv*delta*phix / (q%glb%we*rho)
        
    end do
    end do

    do j = q%loc%js, q%loc%je
    do i = q%loc%is, q%loc%ie
        
        rho = (irho(i,j)+irho(i,j+1))/2.0d0
        mu  = (imu(i,j)+imu(i,j+1))/2.0d0
        delta = (idelta(i,j)+idelta(i,j+1))/2.0d0
        curv = (icurv(i,j)+icurv(i,j+1))/2.0d0
   
        xx = q%loc%tdata%y%ss1(i,j); yy = q%loc%tdata%y%ss2(i,j)
        vx = q%loc%tdata%y%s1(i,j) ; vy = q%loc%tdata%y%s2(i,j) 
        uy = q%loc%tdata%x%l1(i,j)
        
        phix = 0.5d0*( q%loc%normals%x%old(i,j)+q%loc%normals%x%old(i,j+1) )
        phiy = 0.5d0*( q%loc%normals%y%old(i,j)+q%loc%normals%y%old(i,j+1) )
   
        dif_x = mu/rho*xx/q%glb%re + (1.0d0-q%glb%mu_12)*delta*phix*(uy+vx)/(rho*q%glb%re)
        dif_y = mu/rho*yy/q%glb%re + (1.0d0-q%glb%mu_12)*delta*phiy*2.0d0*vy/(rho*q%glb%re)
     
        sy(i,j) = dif_x + dif_y + q%glb%gy*q%glb%btn_g/q%glb%fr - q%glb%btn_sf*curv*delta*phiy / (q%glb%we*rho)
        
    end do
    end do
    
end subroutine


