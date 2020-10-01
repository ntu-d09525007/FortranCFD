subroutine ns_split_solver
use all
!$ use omp_lib
implicit none

p%glb%piter=0

call ppe_sor_init
call ns_split_adv
call ns_split_diff
call ns_check_convergence_vel
!call ppe_sor_solver(p%glb%p_tol)
call ppe_mg_solver(0)
    
end subroutine

subroutine ns_split_adv
use all
!$ use omp_lib
implicit none
integer :: i,j,k,id,iter
integer(8) :: cpustart, cpuend
real(8) :: err


    !$omp parallel do private(i,j,k)
    do id = 0, p%glb%threads-1
        
        call p%of(id)%loc%tdata%init(.true.,p%of(id)%loc%vel%x%old,p%of(id)%loc%vel%y%old,p%of(id)%loc%vel%z%old)

        call p%of(id)%find_stag_vel( p%of(id)%loc%velsrc%x%now, p%of(id)%loc%velsrc%y%now, p%of(id)%loc%velsrc%z%now, &
                                    &p%of(id)%loc%velsrc%x%old, p%of(id)%loc%velsrc%y%old, p%of(id)%loc%velsrc%z%old, &
                                    &p%of(id)%loc%vel%x%old   , p%of(id)%loc%vel%y%old   , p%of(id)%loc%vel%z%old )
    enddo           
    !$omp end parallel do
    
    call pt%velsrc%sync
    call pt%velsrc_old%sync
    
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
    
    !$omp parallel do private(i,j,k)
    do id = 0, p%glb%threads-1
        
        call p%of(id)%loc%tdata%final_srk6()
        
        do k = p%of(id)%loc%ks-p%glb%ghc, p%of(id)%loc%ke+p%glb%ghc
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        
            p%of(id)%loc%vel%x%now(i,j,k) = p%of(id)%loc%tdata%x%target(i,j,k)
            p%of(id)%loc%vel%y%now(i,j,k) = p%of(id)%loc%tdata%y%target(i,j,k)
            p%of(id)%loc%vel%z%now(i,j,k) = p%of(id)%loc%tdata%z%target(i,j,k)
            
        end do
        end do
        end do

        call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now,p%of(id)%loc%vel%z%now)

        enddo
    !$omp end parallel do

    call pt%vel%sync
        
end subroutine

subroutine ns_split_adv_source
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k

    !$omp parallel do 
    do id = 0, p%glb%threads-1

        call p%of(id)%bc(0,p%of(id)%loc%tdata%x%s1);call p%of(id)%bc(0,p%of(id)%loc%tdata%x%s2);call p%of(id)%bc(0,p%of(id)%loc%tdata%x%s3)
        call p%of(id)%bc(0,p%of(id)%loc%tdata%y%s1);call p%of(id)%bc(0,p%of(id)%loc%tdata%y%s2);call p%of(id)%bc(0,p%of(id)%loc%tdata%y%s3)
        call p%of(id)%bc(0,p%of(id)%loc%tdata%z%s1);call p%of(id)%bc(0,p%of(id)%loc%tdata%z%s2);call p%of(id)%bc(0,p%of(id)%loc%tdata%z%s3)
      
    enddo           
    !$omp end parallel do
   
    call pt%tdatax%sync
    call pt%tdatay%sync
    call pt%tdataz%sync
    
    !$omp parallel do private(i,j,k)
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je

            ! ux ------------------------------------------------------------------------
            call p%of(id)%loc%ccd%x%SOLVE_SRKCCD(p%of(id)%loc%vel%x%old(:,j,k),&
                                                &p%of(id)%loc%tdata%x%s1(:,j,k),&
                                                &p%of(id)%loc%tdata%x%ss1(:,j,k))   

            call p%of(id)%loc%ccd%x%SOLVE_SRKCCD(p%of(id)%loc%vel%x%old(:,j,k),&
                                                &p%of(id)%loc%tdata%x%s2(:,j,k),&
                                                &p%of(id)%loc%tdata%x%ss2(:,j,k))
                                                
            call p%of(id)%loc%ccd%x%SOLVE_SRKCCD(p%of(id)%loc%vel%x%old(:,j,k),&
                                                &p%of(id)%loc%tdata%x%s3(:,j,k),&
                                                &p%of(id)%loc%tdata%x%ss3(:,j,k))
            
            ! vx ------------------------------------------------------------------------
            call p%of(id)%loc%ccd%x%SOLVE_SRKCCD(p%of(id)%loc%velsrc%x%now(:,j,k),&
                                                &p%of(id)%loc%tdata%y%s1(:,j,k),&
                                                &p%of(id)%loc%tdata%y%ss1(:,j,k))

            call p%of(id)%loc%ccd%x%SOLVE_SRKCCD(p%of(id)%loc%velsrc%x%now(:,j,k),&
                                                &p%of(id)%loc%tdata%y%s2(:,j,k),&
                                                &p%of(id)%loc%tdata%y%ss2(:,j,k))

            call p%of(id)%loc%ccd%x%SOLVE_SRKCCD(p%of(id)%loc%velsrc%x%now(:,j,k),&
                                                &p%of(id)%loc%tdata%y%s3(:,j,k),&
                                                &p%of(id)%loc%tdata%y%ss3(:,j,k))

            ! wx ------------------------------------------------------------------------        
            call p%of(id)%loc%ccd%x%SOLVE_SRKCCD(p%of(id)%loc%velsrc%x%old(:,j,k),&
                                                &p%of(id)%loc%tdata%z%s1(:,j,k),&
                                                &p%of(id)%loc%tdata%z%ss1(:,j,k))   

            call p%of(id)%loc%ccd%x%SOLVE_SRKCCD(p%of(id)%loc%velsrc%x%old(:,j,k),&
                                                &p%of(id)%loc%tdata%z%s2(:,j,k),&
                                                &p%of(id)%loc%tdata%z%ss2(:,j,k))

            call p%of(id)%loc%ccd%x%SOLVE_SRKCCD(p%of(id)%loc%velsrc%x%old(:,j,k),&
                                                &p%of(id)%loc%tdata%z%s3(:,j,k),&
                                                &p%of(id)%loc%tdata%z%ss3(:,j,k))
            
        end do
        end do
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            p%of(id)%loc%tdata%x%l1(i,j,k) = - p%of(id)%loc%vel%x%old(i,j,k)*p%of(id)%loc%tdata%x%ss1(i,j,k)
            p%of(id)%loc%tdata%x%l2(i,j,k) = - p%of(id)%loc%vel%x%old(i,j,k)*p%of(id)%loc%tdata%x%ss2(i,j,k)
            p%of(id)%loc%tdata%x%l3(i,j,k) = - p%of(id)%loc%vel%x%old(i,j,k)*p%of(id)%loc%tdata%x%ss3(i,j,k)

            p%of(id)%loc%tdata%y%l1(i,j,k) = - p%of(id)%loc%velsrc%x%now(i,j,k)*p%of(id)%loc%tdata%y%ss1(i,j,k)
            p%of(id)%loc%tdata%y%l2(i,j,k) = - p%of(id)%loc%velsrc%x%now(i,j,k)*p%of(id)%loc%tdata%y%ss2(i,j,k)
            p%of(id)%loc%tdata%y%l3(i,j,k) = - p%of(id)%loc%velsrc%x%now(i,j,k)*p%of(id)%loc%tdata%y%ss3(i,j,k)

            p%of(id)%loc%tdata%z%l1(i,j,k) = - p%of(id)%loc%velsrc%x%old(i,j,k)*p%of(id)%loc%tdata%z%ss1(i,j,k)
            p%of(id)%loc%tdata%z%l2(i,j,k) = - p%of(id)%loc%velsrc%x%old(i,j,k)*p%of(id)%loc%tdata%z%ss2(i,j,k)
            p%of(id)%loc%tdata%z%l3(i,j,k) = - p%of(id)%loc%velsrc%x%old(i,j,k)*p%of(id)%loc%tdata%z%ss3(i,j,k)
            
        end do
        end do
        end do
            
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do i = p%of(id)%loc%is, p%of(id)%loc%ie

            ! uy ------------------------------------------------------------------------
            call p%of(id)%loc%ccd%y%SOLVE_SRKCCD(p%of(id)%loc%velsrc%y%now(i,:,k),&
                                                &p%of(id)%loc%tdata%x%s1(i,:,k),&
                                                &p%of(id)%loc%tdata%x%ss1(i,:,k))   

            call p%of(id)%loc%ccd%y%SOLVE_SRKCCD(p%of(id)%loc%velsrc%y%now(i,:,k),&
                                                &p%of(id)%loc%tdata%x%s2(i,:,k),&
                                                &p%of(id)%loc%tdata%x%ss2(i,:,k))
                                                
            call p%of(id)%loc%ccd%y%SOLVE_SRKCCD(p%of(id)%loc%velsrc%y%now(i,:,k),&
                                                &p%of(id)%loc%tdata%x%s3(i,:,k),&
                                                &p%of(id)%loc%tdata%x%ss3(i,:,k))

            ! vy ------------------------------------------------------------------------
            call p%of(id)%loc%ccd%y%SOLVE_SRKCCD(p%of(id)%loc%vel%y%old(i,:,k),&
                                                &p%of(id)%loc%tdata%y%s1(i,:,k),&
                                                &p%of(id)%loc%tdata%y%ss1(i,:,k))   

            call p%of(id)%loc%ccd%y%SOLVE_SRKCCD(p%of(id)%loc%vel%y%old(i,:,k),&
                                                &p%of(id)%loc%tdata%y%s2(i,:,k),&
                                                &p%of(id)%loc%tdata%y%ss2(i,:,k))
                                                
            call p%of(id)%loc%ccd%y%SOLVE_SRKCCD(p%of(id)%loc%vel%y%old(i,:,k),&
                                                &p%of(id)%loc%tdata%y%s3(i,:,k),&
                                                &p%of(id)%loc%tdata%y%ss3(i,:,k))

            ! wy ------------------------------------------------------------------------
            call p%of(id)%loc%ccd%y%SOLVE_SRKCCD(p%of(id)%loc%velsrc%y%old(i,:,k),&
                                                &p%of(id)%loc%tdata%z%s1(i,:,k),&
                                                &p%of(id)%loc%tdata%z%ss1(i,:,k))   

            call p%of(id)%loc%ccd%y%SOLVE_SRKCCD(p%of(id)%loc%velsrc%y%old(i,:,k),&
                                                &p%of(id)%loc%tdata%z%s2(i,:,k),&
                                                &p%of(id)%loc%tdata%z%ss2(i,:,k))
                                                
            call p%of(id)%loc%ccd%y%SOLVE_SRKCCD(p%of(id)%loc%velsrc%y%old(i,:,k),&
                                                &p%of(id)%loc%tdata%z%s3(i,:,k),&
                                                &p%of(id)%loc%tdata%z%ss3(i,:,k))
                                                
        end do
        end do
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            p%of(id)%loc%tdata%x%l1(i,j,k) = p%of(id)%loc%tdata%x%l1(i,j,k) - p%of(id)%loc%velsrc%y%now(i,j,k)*p%of(id)%loc%tdata%x%ss1(i,j,k)
            p%of(id)%loc%tdata%x%l2(i,j,k) = p%of(id)%loc%tdata%x%l2(i,j,k) - p%of(id)%loc%velsrc%y%now(i,j,k)*p%of(id)%loc%tdata%x%ss2(i,j,k)
            p%of(id)%loc%tdata%x%l3(i,j,k) = p%of(id)%loc%tdata%x%l3(i,j,k) - p%of(id)%loc%velsrc%y%now(i,j,k)*p%of(id)%loc%tdata%x%ss3(i,j,k)

            p%of(id)%loc%tdata%y%l1(i,j,k) = p%of(id)%loc%tdata%y%l1(i,j,k) - p%of(id)%loc%vel%y%old(i,j,k)*p%of(id)%loc%tdata%y%ss1(i,j,k)
            p%of(id)%loc%tdata%y%l2(i,j,k) = p%of(id)%loc%tdata%y%l2(i,j,k) - p%of(id)%loc%vel%y%old(i,j,k)*p%of(id)%loc%tdata%y%ss2(i,j,k)
            p%of(id)%loc%tdata%y%l3(i,j,k) = p%of(id)%loc%tdata%y%l3(i,j,k) - p%of(id)%loc%vel%y%old(i,j,k)*p%of(id)%loc%tdata%y%ss3(i,j,k)
            
            p%of(id)%loc%tdata%z%l1(i,j,k) = p%of(id)%loc%tdata%z%l1(i,j,k) - p%of(id)%loc%velsrc%y%old(i,j,k)*p%of(id)%loc%tdata%z%ss1(i,j,k)
            p%of(id)%loc%tdata%z%l2(i,j,k) = p%of(id)%loc%tdata%z%l2(i,j,k) - p%of(id)%loc%velsrc%y%old(i,j,k)*p%of(id)%loc%tdata%z%ss2(i,j,k)
            p%of(id)%loc%tdata%z%l3(i,j,k) = p%of(id)%loc%tdata%z%l3(i,j,k) - p%of(id)%loc%velsrc%y%old(i,j,k)*p%of(id)%loc%tdata%z%ss3(i,j,k)
            
        end do
        end do
        end do

        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie

            ! uz ------------------------------------------------------------------------
            call p%of(id)%loc%ccd%z%SOLVE_SRKCCD(p%of(id)%loc%velsrc%z%now(i,j,:),&
                                                &p%of(id)%loc%tdata%x%s1(i,j,:),&
                                                &p%of(id)%loc%tdata%x%ss1(i,j,:))   

            call p%of(id)%loc%ccd%z%SOLVE_SRKCCD(p%of(id)%loc%velsrc%z%now(i,j,:),&
                                                &p%of(id)%loc%tdata%x%s2(i,j,:),&
                                                &p%of(id)%loc%tdata%x%ss2(i,j,:))
                                                
            call p%of(id)%loc%ccd%z%SOLVE_SRKCCD(p%of(id)%loc%velsrc%z%now(i,j,:),&
                                                &p%of(id)%loc%tdata%x%s3(i,j,:),&
                                                &p%of(id)%loc%tdata%x%ss3(i,j,:))   

            ! vz ------------------------------------------------------------------------
            call p%of(id)%loc%ccd%z%SOLVE_SRKCCD(p%of(id)%loc%velsrc%z%old(i,j,:),&
                                                &p%of(id)%loc%tdata%y%s1(i,j,:),&
                                                &p%of(id)%loc%tdata%y%ss1(i,j,:))   

            call p%of(id)%loc%ccd%z%SOLVE_SRKCCD(p%of(id)%loc%velsrc%z%old(i,j,:),&
                                                &p%of(id)%loc%tdata%y%s2(i,j,:),&
                                                &p%of(id)%loc%tdata%y%ss2(i,j,:))
                                                
            call p%of(id)%loc%ccd%z%SOLVE_SRKCCD(p%of(id)%loc%velsrc%z%old(i,j,:),&
                                                &p%of(id)%loc%tdata%y%s3(i,j,:),&
                                                &p%of(id)%loc%tdata%y%ss3(i,j,:))

            ! wz ------------------------------------------------------------------------
            call p%of(id)%loc%ccd%z%SOLVE_SRKCCD(p%of(id)%loc%vel%z%old(i,j,:),&
                                                &p%of(id)%loc%tdata%z%s1(i,j,:),&
                                                &p%of(id)%loc%tdata%z%ss1(i,j,:))   

            call p%of(id)%loc%ccd%z%SOLVE_SRKCCD(p%of(id)%loc%vel%z%old(i,j,:),&
                                                &p%of(id)%loc%tdata%z%s2(i,j,:),&
                                                &p%of(id)%loc%tdata%z%ss2(i,j,:))
                                                
            call p%of(id)%loc%ccd%z%SOLVE_SRKCCD(p%of(id)%loc%vel%z%old(i,j,:),&
                                                &p%of(id)%loc%tdata%z%s3(i,j,:),&
                                                &p%of(id)%loc%tdata%z%ss3(i,j,:))
                                                
        end do
        end do

        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            p%of(id)%loc%tdata%x%l1(i,j,k) = p%of(id)%loc%tdata%x%l1(i,j,k) - p%of(id)%loc%velsrc%z%now(i,j,k)*p%of(id)%loc%tdata%x%ss1(i,j,k)
            p%of(id)%loc%tdata%x%l2(i,j,k) = p%of(id)%loc%tdata%x%l2(i,j,k) - p%of(id)%loc%velsrc%z%now(i,j,k)*p%of(id)%loc%tdata%x%ss2(i,j,k)
            p%of(id)%loc%tdata%x%l3(i,j,k) = p%of(id)%loc%tdata%x%l3(i,j,k) - p%of(id)%loc%velsrc%z%now(i,j,k)*p%of(id)%loc%tdata%x%ss3(i,j,k)

            p%of(id)%loc%tdata%y%l1(i,j,k) = p%of(id)%loc%tdata%y%l1(i,j,k) - p%of(id)%loc%velsrc%z%old(i,j,k)*p%of(id)%loc%tdata%y%ss1(i,j,k)
            p%of(id)%loc%tdata%y%l2(i,j,k) = p%of(id)%loc%tdata%y%l2(i,j,k) - p%of(id)%loc%velsrc%z%old(i,j,k)*p%of(id)%loc%tdata%y%ss2(i,j,k)
            p%of(id)%loc%tdata%y%l3(i,j,k) = p%of(id)%loc%tdata%y%l3(i,j,k) - p%of(id)%loc%velsrc%z%old(i,j,k)*p%of(id)%loc%tdata%y%ss3(i,j,k)
            
            p%of(id)%loc%tdata%z%l1(i,j,k) = p%of(id)%loc%tdata%z%l1(i,j,k) - p%of(id)%loc%vel%z%old(i,j,k)*p%of(id)%loc%tdata%z%ss1(i,j,k)
            p%of(id)%loc%tdata%z%l2(i,j,k) = p%of(id)%loc%tdata%z%l2(i,j,k) - p%of(id)%loc%vel%z%old(i,j,k)*p%of(id)%loc%tdata%z%ss2(i,j,k)
            p%of(id)%loc%tdata%z%l3(i,j,k) = p%of(id)%loc%tdata%z%l3(i,j,k) - p%of(id)%loc%vel%z%old(i,j,k)*p%of(id)%loc%tdata%z%ss3(i,j,k)
            
        end do
        end do
        end do
     
    enddo           
    !$omp end parallel do
    
end subroutine

subroutine ns_split_diff
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k,iter
real(8) :: err

!$omp parallel do private(i,j,k)
do id = 0, p%glb%threads-1
   
   do k = p%of(id)%loc%ks-p%glb%ghc, p%of(id)%loc%ke+p%glb%ghc 
   do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
   do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        p%of(id)%loc%vel%x%old2(i,j,k) = p%of(id)%loc%vel%x%now(i,j,k)
        p%of(id)%loc%vel%y%old2(i,j,k) = p%of(id)%loc%vel%y%now(i,j,k)
        p%of(id)%loc%vel%z%old2(i,j,k) = p%of(id)%loc%vel%z%now(i,j,k)
   enddo
   enddo
   enddo
        
enddo    
!$omp end parallel do

!==============================================================================

do iter = 1, 5

call ns_linearize

!$omp parallel do private(i,j,k)
do id = 0, p%glb%threads-1
    
   do k = p%of(id)%loc%ks-p%glb%ghc, p%of(id)%loc%ke+p%glb%ghc 
   do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
   do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        p%of(id)%loc%vel%x%now(i,j,k) = 0.5d0*( p%of(id)%loc%vel%x%now(i,j,k) + p%of(id)%loc%vel%x%old2(i,j,k) ) 
        p%of(id)%loc%vel%y%now(i,j,k) = 0.5d0*( p%of(id)%loc%vel%y%now(i,j,k) + p%of(id)%loc%vel%y%old2(i,j,k) )
        p%of(id)%loc%vel%z%now(i,j,k) = 0.5d0*( p%of(id)%loc%vel%z%now(i,j,k) + p%of(id)%loc%vel%z%old2(i,j,k) )
   enddo
   enddo
   enddo

    call p%of(id)%find_stag_vel( p%of(id)%loc%velsrc%x%now, p%of(id)%loc%velsrc%y%now, p%of(id)%loc%velsrc%z%now, &
                                &p%of(id)%loc%velsrc%x%old, p%of(id)%loc%velsrc%y%old, p%of(id)%loc%velsrc%z%old, &
                                &p%of(id)%loc%vel%x%now   , p%of(id)%loc%vel%y%now   , p%of(id)%loc%vel%z%now )

enddo
!$omp end parallel do

call pt%velsrc%sync
call pt%velsrc_old%sync

!$omp parallel do
do id = 0, p%glb%threads-1
    
    call ns_split_diff_source( p%of(id), &
                              &p%of(id)%loc%velsrc%x%now, p%of(id)%loc%velsrc%y%now, p%of(id)%loc%velsrc%z%now, &
                              &p%of(id)%loc%velsrc%x%old, p%of(id)%loc%velsrc%y%old, p%of(id)%loc%velsrc%z%old, &
                              &p%of(id)%loc%vel%x%now   , p%of(id)%loc%vel%y%now   , p%of(id)%loc%vel%z%now   , &
                              &p%of(id)%loc%velsrc%x%tmp, p%of(id)%loc%velsrc%y%tmp, p%of(id)%loc%velsrc%z%tmp, &
                              &p%of(id)%loc%rho%old,p%of(id)%loc%mu%old,p%of(id)%loc%delta%old,p%of(id)%loc%normals%curv%old)
    
enddo    
!$omp end parallel do


!$omp parallel do private(i,j,k)
do id = 0, p%glb%threads-1
    
    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
        p%of(id)%loc%vel%x%now(i,j,k) = p%of(id)%loc%vel%x%old2(i,j,k) + p%of(id)%loc%velsrc%x%tmp(i,j,k) * p%glb%dt                      
        p%of(id)%loc%vel%y%now(i,j,k) = p%of(id)%loc%vel%y%old2(i,j,k) + p%of(id)%loc%velsrc%y%tmp(i,j,k) * p%glb%dt                      
        p%of(id)%loc%vel%z%now(i,j,k) = p%of(id)%loc%vel%z%old2(i,j,k) + p%of(id)%loc%velsrc%z%tmp(i,j,k) * p%glb%dt

    end do
    end do
    end do
        
    call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now,p%of(id)%loc%vel%z%now)
    
enddo
!$omp end parallel do
    
call pt%vel%sync


end do

end subroutine

subroutine ns_split_diff_source(q,u,v,w,uu,vv,ww,us,vs,ws,sx,sy,sz,irho,imu,idelta,icurv)
use all
implicit none
type(job) :: q
real(8), dimension(q%loc%is-q%glb%ghc:q%loc%ie+q%glb%ghc,&
                  &q%loc%js-q%glb%ghc:q%loc%je+q%glb%ghc,&
                  &q%loc%ks-q%glb%ghc:q%loc%ke+q%glb%ghc) :: u,v,w,uu,vv,ww,us,vs,ws,sx,sy,sz,irho,imu,idelta,icurv
real(8) :: rho,mu,delta,curv,xx,yy,zz
real(8) :: ux,uy,uz,vx,vy,vz,wx,wy,wz,phix,phiy,phiz
real(8) :: dif_x,dif_y,dif_z   
real(8) :: cnew,cold          
integer :: i,j,k

! ---------------------------------------------- Calculate derivatives

    do k = q%loc%ks, q%loc%ke
    do j = q%loc%js, q%loc%je
        call q%loc%ccd%x%SOLVE_SRKCCD(us(:,j,k),us(:,j,k),q%loc%tdata%x%s1(:,j,k),q%loc%tdata%x%ss1(:,j,k))
        call q%loc%ccd%x%SOLVE_SRKCCD(vs(:,j,k),vs(:,j,k),q%loc%tdata%y%s1(:,j,k),q%loc%tdata%y%ss1(:,j,k))
        call q%loc%ccd%x%SOLVE_SRKCCD(ws(:,j,k),ws(:,j,k),q%loc%tdata%z%s1(:,j,k),q%loc%tdata%z%ss1(:,j,k))
    end do
    end do

    do k = q%loc%ks, q%loc%ke
    do i = q%loc%is, q%loc%ie
        call q%loc%ccd%y%SOLVE_SRKCCD(us(i,:,k),us(i,:,k),q%loc%tdata%x%s2(i,:,k),q%loc%tdata%x%ss2(i,:,k))
        call q%loc%ccd%y%SOLVE_SRKCCD(vs(i,:,k),vs(i,:,k),q%loc%tdata%y%s2(i,:,k),q%loc%tdata%y%ss2(i,:,k))
        call q%loc%ccd%y%SOLVE_SRKCCD(ws(i,:,k),ws(i,:,k),q%loc%tdata%z%s2(i,:,k),q%loc%tdata%z%ss2(i,:,k))
    end do
    end do

    do j = q%loc%js, q%loc%je
    do i = q%loc%is, q%loc%ie
        call q%loc%ccd%z%SOLVE_SRKCCD(us(i,:,k),us(i,j,:),q%loc%tdata%x%s3(i,j,:),q%loc%tdata%x%ss3(i,j,:))
        call q%loc%ccd%z%SOLVE_SRKCCD(vs(i,:,k),vs(i,j,:),q%loc%tdata%y%s3(i,j,:),q%loc%tdata%y%ss3(i,j,:))
        call q%loc%ccd%z%SOLVE_SRKCCD(ws(i,:,k),ws(i,j,:),q%loc%tdata%z%s3(i,j,:),q%loc%tdata%z%ss3(i,j,:))
    end do
    end do
! ----------------------------------------------
    
    do k = q%loc%ks, q%loc%ke
    do j = q%loc%js, q%loc%je
    do i = q%loc%is, q%loc%ie
        
        xx = q%loc%tdata%x%ss1(i,j,k); yy = q%loc%tdata%x%ss2(i,j,k); zz = q%loc%tdata%x%ss3(i,j,k)
        
        cnew = 0.5d0*( q%loc%c%now(i,j,k) + q%loc%c%now(i+1,j,k) )
        cold = 0.5d0*( q%loc%c%old(i,j,k) + q%loc%c%old(i+1,j,k) )
        
        dif_x = xx/p%glb%re 
        dif_y = yy/p%glb%re 
        dif_z = zz/p%glb%re 
            
        sx(i,j,k) = dif_x + dif_y + dif_z + 0.5d0*(cnew+cold)*q%glb%gx
        
    end do
    end do
    end do

    do k = q%loc%ks, q%loc%ke
    do j = q%loc%js, q%loc%je
    do i = q%loc%is, q%loc%ie

        xx = q%loc%tdata%y%ss1(i,j,k); yy = q%loc%tdata%y%ss2(i,j,k); zz = q%loc%tdata%y%ss3(i,j,k)
            
        cnew = 0.5d0*( q%loc%phi%now(i,j,k) + q%loc%phi%now(i,j+1,k) )
        cold = 0.5d0*( q%loc%phi%old(i,j,k) + q%loc%phi%old(i,j+1,k) )
        
        dif_x = xx/p%glb%re 
        dif_y = yy/p%glb%re 
        dif_z = zz/p%glb%re 
            
        sy(i,j,k) = dif_x + dif_y + dif_z + 0.5d0*(cnew+cold)*q%glb%gy 
        
    end do
    end do
    end do

    do k = q%loc%ks, q%loc%ke
    do j = q%loc%js, q%loc%je
    do i = q%loc%is, q%loc%ie
        
        xx = q%loc%tdata%z%ss1(i,j,k); yy = q%loc%tdata%z%ss2(i,j,k); zz = q%loc%tdata%z%ss3(i,j,k)
           
        cnew = 0.5d0*( q%loc%c%now(i,j,k) + q%loc%c%now(i,j,k+1) )
        cold = 0.5d0*( q%loc%c%old(i,j,k) + q%loc%c%old(i,j,k+1) )
        
        dif_x = xx/p%glb%re 
        dif_y = yy/p%glb%re 
        dif_z = zz/p%glb%re 
            
        sz(i,j,k) = dif_x + dif_y + dif_z + 0.5d0*(cnew+cold)*q%glb%gz 
        
    end do
    end do
    end do
    
end subroutine


