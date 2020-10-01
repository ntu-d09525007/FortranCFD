subroutine rk3_c_solver()
use all
!$ use omp_lib
implicit none
integer :: i,j,k,id
integer(8) :: cpustart, cpuend
real(8) :: src

    call system_clock(cpustart)
    
    call settling_velocity_c

    !$omp parallel do private(i,j,k,src)
    do id = 0, p%glb%threads-1
        
        call c_rk3_source(p%of(id),p%of(id)%loc%tdata%x%l1)
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            src = p%of(id)%loc%tdata%x%l1(i,j,k)
            p%of(id)%loc%c%now(i,j,k) = p%of(id)%loc%c%now(i,j,k) + src * p%glb%dt          
        end do 
        end do 
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%c%now)
    
    enddo
    !$omp end parallel do

    call flux_blance_c

    !$omp parallel do private(i,j,k,src)
    do id = 0, p%glb%threads-1
        
        call c_rk3_source(p%of(id),p%of(id)%loc%tdata%x%l2)
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            src = ( -3.0d0*p%of(id)%loc%tdata%x%l1(i,j,k) + p%of(id)%loc%tdata%x%l2(i,j,k) ) / 4.0d0
            p%of(id)%loc%c%now(i,j,k) = p%of(id)%loc%c%now(i,j,k) + src * p%glb%dt
        end do 
        end do 
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%c%now)
    
    enddo
    !$omp end parallel do

    call flux_blance_c
    
    !$omp parallel do private(i,j,k,src)
    do id = 0, p%glb%threads-1
        
        call c_rk3_source(p%of(id),p%of(id)%loc%tdata%x%l3)
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            src = ( -p%of(id)%loc%tdata%x%l1(i,j,k)-p%of(id)%loc%tdata%x%l2(i,j,k)+8.0d0*p%of(id)%loc%tdata%x%l3(i,j,k) ) / 12.0d0
            p%of(id)%loc%c%now(i,j,k) = p%of(id)%loc%c%now(i,j,k) + src * p%glb%dt
        end do 
        end do 
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%c%now)
    
    enddo
    !$omp end parallel do

    call flux_blance_c
    
    call system_clock(cpuend)
    p%glb%ls_adv = p%glb%ls_adv + real(cpuend-cpustart,kind=8) / real( p%glb%cpurate, kind=8 )

end subroutine

subroutine flux_blance_c()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: src

!$omp parallel do private(i,j,k,src)
do id = 0, p%glb%threads-1

    ! do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    ! do i = p%of(id)%loc%is, p%of(id)%loc%ie
    ! do j = p%of(id)%loc%js, p%of(id)%loc%je
    !     if( p%of(id)%loc%solid%now(i,j,k)>0.5d0 .and. p%of(id)%loc%solid%now(i,j,k+1)<0.5d0 ) then
    !         p%of(id)%loc%c%now(i,j,k) = ( 9.0d0*p%of(id)%loc%c%now(i,j,k+1)-p%of(id)%loc%c%now(i,j,k+2) ) / 8.0d0
    !     endif
    ! enddo
    ! enddo
    ! enddo
    
    ! if (p%of(id)%loc%idz==p%glb%grid_z-1)then
        
    !     do i = p%of(id)%loc%is, p%of(id)%loc%ie
    !     do j = p%of(id)%loc%js, p%of(id)%loc%je
    !         p%of(id)%loc%c%now(i,j,p%glb%node_z) = ( 9.0d0*p%of(id)%loc%c%now(i,j,p%glb%node_z-1)-p%of(id)%loc%c%now(i,j,p%glb%node_z-2) ) / 8.0d0
    !     enddo
    !     enddo
            
    ! endif
        
    if (p%of(id)%loc%idz==p%glb%grid_z-1)then
        
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        do j = p%of(id)%loc%js, p%of(id)%loc%je
            p%of(id)%loc%c%now(i,j,p%glb%node_z) = p%of(id)%loc%c%now(i,j,p%glb%node_z-1) / ( -p%glb%Re*p%glb%Dz*p%glb%us_c + 1.0d0 )
        enddo
        enddo
            
    endif
    
    call p%of(id)%bc(0,p%of(id)%loc%c%now)

enddo    
!$omp end parallel do

call pt%c%sync

!$omp parallel do private(i,j,k)
do id = 0, p%glb%threads-1

    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        p%of(id)%loc%c%now(i,j,k) = max(min(1.0d0,p%of(id)%loc%c%now(i,j,k)),0.0d0)
    enddo
    enddo
    enddo

call p%of(id)%bc(0,p%of(id)%loc%c%now)

enddo  
!$omp end parallel do

call pt%c%sync

end subroutine

subroutine settling_velocity_c()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k

    !$omp parallel do private(i,j,k)
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks-p%glb%ghc, p%of(id)%loc%ke+p%glb%ghc
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc

            p%of(id)%loc%nvel%z%tmp(i,j,k) = p%of(id)%loc%nvel%z%old(i,j,k) - p%glb%us_c * (1.0d0-p%of(id)%loc%solid%now(i,j,k))
            
        end do
        end do
        end do 

    enddo
    !$omp end parallel do

end subroutine 

subroutine c_rk3_source(q,s)
use all
implicit none
type(job) :: q
real(8), dimension(q%loc%is-q%glb%ghc:q%loc%ie+q%glb%ghc,&
                  &q%loc%js-q%glb%ghc:q%loc%je+q%glb%ghc,&
                  &q%loc%ks-q%glb%ghc:q%loc%ke+q%glb%ghc) :: s
integer :: i,j,k

                                  
    do k = q%loc%ks, q%loc%ke
    do j = q%loc%js, q%loc%je
        call q%loc%uccd%x%solve(q%loc%nvel%x%old(:,j,k),q%loc%c%now(:,j,k),q%loc%tdata%x%s1(:,j,k),q%loc%tdata%x%s2(:,j,k))
    end do 
    end do

    do k = q%loc%ks, q%loc%ke
    do i = q%loc%is, q%loc%ie
        call q%loc%uccd%y%solve(q%loc%nvel%y%old(i,:,k),q%loc%c%now(i,:,k),q%loc%tdata%y%s1(i,:,k),q%loc%tdata%y%s2(i,:,k))
    end do 
    end do

    do j = q%loc%js, q%loc%je
    do i = q%loc%is, q%loc%ie
        call q%loc%uccd%z%solve(q%loc%nvel%z%tmp(i,j,:),q%loc%c%now(i,j,:),q%loc%tdata%z%s1(i,j,:),q%loc%tdata%z%s2(i,j,:))
    end do 
    end do
    
    do k = q%loc%ks, q%loc%ke
    do j = q%loc%js, q%loc%je
    do i = q%loc%is, q%loc%ie
        s(i,j,k) = - q%loc%nvel%x%old(i,j,k)*q%loc%tdata%x%s1(i,j,k) + 1.0d0/q%glb%re * q%loc%tdata%x%s2(i,j,k) &
                &  - q%loc%nvel%y%old(i,j,k)*q%loc%tdata%y%s1(i,j,k) + 1.0d0/q%glb%re * q%loc%tdata%y%s2(i,j,k) &
                &  - q%loc%nvel%z%tmp(i,j,k)*q%loc%tdata%z%s1(i,j,k) + 1.0d0/q%glb%re * q%loc%tdata%z%s2(i,j,k)
    end do 
    end do
    end do
    
end subroutine
