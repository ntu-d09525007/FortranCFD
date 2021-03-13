subroutine level_set_symplectic_solver()
use all
!$ use omp_lib
implicit none
integer :: i,j,id,iter
integer(8) :: cpustart, cpuend
real(8) :: err,perr

    call system_clock(cpustart)

    !$omp parallel do
    do id = 0, p%glb%threads-1  
        
        call p%of(id)%loc%tdata%init(.false.,p%of(id)%loc%phi%old)
        
    enddo
    !$omp end parallel do
    
    iter=0
    
    do 
    
        iter = iter + 1

        call level_set_source()
        
        err=0.0_8

        !$omp parallel do reduction(max:err)
        do id = 0, p%glb%threads-1
    
            call p%of(id)%loc%tdata%solve_srk6(err)
            !call p%of(id)%loc%tdata%solve_srk4(err)
            
        enddo
        !$omp end parallel do
        
    
        if( err < p%glb%t_tol )exit
        
        if( mod(iter,500) == 0)write(*,*)"LS solver,",iter,err
        
    end do
    
    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1
        
        call p%of(id)%loc%tdata%final_srk6()
        !call p%of(id)%loc%tdata%final_srk4()
        
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
            p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%tdata%x%target(i,j)
        end do
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%phi%now)
    
    enddo
    !$omp end parallel do

    call pt%phi%sync
    
    call system_clock(cpuend)
    p%glb%ls_adv = p%glb%ls_adv + real(cpuend-cpustart,kind=8)/real(p%glb%cpurate,kind=8)

end subroutine


subroutine level_set_source()
use all 
!$ use omp_lib
implicit none
integer :: i,j,id

    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1
        
        call p%of(id)%bc(0,p%of(id)%loc%tdata%x%s1)
        call p%of(id)%bc(0,p%of(id)%loc%tdata%x%s2)
        call p%of(id)%bc(0,p%of(id)%loc%tdata%x%s3)
    
    enddo
    !$omp end parallel do
    
    call pt%tdatax%sync
    
    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1
        
        ! calculate x derivatives
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        
            call p%of(id)%loc%ccdsolvers%x%solve("uccd",&
                &p%of(id)%loc%tdata%x%s1(:,j),p%of(id)%loc%tdata%x%ss1(:,j),p%of(id)%loc%nvel%x%tmp(:,j),p%of(id)%loc%nvel%x%old(:,j))
                                        
            call p%of(id)%loc%ccdsolvers%x%solve("uccd",&
                &p%of(id)%loc%tdata%x%s2(:,j),p%of(id)%loc%tdata%x%ss2(:,j),p%of(id)%loc%nvel%x%tmp(:,j),p%of(id)%loc%nvel%x%old(:,j))
                                        
            call p%of(id)%loc%ccdsolvers%x%solve("uccd",&
               &p%of(id)%loc%tdata%x%s3(:,j),p%of(id)%loc%tdata%x%ss3(:,j),p%of(id)%loc%nvel%x%tmp(:,j),p%of(id)%loc%nvel%x%old(:,j))
            
        end do
        
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
            
            p%of(id)%loc%tdata%x%l1(i,j) = - p%of(id)%loc%tdata%x%ss1(i,j) * p%of(id)%loc%nvel%x%old(i,j)
            p%of(id)%loc%tdata%x%l2(i,j) = - p%of(id)%loc%tdata%x%ss2(i,j) * p%of(id)%loc%nvel%x%old(i,j)
            p%of(id)%loc%tdata%x%l3(i,j) = - p%of(id)%loc%tdata%x%ss3(i,j) * p%of(id)%loc%nvel%x%old(i,j)
            
        end do
        end do
        
        ! calculate y derivatives
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
            
            call p%of(id)%loc%ccdsolvers%y%solve("uccd",&
                &p%of(id)%loc%tdata%x%s1(i,:),p%of(id)%loc%tdata%x%ss1(i,:),p%of(id)%loc%nvel%y%tmp(i,:),p%of(id)%loc%nvel%y%old(i,:))
                                        
            call p%of(id)%loc%ccdsolvers%y%solve("uccd",&
                &p%of(id)%loc%tdata%x%s2(i,:),p%of(id)%loc%tdata%x%ss2(i,:),p%of(id)%loc%nvel%y%tmp(i,:),p%of(id)%loc%nvel%y%old(i,:))
                                        
            call p%of(id)%loc%ccdsolvers%y%solve("uccd",&
               &p%of(id)%loc%tdata%x%s3(i,:),p%of(id)%loc%tdata%x%ss3(i,:),p%of(id)%loc%nvel%y%tmp(i,:),p%of(id)%loc%nvel%y%old(i,:))
        end do
        
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        
            p%of(id)%loc%tdata%x%l1(i,j) = p%of(id)%loc%tdata%x%l1(i,j) - p%of(id)%loc%tdata%x%ss1(i,j) * p%of(id)%loc%nvel%y%old(i,j)
            p%of(id)%loc%tdata%x%l2(i,j) = p%of(id)%loc%tdata%x%l2(i,j) - p%of(id)%loc%tdata%x%ss2(i,j) * p%of(id)%loc%nvel%y%old(i,j)
            p%of(id)%loc%tdata%x%l3(i,j) = p%of(id)%loc%tdata%x%l3(i,j) - p%of(id)%loc%tdata%x%ss3(i,j) * p%of(id)%loc%nvel%y%old(i,j)
            
        end do
        end do
    
    enddo
    !$omp end parallel do
    
end subroutine
