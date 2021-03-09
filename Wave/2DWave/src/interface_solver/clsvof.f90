subroutine vof_wlic_solver()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,ii,jj,ib
real(8) :: eps,isgn,alpha,beta,w,a1,a3,a4,a5,xc,s
integer(8) :: cpustart, cpuend

    call system_clock(cpustart)

    eps = 1.0d-12

    call p%surface_norms2
    call pt%normals%sync

    ! x direction 
    !$omp parallel do private(i,j,ii,ib,isgn,alpha,beta,w,a1,a3,a4,a5,xc) 
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is-1, p%of(id)%loc%ie
            
            if( p%of(id)%loc%vel%x%old(i,j) > 0.0d0 )then
                ii = i 
                isgn = 1.0d0
            else
                ii = i+1
                isgn = 0.0d0
            endif
            
            beta = 2.3d0
            
            if( p%of(id)%loc%vof%now(ii,j) > 1.0d0-eps .or. p%of(id)%loc%vof%now(ii,j) < eps )then
               
               p%of(id)%loc%tdata%x%s1(i,j) = p%of(id)%loc%vof%now(ii,j) * p%of(id)%loc%vel%x%old(i,j) * p%glb%dt
              
            else
            
                ib = max(1,ii-1)
                
                if( p%of(id)%loc%vof%now(ib,j) < p%of(id)%loc%vof%now(ii+1,j) )then
                    alpha = 1.0d0
                else
                    alpha = -1.0d0
                endif
                
                a1 = dexp( beta*( 2.0d0 * p%of(id)%loc%vof%now(ii,j) - 1.0d0 ) / alpha )
                a3 = dexp( beta )
                xc = 0.5d0 / beta * dlog( (a3*a3-a1*a3)/(a1*a3-1.0d0) )
                a4 = dcosh( beta * ( isgn - p%of(id)%loc%vel%x%old(i,j)*p%glb%dt/p%glb%dx - xc ) )
                a5 = dcosh( beta * ( isgn - xc ) )
                
                w  = abs(p%of(id)%loc%normals%x%now(ii,j)) + abs(p%of(id)%loc%normals%y%now(ii,j))
                w  = abs(p%of(id)%loc%normals%x%now(ii,j)) / w
            
                p%of(id)%loc%tdata%x%s1(i,j) = 0.5d0*( p%of(id)%loc%vel%x%old(i,j)*p%glb%dt - alpha*p%glb%dx/beta*dlog(a4/a5) )
            
                p%of(id)%loc%tdata%x%s1(i,j) = p%of(id)%loc%tdata%x%s1(i,j)*w + ( 1.0d0 - w ) * p%of(id)%loc%vof%now(ii,j) * p%of(id)%loc%vel%x%old(i,j) * p%glb%dt
                                                
            endif
            
        end do
        end do

    enddo
    !$omp end parallel do

    call pt%tdatax%nodes(1)%sync

    !$omp parallel do private(i,j,s)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie

            s = ( p%of(id)%loc%tdata%x%s1(i,j) - p%of(id)%loc%tdata%x%s1(i-1,j) ) / p%glb%dx

            !if(i==1) s = p%of(id)%loc%tdata%x%s1(i,j) / p%glb%dx
            
            p%of(id)%loc%vof%now(i,j) = p%of(id)%loc%vof%now(i,j) - s + &
                                        & p%of(id)%loc%vof%old(i,j) * p%glb%dt * ( p%of(id)%loc%vel%x%old(i,j) - p%of(id)%loc%vel%x%old(i-1,j) ) / p%glb%dx
            
        end do 
        end do 
        
        call p%of(id)%bc(0,p%of(id)%loc%vof%now)
        
    end do
    !$omp end parallel do
    
    call pt%vof%sync

    ! y direction 
    !$omp parallel do private(i,j,jj,ib,isgn,alpha,beta,w,a1,a3,a4,a5,xc)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js-1, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            
            if( p%of(id)%loc%vel%y%old(i,j) > 0.0d0 )then
                jj = j 
                isgn = 1.0d0
            else
                jj = j+1
                isgn = 0.0d0
            endif
            
            beta = 2.3d0
            
            if( p%of(id)%loc%vof%now(i,jj) > 1.0d0-eps .or. p%of(id)%loc%vof%now(i,jj) < eps )then
               
               p%of(id)%loc%tdata%y%s1(i,j) = p%of(id)%loc%vof%now(i,jj) * p%of(id)%loc%vel%y%old(i,j) * p%glb%dt
              
            else
            
                ib = max(1,jj-1)
                
                if( p%of(id)%loc%vof%now(i,ib) < p%of(id)%loc%vof%now(i,jj+1) )then
                    alpha = 1.0d0
                else
                    alpha = -1.0d0
                endif
                
                a1 = dexp( beta*( 2.0d0 * p%of(id)%loc%vof%now(i,jj) - 1.0d0 ) / alpha )
                a3 = dexp( beta )
                xc = 0.5d0 / beta * dlog( ( a3*a3-a1*a3)/(a1*a3-1.0d0) )
                a4 = dcosh( beta * ( isgn - p%of(id)%loc%vel%y%old(i,j)*p%glb%dt/p%glb%dy - xc ) )
                a5 = dcosh( beta * ( isgn - xc ) )
                
                w  = abs(p%of(id)%loc%normals%x%now(i,jj)) + abs(p%of(id)%loc%normals%y%now(i,jj)) 
                w  = abs(p%of(id)%loc%normals%y%now(i,jj)) / w
            
                p%of(id)%loc%tdata%y%s1(i,j) = 0.5d0*( p%of(id)%loc%vel%y%old(i,j)*p%glb%dt - alpha*p%glb%dy/beta*dlog(a4/a5) )
            
                p%of(id)%loc%tdata%y%s1(i,j) = p%of(id)%loc%tdata%y%s1(i,j)*w + ( 1.0d0 - w ) * p%of(id)%loc%vof%now(i,jj) * p%of(id)%loc%vel%y%old(i,j) * p%glb%dt
                                                
            endif
            
        end do
        end do

    enddo
    !$omp end parallel do

    call pt%tdatay%nodes(1)%sync 

    !$omp parallel do private(i,j,s)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie

            s = ( p%of(id)%loc%tdata%y%s1(i,j) - p%of(id)%loc%tdata%y%s1(i,j-1) ) / p%glb%dy

            !if(j==1) s = p%of(id)%loc%tdata%y%s1(i,j) / p%glb%dy
            
            p%of(id)%loc%vof%now(i,j) = p%of(id)%loc%vof%now(i,j) - s + &
                                        & p%of(id)%loc%vof%old(i,j) * p%glb%dt * ( p%of(id)%loc%vel%y%old(i,j) - p%of(id)%loc%vel%y%old(i,j-1) ) / p%glb%dy
            
            if(p%of(id)%loc%vof%now(i,j)<eps)p%of(id)%loc%vof%now(i,j)=0.0d0
            if(p%of(id)%loc%vof%now(i,j)>1.0-eps)p%of(id)%loc%vof%now(i,j)=1.0d0

        end do 
        end do 
        
        call p%of(id)%bc(0,p%of(id)%loc%vof%now)
    
    enddo
    !$omp end parallel do 
    
    call pt%vof%sync
    
    call system_clock(cpuend)
    p%glb%ls_adv = p%glb%ls_adv + real(cpuend-cpustart,kind=8) / real( p%glb%cpurate, kind=8 )
    
end subroutine

subroutine clsvof_recon()
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: time
integer(8) :: cpustart, cpuend

    call system_clock(cpustart)
    
    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
            p%of(id)%loc%vof%tmp(i,j) = p%of(id)%loc%phi%now(i,j)   
            p%of(id)%loc%phi%now(i,j) = 2.0d0*p%of(id)%loc%vof%now(i,j) - 1.0d0
        end do
        end do
    enddo
    !$omp end parallel do

    call level_set_rk3_redis(0,3.0d0*p%glb%ls_wid)

    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1

        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
            if( abs(p%of(id)%loc%vof%tmp(i,j)) > 1.5*p%glb%ls_wid )p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%vof%tmp(i,j)
        end do
        end do

        call p%of(id)%bc(0,p%of(id)%loc%phi%now)
        
    end do
    !$omp end parallel do

    call pt%phi%sync

    call system_clock(cpuend)
    p%glb%ls_red = p%glb%ls_red + real(cpuend-cpustart,kind=8)/real(p%glb%cpurate,kind=8)   
    
    call level_set_rk3_redis(1)

end subroutine





