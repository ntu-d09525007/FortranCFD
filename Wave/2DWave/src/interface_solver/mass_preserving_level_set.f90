subroutine mass_preserving_level_set
use all
!$ use omp_lib
implicit none
integer :: i,j,id,iter
real(8) :: lam, plam

do iter = 1, 3
    
    call p%ls_mv()
    call p%surface_norms2()
    call pt%normals%sync
    
    lam = 0.0_8
    
    !$omp parallel do private(i,j,plam), reduction(+:lam)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            p%of(id)%loc%grad%now(i,j) = dsqrt( p%of(id)%loc%normals%x%now(i,j)**2.0_8 + &
                                                & p%of(id)%loc%normals%y%now(i,j)**2.0_8 )
            
            plam =  p%of(id)%loc%delta%now(i,j)**2.0_8 * p%of(id)%loc%grad%now(i,j) 
            plam = plam * ( 2.0_8*(1.0_8-p%glb%rho_12)*p%of(id)%loc%heavy%now(i,j) + p%glb%rho_12 ) * p%glb%dx * p%glb%dy
            
            lam = lam + plam
            
        end do
        end do
    
    enddo   
    !$omp end parallel do
    
    lam = ( p%glb%imass - p%glb%mass ) / lam 

    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1
            
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%phi%now(i,j) + lam * p%of(id)%loc%delta%now(i,j) * p%of(id)%loc%grad%now(i,j)
            
        end do
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%phi%now)
    
    enddo   
    !$omp end parallel do
    
    call pt%phi%sync
    
end do  

end subroutine

