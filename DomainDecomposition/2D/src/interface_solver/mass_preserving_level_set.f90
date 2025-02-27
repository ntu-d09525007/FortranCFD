subroutine mass_preserving_level_set
use all
!$ use omp_lib
implicit none
integer :: i,j,id,iter
real(8) :: lam, plam, src

call p%curv()
call pt%normals%sync

do iter = 1, 20
    
    call p%ls_mv()
    
    lam = 0.0_8
    
    !$omp parallel do private(i,j,plam), reduction(+:lam)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            p%of(id)%loc%grad%now(i,j) = dsqrt( p%of(id)%loc%normals%x%now(i,j)**2.0_8 + &
                                              & p%of(id)%loc%normals%y%now(i,j)**2.0_8 )
            
            plam =  p%of(id)%loc%delta%now(i,j)**2.0_8 * p%of(id)%loc%grad%now(i,j) 
            plam = plam * ( 2.0_8*(1.0_8-p%glb%rho_12)*p%of(id)%loc%heavy%now(i,j) + p%glb%rho_12 ) * p%glb%dx * p%glb%dy
            if(p%glb%mpls==1)plam = plam * p%of(id)%loc%normals%curv%now(i,j)

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
        
            src = lam * p%of(id)%loc%delta%now(i,j) * p%of(id)%loc%grad%now(i,j) 
            if(p%glb%mpls==1)src = src * p%of(id)%loc%normals%curv%now(i,j)

            p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%phi%now(i,j) + src
            
        end do
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%phi%now)
    
    enddo   
    !$omp end parallel do
    
    call pt%phi%sync

    if( abs(p%glb%imass - p%glb%mass)/p%glb%imass < 1.0d-10 )exit
    
end do  

end subroutine

