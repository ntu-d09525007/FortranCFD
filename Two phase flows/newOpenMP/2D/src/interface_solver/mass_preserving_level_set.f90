subroutine mass_preserving_level_set
use all
!$ use omp_lib
implicit none
integer :: i,j,iter
real(8) :: lam, plam

do iter = 1, 2
    
    call ls_mv()
    call surface_norms()
    
    lam = 0.0_8
    
    !$omp parallel do collapse(2), private(plam), reduction(+:lam)   
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
    
        p%loc%grad%now(i,j) = dsqrt( p%loc%normals%x%now(i,j)**2.0_8 + &
                                   & p%loc%normals%y%now(i,j)**2.0_8)
        
        plam =  p%loc%delta%now(i,j)**2.0_8 * p%loc%grad%now(i,j) 
        plam = plam * ( 2.0_8*(1.0_8-p%glb%rho_12)*p%loc%heavy%now(i,j) + p%glb%rho_12 )*p%glb%dx*p%glb%dy
        
        lam = lam + plam
        
    end do
    end do 
    !$omp end parallel do
    
    lam = ( p%glb%imass - p%glb%mass ) / lam 

    !$omp parallel do collapse(2)  
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
    
        p%loc%phi%now(i,j) = p%loc%phi%now(i,j) + lam * p%loc%delta%now(i,j) * p%loc%grad%now(i,j)
        
    end do
    end do
    !$omp end parallel do

    call bc(p%loc%phi%now)
       
end do  

end subroutine

