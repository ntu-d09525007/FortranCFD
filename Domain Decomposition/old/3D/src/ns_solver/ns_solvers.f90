subroutine ns_solver
use all
implicit none
integer(8) :: cpustart, cpuend

    call system_clock(cpustart)
    
    call ns_init
    
    !call ns_ab_solver_SOR
    !call ns_split_solver
    call ns_ab_solver_mg
    
    call ns_check_convergence
    call p%node_vel
    call pt%nvel%sync
    
    call system_clock(cpuend)
    p%glb%ns = p%glb%ns + real(cpuend-cpustart,kind=8)/real(p%glb%cpurate,kind=8)
    
end subroutine

subroutine ns_init
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k

    call p%rho_mu
    
    call p%surface_norms2
    call pt%normals_1%sync
    call pt%normals_2%sync
    call pt%normals_3%sync
    
    call p%curv
    call pt%curv%sync

end subroutine

subroutine ns_linearize
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k

    !$omp parallel private(id,i,j,k), num_threads(p%glb%threads)
    
        id=0
        !$ id = omp_get_thread_num()
        
        do k = p%of(id)%loc%ks-p%glb%ghc, p%of(id)%loc%ke+p%glb%ghc
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
            p%of(id)%loc%vel%x%tmp(i,j,k) = p%of(id)%loc%vel%x%now(i,j,k) 
            p%of(id)%loc%vel%y%tmp(i,j,k) = p%of(id)%loc%vel%y%now(i,j,k) 
            p%of(id)%loc%vel%z%tmp(i,j,k) = p%of(id)%loc%vel%z%now(i,j,k) 
        end do
        end do
        end do
        
    !$omp end parallel 
        
end subroutine

subroutine ns_check_convergence
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: linf, l2f, div, sumdiv
real(8) :: ux,vy,wz

    linf=0.0d0
    l2f=0.0d0
    div=0.0d0
    sumdiv=0.0d0
    
    !$omp parallel private(id,i,j,k,ux,vy,wz), num_threads(p%glb%threads), reduction(max:linf,div), reduction(+:l2f,sumdiv)
    
        id=0
        !$ id = omp_get_thread_num()
        
        linf=0.0d0
        l2f=0.0d0
        div=0.0d0
        sumdiv=0.0d0
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            linf = max(linf,abs(p%of(id)%loc%vel%x%now(i,j,k)-p%of(id)%loc%vel%x%tmp(i,j,k)))
            linf = max(linf,abs(p%of(id)%loc%vel%y%now(i,j,k)-p%of(id)%loc%vel%y%tmp(i,j,k)))
            linf = max(linf,abs(p%of(id)%loc%vel%z%now(i,j,k)-p%of(id)%loc%vel%z%tmp(i,j,k)))
            
            l2f = l2f + (p%of(id)%loc%vel%x%now(i,j,k)-p%of(id)%loc%vel%x%tmp(i,j,k))**2.0d0
            l2f = l2f + (p%of(id)%loc%vel%y%now(i,j,k)-p%of(id)%loc%vel%y%tmp(i,j,k))**2.0d0
            l2f = l2f + (p%of(id)%loc%vel%z%now(i,j,k)-p%of(id)%loc%vel%z%tmp(i,j,k))**2.0d0
            
            ux = (p%of(id)%loc%vel%x%now(i,j,k)-p%of(id)%loc%vel%x%now(i-1,j,k))/p%glb%dx
            vy = (p%of(id)%loc%vel%y%now(i,j,k)-p%of(id)%loc%vel%y%now(i,j-1,k))/p%glb%dy
            wz = (p%of(id)%loc%vel%z%now(i,j,k)-p%of(id)%loc%vel%z%now(i,j,k-1))/p%glb%dz 
                
            div = max( div, abs(ux+vy+wz) ) 
            sumdiv = sumdiv + abs(ux+vy+wz)
                                            
        end do
        end do
        end do
        
    !$omp end parallel
    
    l2f = dsqrt( l2f / (3.0d0*p%glb%node_x*p%glb%node_y*p%glb%node_z) )
    
    p%glb%vel_sdiv = sumdiv / (p%glb%node_x*p%glb%node_y*p%glb%node_z)
    p%glb%vel_div = div
    p%glb%ns_linf = linf
    p%glb%ns_l2f = l2f
    
end subroutine

