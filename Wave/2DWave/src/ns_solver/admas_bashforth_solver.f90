subroutine ns_ab_setup
implicit none

    call ns_ab_adv_source

end subroutine

subroutine ns_ab_solver
use all
implicit none
integer :: iter

call ns_ab_adv_source
    
p%glb%piter=0

do iter = 1, 5
    
    call ns_linearize
    call ns_ab_diff_source
    call ns_ab_predictor   
    
end do

call ns_check_convergence_vel
call ppe_mg_solver(p%glb%piter)
!call ppe_sor_solver(1.0d-6)
    
end subroutine

subroutine ns_ab_predictor
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: src

!$omp parallel do private(i,j,src)
do id = 0, p%glb%threads-1
        
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
            
        src = 1.50*p%of(id)%loc%velsrc%x%now(i,j) - 0.5d0*p%of(id)%loc%velsrc%x%old(i,j) 
        src = src + p%of(id)%loc%velsrc%x%tmp(i,j)            
        p%of(id)%loc%vel%x%now(i,j) = p%of(id)%loc%vel%x%old(i,j) + p%glb%dt * src
        
        src = 1.50*p%of(id)%loc%velsrc%y%now(i,j) - 0.5d0*p%of(id)%loc%velsrc%y%old(i,j) 
        src = src + p%of(id)%loc%velsrc%y%tmp(i,j)    
        p%of(id)%loc%vel%y%now(i,j) = p%of(id)%loc%vel%y%old(i,j) + p%glb%dt * src

    end do 
    end do
    
enddo        
!$omp end parallel do

call ns_velbc
    
end subroutine