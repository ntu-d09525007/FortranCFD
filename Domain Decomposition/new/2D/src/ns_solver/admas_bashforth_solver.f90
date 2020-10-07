subroutine ns_ab_setup
use all
!$ use omp_lib
implicit none

    call ns_ab_adv_source

end subroutine

subroutine ns_ab_solver_mg
use all
implicit none
integer :: iter,initer,relax_iter

relax_iter = 3

call ns_ab_adv_source
    
iter=0
p%glb%piter=0
    
do 
    
    iter=iter+1 
    
    call ns_linearize
    call ns_ab_diff_source
    call ns_ab_predictor   

    if(iter>5)exit
        
end do

call ns_check_convergence_vel

if( p%glb%iter<relax_iter )then
    call ppe_sor_solver(1.0d-8)
else
    call ppe_mg_solver(iter)
endif
    
end subroutine

subroutine ns_ab_solver_SOR
use all
implicit none
integer :: iter,initer
real(8) :: tol

call ns_ab_adv_source
    
iter=0
p%glb%piter=0
    
do 
    
    iter=iter+1 
    
    call ns_linearize
    call ns_ab_diff_source
    call ns_ab_predictor

    if(iter>5)exit
        
end do

call ns_check_convergence_vel
call ppe_sor_solver(p%glb%p_tol)

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
        
    call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now)
    
enddo        
!$omp end parallel do
    
call pt%vel%sync
    
end subroutine