module all
use ptr_tree
use tree
use vtk_plotter
implicit none
save
type(manager) :: p 
type(ptr_family) :: pt
type(plotter) :: plt
end module all

program main
use all
implicit none
real(8) :: dt

  call problem_init
  dt=p%glb%dt
  
  call ppe_mg_solver

  ! do 
        
        ! !call check_memory
        
        ! p%glb%time = p%glb%time + p%glb%dt 
        ! p%glb%iter = p%glb%iter + 1 
        ! call p%sync 
		
		! ! if(p%glb%iter<50)then
			! ! p%glb%dt=dt*0.01
		! ! else
			! ! p%glb%dt=dt
		! ! endif
      
        ! !call deforming_field_velocity
        ! call p%switch
     
        ! call interface_solver
        ! call ns_solver
        ! call p%plots
	
		! call p%ls_mv
		! if( mod(p%glb%iter,5)==0 )then
			! write(*,*)"========================================="
			! write(*,'(A,",",I8,",",F15.6)')trim(p%glb%name),p%glb%iter,p%glb%time
			! write(*,*)''
			! call print_NS_info
			! call print_LS_info
			! call print_CPU_info
			! call find_momentum
			! write(*,*)"========================================="
		! endif
		! call output
	
	! if( p%glb%time > p%glb%t2s ) exit
	
 ! end do
 
 contains
 
 include 'included.f90'

end program main
