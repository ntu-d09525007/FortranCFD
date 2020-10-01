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
real(8) :: total

 call problem_init

 do 
	
 	p%glb%time = p%glb%time + p%glb%dt 
	p%glb%iter = p%glb%iter + 1 
	call p%sync 
      
	!call deforming_field_velocity
    call p%switch
     
	!call level_set_srk6_solver 
	call level_set_rk3_solver
	call mass_preserving_level_set
	if(mod(p%glb%iter,3).eq.0)call level_set_rk3_redis(1)
	call ns_solver
	call plt%plot
	
	call p%ls_mv
	if( mod(p%glb%iter,5)==0 )then
		total = p%glb%ls_adv + p%glb%ls_red + p%glb%ns
		write(*,*)"========================================="
		write(*,'(A,",",I8,",",F15.6)')trim(p%glb%name),p%glb%iter,p%glb%time
		write(*,*)''
		write(*,'("Divergence :",ES15.4)')p%glb%vel_div
		write(*,'("L2 norm    :",ES15.4)')p%glb%ns_l2f
		write(*,'("Linf norm  :",ES15.4)')p%glb%ns_linf
		!write(*,*)''
		!write(*,'("PPE iters  :",I15)')p%glb%piter
		!write(*,'("PPE error  :",ES15.4)')p%glb%ppe_linf		
		write(*,*)''
		write(*,'("Loss of mass  (%) :",ES15.4)')(p%glb%mass-p%glb%imass)/p%glb%imass
		write(*,'("Loss of volume(%) :",ES15.4)')(p%glb%vol-p%glb%ivol)/p%glb%ivol
		write(*,*)''
		write(*,'("Total CPU time(s) :",F15.6)')total
		write(*,'(4A13)')"LS Adv.","LS Redis.","PPE","NS"
		write(*,'(F12.2,"%",F12.2,"%",F12.2,"%",F12.2,"%")')100.0d0*p%glb%ls_adv/total,100.0d0*p%glb%ls_red/total&
														  &,100.0d0*p%glb%ppe/total,100.0d0*(p%glb%ns-p%glb%ppe)/total
		write(*,*)"========================================="
	endif
	call output
	
	if( p%glb%time > p%glb%t2s ) exit
	
 end do
 
 contains
 
 include 'included.f90'

end program main
