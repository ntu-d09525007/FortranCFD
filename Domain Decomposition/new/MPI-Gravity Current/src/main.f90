module all
use ptr_tree
use tree
implicit none
save
type(manager) :: p 
type(ptr_family) :: pt
end module all

program main
use all
implicit none
include 'mpif.h'
real(8) :: dt
integer :: ierr

  call mpi_init(ierr)

  call mpi_comm_rank(mpi_comm_world,p%glb%mpirank,ierr)
  call mpi_comm_size(mpi_comm_world,p%glb%mpisize,ierr)

  call problem_init
  dt=p%glb%dt

  do 
        
        !call check_memory
        
        p%glb%time = p%glb%time + p%glb%dt 
        p%glb%iter = p%glb%iter + 1 
        call p%sync 
      
        !call deforming_field_velocity
        call p%switch
     
        !call interface_solver
        call rk3_c_solver
        call ns_solver
        call plot
    
        !call p%ls_mv; 
        call p%total_c
        if( mod(p%glb%iter,5)==0 .and. p%glb%mpirank==0 )then
            write(*,*)"========================================="
            write(*,'(A,",",I8,",",F15.6)')trim(p%glb%name),p%glb%iter,p%glb%time
            write(*,*)''
            call print_NS_info
            !call print_LS_info
            call print_c_info
            call print_CPU_info
            call find_momentum
            write(*,*)"========================================="
        endif
        call output
    
    if( p%glb%time > p%glb%t2s ) exit
    
 end do
 
  call mpi_finalize()
 
 contains
 
 include 'included.f90'

end program main
