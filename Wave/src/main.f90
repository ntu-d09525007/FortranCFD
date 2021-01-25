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
integer :: iter,id,i,j,k
real(8) :: norm,err

    call problem_init
  
    iter=0

    do

        iter = iter + 1
        write(*,*)"Set velocity"
        call Stokes_Wave_3rd_velocity

        write(*,*)"Switch"
        !$omp parallel do 
        do id = 0, p%glb%threads-1
            call p%of(id)%loc%vel%switch
            call p%of(id)%loc%nvel%switch
            call p%of(id)%loc%velsrc%switch
            call p%of(id)%loc%p%switch
        enddo
        !$omp end parallel do

        write(*,*)"Solve NS"
        call ns_solver 

        write(*,*)"Calculate error"
        norm=0.0d0
        !$omp parallel do private(i,j,k,err), reduction(+:norm)
        do id = 0, p%glb%threads-1
            do k = p%of(id)%loc%ks, p%of(id)%loc%ke
            do j = p%of(id)%loc%js, p%of(id)%loc%je
            do i = p%of(id)%loc%is, p%of(id)%loc%ie 
                if( p%of(id)%loc%phi%now(i,j,k) < 0.0d0 )then
                    err = ( p%of(id)%loc%vel%x%now(i,j,k)-p%of(id)%loc%vel%x%old(i,j,k) )**2.0d0
                    err = err + ( p%of(id)%loc%vel%y%now(i,j,k)-p%of(id)%loc%vel%y%old(i,j,k) )**2.0d0
                    err = err + ( p%of(id)%loc%vel%z%now(i,j,k)-p%of(id)%loc%vel%z%old(i,j,k) )**2.0d0
                    norm = norm + err
                endif
            enddo
            enddo
            enddo
        enddo
        !$omp end parallel do

        norm = dsqrt( norm / (p%glb%node_x*p%glb%node_y*p%glb%node_z) )

        if( mod(iter,5).eq.0 )then
            write(*,'("Difference :",ES15.4)')norm
            write(*,'("Divergence :",2ES15.4)')p%glb%vel_div,p%glb%vel_sdiv
            write(*,'("PPE iters  :",I15)')p%glb%piter
            write(*,'("PPE error  :",ES15.4)')p%glb%ppe_linf
        endif

        if( norm<1.0d-6 ) exit

    enddo

    stop

  do 
        
        p%glb%time = p%glb%time + p%glb%dt 
        p%glb%iter = p%glb%iter + 1 
        call p%sync 
      
        !call deforming_field_velocity
        call p%switch
     
        call interface_solver
        call ns_solver
        call plot
    
        call p%ls_mv
        if( mod(p%glb%iter,5)==0 )then
            write(*,*)"========================================="
            write(*,'(A,",",I8,",",F15.6)')trim(p%glb%name),p%glb%iter,p%glb%time
            write(*,*)''
            call print_NS_info
            call print_LS_info
            call print_CPU_info
            !call find_momentum
            write(*,*)"========================================="
        endif
        call output
    
    if( p%glb%time > p%glb%t2s ) exit
    
 end do
 
 contains
 
 include 'included.f90'

end program main
