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
integer :: id, i, j
real(8) :: dt

  call problem_init
  dt=p%glb%dt

  do 
        
        p%glb%time = p%glb%time + p%glb%dt 
        p%glb%iter = p%glb%iter + 1 
        call p%sync 
      
        call deforming_field_velocity
        call p%switch
     
        call interface_solver(1.0d0) ! phi finished
        !$omp parallel do private(i,j)
        do id = 0, p%glb%threads-1
          do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
          do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
            p%of(id)%loc%phi1%tmp(i,j) = p%of(id)%loc%phi%now(i,j)  ! store phi%now in phi1%tmp
            p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%phi1%now(i,j)
          enddo
          enddo
        enddo
        !$omp end parallel do

        call interface_solver(0.5d0) ! phi1 finished
        !$omp parallel do private(i,j)
        do id = 0, p%glb%threads-1
          do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
          do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
            p%of(id)%loc%phi1%now(i,j) = p%of(id)%loc%phi%now(i,j)
            p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%phi2%now(i,j)
          enddo
          enddo
        enddo
        !$omp end parallel do 

        call interface_solver(0.5d0) ! phi2 finished
        !$omp parallel do private(i,j)
        do id = 0, p%glb%threads-1
          do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
          do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
            p%of(id)%loc%phi2%now(i,j) = p%of(id)%loc%phi%now(i,j)
            p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%phi1%tmp(i,j)
          enddo
          enddo
        enddo
        !$omp end parallel do 


        ! call ns_solver
        call plot
    
        ! call p%ls_mv
        if( mod(p%glb%iter,5)==0 )then
            write(*,*)"========================================="
            write(*,'(A,",",I8,",",F15.6)')trim(p%glb%name),p%glb%iter,p%glb%time
            write(*,*)''
            ! call print_NS_info
            call print_LS_info
            call print_CPU_info
            ! call find_momentum
            write(*,*)"========================================="
        endif
        call output
    
    if( p%glb%time > p%glb%t2s ) exit
    
 end do
 
 contains
 
 include 'included.f90'

end program main
