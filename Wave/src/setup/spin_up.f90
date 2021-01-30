subroutine spin_up()
use all
!$ use omp_lib
implicit none
integer :: iter,id,i,j,k
real(8) :: norm, norm2

do
    p%glb%iter = p%glb%iter + 1
    call Stokes_Wave_3rd_velocity

    !$omp parallel do 
    do id = 0, p%glb%threads-1
        call p%of(id)%loc%vel%switch
        call p%of(id)%loc%nvel%switch
        call p%of(id)%loc%velsrc%switch
        call p%of(id)%loc%p%switch
    enddo
    !$omp end parallel do

    call ns_solver 

    norm=0.0d0
    norm2=0.0d0
    !$omp parallel do private(i,j,k), reduction(max:norm), reduction(+:norm2)
    do id = 0, p%glb%threads-1
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie 
            if( p%of(id)%loc%phi%now(i,j,k) < 0.0d0 .and. p%glb%z(i,j,k)<1.95 )then
                norm = max( abs( p%of(id)%loc%vel%x%now(i,j,k)-p%of(id)%loc%vel%x%old(i,j,k) ), norm)
                norm = max( abs( p%of(id)%loc%vel%y%now(i,j,k)-p%of(id)%loc%vel%y%old(i,j,k) ), norm)
                norm = max( abs( p%of(id)%loc%vel%z%now(i,j,k)-p%of(id)%loc%vel%z%old(i,j,k) ), norm)
                norm2 = norm2 + ( p%of(id)%loc%vel%x%now(i,j,k)-p%of(id)%loc%vel%x%old(i,j,k) )**2.0d0 &
                              + ( p%of(id)%loc%vel%y%now(i,j,k)-p%of(id)%loc%vel%y%old(i,j,k) )**2.0d0 &
                              + ( p%of(id)%loc%vel%z%now(i,j,k)-p%of(id)%loc%vel%z%old(i,j,k) )**2.0d0
            endif
        enddo
        enddo
        enddo
    enddo
    !$omp end parallel do

    norm2 = dsqrt( norm2 / (p%glb%node_x*p%glb%node_y*p%glb%node_z) )

    if( mod(p%glb%iter,5).eq.0 )then
        write(*,*)"========================================="
        write(*,'("Iter:",I10)')p%glb%iter
        write(*,'("Difference :",2ES15.4)')norm, norm2
        write(*,'("Divergence :",2ES15.4)')p%glb%vel_div,p%glb%vel_sdiv
        write(*,'("PPE iters  :",I15)')p%glb%piter
        write(*,'("PPE error  :",ES15.4)')p%glb%ppe_linf
        write(*,*)"========================================="
    endif

    if( norm<1.0d-5 ) exit

enddo

call Stokes_Wave_3rd_velocity

end subroutine
