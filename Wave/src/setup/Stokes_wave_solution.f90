subroutine Stokes_Wave_3rd_solution()
use all 
!$ use omp_lib
implicit none
integer :: i,j,k,id
real(8) :: a0, e0, c, ka, w
real(8) :: x,y,z

e0=0.55
ka=2.0d0*dacos(-1.0d0)
c=3.7d0
w=c/ka

!========= Geometry ================================================

a0=e0/ka/p%glb%L
!$omp parallel do private(i,j,k,x,y,z)
do id = 0, p%glb%threads-1

    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie 

        x = p%glb%x(i,j,k)
        y = p%glb%y(i,j,k)
        z = p%glb%z(i,j,k)
        
        if( z<a0*dcos(ka**(x-c*p%glb%time))+0.5d0*a0*e0*dcos(2.0*ka**(x-c*p%glb%time))+3.0d0/8.0d0*a0*e0**2.0d0*dcos(3.0d0*ka**(x-c*p%glb%time)) )then
            p%of(id)%loc%phi%now(i,j,k) = 1.0d0
        else
            p%of(id)%loc%phi%now(i,j,k) = -1.0d0
        endif

    enddo
    enddo
    enddo

    call p%of(id)%bc(0,p%of(id)%loc%phi%now)

enddo
!$omp end parallel do

call pt%phi%sync
call level_set_rk3_redis(0,10.0d0*p%glb%dx)

!========= Velocity ================================================

a0=e0/ka/p%glb%U
!$omp parallel do private(i,j,k,x,y,z)
do id = 0, p%glb%threads-1

    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie 

        x = p%glb%x(i,j,k)
        y = p%glb%y(i,j,k)
        z = p%glb%z(i,j,k)

        if( p%of(id)%loc%phi%now(i,j,k) > 0.0d0 )then
            p%of(id)%loc%vel%x%now(i,j,k) = a0*w*dexp(ka*z)*dcos(ka*(x-c*p%glb%time))
            p%of(id)%loc%vel%y%now(i,j,k) = 0.0d0
            p%of(id)%loc%vel%z%now(i,j,k) = a0*w*dexp(ka*z)*dsin(ka*(x-c*p%glb%time))
        endif

    enddo
    enddo
    enddo

enddo
!$omp end parallel do

call pt%vel%sync

end subroutine