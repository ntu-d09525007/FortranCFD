subroutine Stokes_wave_error()
use all 
!$ use omp_lib
implicit none
integer :: i,j,k,id
real(8) :: a0, e0, c, ka, w
real(8) :: x,y,z
real(8) :: xnorm2,xnormax,ynorm2,ynormax,znorm2,znormax

xnorm2=0.0d0;ynorm2=0.0d0;znorm2=0.0d0;
xnormax=0.0d0;ynormax=0.0d0;znormax=0.0d0;

e0=0.55
ka=2.0d0*dacos(-1.0d0)
c=3.7d0
w=c*ka
a0=e0/ka

!$omp parallel do private(i,j,k,x,y,z), reduction(+:xnorm2,ynorm2,znorm2), reduction(max:xnormax,ynormax,znormax)
do id = 0, p%glb%threads-1

    do k = p%of(id)%loc%ks, p%of(id)%loc%ke
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie 

        x = p%glb%x(i,j,k)
        y = p%glb%y(i,j,k)
        z = p%glb%z(i,j,k)

        if( p%of(id)%loc%phi%now(i,j,k) > 0.0d0 )then
            xnorm2 = xnorm2 + (p%of(id)%loc%vel%x%now(i,j,k)-a0*w*dexp(ka*z)*dcos(ka*(x-c*p%glb%time+0.5d0*p%glb%dx)) )**2.0d0
            ynorm2 = ynorm2 + (p%of(id)%loc%vel%y%now(i,j,k))**2.0d0
            znorm2 = znorm2 + (p%of(id)%loc%vel%z%now(i,j,k)-a0*w*dexp(ka*z)*dsin(ka*(x-c*p%glb%time)))**2.0d0

            xnormax = max( xnormax, abs(p%of(id)%loc%vel%x%now(i,j,k)-a0*w*dexp(ka*z)*dcos(ka*(x-c*p%glb%time+0.5d0*p%glb%dx))))
            ynormax = max( ynormax, abs(p%of(id)%loc%vel%y%now(i,j,k)))
            znormax = max( znormax, abs(p%of(id)%loc%vel%z%now(i,j,k)-a0*w*dexp(ka*z)*dsin(ka*(x-c*p%glb%time))))
        endif

    enddo
    enddo
    enddo

enddo
!$omp end parallel do

xnorm2=dsqrt(xnorm2/real(p%glb%node_x*p%glb%node_y*p%glb%node_z,kind=8))
ynorm2=dsqrt(ynorm2/real(p%glb%node_x*p%glb%node_y*p%glb%node_z,kind=8))
znorm2=dsqrt(znorm2/real(p%glb%node_x*p%glb%node_y*p%glb%node_z,kind=8))

write(*,*)""
write(*,'("X error:",2ES15.4)')xnorm2,xnormax
write(*,'("Y error:",2ES15.4)')ynorm2,ynormax
write(*,'("Z error:",2ES15.4)')znorm2,znormax

end subroutine