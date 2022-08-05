subroutine Stokes_wave_error()
use all 
!$ use omp_lib
implicit none
integer :: i,j,id
real(8) :: x,y
real(8) :: xnorm2,xnormax,ynorm2,ynormax
real(8) :: u,v,kx,ky,wt,kdx

xnorm2=0.0d0;ynorm2=0.0d0
xnormax=0.0d0;ynormax=0.0d0

!$omp parallel do private(i,j,x,y,u,v,kx,ky,wt), reduction(+:xnorm2,ynorm2), reduction(max:xnormax,ynormax)
do id = 0, p%glb%threads-1

    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie 

        x = p%glb%x(i,j)
        y = p%glb%y(i,j)

        kx = p%wa%k * x
        ky = p%wa%k * y
        wt = p%wa%w * p%glb%time

        u = p%wa%U * dexp(ky) * dcos(kx-wt+p%wa%k*p%glb%dx*0.5d0)
        v = p%wa%U * dexp(ky) * dsin(kx-wt)

        if( p%of(id)%loc%phi%now(i,j) > 0.0d0 .and. y<0.0d0 )then
            xnorm2 = xnorm2 + (p%of(id)%loc%vel%x%now(i,j)-u)**2.0d0
            ynorm2 = ynorm2 + (p%of(id)%loc%vel%y%now(i,j)-v)**2.0d0

            xnormax = max( xnormax, abs(p%of(id)%loc%vel%x%now(i,j)-u))
            ynormax = max( ynormax, abs(p%of(id)%loc%vel%y%now(i,j)-v)) 
        endif

    enddo
    enddo

enddo
!$omp end parallel do

xnorm2=dsqrt(xnorm2/real(p%glb%node_x*p%glb%node_y,kind=8))
ynorm2=dsqrt(ynorm2/real(p%glb%node_x*p%glb%node_y,kind=8))

write(*,*)""
write(*,'("X error:",2ES15.4)')xnorm2,xnormax
write(*,'("Y error:",2ES15.4)')ynorm2,ynormax

end subroutine

function Stokes_wave_interface(t,e,kh) result(eta)
implicit none
real(8), intent(in) :: t,e,kh
real(8) :: eta
real(8) :: ap

ap = 1.0d0 /  dtanh(kh)

eta = dcos(t)
eta = eta + e    * ap/4.0 * (3.0*ap**2-1.0)*dcos(2.0*t)
eta = eta - e**2 * 3.0/8.0 * (ap**4-3.0*ap**2+3.0)*dcos(t)
eta = eta + e**2 * 3.0/64.0 *(8.0*ap**6+(ap**2-1.0)**2)*dcos(3.0*t)

end  function

function Stokes_wave_u(t,e,kh,ky) result(u)
implicit none
real(8), intent(in) :: t,e,kh,ky
real(8) :: u,s,c

s = dsinh(kh)
c = dcosh(kh)

u = dcosh(kh+ky)*dcos(t) / c
u = u + 0.75d0 * e * dcosh(2.0d0*(kh+ky)) * dcos(2.0*t) /  (s**3.0d0*c)

end  function

function Stokes_wave_v(t,e,kh,ky) result(v)
implicit none
real(8), intent(in) :: t,e,kh,ky
real(8) :: v,s,c

s = dsinh(kh)
c = dcosh(kh)

v = dsinh(kh+ky)*dsin(t) / c
v = v + 0.75d0 * e * dsinh(2.0d0*(kh+ky)) * dsin(2.0*t) /  (s**3.0d0*c)

end  function