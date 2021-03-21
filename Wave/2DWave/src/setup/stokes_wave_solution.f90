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

function Stokes_wave_interface(theta,eps,kh,a) result(eta)
implicit none
real(8), intent(in) :: theta,eps,kh,a
real(8) :: eta
real(8) :: c1,s1
real(8) :: eta1,eta2,eta3

c1=dsinh(kh)
s1=dcosh(kh)

eta3=eps**2.0d0*(3.0d0*(8.0d0*s1**6.0d0+24.0d0*s1**4.0d0+24.0d0*s1**2.0d0+9))/(64.0d0*s1**6.0d0)
eta2=eps*c1/s1*(2.0d0*s1**2.0d0+3.0)/(4.0d0*s1**2.0d0)
eta1=1.0d0-eta3

eta = a * ( eta1*dcos(theta) + eta2*dcos(2.0d0*theta) + eta3*dcos(3.0d0*theta) )

end  function