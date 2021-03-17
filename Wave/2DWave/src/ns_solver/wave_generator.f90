subroutine wave_generator(q,sx,sy,w)
use all
implicit none
type(job) :: q
real(8), dimension(q%loc%is-q%glb%ghc:q%loc%ie+q%glb%ghc,&
                  &q%loc%js-q%glb%ghc:q%loc%je+q%glb%ghc) ::sx,sy
real(8) :: w
real(8) :: width,alpha,beta,I1,kh0,D0,D1,D2,x,theta
integer :: i,j

! ------ Wave source input ------ 
width = 2.0d0*dacos(-1.0d0)*p%wa%width
theta = p%wa%theta / 180.0d0 * dacos(-1.0d0)
!--------------------------------

alpha = 0.5d0*0.53d0**2.0d0 - 0.53d0
beta = 5.0d0/(0.5d0*width)**2.0d0
I1 = dsqrt(dacos(-1.0d0)/beta)*dexp(p%wa%k**2.0d0/(4.0d0*beta))
kh0 = p%wa%k*abs(p%glb%ystart)
D0 = 2.0d0*p%wa%L / ( p%wa%w*I1*p%wa%k*(1.0d0-alpha*kh0**2.0d0 ) )
D1 = p%wa%w**2.0d0
D2 = (alpha+1.0d0/3.0d0) * p%wa%k * kh0**3.0d0 / p%glb%fr
D0 = D0 * 2.0d0 * beta * dsin(-p%wa%w*p%glb%time) / p%wa%w / p%glb%fr
D0 = D0 * (D1+D2)

do j = q%loc%js, q%loc%je
do i = q%loc%is, q%loc%ie

    x = p%glb%x(i,j)

    if( abs(x) <= width/2.0d0 .and. q%loc%phi%now(i,j)>0.0d0 )then
        sx(i,j) = sx(i,j) + w * D0*x*dexp(-beta*x**2.0d0) * dcos(theta)
        sy(i,j) = sy(i,j) + w * D0*x*dexp(-beta*x**2.0d0) * dsin(theta)
    endif

end do
end do

end subroutine