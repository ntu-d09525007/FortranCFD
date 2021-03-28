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