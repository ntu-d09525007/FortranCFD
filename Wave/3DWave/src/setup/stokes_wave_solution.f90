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