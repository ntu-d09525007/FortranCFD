subroutine deforming_field_velocity()
use all 
!$ use omp_lib
implicit none
integer :: i,j,k,id
real(8) :: x, y, z, xx, yy, zz, pi, ct, linf, linf2

    pi = dacos(-1.0_8)
    ct = dcos( pi*p%glb%time / p%glb%t2s )

    !$omp parallel do collapse(3), private(x,y,z,xx,yy,zz)   
    do k = p%loc%ks, p%loc%ke
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
    
        x = p%glb%x(i)
        y = p%glb%y(j)
        z = p%glb%z(k)
        
        xx = 0.5_8 * ( p%glb%x(i) + p%glb%x(i+1) )
        yy = 0.5_8 * ( p%glb%y(j) + p%glb%y(j+1) )
        zz = 0.5_8 * ( p%glb%z(k) + p%glb%z(k+1) )
        
        p%loc%nvel%x%now(i,j,k) =   2.0d0*dsin(pi*x)**2.0d0*dsin(2.0d0*pi*y)*dsin(2.0d0*pi*z) * ct
        p%loc%nvel%y%now(i,j,k) =      -  dsin(2.0d0*pi*x)*dsin(pi*y)**2.0d0*dsin(2.0d0*pi*z) * ct
        p%loc%nvel%z%now(i,j,k) =      -  dsin(2.0d0*pi*x)*dsin(2.0d0*pi*y)*dsin(pi*z)**2.0d0 * ct
        
        p%loc%vel%x%now(i,j,k) =   2.0d0*dsin(pi*xx)**2.0d0*dsin(2.0d0*pi*y)*dsin(2.0d0*pi*z) * ct
        p%loc%vel%y%now(i,j,k) =      -  dsin(2.0d0*pi*x)*dsin(pi*yy)**2.0d0*dsin(2.0d0*pi*z) * ct
        p%loc%vel%z%now(i,j,k) =      -  dsin(2.0d0*pi*x)*dsin(2.0d0*pi*y)*dsin(pi*zz)**2.0d0 * ct
        
    enddo
    enddo
    enddo
    !$omp end parallel do
    
    call bc(p%loc%nvel%x%now)
    call bc(p%loc%nvel%y%now)
    call bc(p%loc%nvel%z%now)
    
    call bc(p%loc%vel%x%now)
    call bc(p%loc%vel%y%now)
    call bc(p%loc%vel%z%now)

end subroutine

