subroutine ns_ab_diff_source
implicit none

    call ns_ab_diff_source_sec

end subroutine

subroutine ns_ab_diff_source_sec
use all
!$ use omp_lib
implicit none
integer :: id

!$omp parallel do 
do id = 0, p%glb%threads-1

    call ns_ab_diff_job_sec(p%of(id),p%of(id)%loc%velsrc%x%tmp,p%of(id)%loc%velsrc%y%tmp,&
                                     p%of(id)%loc%vel%x%old,p%of(id)%loc%vel%y%old, &
                                     p%of(id)%loc%phi%old,p%of(id)%loc%normals%curv%old,p%of(id)%loc%delta%old,&
                                     p%of(id)%loc%rho%old,p%of(id)%loc%mu%old,.true.,0.5d0)

    call ns_ab_diff_job_sec(p%of(id),p%of(id)%loc%velsrc%x%tmp,p%of(id)%loc%velsrc%y%tmp,&
                                     p%of(id)%loc%vel%x%tmp,p%of(id)%loc%vel%y%tmp, &
                                     p%of(id)%loc%phi%now,p%of(id)%loc%normals%curv%now,p%of(id)%loc%delta%now,&
                                     p%of(id)%loc%rho%now,p%of(id)%loc%mu%now,.false.,0.5d0)
enddo
!$omp end parallel do
    
end subroutine

subroutine ns_ab_diff_job_sec(q,sx,sy,u,v,phi,curv,delta,rho,mu,reset,alpha)
use all
implicit none
type(job) :: q
real(8), dimension(q%loc%is-q%glb%ghc:q%loc%ie+q%glb%ghc,&
                  &q%loc%js-q%glb%ghc:q%loc%je+q%glb%ghc) :: sx,sy,u,v,phi,curv,delta,rho,mu
integer :: i,j
logical :: reset
real(8) :: alpha
real(8) :: rho_, mu_, delta_, curv_
real(8) :: xx,yy
real(8) :: ux,uy,vx,vy
real(8) :: phix,phiy
real(8) :: dif_x, dif_y

if(reset)then

    do j = q%loc%js, q%loc%je
    do i = q%loc%is, q%loc%ie
        sx(i,j) = 0.0d0
        sy(i,j) = 0.0d0
    enddo
    enddo

endif

do j = q%loc%js, q%loc%je
do i = q%loc%is, q%loc%ie

    rho_   = 0.5d0*(  rho(i,j)+  rho(i+1,j))
    mu_    = 0.5d0*(   mu(i,j)+   mu(i+1,j))
    delta_ = 0.5d0*(delta(i,j)+delta(i+1,j))
    curv_  = 0.5d0*( curv(i,j)+ curv(i+1,j))

    xx = (u(i+1,j)-2.0d0*u(i,j)+u(i-1,j))/q%glb%dx**2.0d0
    yy = (u(i,j+1)-2.0d0*u(i,j)+u(i,j-1))/q%glb%dy**2.0d0

    ux = 0.5d0*(u(i+1,j)-u(i-1,j))/q%glb%dx
    uy = 0.5d0*(u(i,j+1)-u(i,j-1))/q%glb%dy

    vx = 0.5d0*( v(i+1,j)-v(i,j)+v(i+1,j-1)-v(i,j-1) )/q%glb%dx

    phix = (phi(i+1,j)-phi(i,j))/q%glb%dx
    phiy = 0.25d0*(phi(i+1,j+1)-phi(i+1,j-1)+phi(i,j+1)-phi(i,j-1))/q%glb%dy

    dif_x = mu_ / rho_ * xx / q%glb%re + (1.0d0-q%glb%mu_12) * delta_ * phix * 2.0d0*ux / ( rho_ * q%glb%re)
    dif_y = mu_ / rho_ * yy / q%glb%re + (1.0d0-q%glb%mu_12) * delta_ * phiy * (uy+vx)  / ( rho_ * q%glb%re)

    sx(i,j) = sx(i,j) + alpha * ( dif_x+dif_y + q%glb%gx*q%glb%btn_g / q%glb%fr &
            & - q%glb%btn_sf*curv_*delta_*phix / (q%glb%we*rho_)  )

    ! ==========================================================================================

    rho_   = 0.5d0*(  rho(i,j)+  rho(i,j+1))
    mu_    = 0.5d0*(   mu(i,j)+   mu(i,j+1))
    delta_ = 0.5d0*(delta(i,j)+delta(i,j+1))
    curv_  = 0.5d0*( curv(i,j)+ curv(i,j+1))

    xx = (v(i+1,j)-2.0d0*v(i,j)+v(i-1,j))/q%glb%dx**2.0d0
    yy = (v(i,j+1)-2.0d0*v(i,j)+v(i,j-1))/q%glb%dy**2.0d0

    vx = 0.5d0*(v(i+1,j)-v(i-1,j))/q%glb%dx
    vy = 0.5d0*(v(i,j+1)-v(i,j-1))/q%glb%dy

    uy = 0.5d0*( u(i,j+1)-u(i,j)+u(i-1,j+1)-u(i-1,j) )/q%glb%dy

    phix = 0.25d0*(phi(i+1,j)-phi(i-1,j)+phi(i+1,j+1)-phi(i-1,j+1))/q%glb%dx
    phiy = ( phi(i,j+1)-phi(i,j) )/q%glb%dy

    dif_x = mu_ / rho_ * xx / q%glb%re + (1.0d0-q%glb%mu_12) * delta_ * phix*(uy+vx) / ( rho_ * q%glb%re)
    dif_y = mu_ / rho_ * yy / q%glb%re + (1.0d0-q%glb%mu_12) * delta_ * phiy*2.0d0*vy/ ( rho_ * q%glb%re)

    sy(i,j) = sy(i,j) + alpha * ( dif_x+dif_y + q%glb%gy * q%glb%btn_g / q%glb%fr &
            & - q%glb%btn_sf * curv_*delta_ * phiy / (q%glb%we*rho_)  )

enddo
enddo

end subroutine