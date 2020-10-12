subroutine ns_ab_diff_source
implicit none

    call ns_ab_diff_source_sec

end subroutine

subroutine ns_ab_diff_source_sec
implicit none

    call u_source
    call v_source
    
end subroutine

subroutine u_source()
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: rho,mu,delta,xx,yy
real(8) :: ux,uy,vx,phix,phiy,dif_x,dif_y,curv

    !$omp parallel do private(i,j,rho,mu,delta,xx,yy), &
    !$omp& private(ux,uy,vx,phix,phiy,dif_x,dif_y,curv)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            rho = 0.5d0*(p%of(id)%loc%rho%old(i,j)+p%of(id)%loc%rho%old(i+1,j))
            mu = 0.5d0*(p%of(id)%loc%mu%old(i,j)+p%of(id)%loc%mu%old(i+1,j))
            delta = 0.5d0*(p%of(id)%loc%delta%old(i,j)+p%of(id)%loc%delta%old(i+1,j))
            curv = (p%of(id)%loc%normals%curv%old(i,j)+p%of(id)%loc%normals%curv%old(i+1,j))/2.0d0
            
            xx = (p%of(id)%loc%vel%x%old(I+1,J)-2.0d0*p%of(id)%loc%vel%x%old(I,J)+p%of(id)%loc%vel%x%old(I-1,J))/p%glb%dx**2.0d0
            yy = (p%of(id)%loc%vel%x%old(I,J+1)-2.0d0*p%of(id)%loc%vel%x%old(I,J)+p%of(id)%loc%vel%x%old(I,J-1))/p%glb%dy**2.0d0

            ux = 0.5d0*( p%of(id)%loc%vel%x%old(i+1,j)-p%of(id)%loc%vel%x%old(i-1,j) )/p%glb%dx
            uy = 0.5d0*( p%of(id)%loc%vel%x%old(i,j+1)-p%of(id)%loc%vel%x%old(i,j-1) )/p%glb%dy

            vx = 0.5d0*( p%of(id)%loc%vel%y%old(i+1,j)-p%of(id)%loc%vel%y%old(i,j)+p%of(id)%loc%vel%y%old(i+1,j-1)-p%of(id)%loc%vel%y%old(i,j-1) )/p%glb%dx

            phix = (p%of(id)%loc%phi%old(i+1,j)-p%of(id)%loc%phi%old(i,j))/p%glb%dx
            phiy = 0.25d0*(p%of(id)%loc%phi%old(i+1,j+1)-p%of(id)%loc%phi%old(i+1,j-1)+p%of(id)%loc%phi%old(i,j+1)-p%of(id)%loc%phi%old(i,j-1))/p%glb%dy

            dif_x = 0.0d0!mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*2.0d0*ux/(rho*p%glb%re)
            dif_y = 0.0d0!mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*(uy+vx)/(rho*p%glb%re)

            p%of(id)%loc%velsrc%x%tmp(i,j) = dif_x + dif_y + p%glb%gx*p%glb%btn_g / p%glb%fr &
            & - p%glb%btn_sf*curv*delta*phix / (p%glb%we*rho) 
            
            ! =================================================
            
            rho = 0.5d0*(p%of(id)%loc%rho%old(i,j)+p%of(id)%loc%rho%old(i+1,j))
            mu = 0.5d0*(p%of(id)%loc%mu%old(i,j)+p%of(id)%loc%mu%old(i+1,j))
            delta = 0.5d0*(p%of(id)%loc%delta%now(i,j)+p%of(id)%loc%delta%now(i+1,j))
            curv = (p%of(id)%loc%normals%curv%now(i,j)+p%of(id)%loc%normals%curv%now(i+1,j))/2.0d0
            
            xx = (p%of(id)%loc%vel%x%tmp(I+1,J)-2.0d0*p%of(id)%loc%vel%x%tmp(I,J)+p%of(id)%loc%vel%x%tmp(I-1,J))/p%glb%dx**2.0d0
            yy = (p%of(id)%loc%vel%x%tmp(I,J+1)-2.0d0*p%of(id)%loc%vel%x%tmp(I,J)+p%of(id)%loc%vel%x%tmp(I,J-1))/p%glb%dy**2.0d0

            ux = 0.5d0*( p%of(id)%loc%vel%x%tmp(i+1,j)-p%of(id)%loc%vel%x%tmp(i-1,j) )/p%glb%dx
            uy = 0.5d0*( p%of(id)%loc%vel%x%tmp(i,j+1)-p%of(id)%loc%vel%x%tmp(i,j-1) )/p%glb%dy

            vx = 0.5d0*( p%of(id)%loc%vel%y%tmp(i+1,j)-p%of(id)%loc%vel%y%tmp(i,j)+p%of(id)%loc%vel%y%tmp(i+1,j-1)-p%of(id)%loc%vel%y%tmp(i,j-1) )/p%glb%dx

            phix = (p%of(id)%loc%phi%now(i+1,j)-p%of(id)%loc%phi%now(i,j))/p%glb%dx
            phiy = 0.25d0*(p%of(id)%loc%phi%now(i+1,j+1)-p%of(id)%loc%phi%now(i+1,j-1)+p%of(id)%loc%phi%now(i,j+1)-p%of(id)%loc%phi%now(i,j-1))/p%glb%dy

            dif_x = 0.0d0!mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*2.0d0*ux/(rho*p%glb%re)
            dif_y = 0.0d0!mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*(uy+vx)/(rho*p%glb%re)
    
            p%of(id)%loc%velsrc%x%tmp(i,j) = ( p%of(id)%loc%velsrc%x%tmp(i,j) + &
            & dif_x + dif_y + p%glb%gx*p%glb%btn_g / p%glb%fr &
            & - p%glb%btn_sf*curv*delta*phix / (p%glb%we*rho) )/2.0d0
            
        end do
        end do 
     
    enddo
    !$omp end parallel do

end subroutine

subroutine v_source()
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: rho,mu,delta,xx,yy
real(8) :: vx,vy,uy,phix,phiy,dif_x,dif_y,curv

    !$omp parallel do private(i,j,rho,mu,delta,xx,yy), &
    !$omp& private(vx,vy,uy,phix,phiy,dif_x,dif_y,curv)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            rho = 0.5d0*(p%of(id)%loc%rho%old(i,j)+p%of(id)%loc%rho%old(i,j+1))
            mu = 0.5d0*(p%of(id)%loc%mu%old(i,j)+p%of(id)%loc%mu%old(i,j+1))
            delta = 0.5d0*(p%of(id)%loc%delta%old(i,j)+p%of(id)%loc%delta%old(i,j+1))
            curv = (p%of(id)%loc%normals%curv%old(i,j)+p%of(id)%loc%normals%curv%old(i,j+1))/2.0d0
            
            xx = (p%of(id)%loc%vel%y%old(I+1,J)-2.0d0*p%of(id)%loc%vel%y%old(I,J)+p%of(id)%loc%vel%y%old(I-1,J))/p%glb%dx**2.0d0
            yy = (p%of(id)%loc%vel%y%old(I,J+1)-2.0d0*p%of(id)%loc%vel%y%old(I,J)+p%of(id)%loc%vel%y%old(I,J-1))/p%glb%dy**2.0d0

            vx = 0.5d0*( p%of(id)%loc%vel%y%old(i+1,j)-p%of(id)%loc%vel%y%old(i-1,j) )/p%glb%dx
            vy = 0.5d0*( p%of(id)%loc%vel%y%old(i,j+1)-p%of(id)%loc%vel%y%old(i,j-1) )/p%glb%dy

            uy = 0.5d0*( p%of(id)%loc%vel%x%old(i,j+1)-p%of(id)%loc%vel%x%old(i,j)+p%of(id)%loc%vel%x%old(i-1,j+1)-p%of(id)%loc%vel%x%old(i-1,j) )/p%glb%dy

            phix = 0.25d0*(p%of(id)%loc%phi%old(i+1,j)-p%of(id)%loc%phi%old(i-1,j)+p%of(id)%loc%phi%old(i+1,j+1)-p%of(id)%loc%phi%old(i-1,j+1))/p%glb%dx
            phiy = ( p%of(id)%loc%phi%old(i,j+1)-p%of(id)%loc%phi%old(i,j) )/p%glb%dy

            dif_x = 0.0d0!mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*(uy+vx)/(rho*p%glb%re)
            dif_y = 0.0d0!mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*2.0d0*vy/(rho*p%glb%re)

            p%of(id)%loc%velsrc%y%tmp(i,j) = dif_x + dif_y + p%glb%gy*p%glb%btn_g / p%glb%fr &
            & - p%glb%btn_sf*curv*delta*phiy / (p%glb%we*rho)
            
            ! =========================================================================
            
            rho = 0.5d0*(p%of(id)%loc%rho%old(i,j)+p%of(id)%loc%rho%old(i,j+1))
            mu = 0.5d0*(p%of(id)%loc%mu%old(i,j)+p%of(id)%loc%mu%old(i,j+1))
            delta = 0.5d0*(p%of(id)%loc%delta%now(i,j)+p%of(id)%loc%delta%now(i,j+1))
            curv = (p%of(id)%loc%normals%curv%now(i,j)+p%of(id)%loc%normals%curv%now(i,j+1))/2.0d0
            
            xx = (p%of(id)%loc%vel%y%tmp(I+1,J)-2.0d0*p%of(id)%loc%vel%y%tmp(I,J)+p%of(id)%loc%vel%y%tmp(I-1,J))/p%glb%dx**2.0d0
            yy = (p%of(id)%loc%vel%y%tmp(I,J+1)-2.0d0*p%of(id)%loc%vel%y%tmp(I,J)+p%of(id)%loc%vel%y%tmp(I,J-1))/p%glb%dy**2.0d0

            vx = 0.5d0*( p%of(id)%loc%vel%y%tmp(i+1,j)-p%of(id)%loc%vel%y%tmp(i-1,j) )/p%glb%dx
            vy = 0.5d0*( p%of(id)%loc%vel%y%tmp(i,j+1)-p%of(id)%loc%vel%y%tmp(i,j-1) )/p%glb%dy

            uy = 0.5d0*( p%of(id)%loc%vel%x%tmp(i,j+1)-p%of(id)%loc%vel%x%tmp(i,j)+p%of(id)%loc%vel%x%tmp(i-1,j+1)-p%of(id)%loc%vel%x%tmp(i-1,j) )/p%glb%dy

            phix = 0.25d0*(p%of(id)%loc%phi%now(i+1,j)-p%of(id)%loc%phi%now(i-1,j)+p%of(id)%loc%phi%now(i+1,j+1)-p%of(id)%loc%phi%now(i-1,j+1))/p%glb%dx
            phiy = ( p%of(id)%loc%phi%now(i,j+1)-p%of(id)%loc%phi%now(i,j) )/p%glb%dy
 
            dif_x = 0.0d0!mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*(uy+vx)/(rho*p%glb%re)
            dif_y = 0.0d0!mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*2.0d0*vy/(rho*p%glb%re)

            p%of(id)%loc%velsrc%y%tmp(i,j) = ( p%of(id)%loc%velsrc%y%tmp(i,j) + &
            & dif_x + dif_y + p%glb%gy*p%glb%btn_g / p%glb%fr &
            & - p%glb%btn_sf*curv*delta*phiy / (p%glb%we*rho) )/2.0d0
            
        end do
        end do
     
    enddo       
    !$omp end parallel do

end subroutine
