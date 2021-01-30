subroutine ns_ab_diff_source
implicit none

    call ns_ab_diff_source_sec

end subroutine

subroutine ns_ab_diff_source_sec
implicit none

    call u_source
    call v_source
    call w_source 
    
end subroutine

subroutine u_source()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: rho,mu,delta,xx,yy,zz
real(8) :: ux,uy,uz,vx,wx,phix,phiy,phiz,dif_x,dif_y,dif_z,curv

    !$omp parallel do private(i,j,k,rho,mu,delta,xx,yy,zz), &
    !$omp& private(ux,uy,uz,vx,wx,phix,phiy,phiz,dif_x,dif_y,dif_z,curv)
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            rho = 0.5d0*(p%of(id)%loc%rho%old(i,j,k)+p%of(id)%loc%rho%old(i+1,j,k))
            mu = 0.5d0*(p%of(id)%loc%mu%old(i,j,k)+p%of(id)%loc%mu%old(i+1,j,k))
            delta = 0.5d0*(p%of(id)%loc%delta%old(i,j,k)+p%of(id)%loc%delta%old(i+1,j,k))
            curv = (p%of(id)%loc%normals%curv%old(i,j,k)+p%of(id)%loc%normals%curv%old(i+1,j,k))/2.0d0
            
            xx = (p%of(id)%loc%vel%x%old(I+1,J,K)-2.0d0*p%of(id)%loc%vel%x%old(I,J,K)+p%of(id)%loc%vel%x%old(I-1,J,K))/p%glb%dx**2.0d0
            yy = (p%of(id)%loc%vel%x%old(I,J+1,K)-2.0d0*p%of(id)%loc%vel%x%old(I,J,K)+p%of(id)%loc%vel%x%old(I,J-1,K))/p%glb%dy**2.0d0
            zz = (p%of(id)%loc%vel%x%old(I,J,K+1)-2.0d0*p%of(id)%loc%vel%x%old(I,J,K)+p%of(id)%loc%vel%x%old(I,J,K-1))/p%glb%dz**2.0d0
            
            ux = 0.5d0*( p%of(id)%loc%vel%x%old(i+1,j,k)-p%of(id)%loc%vel%x%old(i-1,j,k) )/p%glb%dx
            uy = 0.5d0*( p%of(id)%loc%vel%x%old(i,j+1,k)-p%of(id)%loc%vel%x%old(i,j-1,k) )/p%glb%dy
            uz = 0.5d0*( p%of(id)%loc%vel%x%old(i,j,k+1)-p%of(id)%loc%vel%x%old(i,j,k-1) )/p%glb%dz

            if(k==p%glb%node_z)then
                uz=1.0d0
            endif
    
            vx = 0.5d0*( p%of(id)%loc%vel%y%old(i+1,j,k)-p%of(id)%loc%vel%y%old(i,j,k)+p%of(id)%loc%vel%y%old(i+1,j-1,k)-p%of(id)%loc%vel%y%old(i,j-1,k) )/p%glb%dx
            wx = 0.5d0*( p%of(id)%loc%vel%z%old(i+1,j,k)-p%of(id)%loc%vel%z%old(i,j,k)+p%of(id)%loc%vel%z%old(i+1,j,k-1)-p%of(id)%loc%vel%z%old(i,j,k-1) )/p%glb%dx
    
            phix = (p%of(id)%loc%phi%old(i+1,j,k)-p%of(id)%loc%phi%old(i,j,k))/p%glb%dx
            phiy = 0.25d0*(p%of(id)%loc%phi%old(i+1,j+1,k)-p%of(id)%loc%phi%old(i+1,j-1,k)+p%of(id)%loc%phi%old(i,j+1,k)-p%of(id)%loc%phi%old(i,j-1,k))/p%glb%dy
            phiz = 0.25d0*(p%of(id)%loc%phi%old(i+1,j,k+1)-p%of(id)%loc%phi%old(i+1,j,k-1)+p%of(id)%loc%phi%old(i,j,k+1)-p%of(id)%loc%phi%old(i,j,k-1))/p%glb%dz
            
            dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*2.0d0*ux/(rho*p%glb%re)
            dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*(uy+vx)/(rho*p%glb%re)
            dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*(uz+wx)/(rho*p%glb%re)
            
            p%of(id)%loc%velsrc%x%tmp(i,j,k) = dif_x + dif_y + dif_z &
            & + p%glb%gx*p%glb%btn_g / p%glb%fr &
            & - p%glb%btn_sf*curv*delta*phix / (p%glb%we*rho) 
            
            ! =================================================
            
            rho = 0.5d0*(p%of(id)%loc%rho%old(i,j,k)+p%of(id)%loc%rho%old(i+1,j,k))
            mu = 0.5d0*(p%of(id)%loc%mu%old(i,j,k)+p%of(id)%loc%mu%old(i+1,j,k))
            delta = 0.5d0*(p%of(id)%loc%delta%now(i,j,k)+p%of(id)%loc%delta%now(i+1,j,k))
            curv = (p%of(id)%loc%normals%curv%now(i,j,k)+p%of(id)%loc%normals%curv%now(i+1,j,k))/2.0d0
            
            xx = (p%of(id)%loc%vel%x%tmp(I+1,J,K)-2.0d0*p%of(id)%loc%vel%x%tmp(I,J,K)+p%of(id)%loc%vel%x%tmp(I-1,J,K))/p%glb%dx**2.0d0
            yy = (p%of(id)%loc%vel%x%tmp(I,J+1,K)-2.0d0*p%of(id)%loc%vel%x%tmp(I,J,K)+p%of(id)%loc%vel%x%tmp(I,J-1,K))/p%glb%dy**2.0d0
            zz = (p%of(id)%loc%vel%x%tmp(I,J,K+1)-2.0d0*p%of(id)%loc%vel%x%tmp(I,J,K)+p%of(id)%loc%vel%x%tmp(I,J,K-1))/p%glb%dz**2.0d0
            
            ux = 0.5d0*( p%of(id)%loc%vel%x%tmp(i+1,j,k)-p%of(id)%loc%vel%x%tmp(i-1,j,k) )/p%glb%dx
            uy = 0.5d0*( p%of(id)%loc%vel%x%tmp(i,j+1,k)-p%of(id)%loc%vel%x%tmp(i,j-1,k) )/p%glb%dy
            uz = 0.5d0*( p%of(id)%loc%vel%x%tmp(i,j,k+1)-p%of(id)%loc%vel%x%tmp(i,j,k-1) )/p%glb%dz
    
            if(k==p%glb%node_z)then
                uz=1.0d0
            endif

            vx = 0.5d0*( p%of(id)%loc%vel%y%tmp(i+1,j,k)-p%of(id)%loc%vel%y%tmp(i,j,k)+p%of(id)%loc%vel%y%tmp(i+1,j-1,k)-p%of(id)%loc%vel%y%tmp(i,j-1,k) )/p%glb%dx
            wx = 0.5d0*( p%of(id)%loc%vel%z%tmp(i+1,j,k)-p%of(id)%loc%vel%z%tmp(i,j,k)+p%of(id)%loc%vel%z%tmp(i+1,j,k-1)-p%of(id)%loc%vel%z%tmp(i,j,k-1) )/p%glb%dx
    
            phix = (p%of(id)%loc%phi%now(i+1,j,k)-p%of(id)%loc%phi%now(i,j,k))/p%glb%dx
            phiy = 0.25d0*(p%of(id)%loc%phi%now(i+1,j+1,k)-p%of(id)%loc%phi%now(i+1,j-1,k)+p%of(id)%loc%phi%now(i,j+1,k)-p%of(id)%loc%phi%now(i,j-1,k))/p%glb%dy
            phiz = 0.25d0*(p%of(id)%loc%phi%now(i+1,j,k+1)-p%of(id)%loc%phi%now(i+1,j,k-1)+p%of(id)%loc%phi%now(i,j,k+1)-p%of(id)%loc%phi%now(i,j,k-1))/p%glb%dz
            
            dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*2.0d0*ux/(rho*p%glb%re)
            dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*(uy+vx)/(rho*p%glb%re)
            dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*(uz+wx)/(rho*p%glb%re)
            
            p%of(id)%loc%velsrc%x%tmp(i,j,k) = ( p%of(id)%loc%velsrc%x%tmp(i,j,k) + &
            & dif_x + dif_y + dif_z &
            & + p%glb%gx*p%glb%btn_g / p%glb%fr &
            & - p%glb%btn_sf*curv*delta*phix / (p%glb%we*rho) )/2.0d0
            
        end do
        end do
        end do 
     
    enddo
    !$omp end parallel do

end subroutine

subroutine v_source()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: rho,mu,delta,xx,yy,zz
real(8) :: vx,vy,vz,wy,uy,phix,phiy,phiz,dif_x,dif_y,dif_z,curv

    !$omp parallel do private(i,j,k,rho,mu,delta,xx,yy,zz), &
    !$omp& private(vx,vy,vz,wy,uy,phix,phiy,phiz,dif_x,dif_y,dif_z,curv)
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            rho = 0.5d0*(p%of(id)%loc%rho%old(i,j,k)+p%of(id)%loc%rho%old(i,j+1,k))
            mu = 0.5d0*(p%of(id)%loc%mu%old(i,j,k)+p%of(id)%loc%mu%old(i,j+1,k))
            delta = 0.5d0*(p%of(id)%loc%delta%old(i,j,k)+p%of(id)%loc%delta%old(i,j+1,k))
            curv = (p%of(id)%loc%normals%curv%old(i,j,k)+p%of(id)%loc%normals%curv%old(i,j+1,k))/2.0d0
            
            xx = (p%of(id)%loc%vel%y%old(I+1,J,K)-2.0d0*p%of(id)%loc%vel%y%old(I,J,K)+p%of(id)%loc%vel%y%old(I-1,J,K))/p%glb%dx**2.0d0
            yy = (p%of(id)%loc%vel%y%old(I,J+1,K)-2.0d0*p%of(id)%loc%vel%y%old(I,J,K)+p%of(id)%loc%vel%y%old(I,J-1,K))/p%glb%dy**2.0d0
            zz = (p%of(id)%loc%vel%y%old(I,J,K+1)-2.0d0*p%of(id)%loc%vel%y%old(I,J,K)+p%of(id)%loc%vel%y%old(I,J,K-1))/p%glb%dz**2.0d0
            
            vx = 0.5d0*( p%of(id)%loc%vel%y%old(i+1,j,k)-p%of(id)%loc%vel%y%old(i-1,j,k) )/p%glb%dx
            vy = 0.5d0*( p%of(id)%loc%vel%y%old(i,j+1,k)-p%of(id)%loc%vel%y%old(i,j-1,k) )/p%glb%dy
            vz = 0.5d0*( p%of(id)%loc%vel%y%old(i,j,k+1)-p%of(id)%loc%vel%y%old(i,j,k-1) )/p%glb%dz
    
            if(k==p%glb%node_z)then
                vz=0.0d0
            endif

            uy = 0.5d0*( p%of(id)%loc%vel%x%old(i,j+1,k)-p%of(id)%loc%vel%x%old(i,j,k)+p%of(id)%loc%vel%x%old(i-1,j+1,k)-p%of(id)%loc%vel%x%old(i-1,j,k) )/p%glb%dy
            wy = 0.5d0*( p%of(id)%loc%vel%z%old(i,j+1,k)-p%of(id)%loc%vel%z%old(i,j,k)+p%of(id)%loc%vel%z%old(i,j+1,k-1)-p%of(id)%loc%vel%z%old(i,j,k-1) )/p%glb%dy
    
            phix = 0.25d0*(p%of(id)%loc%phi%old(i+1,j,k)-p%of(id)%loc%phi%old(i-1,j,k)+p%of(id)%loc%phi%old(i+1,j+1,k)-p%of(id)%loc%phi%old(i-1,j+1,k))/p%glb%dx
            phiy = ( p%of(id)%loc%phi%old(i,j+1,k)-p%of(id)%loc%phi%old(i,j,k) )/p%glb%dy
            phiz = 0.25d0*(p%of(id)%loc%phi%old(i,j+1,k+1)-p%of(id)%loc%phi%old(i,j+1,k-1)+p%of(id)%loc%phi%old(i,j,k+1)-p%of(id)%loc%phi%old(i,j,k-1))/p%glb%dz
            
            dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*(uy+vx)/(rho*p%glb%re)
            dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*2.0d0*vy/(rho*p%glb%re)
            dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*(wy+vz)/(rho*p%glb%re)
            
            p%of(id)%loc%velsrc%y%tmp(i,j,k) = dif_x + dif_y + dif_z &
            & + p%glb%gy*p%glb%btn_g / p%glb%fr &
            & - p%glb%btn_sf*curv*delta*phiy / (p%glb%we*rho)
            
            ! =========================================================================
            
            rho = 0.5d0*(p%of(id)%loc%rho%old(i,j,k)+p%of(id)%loc%rho%old(i,j+1,k))
            mu = 0.5d0*(p%of(id)%loc%mu%old(i,j,k)+p%of(id)%loc%mu%old(i,j+1,k))
            delta = 0.5d0*(p%of(id)%loc%delta%now(i,j,k)+p%of(id)%loc%delta%now(i,j+1,k))
            curv = (p%of(id)%loc%normals%curv%now(i,j,k)+p%of(id)%loc%normals%curv%now(i,j+1,k))/2.0d0
            
            xx = (p%of(id)%loc%vel%y%tmp(I+1,J,K)-2.0d0*p%of(id)%loc%vel%y%tmp(I,J,K)+p%of(id)%loc%vel%y%tmp(I-1,J,K))/p%glb%dx**2.0d0
            yy = (p%of(id)%loc%vel%y%tmp(I,J+1,K)-2.0d0*p%of(id)%loc%vel%y%tmp(I,J,K)+p%of(id)%loc%vel%y%tmp(I,J-1,K))/p%glb%dy**2.0d0
            zz = (p%of(id)%loc%vel%y%tmp(I,J,K+1)-2.0d0*p%of(id)%loc%vel%y%tmp(I,J,K)+p%of(id)%loc%vel%y%tmp(I,J,K-1))/p%glb%dz**2.0d0
            
            vx = 0.5d0*( p%of(id)%loc%vel%y%tmp(i+1,j,k)-p%of(id)%loc%vel%y%tmp(i-1,j,k) )/p%glb%dx
            vy = 0.5d0*( p%of(id)%loc%vel%y%tmp(i,j+1,k)-p%of(id)%loc%vel%y%tmp(i,j-1,k) )/p%glb%dy
            vz = 0.5d0*( p%of(id)%loc%vel%y%tmp(i,j,k+1)-p%of(id)%loc%vel%y%tmp(i,j,k-1) )/p%glb%dz
    
            if(k==p%glb%node_z)then
                vz=0.0d0
            endif

            uy = 0.5d0*( p%of(id)%loc%vel%x%tmp(i,j+1,k)-p%of(id)%loc%vel%x%tmp(i,j,k)+p%of(id)%loc%vel%x%tmp(i-1,j+1,k)-p%of(id)%loc%vel%x%tmp(i-1,j,k) )/p%glb%dy
            wy = 0.5d0*( p%of(id)%loc%vel%z%tmp(i,j+1,k)-p%of(id)%loc%vel%z%tmp(i,j,k)+p%of(id)%loc%vel%z%tmp(i,j+1,k-1)-p%of(id)%loc%vel%z%tmp(i,j,k-1) )/p%glb%dy
    
            phix = 0.25d0*(p%of(id)%loc%phi%now(i+1,j,k)-p%of(id)%loc%phi%now(i-1,j,k)+p%of(id)%loc%phi%now(i+1,j+1,k)-p%of(id)%loc%phi%now(i-1,j+1,k))/p%glb%dx
            phiy = ( p%of(id)%loc%phi%now(i,j+1,k)-p%of(id)%loc%phi%now(i,j,k) )/p%glb%dy
            phiz = 0.25d0*(p%of(id)%loc%phi%now(i,j+1,k+1)-p%of(id)%loc%phi%now(i,j+1,k-1)+p%of(id)%loc%phi%now(i,j,k+1)-p%of(id)%loc%phi%now(i,j,k-1))/p%glb%dz
            
            dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*(uy+vx)/(rho*p%glb%re)
            dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*2.0d0*vy/(rho*p%glb%re)
            dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*(wy+vz)/(rho*p%glb%re)
            
            p%of(id)%loc%velsrc%y%tmp(i,j,k) = ( p%of(id)%loc%velsrc%y%tmp(i,j,k) + &
            & dif_x + dif_y + dif_z + &
            & p%glb%gy*p%glb%btn_g / p%glb%fr &
            & - p%glb%btn_sf*curv*delta*phiy / (p%glb%we*rho) )/2.0d0
            
        end do
        end do
        end do 
     
    enddo       
    !$omp end parallel do

end subroutine

subroutine w_source()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: rho,mu,delta,xx,yy,zz
real(8) :: wx,wy,wz,uz,vz,phix,phiy,phiz,dif_x,dif_y,dif_z,curv

    !$omp parallel do private(i,j,k,rho,mu,delta,xx,yy,zz), &
    !$omp& private(wx,wy,wz,uz,vz,phix,phiy,phiz,dif_x,dif_y,dif_z,curv)
    do id = 0, p%glb%threads-1
        
        do k = p%of(id)%loc%ks, p%of(id)%loc%ke
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            rho = 0.5d0*(p%of(id)%loc%rho%old(i,j,k)+p%of(id)%loc%rho%old(i,j,k+1))
            mu = 0.5d0*(p%of(id)%loc%mu%old(i,j,k)+p%of(id)%loc%mu%old(i,j,k+1))
            delta = 0.5d0*(p%of(id)%loc%delta%old(i,j,k)+p%of(id)%loc%delta%old(i,j,k+1))
            curv = (p%of(id)%loc%normals%curv%old(i,j,k)+p%of(id)%loc%normals%curv%old(i,j+1,k))/2.0d0
            
            xx = (p%of(id)%loc%vel%z%old(I+1,J,K)-2.0d0*p%of(id)%loc%vel%z%old(I,J,K)+p%of(id)%loc%vel%z%old(I-1,J,K))/p%glb%dx**2.0d0
            yy = (p%of(id)%loc%vel%z%old(I,J+1,K)-2.0d0*p%of(id)%loc%vel%z%old(I,J,K)+p%of(id)%loc%vel%z%old(I,J-1,K))/p%glb%dy**2.0d0
            zz = (p%of(id)%loc%vel%z%old(I,J,K+1)-2.0d0*p%of(id)%loc%vel%z%old(I,J,K)+p%of(id)%loc%vel%z%old(I,J,K-1))/p%glb%dz**2.0d0
            
            wx = 0.5d0*( p%of(id)%loc%vel%z%old(i+1,j,k)-p%of(id)%loc%vel%z%old(i-1,j,k) )/p%glb%dx
            wy = 0.5d0*( p%of(id)%loc%vel%z%old(i,j+1,k)-p%of(id)%loc%vel%z%old(i,j-1,k) )/p%glb%dy
            wz = 0.5d0*( p%of(id)%loc%vel%z%old(i,j,k+1)-p%of(id)%loc%vel%z%old(i,j,k-1) )/p%glb%dz
    
            uz = 0.5d0*( p%of(id)%loc%vel%x%old(i,j,k+1)-p%of(id)%loc%vel%x%old(i,j,k)+p%of(id)%loc%vel%x%old(i-1,j,k+1)-p%of(id)%loc%vel%x%old(i-1,j,k) )/p%glb%dz
            vz = 0.5d0*( p%of(id)%loc%vel%y%old(i,j,k+1)-p%of(id)%loc%vel%y%old(i,j,k)+p%of(id)%loc%vel%y%old(i,j-1,k+1)-p%of(id)%loc%vel%y%old(i,j-1,k) )/p%glb%dz
    
            phix = 0.25d0*( p%of(id)%loc%phi%old(i+1,j,k+1)-p%of(id)%loc%phi%old(i-1,j,k+1) + p%of(id)%loc%phi%old(i+1,j,k)-p%of(id)%loc%phi%old(i-1,j,k) )/p%glb%dx
            phiy = 0.25d0*( p%of(id)%loc%phi%old(i,j+1,k+1)-p%of(id)%loc%phi%old(i,j-1,k+1) + p%of(id)%loc%phi%old(i,j+1,k)-p%of(id)%loc%phi%old(i,j-1,k) )/p%glb%dy
            phiz = ( p%of(id)%loc%phi%old(i,j,k+1)-p%of(id)%loc%phi%old(i,j,k) )/p%glb%dz
            
            dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*(uz+wx)/(rho*p%glb%re)
            dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*(vz+wy)/(rho*p%glb%re)
            dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*2.0d0*wz/(rho*p%glb%re)
            
            p%of(id)%loc%velsrc%z%tmp(i,j,k) = dif_x + dif_y + dif_z &
            & + p%glb%gz*p%glb%btn_g / p%glb%fr &
            & - p%glb%btn_sf*curv*delta*phiz / (p%glb%we*rho) 
            
            ! ==========================================================================
            
            rho = 0.5d0*(p%of(id)%loc%rho%old(i,j,k)+p%of(id)%loc%rho%old(i,j,k+1))
            mu = 0.5d0*(p%of(id)%loc%mu%old(i,j,k)+p%of(id)%loc%mu%old(i,j,k+1))
            delta = 0.5d0*(p%of(id)%loc%delta%now(i,j,k)+p%of(id)%loc%delta%now(i,j,k+1))
            curv = (p%of(id)%loc%normals%curv%now(i,j,k)+p%of(id)%loc%normals%curv%now(i,j+1,k))/2.0d0
            
            xx = (p%of(id)%loc%vel%z%tmp(I+1,J,K)-2.0d0*p%of(id)%loc%vel%z%tmp(I,J,K)+p%of(id)%loc%vel%z%tmp(I-1,J,K))/p%glb%dx**2.0d0
            yy = (p%of(id)%loc%vel%z%tmp(I,J+1,K)-2.0d0*p%of(id)%loc%vel%z%tmp(I,J,K)+p%of(id)%loc%vel%z%tmp(I,J-1,K))/p%glb%dy**2.0d0
            zz = (p%of(id)%loc%vel%z%tmp(I,J,K+1)-2.0d0*p%of(id)%loc%vel%z%tmp(I,J,K)+p%of(id)%loc%vel%z%tmp(I,J,K-1))/p%glb%dz**2.0d0
            
            wx = 0.5d0*( p%of(id)%loc%vel%z%tmp(i+1,j,k)-p%of(id)%loc%vel%z%tmp(i-1,j,k) )/p%glb%dx
            wy = 0.5d0*( p%of(id)%loc%vel%z%tmp(i,j+1,k)-p%of(id)%loc%vel%z%tmp(i,j-1,k) )/p%glb%dy
            wz = 0.5d0*( p%of(id)%loc%vel%z%tmp(i,j,k+1)-p%of(id)%loc%vel%z%tmp(i,j,k-1) )/p%glb%dz
    
            uz = 0.5d0*( p%of(id)%loc%vel%x%tmp(i,j,k+1)-p%of(id)%loc%vel%x%tmp(i,j,k)+p%of(id)%loc%vel%x%tmp(i-1,j,k+1)-p%of(id)%loc%vel%x%tmp(i-1,j,k) )/p%glb%dz
            vz = 0.5d0*( p%of(id)%loc%vel%y%tmp(i,j,k+1)-p%of(id)%loc%vel%y%tmp(i,j,k)+p%of(id)%loc%vel%y%tmp(i,j-1,k+1)-p%of(id)%loc%vel%y%tmp(i,j-1,k) )/p%glb%dz
    
            phix = 0.25d0*( p%of(id)%loc%phi%now(i+1,j,k+1)-p%of(id)%loc%phi%now(i-1,j,k+1) + p%of(id)%loc%phi%now(i+1,j,k)-p%of(id)%loc%phi%now(i-1,j,k) )/p%glb%dx
            phiy = 0.25d0*( p%of(id)%loc%phi%now(i,j+1,k+1)-p%of(id)%loc%phi%now(i,j-1,k+1) + p%of(id)%loc%phi%now(i,j+1,k)-p%of(id)%loc%phi%now(i,j-1,k) )/p%glb%dy
            phiz = ( p%of(id)%loc%phi%now(i,j,k+1)-p%of(id)%loc%phi%now(i,j,k) )/p%glb%dz
            
            dif_x = mu/rho*xx/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phix*(uz+wx)/(rho*p%glb%re)
            dif_y = mu/rho*yy/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiy*(vz+wy)/(rho*p%glb%re)
            dif_z = mu/rho*zz/p%glb%re + (1.0d0-p%glb%mu_12)*delta*phiz*2.0d0*wz/(rho*p%glb%re)
            
            p%of(id)%loc%velsrc%z%tmp(i,j,k) = ( p%of(id)%loc%velsrc%z%tmp(i,j,k) + &
            & dif_x + dif_y + dif_z + &
            & p%glb%gz*p%glb%btn_g / p%glb%fr &
            & - p%glb%btn_sf*curv*delta*phiz / (p%glb%we*rho) ) /2.0d0
            
        end do
        end do
        end do 
        
    enddo
    !$omp end parallel do

end subroutine

subroutine ns_ab_diff_source_uccd


end subroutine
