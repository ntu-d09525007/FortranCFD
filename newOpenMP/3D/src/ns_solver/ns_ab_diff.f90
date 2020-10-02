subroutine ns_ab_diff_sec
implicit none

    call sec_part1
    call sec_part2
    call sec_part3

end subroutine

subroutine ns_ab_diff_uccd
use all
implicit none
integer :: i,j,k

call bc(p%loc%vel%x%old)
call bc(p%loc%vel%y%old)
call bc(p%loc%vel%z%old)

!$omp parallel do collapse(2)        
do k = p%loc%ks-p%glb%ghc, p%loc%ke+p%glb%ghc
do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc

    call p%loc%ccdsolvers%x%solve("uccd",p%loc%vel%x%old(:,j,k),&
        p%loc%vel_ten%xx(:,j,k),p%loc%vel_ten%xxx(:,j,k),p%loc%vel%x%tmp(:,j,k))
    call p%loc%ccdsolvers%x%solve("uccd",p%loc%vel%y%old(:,j,k),&
        p%loc%vel_ten%yx(:,j,k),p%loc%vel_ten%yxx(:,j,k),p%loc%tdata%x%s1(:,j,k))
    call p%loc%ccdsolvers%x%solve("uccd",p%loc%vel%z%old(:,j,k),&
        p%loc%vel_ten%zx(:,j,k),p%loc%vel_ten%zxx(:,j,k),p%loc%tdata%x%s2(:,j,k))

end do
end do
!$omp end parallel do

!$omp parallel do collapse(2)   
do k = p%loc%ks-p%glb%ghc, p%loc%ke+p%glb%ghc
do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc

    call p%loc%ccdsolvers%y%solve("uccd",p%loc%vel%x%old(i,:,k),&
        p%loc%vel_ten%xy(i,:,k),p%loc%vel_ten%xyy(i,:,k),p%loc%tdata%y%s1(i,:,k))
    call p%loc%ccdsolvers%y%solve("uccd",p%loc%vel%y%old(i,:,k),&
        p%loc%vel_ten%yy(i,:,k),p%loc%vel_ten%yyy(i,:,k), p%loc%vel%y%tmp(i,:,k))
    call p%loc%ccdsolvers%y%solve("uccd",p%loc%vel%z%old(i,:,k),&
        p%loc%vel_ten%zy(i,:,k),p%loc%vel_ten%zyy(i,:,k),p%loc%tdata%y%s2(i,:,k))

end do
end do
!$omp end parallel do

!$omp parallel do collapse(2)
do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc

    call p%loc%ccdsolvers%z%solve("uccd",p%loc%vel%x%old(i,j,:),&
        p%loc%vel_ten%xz(i,j,:),p%loc%vel_ten%xzz(i,j,:),p%loc%tdata%z%s1(i,j,:))
    call p%loc%ccdsolvers%z%solve("uccd",p%loc%vel%y%old(i,j,:),&
        p%loc%vel_ten%yz(i,j,:),p%loc%vel_ten%yzz(i,j,:),p%loc%tdata%z%s2(i,j,:))
    call p%loc%ccdsolvers%z%solve("uccd",p%loc%vel%z%old(i,j,:),&
        p%loc%vel_ten%zz(i,j,:),p%loc%vel_ten%zzz(i,j,:), p%loc%vel%z%tmp(i,j,:))
                    
end do
end do
!$omp end parallel do

call velbc(p%loc%vel%x%old,p%loc%vel%y%old,p%loc%vel%z%old)

call uccd_part1
call uccd_part2
call uccd_part3

end subroutine

subroutine sec_part1
use all 
implicit none
integer :: i,j,k
real(8) :: rho,mu,xx,yy,zz,dif_x,dif_y,dif_z

    !$omp parallel do collapse(3), private(rho,mu,xx,yy,zz,dif_x,dif_y,dif_z)
    do k = p%loc%ks, p%loc%ke
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
    
        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i+1,j,k))
         mu = 0.5d0*(p%loc%mu%old(i,j,k)+p%loc%mu%old(i+1,j,k))

        xx = (p%loc%vel%x%old(I+1,J,K)-2.0d0*p%loc%vel%x%old(I,J,K)+p%loc%vel%x%old(I-1,J,K))/p%glb%dx**2.0d0
        yy = (p%loc%vel%x%old(I,J+1,K)-2.0d0*p%loc%vel%x%old(I,J,K)+p%loc%vel%x%old(I,J-1,K))/p%glb%dy**2.0d0
        zz = (p%loc%vel%x%old(I,J,K+1)-2.0d0*p%loc%vel%x%old(I,J,K)+p%loc%vel%x%old(I,J,K-1))/p%glb%dz**2.0d0
        
        dif_x = mu/rho*xx/p%glb%re 
        dif_y = mu/rho*yy/p%glb%re
        dif_z = mu/rho*zz/p%glb%re 

        p%loc%velsrc%x%tmp(i,j,k) = dif_x + dif_y + dif_z 

        !----------------------------------------------------

        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i,j+1,k))
         mu = 0.5d0*(p%loc%mu%old(i,j,k)+p%loc%mu%old(i,j+1,k))

        xx = (p%loc%vel%y%old(I+1,J,K)-2.0d0*p%loc%vel%y%old(I,J,K)+p%loc%vel%y%old(I-1,J,K))/p%glb%dx**2.0d0
        yy = (p%loc%vel%y%old(I,J+1,K)-2.0d0*p%loc%vel%y%old(I,J,K)+p%loc%vel%y%old(I,J-1,K))/p%glb%dy**2.0d0
        zz = (p%loc%vel%y%old(I,J,K+1)-2.0d0*p%loc%vel%y%old(I,J,K)+p%loc%vel%y%old(I,J,K-1))/p%glb%dz**2.0d0
        
        dif_x = mu/rho*xx/p%glb%re 
        dif_y = mu/rho*yy/p%glb%re
        dif_z = mu/rho*zz/p%glb%re 

        p%loc%velsrc%y%tmp(i,j,k) = dif_x + dif_y + dif_z 

        !----------------------------------------------------

        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i,j,k+1))
         mu = 0.5d0*(p%loc%mu%old(i,j,k)+p%loc%mu%old(i,j,k+1))

        xx = (p%loc%vel%z%old(I+1,J,K)-2.0d0*p%loc%vel%z%old(I,J,K)+p%loc%vel%z%old(I-1,J,K))/p%glb%dx**2.0d0
        yy = (p%loc%vel%z%old(I,J+1,K)-2.0d0*p%loc%vel%z%old(I,J,K)+p%loc%vel%z%old(I,J-1,K))/p%glb%dy**2.0d0
        zz = (p%loc%vel%z%old(I,J,K+1)-2.0d0*p%loc%vel%z%old(I,J,K)+p%loc%vel%z%old(I,J,K-1))/p%glb%dz**2.0d0
        
        dif_x = mu/rho*xx/p%glb%re 
        dif_y = mu/rho*yy/p%glb%re
        dif_z = mu/rho*zz/p%glb%re 

        p%loc%velsrc%z%tmp(i,j,k) = dif_x + dif_y + dif_z 
           
    end do
    end do
    end do 
    !$omp end parallel do

end subroutine

subroutine sec_part2
use all 
implicit none
integer :: i,j,k
real(8) :: rho,ux,vx,wx,uy,vy,wy,uz,vz,wz
real(8) :: dif_x,dif_y,dif_z
real(8) :: phix,phiy,phiz

    !$omp parallel do collapse(3), private(rho,ux,vx,wx,uy,vy,wy,uz,vz,wz,dif_x,dif_y,dif_z,phix,phiy,phiz)
    do k = p%loc%ks, p%loc%ke
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie

        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i+1,j,k))

        ux = 0.5d0*( p%loc%vel%x%old(i+1,j,k)-p%loc%vel%x%old(i-1,j,k) )/p%glb%dx
        uy = 0.5d0*( p%loc%vel%x%old(i,j+1,k)-p%loc%vel%x%old(i,j-1,k) )/p%glb%dy
        uz = 0.5d0*( p%loc%vel%x%old(i,j,k+1)-p%loc%vel%x%old(i,j,k-1) )/p%glb%dz

        vx = 0.5d0*( p%loc%vel%y%old(i+1,j,k)-p%loc%vel%y%old(i,j,k)+p%loc%vel%y%old(i+1,j-1,k)-p%loc%vel%y%old(i,j-1,k) )/p%glb%dx
        wx = 0.5d0*( p%loc%vel%z%old(i+1,j,k)-p%loc%vel%z%old(i,j,k)+p%loc%vel%z%old(i+1,j,k-1)-p%loc%vel%z%old(i,j,k-1) )/p%glb%dx

        phix = (p%loc%phi%old(i+1,j,k)-p%loc%phi%old(i,j,k))/p%glb%dx
        phiy = 0.25d0*(p%loc%phi%old(i+1,j+1,k)-p%loc%phi%old(i+1,j-1,k)+p%loc%phi%old(i,j+1,k)-p%loc%phi%old(i,j-1,k))/p%glb%dy
        phiz = 0.25d0*(p%loc%phi%old(i+1,j,k+1)-p%loc%phi%old(i+1,j,k-1)+p%loc%phi%old(i,j,k+1)-p%loc%phi%old(i,j,k-1))/p%glb%dz
        
        dif_x = phix*2.0d0*ux/(rho*p%glb%re)
        dif_y = phiy*(uy+vx)/(rho*p%glb%re)
        dif_z = phiz*(uz+wx)/(rho*p%glb%re)

        p%loc%velsrc%x%tmp(i,j,k) = p%loc%velsrc%x%tmp(i,j,k) + dif_x + dif_y + dif_z 

        !-------------------------------------------------------------

        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i,j+1,k))

        vx = 0.5d0*( p%loc%vel%y%old(i+1,j,k)-p%loc%vel%y%old(i-1,j,k) )/p%glb%dx
        vy = 0.5d0*( p%loc%vel%y%old(i,j+1,k)-p%loc%vel%y%old(i,j-1,k) )/p%glb%dy
        vz = 0.5d0*( p%loc%vel%y%old(i,j,k+1)-p%loc%vel%y%old(i,j,k-1) )/p%glb%dz

        uy = 0.5d0*( p%loc%vel%x%old(i,j+1,k)-p%loc%vel%x%old(i,j,k)+p%loc%vel%x%old(i-1,j+1,k)-p%loc%vel%x%old(i-1,j,k) )/p%glb%dy
        wy = 0.5d0*( p%loc%vel%z%old(i,j+1,k)-p%loc%vel%z%old(i,j,k)+p%loc%vel%z%old(i,j+1,k-1)-p%loc%vel%z%old(i,j,k-1) )/p%glb%dy

        phix = 0.25d0*(p%loc%phi%old(i+1,j,k)-p%loc%phi%old(i-1,j,k)+p%loc%phi%old(i+1,j+1,k)-p%loc%phi%old(i-1,j+1,k))/p%glb%dx
        phiy = ( p%loc%phi%old(i,j+1,k)-p%loc%phi%old(i,j,k) )/p%glb%dy
        phiz = 0.25d0*(p%loc%phi%old(i,j+1,k+1)-p%loc%phi%old(i,j+1,k-1)+p%loc%phi%old(i,j,k+1)-p%loc%phi%old(i,j,k-1))/p%glb%dz

        dif_x = phix*(uy+vx)/(rho*p%glb%re)
        dif_y = phiy*2.0d0*vy/(rho*p%glb%re)
        dif_z = phiz*(wy+vz)/(rho*p%glb%re)

        p%loc%velsrc%y%tmp(i,j,k) = p%loc%velsrc%y%tmp(i,j,k) + dif_x + dif_y + dif_z 

        !-------------------------------------------------------------

        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i,j,k+1))
        
        wx = 0.5d0*( p%loc%vel%z%old(i+1,j,k)-p%loc%vel%z%old(i-1,j,k) )/p%glb%dx
        wy = 0.5d0*( p%loc%vel%z%old(i,j+1,k)-p%loc%vel%z%old(i,j-1,k) )/p%glb%dy
        wz = 0.5d0*( p%loc%vel%z%old(i,j,k+1)-p%loc%vel%z%old(i,j,k-1) )/p%glb%dz

        uz = 0.5d0*( p%loc%vel%x%old(i,j,k+1)-p%loc%vel%x%old(i,j,k)+p%loc%vel%x%old(i-1,j,k+1)-p%loc%vel%x%old(i-1,j,k) )/p%glb%dz
        vz = 0.5d0*( p%loc%vel%y%old(i,j,k+1)-p%loc%vel%y%old(i,j,k)+p%loc%vel%y%old(i,j-1,k+1)-p%loc%vel%y%old(i,j-1,k) )/p%glb%dz

        phix = 0.25d0*( p%loc%phi%old(i+1,j,k+1)-p%loc%phi%old(i-1,j,k+1) + p%loc%phi%old(i+1,j,k)-p%loc%phi%old(i-1,j,k) )/p%glb%dx
        phiy = 0.25d0*( p%loc%phi%old(i,j+1,k+1)-p%loc%phi%old(i,j-1,k+1) + p%loc%phi%old(i,j+1,k)-p%loc%phi%old(i,j-1,k) )/p%glb%dy
        phiz = ( p%loc%phi%old(i,j,k+1)-p%loc%phi%old(i,j,k) )/p%glb%dz

        dif_x = phix*(uz+wx)/(rho*p%glb%re)
        dif_y = phiy*(vz+wy)/(rho*p%glb%re)
        dif_z = phiz*2.0d0*wz/(rho*p%glb%re)

        p%loc%velsrc%z%tmp(i,j,k) = p%loc%velsrc%z%tmp(i,j,k) + dif_x + dif_y + dif_z 

    end do
    end do
    end do 
    !$omp end parallel do

end subroutine

subroutine sec_part3
use all 
implicit none
integer :: i,j,k
real(8) :: rho,delta,curv,phix,phiy,phiz

!$omp parallel do collapse(3), private(rho,delta,curv,phix,phiy,phiz)
do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie

    rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i+1,j,k))
    delta = 0.5d0*(p%loc%delta%old(i,j,k)+p%loc%delta%old(i+1,j,k))
    curv = (p%loc%normals%curv%old(i,j,k)+p%loc%normals%curv%old(i+1,j,k))/2.0d0

    phix = 0.5d0*( p%loc%normals%x%old(i,j,k)+p%loc%normals%x%old(i+1,j,k) )
    p%loc%velsrc%x%tmp(i,j,k) = p%loc%velsrc%x%tmp(i,j,k) &
        + p%glb%gx*p%glb%btn_g / p%glb%fr - p%glb%btn_sf*curv*delta*phix / (p%glb%we*rho) 

    !--------------------------------------------------

    rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i,j+1,k))
    delta = 0.5d0*(p%loc%delta%old(i,j,k)+p%loc%delta%old(i,j+1,k))
    curv = (p%loc%normals%curv%old(i,j,k)+p%loc%normals%curv%old(i,j+1,k))/2.0d0

    phiy = 0.5d0*( p%loc%normals%y%old(i,j,k)+p%loc%normals%y%old(i,j+1,k) )
    p%loc%velsrc%y%tmp(i,j,k) = p%loc%velsrc%y%tmp(i,j,k) &
        + p%glb%gy*p%glb%btn_g / p%glb%fr - p%glb%btn_sf*curv*delta*phiy / (p%glb%we*rho) 

    !--------------------------------------------------

    rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i,j,k+1))
    delta = 0.5d0*(p%loc%delta%old(i,j,k)+p%loc%delta%old(i,j,k+1))
    curv = (p%loc%normals%curv%old(i,j,k)+p%loc%normals%curv%old(i,j,k+1))/2.0d0

    phiz = 0.5d0*( p%loc%normals%z%old(i,j,k)+p%loc%normals%z%old(i,j,k+1) )
    p%loc%velsrc%z%tmp(i,j,k) = p%loc%velsrc%z%tmp(i,j,k) &
        + p%glb%gz*p%glb%btn_g / p%glb%fr - p%glb%btn_sf*curv*delta*phiz / (p%glb%we*rho) 

enddo
enddo
enddo
!$omp end parallel do

end subroutine

subroutine uccd_part1
use all 
implicit none
integer :: i,j,k
real(8) :: rho,mu,xx,yy,zz,dif_x,dif_y,dif_z

    !$omp parallel do collapse(3), private(rho,mu,xx,yy,zz,dif_x,dif_y,dif_z)
    do k = p%loc%ks, p%loc%ke
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
    
        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i+1,j,k))
         mu = 0.5d0*(p%loc%mu%old(i,j,k)+p%loc%mu%old(i+1,j,k))

        xx = p%loc%vel_ten%xxx(i,j,k);yy = p%loc%vel_ten%xyy(i,j,k);zz = p%loc%vel_ten%xzz(i,j,k)
        
        dif_x = mu/rho*xx/p%glb%re 
        dif_y = mu/rho*yy/p%glb%re
        dif_z = mu/rho*zz/p%glb%re 

        p%loc%velsrc%x%tmp(i,j,k) = dif_x + dif_y + dif_z 

        !----------------------------------------------------

        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i,j+1,k))
         mu = 0.5d0*(p%loc%mu%old(i,j,k)+p%loc%mu%old(i,j+1,k))

        xx = p%loc%vel_ten%yxx(i,j,k);yy = p%loc%vel_ten%yyy(i,j,k);zz = p%loc%vel_ten%yzz(i,j,k)
        
        dif_x = mu/rho*xx/p%glb%re 
        dif_y = mu/rho*yy/p%glb%re
        dif_z = mu/rho*zz/p%glb%re 

        p%loc%velsrc%y%tmp(i,j,k) = dif_x + dif_y + dif_z 

        !----------------------------------------------------

        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i,j,k+1))
         mu = 0.5d0*(p%loc%mu%old(i,j,k)+p%loc%mu%old(i,j,k+1))

        xx = p%loc%vel_ten%zxx(i,j,k);yy = p%loc%vel_ten%zyy(i,j,k);zz = p%loc%vel_ten%zzz(i,j,k)
        
        dif_x = mu/rho*xx/p%glb%re 
        dif_y = mu/rho*yy/p%glb%re
        dif_z = mu/rho*zz/p%glb%re 

        p%loc%velsrc%z%tmp(i,j,k) = dif_x + dif_y + dif_z 
           
    end do
    end do
    end do 
    !$omp end parallel do

end subroutine

subroutine uccd_part2
use all 
implicit none
integer :: i,j,k
real(8) :: rho,ux,vx,wx,uy,vy,wy,uz,vz,wz
real(8) :: dif_x,dif_y,dif_z
real(8) :: phix,phiy,phiz

    !$omp parallel do collapse(3), private(rho,ux,vx,wx,uy,vy,wy,uz,vz,wz,dif_x,dif_y,dif_z,phix,phiy,phiz)
    do k = p%loc%ks, p%loc%ke
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie

        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i+1,j,k))

        ux = p%loc%vel_ten%xx(i,j,k)
        uy = p%loc%vel_ten%xy(i,j,k)
        uz = p%loc%vel_ten%xz(i,j,k)

        vx = 0.25d0*( p%loc%vel_ten%yx(i,j,k)+p%loc%vel_ten%yx(i+1,j,k)&
                    &+p%loc%vel_ten%yx(i,j-1,k)+p%loc%vel_ten%yx(i+1,j-1,k) )
        wx = 0.25d0*( p%loc%vel_ten%zx(i,j,k)+p%loc%vel_ten%zx(i+1,j,k)&
                    &+p%loc%vel_ten%zx(i,j,k-1)+p%loc%vel_ten%zx(i+1,j,k-1) )

        phix = 0.5d0*( p%loc%mu_ten%x(i,j,k)+p%loc%mu_ten%x(i+1,j,k) )
        phiy = 0.5d0*( p%loc%mu_ten%y(i,j,k)+p%loc%mu_ten%y(i+1,j,k) )
        phiz = 0.5d0*( p%loc%mu_ten%z(i,j,k)+p%loc%mu_ten%z(i+1,j,k) )
        
        dif_x = phix*2.0d0*ux/(rho*p%glb%re)
        dif_y = phiy*(uy+vx)/(rho*p%glb%re)
        dif_z = phiz*(uz+wx)/(rho*p%glb%re)

        p%loc%velsrc%x%tmp(i,j,k) = p%loc%velsrc%x%tmp(i,j,k) + dif_x + dif_y + dif_z 

        !-------------------------------------------------------------

        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i,j+1,k))

        vx = p%loc%vel_ten%yx(i,j,k)
        vy = p%loc%vel_ten%yy(i,j,k)
        vz = p%loc%vel_ten%yz(i,j,k)

        uy = 0.25d0*( p%loc%vel_ten%xy(i,j,k)+p%loc%vel_ten%xy(i,j+1,k)&
                    &+p%loc%vel_ten%xy(i-1,j,k)+p%loc%vel_ten%xy(i-1,j+1,k) )
        wy = 0.25d0*( p%loc%vel_ten%zy(i,j,k)+p%loc%vel_ten%zy(i,j+1,k)&
                    &+p%loc%vel_ten%zy(i,j,k-1)+p%loc%vel_ten%zy(i,j+1,k-1) )

        phix = 0.5d0*( p%loc%mu_ten%x(i,j,k)+p%loc%mu_ten%x(i,j+1,k) )
        phiy = 0.5d0*( p%loc%mu_ten%y(i,j,k)+p%loc%mu_ten%y(i,j+1,k) )
        phiz = 0.5d0*( p%loc%mu_ten%z(i,j,k)+p%loc%mu_ten%z(i,j+1,k) )

        dif_x = phix*(uy+vx)/(rho*p%glb%re)
        dif_y = phiy*2.0d0*vy/(rho*p%glb%re)
        dif_z = phiz*(wy+vz)/(rho*p%glb%re)

        p%loc%velsrc%y%tmp(i,j,k) = p%loc%velsrc%y%tmp(i,j,k) + dif_x + dif_y + dif_z 

        !-------------------------------------------------------------

        rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i,j,k+1))

        wx = p%loc%vel_ten%zx(i,j,k)
        wy = p%loc%vel_ten%zy(i,j,k)
        wz = p%loc%vel_ten%zz(i,j,k)

        uz = 0.25d0*( p%loc%vel_ten%xz(i,j,k)+p%loc%vel_ten%xz(i,j,k+1)&
                    &+p%loc%vel_ten%xz(i-1,j,k)+p%loc%vel_ten%xz(i-1,j,k+1) )
        vz = 0.25d0*( p%loc%vel_ten%yz(i,j,k)+p%loc%vel_ten%yz(i,j,k+1)&
                    &+p%loc%vel_ten%yz(i,j-1,k)+p%loc%vel_ten%yz(i,j-1,k+1) )

        phix = 0.5d0*( p%loc%mu_ten%x(i,j,k)+p%loc%mu_ten%x(i,j,k+1) )
        phiy = 0.5d0*( p%loc%mu_ten%y(i,j,k)+p%loc%mu_ten%y(i,j,k+1) )
        phiz = 0.5d0*( p%loc%mu_ten%z(i,j,k)+p%loc%mu_ten%z(i,j,k+1) )

        dif_x = phix*(uz+wx)/(rho*p%glb%re)
        dif_y = phiy*(vz+wy)/(rho*p%glb%re)
        dif_z = phiz*2.0d0*wz/(rho*p%glb%re)

        p%loc%velsrc%z%tmp(i,j,k) = p%loc%velsrc%z%tmp(i,j,k) + dif_x + dif_y + dif_z 

    end do
    end do
    end do 
    !$omp end parallel do

end subroutine

subroutine uccd_part3
use all 
implicit none
integer :: i,j,k
real(8) :: rho,delta,curv,phix,phiy,phiz

!$omp parallel do collapse(3), private(rho,delta,curv,phix,phiy,phiz)
do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie

    rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i+1,j,k))
    delta = 0.5d0*(p%loc%delta%old(i,j,k)+p%loc%delta%old(i+1,j,k))
    curv = (p%loc%normals%curv%old(i,j,k)+p%loc%normals%curv%old(i+1,j,k))/2.0d0

    phix = 0.5d0*( p%loc%normals%x%old(i,j,k)+p%loc%normals%x%old(i+1,j,k) )
    p%loc%velsrc%x%tmp(i,j,k) = p%loc%velsrc%x%tmp(i,j,k) &
        + p%glb%gx*p%glb%btn_g / p%glb%fr - p%glb%btn_sf*curv*delta*phix / (p%glb%we*rho) 

    !--------------------------------------------------

    rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i,j+1,k))
    delta = 0.5d0*(p%loc%delta%old(i,j,k)+p%loc%delta%old(i,j+1,k))
    curv = (p%loc%normals%curv%old(i,j,k)+p%loc%normals%curv%old(i,j+1,k))/2.0d0

    phiy = 0.5d0*( p%loc%normals%y%old(i,j,k)+p%loc%normals%y%old(i,j+1,k) )
    p%loc%velsrc%y%tmp(i,j,k) = p%loc%velsrc%y%tmp(i,j,k) &
        + p%glb%gy*p%glb%btn_g / p%glb%fr - p%glb%btn_sf*curv*delta*phiy / (p%glb%we*rho) 

    !--------------------------------------------------

    rho = 0.5d0*(p%loc%rho%old(i,j,k)+p%loc%rho%old(i,j,k+1))
    delta = 0.5d0*(p%loc%delta%old(i,j,k)+p%loc%delta%old(i,j,k+1))
    curv = (p%loc%normals%curv%old(i,j,k)+p%loc%normals%curv%old(i,j,k+1))/2.0d0

    phiz = 0.5d0*( p%loc%normals%z%old(i,j,k)+p%loc%normals%z%old(i,j,k+1) )
    p%loc%velsrc%z%tmp(i,j,k) = p%loc%velsrc%z%tmp(i,j,k) &
        + p%glb%gz*p%glb%btn_g / p%glb%fr - p%glb%btn_sf*curv*delta*phiz / (p%glb%we*rho) 

enddo
enddo
enddo
!$omp end parallel do

end subroutine
