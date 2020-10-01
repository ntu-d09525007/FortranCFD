subroutine ns_ab_diff_sec
implicit none

    call sec_part1
    !call sec_part2
    call sec_part3

end subroutine

subroutine ns_ab_diff_uccd
use all
implicit none
integer :: i,j

call bc(p%loc%vel%x%old)
call bc(p%loc%vel%y%old)

!$omp parallel do         
do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc

    call p%loc%ccdsolvers%x%solve("uccd",p%loc%vel%x%old(:,j),&
        p%loc%vel_ten%xx(:,j),p%loc%vel_ten%xxx(:,j),p%loc%vel%x%tmp(:,j))
    call p%loc%ccdsolvers%x%solve("uccd",p%loc%vel%y%old(:,j),&
        p%loc%vel_ten%yx(:,j),p%loc%vel_ten%yxx(:,j),p%loc%tdata%x%s1(:,j))
    call p%loc%ccdsolvers%x%solve("uccd",p%loc%vel%z%old(:,j),&
        p%loc%vel_ten%zx(:,j),p%loc%vel_ten%zxx(:,j),p%loc%tdata%x%s2(:,j))

end do
!$omp end parallel do

!$omp parallel do
do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc

    call p%loc%ccdsolvers%y%solve("uccd",p%loc%vel%x%old(i,:),&
        p%loc%vel_ten%xy(i,:),p%loc%vel_ten%xyy(i,:),p%loc%tdata%y%s1(i,:))
    call p%loc%ccdsolvers%y%solve("uccd",p%loc%vel%y%old(i,:),&
        p%loc%vel_ten%yy(i,:),p%loc%vel_ten%yyy(i,:), p%loc%vel%y%tmp(i,:))
    call p%loc%ccdsolvers%y%solve("uccd",p%loc%vel%z%old(i,:),&
        p%loc%vel_ten%zy(i,:),p%loc%vel_ten%zyy(i,:),p%loc%tdata%y%s2(i,:))

end do
!$omp end parallel do

call velbc(p%loc%vel%x%old,p%loc%vel%y%old)

call uccd_part1
!call uccd_part2
call uccd_part3

end subroutine

subroutine sec_part1
use all 
implicit none
integer :: i,j
real(8) :: rho,mu,xx,yy,dif_x,dif_y

    !$omp parallel do collapse(2), private(rho,mu,xx,yy,dif_x,dif_y)
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
    
        rho = 0.5d0*(p%loc%rho%old(i,j)+p%loc%rho%old(i+1,j))
         mu = 0.5d0*(p%loc%mu%old(i,j)+p%loc%mu%old(i+1,j))

        xx = (p%loc%vel%x%old(I+1,J)-2.0d0*p%loc%vel%x%old(I,J)+p%loc%vel%x%old(I-1,J))/p%glb%dx**2.0d0
        yy = (p%loc%vel%x%old(I,J+1)-2.0d0*p%loc%vel%x%old(I,J)+p%loc%vel%x%old(I,J-1))/p%glb%dy**2.0d0

        dif_x = mu/rho*xx/p%glb%re 
        dif_y = mu/rho*yy/p%glb%re

        p%loc%velsrc%x%tmp(i,j) = dif_x + dif_y 

        !----------------------------------------------------

        rho = 0.5d0*(p%loc%rho%old(i,j)+p%loc%rho%old(i,j+1))
         mu = 0.5d0*(p%loc%mu%old(i,j)+p%loc%mu%old(i,j+1))

        xx = (p%loc%vel%y%old(I+1,J)-2.0d0*p%loc%vel%y%old(I,J)+p%loc%vel%y%old(I-1,J))/p%glb%dx**2.0d0
        yy = (p%loc%vel%y%old(I,J+1)-2.0d0*p%loc%vel%y%old(I,J)+p%loc%vel%y%old(I,J-1))/p%glb%dy**2.0d0

        dif_x = mu/rho*xx/p%glb%re 
        dif_y = mu/rho*yy/p%glb%re

        p%loc%velsrc%y%tmp(i,j) = dif_x + dif_y

    end do
    end do 
    !$omp end parallel do

end subroutine

subroutine sec_part2
use all 
implicit none
integer :: i,j
real(8) :: rho,ux,vx,uy,vy
real(8) :: dif_x,dif_y
real(8) :: phix,phiy

    !$omp parallel do collapse(2), private(rho,ux,vx,uy,vy,dif_x,dif_y,phix,phiy)
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie

        rho = 0.5d0*(p%loc%rho%old(i,j)+p%loc%rho%old(i+1,j))

        ux = 0.5d0*( p%loc%vel%x%old(i+1,j)-p%loc%vel%x%old(i-1,j) )/p%glb%dx
        uy = 0.5d0*( p%loc%vel%x%old(i,j+1)-p%loc%vel%x%old(i,j-1) )/p%glb%dy

        vx = 0.5d0*( p%loc%vel%y%old(i+1,j)-p%loc%vel%y%old(i,j)+p%loc%vel%y%old(i+1,j-1)-p%loc%vel%y%old(i,j-1) )/p%glb%dx

        phix = (p%loc%phi%old(i+1,j)-p%loc%phi%old(i,j))/p%glb%dx
        phiy = 0.25d0*(p%loc%phi%old(i+1,j+1)-p%loc%phi%old(i+1,j-1)+p%loc%phi%old(i,j+1)-p%loc%phi%old(i,j-1))/p%glb%dy
    
        dif_x = phix*2.0d0*ux/(rho*p%glb%re)
        dif_y = phiy*(uy+vx)/(rho*p%glb%re)

        p%loc%velsrc%x%tmp(i,j) = p%loc%velsrc%x%tmp(i,j) + dif_x + dif_y 

        !-------------------------------------------------------------

        rho = 0.5d0*(p%loc%rho%old(i,j)+p%loc%rho%old(i,j+1))

        vx = 0.5d0*( p%loc%vel%y%old(i+1,j)-p%loc%vel%y%old(i-1,j) )/p%glb%dx
        vy = 0.5d0*( p%loc%vel%y%old(i,j+1)-p%loc%vel%y%old(i,j-1) )/p%glb%dy

        uy = 0.5d0*( p%loc%vel%x%old(i,j+1)-p%loc%vel%x%old(i,j)+p%loc%vel%x%old(i-1,j+1)-p%loc%vel%x%old(i-1,j) )/p%glb%dy

        phix = 0.25d0*(p%loc%phi%old(i+1,j)-p%loc%phi%old(i-1,j)+p%loc%phi%old(i+1,j+1)-p%loc%phi%old(i-1,j+1))/p%glb%dx
        phiy = ( p%loc%phi%old(i,j+1)-p%loc%phi%old(i,j) )/p%glb%dy

        dif_x = phix*(uy+vx)/(rho*p%glb%re)
        dif_y = phiy*2.0d0*vy/(rho*p%glb%re)

        p%loc%velsrc%y%tmp(i,j) = p%loc%velsrc%y%tmp(i,j) + dif_x + dif_y

    end do
    end do 
    !$omp end parallel do

end subroutine

subroutine sec_part3
use all 
implicit none
integer :: i,j
real(8) :: rho,delta,curv,phix,phiy

!$omp parallel do collapse(2), private(rho,delta,curv,phix,phiy)
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie

    rho = 0.5d0*(p%loc%rho%old(i,j)+p%loc%rho%old(i+1,j))
    delta = 0.5d0*(p%loc%delta%old(i,j)+p%loc%delta%old(i+1,j))
    curv = (p%loc%normals%curv%old(i,j)+p%loc%normals%curv%old(i+1,j))/2.0d0

    phix = 0.5d0*( p%loc%normals%x%old(i,j)+p%loc%normals%x%old(i+1,j) )
    p%loc%velsrc%x%tmp(i,j) = p%loc%velsrc%x%tmp(i,j) &
        + p%glb%gx*p%glb%btn_g / p%glb%fr - p%glb%btn_sf*curv*delta*phix / (p%glb%we*rho) 

    !--------------------------------------------------

    rho = 0.5d0*(p%loc%rho%old(i,j)+p%loc%rho%old(i,j+1))
    delta = 0.5d0*(p%loc%delta%old(i,j)+p%loc%delta%old(i,j+1))
    curv = (p%loc%normals%curv%old(i,j)+p%loc%normals%curv%old(i,j+1))/2.0d0

    phiy = 0.5d0*( p%loc%normals%y%old(i,j)+p%loc%normals%y%old(i,j+1) )
    p%loc%velsrc%y%tmp(i,j) = p%loc%velsrc%y%tmp(i,j) &
        + p%glb%gy*p%glb%btn_g / p%glb%fr - p%glb%btn_sf*curv*delta*phiy / (p%glb%we*rho) 


enddo
enddo
!$omp end parallel do

end subroutine

subroutine uccd_part1
use all 
implicit none
integer :: i,j
real(8) :: rho,mu,xx,yy,dif_x,dif_y

    !$omp parallel do collapse(2), private(rho,mu,xx,yy,dif_x,dif_y)
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
    
        rho = 0.5d0*(p%loc%rho%old(i,j)+p%loc%rho%old(i+1,j))
         mu = 0.5d0*(p%loc%mu%old(i,j)+p%loc%mu%old(i+1,j))

        xx = p%loc%vel_ten%xxx(i,j)
        yy = p%loc%vel_ten%xyy(i,j)
        
        dif_x = mu/rho*xx/p%glb%re 
        dif_y = mu/rho*yy/p%glb%re

        p%loc%velsrc%x%tmp(i,j) = dif_x + dif_y 

        !----------------------------------------------------

        rho = 0.5d0*(p%loc%rho%old(i,j)+p%loc%rho%old(i,j+1))
         mu = 0.5d0*(p%loc%mu%old(i,j)+p%loc%mu%old(i,j+1))

        xx = p%loc%vel_ten%yxx(i,j)
        yy = p%loc%vel_ten%yyy(i,j)
        
        dif_x = mu/rho*xx/p%glb%re 
        dif_y = mu/rho*yy/p%glb%re

        p%loc%velsrc%y%tmp(i,j) = dif_x + dif_y 
           
    end do
    end do 
    !$omp end parallel do

end subroutine

subroutine uccd_part2
use all 
implicit none
integer :: i,j
real(8) :: rho,ux,vx,uy,vy
real(8) :: dif_x,dif_y
real(8) :: phix,phiy

    !$omp parallel do collapse(2), private(rho,ux,vx,uy,vy,dif_x,dif_y,phix,phiy)
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie

        rho = 0.5d0*(p%loc%rho%old(i,j)+p%loc%rho%old(i+1,j))

        ux = p%loc%vel_ten%xx(i,j)
        uy = p%loc%vel_ten%xy(i,j)

        vx = 0.25d0*( p%loc%vel_ten%yx(i,j)+p%loc%vel_ten%yx(i+1,j)&
                    &+p%loc%vel_ten%yx(i,j-1)+p%loc%vel_ten%yx(i+1,j-1) )

        phix = 0.5d0*( p%loc%mu_ten%x(i,j)+p%loc%mu_ten%x(i+1,j) )
        phiy = 0.5d0*( p%loc%mu_ten%y(i,j)+p%loc%mu_ten%y(i+1,j) )

        dif_x = phix*2.0d0*ux/(rho*p%glb%re)
        dif_y = phiy*(uy+vx)/(rho*p%glb%re)

        p%loc%velsrc%x%tmp(i,j) = p%loc%velsrc%x%tmp(i,j) + dif_x + dif_y 

        !-------------------------------------------------------------

        rho = 0.5d0*(p%loc%rho%old(i,j)+p%loc%rho%old(i,j+1))

        vx = p%loc%vel_ten%yx(i,j)
        vy = p%loc%vel_ten%yy(i,j)

        uy = 0.25d0*( p%loc%vel_ten%xy(i,j)+p%loc%vel_ten%xy(i,j+1)&
                    &+p%loc%vel_ten%xy(i-1,j)+p%loc%vel_ten%xy(i-1,j+1) )

        phix = 0.5d0*( p%loc%mu_ten%x(i,j)+p%loc%mu_ten%x(i,j+1) )
        phiy = 0.5d0*( p%loc%mu_ten%y(i,j)+p%loc%mu_ten%y(i,j+1) )

        dif_x = phix*(uy+vx)/(rho*p%glb%re)
        dif_y = phiy*2.0d0*vy/(rho*p%glb%re)

        p%loc%velsrc%y%tmp(i,j) = p%loc%velsrc%y%tmp(i,j) + dif_x + dif_y 

    end do
    end do
    !$omp end parallel do

end subroutine

subroutine uccd_part3
use all 
implicit none
integer :: i,j
real(8) :: rho,delta,curv,phix,phiy

!$omp parallel do collapse(2), private(rho,delta,curv,phix,phiy)
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie

    rho = 0.5d0*(p%loc%rho%old(i,j)+p%loc%rho%old(i+1,j))
    delta = 0.5d0*(p%loc%delta%old(i,j)+p%loc%delta%old(i+1,j))
    curv = (p%loc%normals%curv%old(i,j)+p%loc%normals%curv%old(i+1,j))/2.0d0

    phix = 0.5d0*( p%loc%normals%x%old(i,j)+p%loc%normals%x%old(i+1,j) )
    p%loc%velsrc%x%tmp(i,j) = p%loc%velsrc%x%tmp(i,j) &
        + p%glb%gx*p%glb%btn_g / p%glb%fr - p%glb%btn_sf*curv*delta*phix / (p%glb%we*rho) 

    !--------------------------------------------------

    rho = 0.5d0*(p%loc%rho%old(i,j)+p%loc%rho%old(i,j+1))
    delta = 0.5d0*(p%loc%delta%old(i,j)+p%loc%delta%old(i,j+1))
    curv = (p%loc%normals%curv%old(i,j)+p%loc%normals%curv%old(i,j+1))/2.0d0

    phiy = 0.5d0*( p%loc%normals%y%old(i,j)+p%loc%normals%y%old(i,j+1) )
    p%loc%velsrc%y%tmp(i,j) = p%loc%velsrc%y%tmp(i,j) &
        + p%glb%gy*p%glb%btn_g / p%glb%fr - p%glb%btn_sf*curv*delta*phiy / (p%glb%we*rho) 

enddo
enddo
!$omp end parallel do

end subroutine
