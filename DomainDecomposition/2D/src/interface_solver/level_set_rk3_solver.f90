subroutine level_set_rk3_solver()
use all
!$ use omp_lib
implicit none
integer :: i,j,id
integer(8) :: cpustart, cpuend

    call system_clock(cpustart)
    
    !---------------------------
    ! RK3 step 1
    !---------------------------
    call level_set_rk3_source

    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1

        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%phi%now(i,j) + p%of(id)%loc%phi%tmp(i,j) * p%glb%dt
        end do 
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%phi%now)
    
    enddo
    !$omp end parallel do
    
    call pt%phi%sync

    !---------------------------
    ! RK3 step 2
    !---------------------------
    call level_set_rk3_source

    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1

        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%phi%now(i,j) + p%of(id)%loc%phi%tmp(i,j) * p%glb%dt
            p%of(id)%loc%phi%now(i,j) = ( p%of(id)%loc%phi%now(i,j) + 3.0d0*p%of(id)%loc%phi%old(i,j) ) / 4.0d0
        end do 
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%phi%now)
    
    enddo
    !$omp end parallel do
    
    call pt%phi%sync

    !---------------------------
    ! RK3 step 3
    !---------------------------
    call level_set_rk3_source

    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1

        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%phi%now(i,j) + p%of(id)%loc%phi%tmp(i,j) * p%glb%dt
            p%of(id)%loc%phi%now(i,j) = ( p%of(id)%loc%phi%now(i,j) + 2.0d0*p%of(id)%loc%phi%old(i,j) ) / 3.0d0
        end do 
        end do
        
        call p%of(id)%bc(0,p%of(id)%loc%phi%now)
    
    enddo
    !$omp end parallel do
    
    call pt%phi%sync    
    call system_clock(cpuend)
    p%glb%ls_adv = p%glb%ls_adv + real(cpuend-cpustart,kind=8) / real( p%glb%cpurate, kind=8 )

end subroutine

subroutine level_set_rk3_source()
implicit none

!call level_set_rk3_source_weno
call level_set_rk3_source_ccd

end subroutine

subroutine level_set_rk3_source_weno
use all
!$ use omp_lib
implicit none
integer :: id,i,j

!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
    
        p%of(id)%loc%tdata%x%s1(i,j) = 0.5d0*(p%of(id)%loc%nvel%x%old(i,j)+abs(p%of(id)%loc%nvel%x%old(i,j)))*p%of(id)%loc%phi%now(i,j)
        p%of(id)%loc%tdata%x%s2(i,j) = 0.5d0*(p%of(id)%loc%nvel%x%old(i,j)-abs(p%of(id)%loc%nvel%x%old(i,j)))*p%of(id)%loc%phi%now(i,j)
    
        p%of(id)%loc%tdata%y%s1(i,j) = 0.5d0*(p%of(id)%loc%nvel%y%old(i,j)+abs(p%of(id)%loc%nvel%y%old(i,j)))*p%of(id)%loc%phi%now(i,j)
        p%of(id)%loc%tdata%y%s2(i,j) = 0.5d0*(p%of(id)%loc%nvel%y%old(i,j)-abs(p%of(id)%loc%nvel%y%old(i,j)))*p%of(id)%loc%phi%now(i,j)

    end do
    end do 
    
    call p%of(id)%bc(0,p%of(id)%loc%tdata%x%s1);call p%of(id)%bc(0,p%of(id)%loc%tdata%x%s2)
    call p%of(id)%bc(0,p%of(id)%loc%tdata%y%s1);call p%of(id)%bc(0,p%of(id)%loc%tdata%y%s2)

enddo
!$omp end parallel do

call pt%tdatax%sync
call pt%tdatay%sync

!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1

    do j = p%of(id)%loc%js, p%of(id)%loc%je
        call wenojs_flux_split(p%of(id)%loc%tdata%x%s2(:,j),p%of(id)%loc%tdata%x%s1(:,j),&
                               p%of(id)%loc%tdata%x%ss2(:,j),p%of(id)%loc%tdata%x%ss1(:,j),&
                               p%of(id)%loc%is,p%of(id)%loc%ie,p%of(id)%glb%ghc)

        ! call crweno_flux_split(p%of(id)%loc%tdata%x%s2(:,j),p%of(id)%loc%tdata%x%s1(:,j),&
        !                        p%of(id)%loc%tdata%x%ss2(:,j),p%of(id)%loc%tdata%x%ss1(:,j),&
        !                        p%of(id)%loc%is,p%of(id)%loc%ie,p%of(id)%glb%ghc)
    end do 

    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        call wenojs_flux_split(p%of(id)%loc%tdata%y%s2(i,:),p%of(id)%loc%tdata%y%s1(i,:),&
                               p%of(id)%loc%tdata%y%ss2(i,:),p%of(id)%loc%tdata%y%ss1(i,:),&
                               p%of(id)%loc%js,p%of(id)%loc%je,p%of(id)%glb%ghc)

        ! call crweno_flux_split(p%of(id)%loc%tdata%y%s2(i,:),p%of(id)%loc%tdata%y%s1(i,:),&
        !                        p%of(id)%loc%tdata%y%ss2(i,:),p%of(id)%loc%tdata%y%ss1(i,:),&
        !                        p%of(id)%loc%js,p%of(id)%loc%je,p%of(id)%glb%ghc)
    end do

    do j = p%of(id)%loc%js-p%glb%ghc+2, p%of(id)%loc%je+p%glb%ghc-3
    do i = p%of(id)%loc%is-p%glb%ghc+2, p%of(id)%loc%ie+p%glb%ghc-3
        p%of(id)%loc%tdata%x%s1(i,j) = p%of(id)%loc%tdata%x%ss1(i,j)
        p%of(id)%loc%tdata%x%s2(i,j) = p%of(id)%loc%tdata%x%ss2(i,j)

        p%of(id)%loc%tdata%y%s1(i,j) = p%of(id)%loc%tdata%y%ss1(i,j)
        p%of(id)%loc%tdata%y%s2(i,j) = p%of(id)%loc%tdata%y%ss2(i,j)
    enddo
    enddo

enddo
!$omp end parallel do

call pt%tdatax%sync
call pt%tdatay%sync

!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        p%of(id)%loc%phi%tmp(i,j) = ( p%of(id)%loc%tdata%x%s1(i-1,j)-p%of(id)%loc%tdata%x%s1(i,j) ) / p%glb%dx + &
                                &   ( p%of(id)%loc%tdata%x%s2(i-1,j)-p%of(id)%loc%tdata%x%s2(i,j) ) / p%glb%dx + &
                                &   ( p%of(id)%loc%tdata%y%s1(i,j-1)-p%of(id)%loc%tdata%y%s1(i,j) ) / p%glb%dy + &
                                &   ( p%of(id)%loc%tdata%y%s2(i,j-1)-p%of(id)%loc%tdata%y%s2(i,j) ) / p%glb%dy 
    enddo
    enddo
enddo
!$omp end parallel do 

end subroutine

subroutine level_set_rk3_source_ccd
use all
!$ use omp_lib
implicit none
integer :: id,i,j

!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1

    do j = p%of(id)%loc%js, p%of(id)%loc%je
        call p%of(id)%loc%ccdsolvers%x%solve("uccd",p%of(id)%loc%phi%now(:,j),&
            &p%of(id)%loc%tdata%x%s1(:,j),p%of(id)%loc%tdata%x%s2(:,j),p%of(id)%loc%nvel%x%old(:,j))
    end do 

    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        call p%of(id)%loc%ccdsolvers%y%solve("uccd",p%of(id)%loc%phi%now(i,:),&
            p%of(id)%loc%tdata%y%s1(i,:),p%of(id)%loc%tdata%y%s2(i,:),p%of(id)%loc%nvel%y%old(i,:))
    end do

    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        p%of(id)%loc%phi%tmp(i,j) = - p%of(id)%loc%nvel%x%old(i,j) * p%of(id)%loc%tdata%x%s1(i,j) &
                                  & - p%of(id)%loc%nvel%y%old(i,j) * p%of(id)%loc%tdata%y%s1(i,j) 
    enddo
    enddo

enddo
!$omp end parallel do 
    
end subroutine

subroutine wenojs_flux_split(f,g,fp,gm,is,ie,ghc)
implicit none
integer :: i, is, ie, ghc
real(8),dimension(is-ghc:ie+ghc) :: f, g, fp, gm
real(8) :: a1,a2,a3,b1,b2,b3,w1,w2,w3,eps

EPS = 1.0D-13

do i = is-ghc+2, ie+ghc-3
    
    b1 = 13.0d0*(g(i-2)-2.0d0*g(i-1)+g(i))**2.0d0 + 3.0d0*(g(i-2)-4.0d0*g(i-1)+3.0d0*g(i))**2.0d0
    b2 = 13.0d0*(g(i-1)-2.0d0*g(i)+g(i+1))**2.0d0 + 3.0d0*(g(i-1)-g(i+1))**2.0d0
    b3 = 13.0d0*(g(i)-2.0d0*g(i+1)+g(i+2))**2.0d0 + 3.0d0*(3.0d0*g(i)-4.0d0*g(i+1)+g(i+2))**2.0d0
    
    a1 = 1.0d0/(EPS+b1)**2.0d0
    a2 = 6.0d0/(EPS+b2)**2.0d0
    a3 = 3.0d0/(EPS+b3)**2.0d0

    !a1 = 1.0d0*(1.0d0+abs(b3-b1)/(EPS+b1))
    !a2 = 6.0d0*(1.0d0+abs(b3-b1)/(EPS+b2))
    !a3 = 3.0d0*(1.0d0+abs(b3-b1)/(EPS+b3))
    
    w1 = a1/(a1+a2+a3)
    w2 = a2/(a1+a2+a3)
    w3 = a3/(a1+a2+a3)
    
    gm(i) = w1/3.0d0*g(i-2) - (7.0d0*w1+w2)/6.0d0*g(i-1) + (11.0d0*w1+5.0d0*w2+2.0d0*w3)/6.0d0*g(i) &
            + (2.0d0*w2+5.0d0*w3)/6.0d0*g(i+1) - w3/6.0d0*g(i+2)
            
end do
    
do i = is-ghc+2, ie+ghc-3
    
    b3 = 13.0d0*(f(i-1)-2.0d0*f(i)  +f(i+1))**2.0d0 + 3.0d0*(f(i-1)-4.0d0*f(i)+3.0d0*f(i+1))**2.0d0
    b2 = 13.0d0*(f(i)  -2.0d0*f(i+1)+f(i+2))**2.0d0 + 3.0d0*(f(i)-f(i+2))**2.0d0
    b1 = 13.0d0*(f(i+1)-2.0d0*f(i+2)+f(i+3))**2.0d0 + 3.0d0*(3.0d0*f(i+1)-4.0d0*f(i+2)+f(i+3))**2.0d0
    
    a1 = 1.0d0/(EPS+b1)**2.0d0
    a2 = 6.0d0/(EPS+b2)**2.0d0
    a3 = 3.0d0/(EPS+b3)**2.0d0

    !a1 = 1.0d0*(1.0d0+abs(b3-b1)/(EPS+b1))
    !a2 = 6.0d0*(1.0d0+abs(b3-b1)/(EPS+b2))
    !a3 = 3.0d0*(1.0d0+abs(b3-b1)/(EPS+b3))
    
    w1 = a1 / (a1+a2+a3)
    w2 = a2 / (a1+a2+a3)
    w3 = a3 / (a1+a2+a3)
    
    fp(i) =  w3*(-f(i-1)+5.0d0*f(i)+2.0d0*f(i+1))/6.0d0 &
            +w2*(2.0d0*f(i)+5.0d0*f(i+1)-f(i+2))/6.0d0 &
            +w1*(11.0d0*f(i+1)-7.0d0*f(i+2)+2.0d0*f(i+3))/6.0d0 
    
    
end do


end subroutine

subroutine crweno_flux_split(f,g,fp,gm,is,ie,ghc)
implicit none
integer :: i, is, ie, ghc
real(8),dimension(is-ghc:ie+ghc) :: f, g, fp, gm
real(8),dimension(is-ghc+3:ie+ghc-4) :: A,B,C,S
real(8) :: a1,a2,a3,b1,b2,b3,w1,w2,w3,eps,c1,c2,c3

EPS = 1.0D-13

c1 = 0.2089141306d0
c2 = 0.4999999998d0
c3 = 0.2910858692d0
    
do i = is-ghc+3, ie+ghc-4
    
    b1 = 13.0d0*(g(i-2)-2.0d0*g(i-1)+g(i))**2.0d0   + 3.0d0*(    g(i-2)-4.0d0*g(i-1)+3.0d0*g(i))**2.0d0
    b2 = 13.0d0*(g(i-1)-2.0d0*g(i)  +g(i+1))**2.0d0 + 3.0d0*(    g(i-1)  -g(i+1))**2.0d0
    b3 = 13.0d0*(g(i)  -2.0d0*g(i+1)+g(i+2))**2.0d0 + 3.0d0*(3.0d0*g(i)-4.0d0*g(i+1)+g(i+2))**2.0d0
    
    a1 = c1*(1.0d0+abs(b3-b1)/(EPS+b1))
    a2 = c2*(1.0d0+abs(b3-b1)/(EPS+b2))
    a3 = c3*(1.0d0+abs(b3-b1)/(EPS+b3))
    
    w1 = a1/(a1+a2+a3)
    w2 = a2/(a1+a2+a3)
    w3 = a3/(a1+a2+a3)  
    
    A(i) = (2.0d0*w1+w2)/3.0d0
    B(i) = (w1+2.0d0*(w2+w3))/3.0d0
    C(i) =  w3/3.0d0
    S(i) = w1/6.0d0*g(i-1) + (5.0d0*(w1+w2)+w3)/6.0d0*g(i) + (w2+5.0d0*w3)/6.0d0*g(i+1)
    
end do  

    S(is-ghc+3) = S(is-ghc+3) - A(is-ghc+3)*gm(is-ghc+2)
    S(ie+ghc-4) = S(ie+ghc-4) - C(ie+ghc-4)*gm(ie+ghc-3)

    call solve_tridiagonal(A,B,C,S,gm(is-ghc+3:ie+ghc-4),is-ghc+3,ie+ghc-4)

do i = is-ghc+3, ie+ghc-4
    
    b3 = 13.0d0*(f(i-1)-2.0d0*f(i)  +f(i+1))**2.0d0 + 3.0d0*(f(i-1)-4.0d0*f(i)+3.0d0*f(i+1))**2.0d0
    b2 = 13.0d0*(f(i)  -2.0d0*f(i+1)+f(i+2))**2.0d0 + 3.0d0*(f(i)-f(i+2))**2.0d0
    b1 = 13.0d0*(f(i+1)-2.0d0*f(i+2)+f(i+3))**2.0d0 + 3.0d0*(3.0d0*f(i+1)-4.0d0*f(i+2)+f(i+3))**2.0d0
    
    a1 = c1*(1.0d0+abs(b3-b1)/(EPS+b1))
    a2 = c2*(1.0d0+abs(b3-b1)/(EPS+b2))
    a3 = c3*(1.0d0+abs(b3-b1)/(EPS+b3))
    
    w1 = a1 / (a1+a2+a3)
    w2 = a2 / (a1+a2+a3)
    w3 = a3 / (a1+a2+a3)
    
    A(i) = (w3)/3.0d0
    B(i) = (w1+2.0d0*(w2+w3))/3.0d0
    C(i) = (w2+2.0d0*w1)/3.0d0
    S(i) = (5.0d0*w3+w2)/6.0d0*f(i) + (w3+5.0d0*(w2+w1))/6.0d0*f(i+1) + w1/6.0d0*f(i+2)
    
end do
    
    S(is-ghc+3) = S(is-ghc+3) - A(is-ghc+3)*fp(is-ghc+2)
    S(ie+ghc-4) = S(ie+ghc-4) - C(ie+ghc-4)*fp(ie+ghc-3)

    call solve_tridiagonal(A,B,C,S,fp(is-ghc+3:ie+ghc-4),is-ghc+3,ie+ghc-4)
    
end subroutine
