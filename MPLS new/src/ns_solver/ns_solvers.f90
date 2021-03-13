subroutine ns_solver
use all
implicit none
integer(8) :: cpustart, cpuend

call system_clock(cpustart)

call ns_init

call ns_ab_solver
call ns_split_solver

call ns_check_convergence_div

call p%node_vel
call pt%nvel%sync

call system_clock(cpuend)
p%glb%ns = p%glb%ns + real(cpuend-cpustart,kind=8)/real(p%glb%cpurate,kind=8)
    
end subroutine

subroutine ns_init
use all
!$ use omp_lib
implicit none

call p%rho_mu
call p%curv
call pt%normals%sync
    
end subroutine

subroutine ns_linearize
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k

!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
    do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        p%of(id)%loc%vel%x%tmp(i,j) = p%of(id)%loc%vel%x%now(i,j) 
        p%of(id)%loc%vel%y%tmp(i,j) = p%of(id)%loc%vel%y%now(i,j) 
    end do
    end do
  
enddo       
!$omp end parallel do
        
end subroutine

subroutine ns_check_convergence_div
use all
!$ use omp_lib
implicit none
integer :: id,i,j,n
real(8) :: div, sumdiv
real(8) :: ux,vy

div=0.0d0
sumdiv=0.0d0
n=0
!$omp parallel do private(i,j,ux,vy), reduction(max:div), reduction(+:sumdiv,n)
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie

        n=n+1
        ux = (p%of(id)%loc%vel%x%now(i,j)-p%of(id)%loc%vel%x%now(i-1,j))/p%glb%dx
        vy = (p%of(id)%loc%vel%y%now(i,j)-p%of(id)%loc%vel%y%now(i,j-1))/p%glb%dy
   
        div = max( div, abs(ux+vy) ) 
        sumdiv = sumdiv + abs(ux+vy)

    end do
    end do
  
enddo   
!$omp end parallel do
  
p%glb%vel_sdiv = sumdiv / (n)
p%glb%vel_div = div
    
end subroutine

subroutine ns_check_convergence_vel
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: linf, l2f

linf=0.0d0
l2f=0.0d0

!$omp parallel do private(i,j), reduction(max:linf), reduction(+:l2f)
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
    
        linf = max(linf,abs(p%of(id)%loc%vel%x%now(i,j)-p%of(id)%loc%vel%x%tmp(i,j)))
        linf = max(linf,abs(p%of(id)%loc%vel%y%now(i,j)-p%of(id)%loc%vel%y%tmp(i,j)))
   
        l2f = l2f + (p%of(id)%loc%vel%x%now(i,j)-p%of(id)%loc%vel%x%tmp(i,j))**2.0d0
        l2f = l2f + (p%of(id)%loc%vel%y%now(i,j)-p%of(id)%loc%vel%y%tmp(i,j))**2.0d0

    end do
    end do
  
enddo   
!$omp end parallel do

l2f = dsqrt( l2f / (3.0d0*p%glb%node_x*p%glb%node_y) )

p%glb%ns_linf = linf
p%glb%ns_l2f = l2f
    
end subroutine


