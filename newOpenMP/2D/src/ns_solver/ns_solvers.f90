subroutine ns_solver
use all
implicit none
integer(8) :: cpustart, cpuend

    call system_clock(cpustart)
    
    call ns_init
    
    call ns_ab_solver_SOR
    !call ns_split_solver
    
    call ns_check_convergence_div
    call node_vel
    
    call system_clock(cpuend)
    p%glb%ns = p%glb%ns + real(cpuend-cpustart,kind=8)/real(p%glb%cpurate,kind=8)
    
end subroutine

subroutine ns_init
use all
implicit none

call rho_mu
call find_gradient(p%loc%mu_ten,p%loc%mu%old)

if( p%glb%btn_sf > 0 )then
    call curv
endif
    
end subroutine

subroutine ns_linearize
use all
!$ use omp_lib
implicit none
integer :: i,j

!$omp parallel do collapse(2)  
do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
    p%loc%vel%x%tmp(i,j) = p%loc%vel%x%now(i,j) 
    p%loc%vel%y%tmp(i,j) = p%loc%vel%y%now(i,j) 
end do
end do   
!$omp end parallel do

   
call find_stag_vel(  p%loc%tdata%x%s1, p%loc%tdata%y%s1, &
                     p%loc%vel%x%tmp, p%loc%vel%y%tmp )


end subroutine

subroutine ns_check_convergence_div
use all
!$ use omp_lib
implicit none
integer :: i,j
real(8) :: div, sumdiv
real(8) :: ux,vy

div=0.0d0
sumdiv=0.0d0

!$omp parallel do collapse(2), private(ux,vy), reduction(max:div), reduction(+:sumdiv) 
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie

    ux = (p%loc%vel%x%now(i,j)-p%loc%vel%x%now(i-1,j))/p%glb%dx
    vy = (p%loc%vel%y%now(i,j)-p%loc%vel%y%now(i,j-1))/p%glb%dy

    div = max( div, abs(ux+vy) ) 
    sumdiv = sumdiv + abs(ux+vy)
                                    
end do
end do
!$omp end parallel do
  
p%glb%vel_sdiv = sumdiv / (p%glb%node_x*p%glb%node_y)
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
id=0

!$omp parallel do collapse(2), reduction(max:linf), reduction(+:l2f)
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie

    linf = max(linf,abs(p%loc%vel%x%now(i,j)-p%loc%vel%x%tmp(i,j)))
    linf = max(linf,abs(p%loc%vel%y%now(i,j)-p%loc%vel%y%tmp(i,j)))
  
    l2f = l2f + (p%loc%vel%x%now(i,j)-p%loc%vel%x%tmp(i,j))**2.0d0
    l2f = l2f + (p%loc%vel%y%now(i,j)-p%loc%vel%y%tmp(i,j))**2.0d0
                              
end do
end do
!$omp end parallel do

l2f = dsqrt( l2f / (2.0d0*p%glb%node_x*p%glb%node_y) )

p%glb%ns_linf = linf
p%glb%ns_l2f = l2f
    
end subroutine

subroutine node_vel()
use all
implicit none
integer :: i,j

!$omp parallel do collapse(2)
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie
    p%loc%nvel%x%now(i,j) = 0.5d0 * ( p%loc%vel%x%now(i-1,j) + p%loc%vel%x%now(i,j) )
    p%loc%nvel%y%now(i,j) = 0.5d0 * ( p%loc%vel%y%now(i,j-1) + p%loc%vel%y%now(i,j) )
end do
end do
!$omp end parallel do

call nvelbc(p%loc%nvel%x%now,p%loc%nvel%y%now)

end subroutine
