subroutine ns_velbc()
use all
!$ use omp_lib
implicit none
integer :: id,i,j

!$omp parallel do 
do id = 0, p%glb%threads-1
    call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now)
enddo
!$omp end parallel do

call pt%vel%sync

end subroutine