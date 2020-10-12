subroutine find_momentum()
use all
!$ use omp_lib
implicit none
integer :: id,i,j
real(8) :: momx,momy,dv
    
    momx=0.0;momy=0.0
    dv = p%glb%dx*p%glb%dy
    
    !$omp parallel do private(i,j), reduction(+:momx,momy)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
            momx = momx + p%of(id)%loc%heavy%now(i,j)*p%of(id)%loc%vel%x%now(i,j)*dv
            momy = momy + p%of(id)%loc%heavy%now(i,j)*p%of(id)%loc%vel%y%now(i,j)*dv
  
        end do
        end do

    enddo
    !$omp end parallel do
    
    momx = momx / p%glb%vol
    momy = momy / p%glb%vol
    
    write(*,*)''
    write(*,'("X momentum  :",F12.5)')momx
    write(*,'("Y momentum  :",F12.5)')momy

end subroutine