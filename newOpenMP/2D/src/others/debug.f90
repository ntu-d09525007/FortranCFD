subroutine find_momentum()
use all
implicit none
integer :: i,j
real(8) :: momx,momy,dv
    
    momx=0.0;momy=0.0
    dv = p%glb%dx * p%glb%dy
    
    !$omp parallel do collapse(2), reduction(+:momx,momy)
    do j = p%loc%js, p%loc%je
    do i = p%loc%is, p%loc%ie
    
        momx = momx + p%loc%heavy%now(i,j)*p%loc%vel%x%now(i,j)*dv
        momy = momy + p%loc%heavy%now(i,j)*p%loc%vel%y%now(i,j)*dv
        
    end do
    end do
    !$omp end parallel do
    
    momx = momx / p%glb%vol
    momy = momy / p%glb%vol
    
    write(*,*)''
    write(*,'("X momentum  :",F12.5)')momx
    write(*,'("Y momentum  :",F12.5)')momy

end subroutine