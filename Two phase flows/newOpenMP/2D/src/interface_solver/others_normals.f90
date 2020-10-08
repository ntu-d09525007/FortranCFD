subroutine surface_norms()
use all
implicit none
integer :: id,I,J 

!$omp parallel do
do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc

    call p%loc%ccdsolvers%x%solve("ccd",p%loc%phi%now(:,j),p%loc%normals%x%now(:,j),p%loc%normals%xx%now(:,j))
                                                            
enddo
!$omp end parallel do

!$omp parallel do 
do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc

    call p%loc%ccdsolvers%y%solve("ccd",p%loc%phi%now(i,:),p%loc%normals%y%now(i,:),p%loc%normals%yy%now(i,:))
                                                                    
enddo
!$omp end parallel do  
        
!===========================================

!$omp parallel do
do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc

    call p%loc%ccdsolvers%y%solve("ccd",p%loc%normals%x%now(i,:),p%loc%normals%xy%now(i,:),p%loc%normals%curv%now(i,:))
                                                                    
enddo
!$omp end parallel do

call bc(p%loc%normals%x%now)
call bc(p%loc%normals%xx%now)

call bc(p%loc%normals%y%now)
call bc(p%loc%normals%yy%now)
 
call bc(p%loc%normals%xy%now)

end subroutine

subroutine surface_norms2()
use all
implicit none
integer :: id,I,J

!$omp parallel do collapse(2)
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie

    p%loc%normals%x%now(i,j) = 0.5d0*( p%loc%phi%now(i+1,j) - p%loc%phi%now(i-1,j) )/p%glb%dx
    p%loc%normals%y%now(i,j) = 0.5d0*( p%loc%phi%now(i,j+1) - p%loc%phi%now(i,j-1) )/p%glb%dy

    p%loc%normals%xx%now(i,j) = ( p%loc%phi%now(i+1,j) - 2.0d0*p%loc%phi%now(i,j) - p%loc%phi%now(i-1,j) )/p%glb%dx**2.0d0
    p%loc%normals%yy%now(i,j) = ( p%loc%phi%now(i,j+1) - 2.0d0*p%loc%phi%now(i,j) - p%loc%phi%now(i,j-1) )/p%glb%dy**2.0d0

    p%loc%normals%xy%now(i,j) = ( p%loc%phi%now(i+1,j+1)+p%loc%phi%now(i-1,j-1) &
                                & - p%loc%phi%now(i-1,j+1)-p%loc%phi%now(i+1,j-1) )/(4.0d0*p%glb%dx*p%glb%dy)

end do
end do
!$omp end parallel do

!===========================================

call bc(p%loc%normals%x%now)
call bc(p%loc%normals%xx%now)

call bc(p%loc%normals%y%now)
call bc(p%loc%normals%yy%now)
 
call bc(p%loc%normals%xy%now)
    
end subroutine

subroutine curv()
use all
implicit none
integer :: id,i,j,k
real(8) :: fx,fxx,fy,fyy,fxy

call surface_norms

!$omp parallel do collapse(2), private(fx,fxx,fy,fyy,fxy)
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie

    fx  = p%loc%normals%x%now(i,j)
    fy  = p%loc%normals%y%now(i,j)
    
    fxx = p%loc%normals%xx%now(i,j)
    fyy = p%loc%normals%yy%now(i,j)
    
    fxy = p%loc%normals%xy%now(i,j)
    
    p%loc%normals%curv%now(i,j) = ( fyy*fx**2.0d0+fxx*fy**2.0d0-2.0d0*fxy*fx*fy ) / (fx**2.0d0+fy**2.0d0+1.0d-12)**1.5d0 
    
end do
end do
!$omp end parallel do

call bc(p%loc%normals%curv%now)
    
end subroutine