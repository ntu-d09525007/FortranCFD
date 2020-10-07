subroutine job_find_tensor(p,tens,u,v)
implicit none
class(job) :: p
type(tensor) :: tens
real(8), dimension(p%loc%is-p%glb%ghc:p%loc%ie+p%glb%ghc,&
                  &p%loc%js-p%glb%ghc:p%loc%je+p%glb%ghc) :: u,v
integer :: i, j

do j = p%loc%js, p%loc%je
    call p%loc%ccdsolvers%x%solve("ccd",u(:,j),tens%xx(:,j),tens%xxx(:,j))
    call p%loc%ccdsolvers%x%solve("ccd",v(:,j),tens%yx(:,j),tens%yxx(:,j))
enddo

do i = p%loc%is, p%loc%ie
    call p%loc%ccdsolvers%y%solve("ccd",u(i,:),tens%xy(i,:),tens%xyy(i,:))
    call p%loc%ccdsolvers%y%solve("ccd",v(i,:),tens%yy(i,:),tens%yyy(i,:))
enddo


end subroutine

subroutine job_find_gradient(p,tens,phi)
implicit none
class(job) :: p
type(tensor) :: tens
real(8), dimension(p%loc%is-p%glb%ghc:p%loc%ie+p%glb%ghc,&
                  &p%loc%js-p%glb%ghc:p%loc%je+p%glb%ghc) :: phi
integer :: i,j

do j = p%loc%js, p%loc%je
    call p%loc%ccdsolvers%x%solve("ccd",phi(:,j),tens%x(:,j),tens%tmp(:,j) )
enddo

do i = p%loc%is, p%loc%ie
    call p%loc%ccdsolvers%y%solve("ccd",phi(i,:),tens%y(i,:),tens%tmp(i,:) )
enddo


end subroutine