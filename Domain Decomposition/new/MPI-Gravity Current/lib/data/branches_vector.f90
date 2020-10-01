subroutine job_find_tensor(p,tens,u,v,w)
implicit none
class(job) :: p
type(tensor) :: tens
real(8), dimension(p%loc%is-p%glb%ghc:p%loc%ie+p%glb%ghc,&
                  &p%loc%js-p%glb%ghc:p%loc%je+p%glb%ghc,&
                  &p%loc%ks-p%glb%ghc:p%loc%ke+p%glb%ghc) :: u,v,w
integer :: i, j, k

do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
    call p%loc%ccd%x%solve_fixed_central(15.0_8/16.0_8,u(:,j,k),tens%xx(:,j,k),tens%xxx(:,j,k))
    call p%loc%ccd%x%solve_fixed_central(15.0_8/16.0_8,v(:,j,k),tens%yx(:,j,k),tens%yxx(:,j,k))
    call p%loc%ccd%x%solve_fixed_central(15.0_8/16.0_8,w(:,j,k),tens%zx(:,j,k),tens%zxx(:,j,k))
enddo
enddo

do k = p%loc%ks, p%loc%ke
do i = p%loc%is, p%loc%ie
    call p%loc%ccd%y%solve_fixed_central(15.0_8/16.0_8,u(i,:,k),tens%xy(i,:,k),tens%xyy(i,:,k))
    call p%loc%ccd%y%solve_fixed_central(15.0_8/16.0_8,v(i,:,k),tens%yy(i,:,k),tens%yyy(i,:,k))
    call p%loc%ccd%y%solve_fixed_central(15.0_8/16.0_8,w(i,:,k),tens%zy(i,:,k),tens%zyy(i,:,k))
enddo
enddo

do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie
    call p%loc%ccd%z%solve_fixed_central(15.0_8/16.0_8,u(i,j,:),tens%xz(i,j,:),tens%xzz(i,j,:))
    call p%loc%ccd%z%solve_fixed_central(15.0_8/16.0_8,v(i,j,:),tens%yz(i,j,:),tens%yzz(i,j,:))
    call p%loc%ccd%z%solve_fixed_central(15.0_8/16.0_8,w(i,j,:),tens%zz(i,j,:),tens%zzz(i,j,:))
enddo
enddo

end subroutine

subroutine job_find_gradient(p,tens,phi)
implicit none
class(job) :: p
type(tensor) :: tens
real(8), dimension(p%loc%is-p%glb%ghc:p%loc%ie+p%glb%ghc,&
                  &p%loc%js-p%glb%ghc:p%loc%je+p%glb%ghc,&
                  &p%loc%ks-p%glb%ghc:p%loc%ke+p%glb%ghc) :: phi
integer :: i,j,k

do k = p%loc%ks, p%loc%ke
do j = p%loc%js, p%loc%je
    call p%loc%ccd%x%solve_fixed_central(15.0_8/16.0_8,phi(:,j,k),tens%x(:,j,k) )
enddo
enddo

do k = p%loc%ks, p%loc%ke
do i = p%loc%is, p%loc%ie
    call p%loc%ccd%y%solve_fixed_central(15.0_8/16.0_8,phi(i,:,k),tens%y(i,:,k) )
enddo
enddo

do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie
    call p%loc%ccd%z%solve_fixed_central(15.0_8/16.0_8,phi(i,j,:),tens%z(i,j,:) )
enddo
enddo

end subroutine