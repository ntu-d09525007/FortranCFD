<<<<<<< HEAD
subroutine ns_velbc()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: tau

tau = p%glb%rho_12

!$omp parallel do private(i,j,k)
do id = 0, p%glb%threads-1

    call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now,p%of(id)%loc%vel%z%now)

    if(p%of(id)%loc%idz==p%glb%grid_z-1)then
        do k = 1, p%glb%ghc
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            p%of(id)%loc%vel%z%now(i,j,p%of(id)%loc%ke+k) = 0.0d0
            p%of(id)%loc%vel%y%now(i,j,p%of(id)%loc%ke+k) = p%of(id)%loc%vel%y%now(i,j,p%of(id)%loc%ke-1+k)
            p%of(id)%loc%vel%x%now(i,j,p%of(id)%loc%ke+k) = p%of(id)%loc%vel%x%now(i,j,p%of(id)%loc%ke-1+k) + tau*p%glb%dz
        enddo
        enddo
        enddo
    endif

enddo
!$omp end parallel do

call pt%vel%sync

end subroutine

subroutine ns_nvel()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: tau

call p%node_vel

tau = p%glb%rho_12

!$omp parallel do private(i,j,k)
do id = 0, p%glb%threads-1

    if(p%of(id)%loc%idz==p%glb%grid_z-1)then
        do k = 1, p%glb%ghc
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            p%of(id)%loc%nvel%z%now(i,j,p%of(id)%loc%ke+k) = 0.0d0
            p%of(id)%loc%nvel%y%now(i,j,p%of(id)%loc%ke+k) = p%of(id)%loc%nvel%y%now(i,j,p%of(id)%loc%ke-1+k)
            p%of(id)%loc%nvel%x%now(i,j,p%of(id)%loc%ke+k) = p%of(id)%loc%nvel%x%now(i,j,p%of(id)%loc%ke-1+k) + tau*p%glb%dz
        enddo
        enddo
        enddo
    endif

enddo
!$omp end parallel do

call pt%nvel%sync

end subroutine
=======
subroutine ns_velbc()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k

!$omp parallel do private(i,j,k)
do id = 0, p%glb%threads-1

    call p%of(id)%velbc(p%of(id)%loc%vel%x%now,p%of(id)%loc%vel%y%now,p%of(id)%loc%vel%z%now)

    if(p%of(id)%loc%idz==p%glb%grid_z-1)then
        do k = 1, p%glb%ghc
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            p%of(id)%loc%vel%z%now(i,j,p%of(id)%loc%ke+k) = 0.0d0
            p%of(id)%loc%vel%y%now(i,j,p%of(id)%loc%ke+k) = p%of(id)%loc%vel%y%now(i,j,p%of(id)%loc%ke)
            p%of(id)%loc%vel%x%now(i,j,p%of(id)%loc%ke+k) = p%of(id)%loc%vel%x%now(i,j,p%of(id)%loc%ke-1+k) + p%glb%dz
        enddo
        enddo
        enddo
    endif

enddo
!$omp end parallel do

call pt%vel%sync

end subroutine

subroutine ns_nvel()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k

call p%node_vel

!$omp parallel do private(i,j,k)
do id = 0, p%glb%threads-1

    if(p%of(id)%loc%idz==p%glb%grid_z-1)then
        do k = 1, p%glb%ghc
        do j = p%of(id)%loc%js, p%of(id)%loc%je
        do i = p%of(id)%loc%is, p%of(id)%loc%ie
            p%of(id)%loc%nvel%z%now(i,j,p%of(id)%loc%ke+k) = 0.0d0
            p%of(id)%loc%nvel%y%now(i,j,p%of(id)%loc%ke+k) = p%of(id)%loc%nvel%y%now(i,j,p%of(id)%loc%ke)
            p%of(id)%loc%nvel%x%now(i,j,p%of(id)%loc%ke+k) = p%of(id)%loc%nvel%x%now(i,j,p%of(id)%loc%ke-1+k) + p%glb%dz
        enddo
        enddo
        enddo
    endif

enddo
!$omp end parallel do

call pt%nvel%sync

end subroutine
>>>>>>> fc06937481dd4e7c913c59d569e8c59aa8fdf971
