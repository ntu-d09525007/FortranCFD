subroutine check_memory
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k

	call check_velocity
	call check_phi
	call check_pressure
	
end subroutine

subroutine check_velocity
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k

	do id = 0, p%glb%threads-2
		do i = 1, 3
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do k = p%of(id)%loc%ks, p%of(id)%loc%ke
		
			 if( abs(p%of(id)%loc%vel%x%now(p%of(id)%loc%ie+i,j,k)-p%of(id+1)%loc%vel%x%now(p%of(id)%loc%ie+i,j,k))>1.0d-12 )then
				 write(*,*)"vel.x",id,i,j,k
			 endif
			
			 if( abs(p%of(id)%loc%vel%y%now(p%of(id)%loc%ie+i,j,k)-p%of(id+1)%loc%vel%y%now(p%of(id)%loc%ie+i,j,k))>1.0d-12 )then
				 write(*,*)"vel.y",id,i,j,k
			 endif
			
			 if( abs(p%of(id)%loc%vel%z%now(p%of(id)%loc%ie+i,j,k)-p%of(id+1)%loc%vel%z%now(p%of(id)%loc%ie+i,j,k))>1.0d-12 )then
				 write(*,*)"vel.z",id,i,j,k
			 endif
			
		end do
		end do
		end do
	end do

	do id = 1, p%glb%threads-1
		do i = 1, 3
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do k = p%of(id)%loc%ks, p%of(id)%loc%ke
		
			 if( abs(p%of(id)%loc%vel%x%now(p%of(id)%loc%is-i,j,k)-p%of(id-1)%loc%vel%x%now(p%of(id)%loc%is-i,j,k))>1.0d-12 )then
				 write(*,*)"vel.x",id,i,j,k
			 endif
			
			 if( abs(p%of(id)%loc%vel%y%now(p%of(id)%loc%is-i,j,k)-p%of(id-1)%loc%vel%y%now(p%of(id)%loc%is-i,j,k))>1.0d-12 )then
				 write(*,*)"vel.y",id,i,j,k
			 endif
			
			 if( abs(p%of(id)%loc%vel%z%now(p%of(id)%loc%is-i,j,k)-p%of(id-1)%loc%vel%z%now(p%of(id)%loc%is-i,j,k))>1.0d-12 )then
				 write(*,*)"vel.z",id,i,j,k
			 endif
			
		end do
		end do
		end do
	end do
	
end subroutine

subroutine check_pressure
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k

	do id = 0, p%glb%threads-2
		do i = 1, 3
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do k = p%of(id)%loc%ks, p%of(id)%loc%ke
		
			 if( abs(p%of(id)%loc%p%now(p%of(id)%loc%ie+i,j,k)-p%of(id+1)%loc%p%now(p%of(id)%loc%ie+i,j,k))>1.0d-12 )then
				 write(*,*)"pressure",id,i,j,k
			 endif
			
		end do
		end do
		end do
	end do

	do id = 1, p%glb%threads-1
		do i = 1, 3
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do k = p%of(id)%loc%ks, p%of(id)%loc%ke
		
			 if( abs(p%of(id)%loc%p%now(p%of(id)%loc%is-i,j,k)-p%of(id-1)%loc%p%now(p%of(id)%loc%is-i,j,k))>1.0d-12 )then
				 write(*,*)"pressure",id,i,j,k
			 endif
			
		end do
		end do
		end do
	end do
	
end subroutine

subroutine check_phi
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k

	do id = 0, p%glb%threads-2
		do i = 1, 3
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do k = p%of(id)%loc%ks, p%of(id)%loc%ke
		
			 if( abs(p%of(id)%loc%phi%now(p%of(id)%loc%ie+i,j,k)-p%of(id+1)%loc%phi%now(p%of(id)%loc%ie+i,j,k))>1.0d-12 )then
				 write(*,*)"phi",id,i,j,k
			 endif
			 if( abs(p%of(id)%loc%rho%now(p%of(id)%loc%ie+i,j,k)-p%of(id+1)%loc%rho%now(p%of(id)%loc%ie+i,j,k))>1.0d-12 )then
				 write(*,*)"rho",id,i,j,k
			 endif
			 if( abs(p%of(id)%loc%mu%now(p%of(id)%loc%ie+i,j,k)-p%of(id+1)%loc%mu%now(p%of(id)%loc%ie+i,j,k))>1.0d-12 )then
				 write(*,*)"mu",id,i,j,k
			 endif
			
		end do
		end do
		end do
	end do

	do id = 1, p%glb%threads-1
		do i = 1, 3
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do k = p%of(id)%loc%ks, p%of(id)%loc%ke
		
			 if( abs(p%of(id)%loc%phi%now(p%of(id)%loc%is-i,j,k)-p%of(id-1)%loc%phi%now(p%of(id)%loc%is-i,j,k))>1.0d-12 )then
				 write(*,*)"phi",id,i,j,k
			 endif
			 if( abs(p%of(id)%loc%rho%now(p%of(id)%loc%is-i,j,k)-p%of(id-1)%loc%rho%now(p%of(id)%loc%is-i,j,k))>1.0d-12 )then
				 write(*,*)"rho",id,i,j,k
			 endif
			 if( abs(p%of(id)%loc%mu%now(p%of(id)%loc%is-i,j,k)-p%of(id-1)%loc%mu%now(p%of(id)%loc%is-i,j,k))>1.0d-12 )then
				 write(*,*)"mu",id,i,j,k
			 endif
			
		end do
		end do
		end do
	end do
	
end subroutine

subroutine find_momentum()
use all
!$ use omp_lib
implicit none
integer :: id,i,j,k
real(8) :: momx,momy,momz,dv
	
	momx=0.0;momy=0.0;momz=0.0
	dv = p%glb%dx*p%glb%dy*p%glb%dz
	
	!$omp parallel private(id,i,j,k), num_threads(p%glb%threads), reduction(+:momx,momy,momz)
		
		id=0
		!$ id = omp_get_thread_num()
		
		do k = p%of(id)%loc%ks, p%of(id)%loc%ke
		do j = p%of(id)%loc%js, p%of(id)%loc%je
		do i = p%of(id)%loc%is, p%of(id)%loc%ie
		
			momx = momx + p%of(id)%loc%heavy%now(i,j,k)*p%of(id)%loc%vel%x%now(i,j,k)*dv
			momy = momy + p%of(id)%loc%heavy%now(i,j,k)*p%of(id)%loc%vel%y%now(i,j,k)*dv
			momz = momz + p%of(id)%loc%heavy%now(i,j,k)*p%of(id)%loc%vel%z%now(i,j,k)*dv
			
		end do
		end do
		end do
		
	!$omp end parallel
	
	momx = momx / p%glb%vol
	momy = momy / p%glb%vol
	momz = momz / p%glb%vol
	
	write(*,*)''
	write(*,'("X momentum  :",F12.5)')momx
	write(*,'("Y momentum  :",F12.5)')momy
	write(*,'("Z momentum  :",F12.5)')momz

end subroutine