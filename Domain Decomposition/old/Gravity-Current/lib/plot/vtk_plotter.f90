module vtk_plotter
use tree
!$ use omp_lib
implicit none

type plotter_data
integer :: is, ie
real(8), dimension(:,:,:), pointer :: dat, x, y, z
end type

type plotter_unit
integer :: s, e, ny, nz
type(plotter_data), allocatable :: from(:)
contains
procedure plot_scalar => unit_plot_scalar
procedure clean_scalar => unit_clean_scalar
procedure plot_vector => unit_plot_vector
procedure clean_vector => unit_clean_vector
end type plotter_unit

type plotter
integer :: pthreads
type(plotter_unit), allocatable :: of(:)
type(manager), pointer :: src
contains
procedure init => plotter_init
procedure plot => plotter_plot
end type plotter

contains

subroutine plotter_init(this,src)
implicit none
class(plotter) :: this
type(manager), target :: src
integer :: id, i

	this%pthreads = src%glb%gz
	
	this%src => src
	
	allocate( this%of(0:this%pthreads-1) )
	
	!$omp parallel private(id,i), num_threads(this%pthreads)
	
		id=0
		!$ id = omp_get_thread_num()
		
		this%of(id)%s = int( real(id*src%glb%threads,kind=8) / real(src%glb%pthreads,kind=8) ) 
		this%of(id)%e = int( real((id+1)*src%glb%threads,kind=8) / real(src%glb%pthreads,kind=8) ) -1
		
		allocate( this%of(id)%from( this%of(id)%s:this%of(id)%e ) )
		
		do i = this%of(id)%s, this%of(id)%e
			this%of(id)%from(i)%is = src%of(i)%loc%is
			this%of(id)%from(i)%ie = src%of(i)%loc%ie
			
			this%of(id)%from(i)%dat => null()
			this%of(id)%from(i)%x => null()
			this%of(id)%from(i)%y => null()
			this%of(id)%from(i)%z => null()
		enddo
		
		this%of(id)%ny = src%glb%node_y
		this%of(id)%nz = src%glb%node_z
		
	!$omp end parallel
	
	!do id = 0, this%pthreads-1
	!	write(*,*)id,this%of(id)%s,this%of(id)%e,this%of(id)%from(this%of(id)%s)%is,this%of(id)%from(this%of(id)%e)%ie
	!enddo
	
end subroutine

subroutine plotter_plot(p)
implicit none
class(plotter) :: p
integer :: id,s,e,nx,ny,nz,i,j,k,m
character(6) :: name

	if ( abs(p%src%glb%time - p%src%glb%pid * p%src%glb%t2p) > p%src%glb%dt ) return

	!call p%src%vortex
	
	!$omp parallel private(id,name,s,e,nx,ny,nz,i,j,k,m), num_threads(p%pthreads)
	
		id=0
		!$ id = omp_get_thread_num()

		write(name,'(I2.2,"_",I3.3)')id,p%src%glb%pid
		
		open(unit=777+id,file="./out/"//trim(p%src%glb%name)//'_'//name//".vtk")
		
		s = p%of(id)%s; e = p%of(id)%e; 
		
		nx = p%of(id)%from(e)%ie-p%of(id)%from(s)%is+2
		ny = p%src%glb%node_y+1
		nz = p%src%glb%node_z+1
		
		WRITE(id+777,'(A)')"# vtk DataFile Version 3.0"
		write(id+777,'(A)')"vtk TEST"
		WRITE(id+777,'(A)')"ASCII"
		WRITE(id+777,'(A)')"DATASET STRUCTURED_POINTS"
	
		WRITE(id+777,'(A,3I6)')"DIMENSIONS ",nx,ny,nz
		
		WRITE(id+777,'(A,3ES15.4)')"SPACING ",p%src%glb%dx, p%src%glb%dy, p%src%glb%dz
		
		IF( ID/=0 )then
			WRITE(id+777,'(A,3ES15.4)')"ORIGIN ",p%src%glb%x(p%of(id)%from(s)%is-1),p%src%glb%ystart,p%src%glb%zstart
		ELSE
			WRITE(id+777,'(A,3ES15.4)')"ORIGIN ",p%src%glb%x(p%of(id)%from(s)%is-1)-0.5D0*p%src%glb%DX,p%src%glb%ystart,p%src%glb%zstart
		ENDIF
		
		WRITE(id+777,'(A,I12)')"POINT_DATA ",nx*ny*nz
		
		! ============= Plot a scalar =====================
		do i = s, e
			p%of(id)%from(i)%dat => p%src%of(i)%loc%phi%now
		enddo
		
		call p%of(id)%plot_scalar(id+777,"PHI")
			
		call p%of(id)%clean_scalar
		
		! -------------------------------------------------
		
		if( p%src%glb%method == 3)then
		
		do i = s, e
			p%of(id)%from(i)%dat => p%src%of(i)%loc%vof%now
		enddo
		
		call p%of(id)%plot_scalar(id+777,"VOF")
			
		call p%of(id)%clean_scalar
		
		endif
		
		! -------------------------------------------------
		
		do i = s, e
			p%of(id)%from(i)%dat => p%src%of(i)%loc%q_cri%now
		enddo
		
		call p%of(id)%plot_scalar(id+777,"Q")
			
		call p%of(id)%clean_scalar
		
		! -------------------------------------------------
		
		do i = s, e
			p%of(id)%from(i)%dat => p%src%of(i)%loc%helicity%now
		enddo
		
		call p%of(id)%plot_scalar(id+777,"Helicity")
			
		call p%of(id)%clean_scalar
		
		! ============= Plot a vector =====================
		do i = s, e
			p%of(id)%from(i)%x => p%src%of(i)%loc%nvel%x%now
			p%of(id)%from(i)%y => p%src%of(i)%loc%nvel%y%now
			p%of(id)%from(i)%z => p%src%of(i)%loc%nvel%z%now
		enddo
		
		call p%of(id)%plot_vector(id+777,"Velocity")
		
		call p%of(id)%clean_vector
		! -------------------------------------------------
		do i = s, e
			p%of(id)%from(i)%x => p%src%of(i)%loc%vort%x%now
			p%of(id)%from(i)%y => p%src%of(i)%loc%vort%y%now
			p%of(id)%from(i)%z => p%src%of(i)%loc%vort%z%now
		enddo
		
		call p%of(id)%plot_vector(id+777,"Vorticity")
		
		call p%of(id)%clean_vector
		! -------------------------------------------------
		do i = s, e
			p%of(id)%from(i)%x => p%src%of(i)%loc%lamb%x%now
			p%of(id)%from(i)%y => p%src%of(i)%loc%lamb%y%now
			p%of(id)%from(i)%z => p%src%of(i)%loc%lamb%z%now
		enddo
		! -------------------------------------------------
		call p%of(id)%plot_vector(id+777,"Lamb")
		
		call p%of(id)%clean_vector
		
	!$omp end parallel
	
	p%src%glb%pid = p%src%glb%pid + 1
	call p%src%sync

end subroutine

subroutine unit_plot_scalar(p,fid,name)
implicit none
class(plotter_unit) :: p
integer :: i, j, k, m, fid
character(*) :: name

	write(fid,'(A)')"SCALARS "//NAME//" FLOAT"
	write(fid,'(A)')"LOOKUP_TABLE DEFAULT"
 	
 	do k = 0, p%nz
	do j = 0, p%ny
		
		write(fid,*)p%from(p%s)%dat(p%from(p%s)%is-1,j,k)
		
		do m = p%s, p%e
		do i = p%from(m)%is, p%from(m)%ie
			write(fid,*)p%from(m)%dat(i,j,k)
		enddo
		enddo
		
		!write(fid,*)p%from(p%e)%dat(p%from(p%e)%ie+1,j,k)
		
	enddo
	enddo
		

end subroutine

subroutine unit_clean_scalar(p)
implicit none
class(plotter_unit) :: p
integer :: i

 do i = p%s, p%e
	  nullify(p%from(i)%dat)
 enddo
 
end subroutine

subroutine unit_plot_vector(p,fid,name)
implicit none
class(plotter_unit) :: p
integer :: i, j, k, m, fid
character(*) :: name

	write(fid,'(A)')"VECTORS "//NAME//" FLOAT"
  
    do k = 0, p%nz
	do j = 0, p%ny
	
		write(fid,*)p%from(p%s)%x(p%from(p%s)%is-1,j,k),p%from(p%s)%y(p%from(p%s)%is-1,j,k),p%from(p%s)%z(p%from(p%s)%is-1,j,k)
		
		do m = p%s, p%e
		do i = p%from(m)%is, p%from(m)%ie
			write(fid,*)p%from(m)%x(i,j,k),p%from(m)%y(i,j,k),p%from(m)%z(i,j,k)
		enddo
		enddo
  
		!write(fid,*)p%from(p%e)%x(p%from(p%e)%ie+1,j,k),p%from(p%e)%y(p%from(p%e)%ie+1,j,k),p%from(p%e)%z(p%from(p%e)%ie+1,j,k)
		
	enddo
	enddo

end subroutine

subroutine unit_clean_vector(p)
implicit none
class(plotter_unit) :: p
integer :: i

 do i = p%s, p%e
	  nullify(p%from(i)%x,p%from(i)%y,p%from(i)%z)
 enddo
 
end subroutine


end module vtk_plotter