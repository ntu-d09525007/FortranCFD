subroutine plot
USE Lib_VTK_IO
use all
implicit none 
integer :: E_IO, id
integer :: nn, nx1,nx2,ny1,ny2,nz1,nz2
integer :: i,j,k
real(8),dimension(1:p%glb%node_x,1:p%glb%node_y,1:p%glb%node_z) :: x,y,z
character(3) :: name

if ( abs(p%glb%time - p%glb%pid * p%glb%t2p) > p%glb%dt ) return

!$omp parallel do collapse(3)
do k = 1, p%glb%node_z
do j = 1, p%glb%node_y
do i = 1, p%glb%node_x
    x(i,j,k) = p%glb%x(i)
    y(i,j,k) = p%glb%y(j)
    z(i,j,k) = p%glb%z(k)
enddo
enddo
enddo
!$omp end parallel do

write(name,'(I3.3)')p%glb%pid

E_IO = VTK_INI_XML(output_format='raw', filename='./out/'//trim(p%glb%name)//'_'//name//'.vts', &
                   mesh_topology='StructuredGrid', nx1=1, nx2=p%glb%node_x, ny1=1, ny2=p%glb%node_y, nz1=1, nz2=p%glb%node_z) 

do id = 0, p%glb%threads-1
nx1 = p%of(id)%loc%is; nx2=p%of(id)%loc%ie
ny1 = p%of(id)%loc%js; ny2=p%of(id)%loc%je
nz1 = p%of(id)%loc%ks; nz2=p%of(id)%loc%ke

nn=(nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1)

E_IO = VTK_GEO_XML(nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2, nz1=nz1, nz2=nz2, NN=nn, &
                   X=x(nx1:nx2,ny1:ny2,nz1:nz2),Y=y(nx1:nx2,ny1:ny2,nz1:nz2),Z=z(nx1:nx2,ny1:ny2,nz1:nz2))

E_IO = VTK_DAT_XML(var_location = 'node', var_block_action = 'open')
E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'concentration', var = p%of(id)%loc%c%now(nx1:nx2,ny1:ny2,nz1:nz2) )
E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'salinity', var = p%of(id)%loc%s%now(nx1:nx2,ny1:ny2,nz1:nz2) )
E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'solid', var = p%of(id)%loc%solid%now(nx1:nx2,ny1:ny2,nz1:nz2) )
E_IO = VTK_DAT_XML(var_location = 'node', var_block_action = 'close')
E_IO = VTK_GEO_XML()
enddo

E_IO = VTK_END_XML()

p%glb%pid = p%glb%pid + 1
call p%sync

end subroutine