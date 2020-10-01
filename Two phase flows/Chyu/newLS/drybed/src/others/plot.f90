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

!call vortex_dynamics
call ls_funs

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


nx1 = p%loc%is; nx2=p%loc%ie
ny1 = p%loc%js; ny2=p%loc%je
nz1 = p%loc%ks; nz2=p%loc%ke

nn=(nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1)

E_IO = VTK_GEO_XML(nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2, nz1=nz1, nz2=nz2, NN=nn, &
                   X=x(nx1:nx2,ny1:ny2,nz1:nz2),Y=y(nx1:nx2,ny1:ny2,nz1:nz2),Z=z(nx1:nx2,ny1:ny2,nz1:nz2))
                   
E_IO = VTK_DAT_XML(var_location = 'node', var_block_action = 'open')
E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Phi', var = p%loc%phi%now(nx1:nx2,ny1:ny2,nz1:nz2) )
E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'H', var = p%loc%heavy%now(nx1:nx2,ny1:ny2,nz1:nz2) )
E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'S', var = p%loc%sign%now(nx1:nx2,ny1:ny2,nz1:nz2) )
! E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Pressure', var = p%loc%p%now(nx1:nx2,ny1:ny2,nz1:nz2) )
! E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Density', var = p%loc%rho%now(nx1:nx2,ny1:ny2,nz1:nz2) )
! E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Q', var = p%loc%q_cri%now(nx1:nx2,ny1:ny2,nz1:nz2) )
! E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Omega', var = p%loc%omega_cri%now(nx1:nx2,ny1:ny2,nz1:nz2) )
! E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Lamb_div', var = p%loc%lamb_div%now(nx1:nx2,ny1:ny2,nz1:nz2) )

E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Velocity', varX = p%loc%nvel%x%now(nx1:nx2,ny1:ny2,nz1:nz2),&
                                                    &varY = p%loc%nvel%y%now(nx1:nx2,ny1:ny2,nz1:nz2),&
                                                    &varZ = p%loc%nvel%z%now(nx1:nx2,ny1:ny2,nz1:nz2))
                                                    
! E_IO = VTK_VAR_XML(NC_NN = nn, varname= 'Vorticity', varX = p%loc%vort%x%now(nx1:nx2,ny1:ny2,nz1:nz2),&
!                                                     &varY = p%loc%vort%y%now(nx1:nx2,ny1:ny2,nz1:nz2),&
!                                                     &varZ = p%loc%vort%z%now(nx1:nx2,ny1:ny2,nz1:nz2))
                                                    
! E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Lamb',     varX = p%loc%lamb%x%now(nx1:nx2,ny1:ny2,nz1:nz2),&
!                                                     &varY = p%loc%lamb%y%now(nx1:nx2,ny1:ny2,nz1:nz2),&
!                                                     &varZ = p%loc%lamb%z%now(nx1:nx2,ny1:ny2,nz1:nz2))   

! E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Vort_adv', varX = p%loc%Vort_adv%x%now(nx1:nx2,ny1:ny2,nz1:nz2),&
!                                                     &varY = p%loc%Vort_adv%y%now(nx1:nx2,ny1:ny2,nz1:nz2),&
!                                                     &varZ = p%loc%Vort_adv%z%now(nx1:nx2,ny1:ny2,nz1:nz2)) 

! E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Vort_tws', varX = p%loc%Vort_tws%x%now(nx1:nx2,ny1:ny2,nz1:nz2),&
!                                                     &varY = p%loc%Vort_tws%y%now(nx1:nx2,ny1:ny2,nz1:nz2),&
!                                                     &varZ = p%loc%Vort_tws%z%now(nx1:nx2,ny1:ny2,nz1:nz2))  

! E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Vort_baro',varX = p%loc%Vort_baro%x%now(nx1:nx2,ny1:ny2,nz1:nz2),&
!                                                     &varY = p%loc%Vort_baro%y%now(nx1:nx2,ny1:ny2,nz1:nz2),&
!                                                     &varZ = p%loc%Vort_baro%z%now(nx1:nx2,ny1:ny2,nz1:nz2)) 

! E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Vort_visc',varX = p%loc%Vort_visc%x%now(nx1:nx2,ny1:ny2,nz1:nz2),&
!                                                     &varY = p%loc%Vort_visc%y%now(nx1:nx2,ny1:ny2,nz1:nz2),&
!                                                     &varZ = p%loc%Vort_visc%z%now(nx1:nx2,ny1:ny2,nz1:nz2)) 

E_IO = VTK_DAT_XML(var_location = 'node', var_block_action = 'close')
E_IO = VTK_GEO_XML()

E_IO = VTK_END_XML()

p%glb%pid = p%glb%pid + 1

end subroutine