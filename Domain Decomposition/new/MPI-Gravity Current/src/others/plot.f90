subroutine plot
USE Lib_VTK_IO
use all
implicit none 
integer :: E_IO, id, num
integer :: nn, nx1,nx2,ny1,ny2,nz1,nz2
integer :: i,j,k
real(8),dimension(1:p%glb%node_x,1:p%glb%node_y,1:p%glb%node_z) :: x,y,z
character(3) :: name
character(2) :: rank

if ( abs(p%glb%time - p%glb%pid * p%glb%t2p) > p%glb%dt ) return

call vortex_dynamics

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
write(rank,'(I2.2)')p%glb%mpirank


if( p%glb%split==1 )then

    if( p%glb%mpirank==0 )then
        E_IO = VTK_INI_XML(output_format='raw', filename='./out/'//rank//'_'//name//'.vts', &
                           mesh_topology='StructuredGrid', nx1=p%glb%mpis(p%glb%mpirank), nx2=p%glb%mpie(p%glb%mpirank), ny1=1, ny2=p%glb%node_y, nz1=1, nz2=p%glb%node_z)
    else
        E_IO = VTK_INI_XML(output_format='raw', filename='./out/'//rank//'_'//name//'.vts', &
                           mesh_topology='StructuredGrid', nx1=p%glb%mpis(p%glb%mpirank)-1, nx2=p%glb%mpie(p%glb%mpirank), ny1=1, ny2=p%glb%node_y, nz1=1, nz2=p%glb%node_z)
    endif

else if (p%glb%split==2 )then   

    if( p%glb%mpirank==0 )then
        E_IO = VTK_INI_XML(output_format='raw', filename='./out/'//rank//'_'//name//'.vts', &
                           mesh_topology='StructuredGrid', nx1=1, nx2=p%glb%node_x, ny1=p%glb%mpis(p%glb%mpirank), ny2=p%glb%mpie(p%glb%mpirank), nz1=1, nz2=p%glb%node_z)
    else
        E_IO = VTK_INI_XML(output_format='raw', filename='./out/'//rank//'_'//name//'.vts', &
                           mesh_topology='StructuredGrid', nx1=1, nx2=p%glb%node_x, ny1=p%glb%mpis(p%glb%mpirank)-1, ny2=p%glb%mpie(p%glb%mpirank), nz1=1, nz2=p%glb%node_z)
    endif

else

    if( p%glb%mpirank==0 )then
        E_IO = VTK_INI_XML(output_format='raw', filename='./out/'//rank//'_'//name//'.vts', &
                           mesh_topology='StructuredGrid', nx1=1, nx2=p%glb%node_x, ny1=1, ny2=p%glb%node_y, nz1=p%glb%mpis(p%glb%mpirank), nz2=p%glb%mpie(p%glb%mpirank))
    else
        E_IO = VTK_INI_XML(output_format='raw', filename='./out/'//rank//'_'//name//'.vts', &
                           mesh_topology='StructuredGrid', nx1=1, nx2=p%glb%node_x, ny1=1, ny2=p%glb%node_y, nz1=p%glb%mpis(p%glb%mpirank)-1, nz2=p%glb%mpie(p%glb%mpirank))
    endif

endif

do id = 0, p%glb%threads-1

nx1 = p%of(id)%loc%is; nx2=p%of(id)%loc%ie
ny1 = p%of(id)%loc%js; ny2=p%of(id)%loc%je
nz1 = p%of(id)%loc%ks; nz2=p%of(id)%loc%ke

if( id<p%glb%grid_A .and. p%glb%mpirank>0 )then

    if(p%glb%split==1)then
        nx1=nx1-1
    else if (p%glb%split==2)then
        ny1=ny1-1
    else
        nz1=nz1-1
    endif

endif

nn=(nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1)

E_IO = VTK_GEO_XML(nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2, nz1=nz1, nz2=nz2, NN=nn, &
                   X=x(nx1:nx2,ny1:ny2,nz1:nz2),Y=y(nx1:nx2,ny1:ny2,nz1:nz2),Z=z(nx1:nx2,ny1:ny2,nz1:nz2))
                   
E_IO = VTK_DAT_XML(var_location = 'node', var_block_action = 'open')
E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Concentration', var = p%of(id)%loc%c%now(nx1:nx2,ny1:ny2,nz1:nz2) )
E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Pressure', var = p%of(id)%loc%p%now(nx1:nx2,ny1:ny2,nz1:nz2) )
E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Density', var = p%of(id)%loc%rho%now(nx1:nx2,ny1:ny2,nz1:nz2) )
E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Q', var = p%of(id)%loc%q_cri%now(nx1:nx2,ny1:ny2,nz1:nz2) )
E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Omega', var = p%of(id)%loc%omega_cri%now(nx1:nx2,ny1:ny2,nz1:nz2) )
E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Lamb_div', var = p%of(id)%loc%lamb_div%now(nx1:nx2,ny1:ny2,nz1:nz2) )

E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Velocity', varX = p%of(id)%loc%nvel%x%now(nx1:nx2,ny1:ny2,nz1:nz2),&
                                                    &varY = p%of(id)%loc%nvel%y%now(nx1:nx2,ny1:ny2,nz1:nz2),&
                                                    &varZ = p%of(id)%loc%nvel%z%now(nx1:nx2,ny1:ny2,nz1:nz2))
                                                    
E_IO = VTK_VAR_XML(NC_NN = nn, varname= 'Vorticity', varX = p%of(id)%loc%vort%x%now(nx1:nx2,ny1:ny2,nz1:nz2),&
                                                    &varY = p%of(id)%loc%vort%y%now(nx1:nx2,ny1:ny2,nz1:nz2),&
                                                    &varZ = p%of(id)%loc%vort%z%now(nx1:nx2,ny1:ny2,nz1:nz2))
                                                    
E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Lamb',     varX = p%of(id)%loc%lamb%x%now(nx1:nx2,ny1:ny2,nz1:nz2),&
                                                    &varY = p%of(id)%loc%lamb%y%now(nx1:nx2,ny1:ny2,nz1:nz2),&
                                                    &varZ = p%of(id)%loc%lamb%z%now(nx1:nx2,ny1:ny2,nz1:nz2))   

E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Vort_adv', varX = p%of(id)%loc%Vort_adv%x%now(nx1:nx2,ny1:ny2,nz1:nz2),&
                                                    &varY = p%of(id)%loc%Vort_adv%y%now(nx1:nx2,ny1:ny2,nz1:nz2),&
                                                    &varZ = p%of(id)%loc%Vort_adv%z%now(nx1:nx2,ny1:ny2,nz1:nz2)) 

E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Vort_tws', varX = p%of(id)%loc%Vort_tws%x%now(nx1:nx2,ny1:ny2,nz1:nz2),&
                                                    &varY = p%of(id)%loc%Vort_tws%y%now(nx1:nx2,ny1:ny2,nz1:nz2),&
                                                    &varZ = p%of(id)%loc%Vort_tws%z%now(nx1:nx2,ny1:ny2,nz1:nz2))  

E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Vort_baro',varX = p%of(id)%loc%Vort_baro%x%now(nx1:nx2,ny1:ny2,nz1:nz2),&
                                                    &varY = p%of(id)%loc%Vort_baro%y%now(nx1:nx2,ny1:ny2,nz1:nz2),&
                                                    &varZ = p%of(id)%loc%Vort_baro%z%now(nx1:nx2,ny1:ny2,nz1:nz2)) 

E_IO = VTK_VAR_XML(NC_NN = nn, varname = 'Vort_visc',varX = p%of(id)%loc%Vort_visc%x%now(nx1:nx2,ny1:ny2,nz1:nz2),&
                                                    &varY = p%of(id)%loc%Vort_visc%y%now(nx1:nx2,ny1:ny2,nz1:nz2),&
                                                    &varZ = p%of(id)%loc%Vort_visc%z%now(nx1:nx2,ny1:ny2,nz1:nz2)) 

E_IO = VTK_DAT_XML(var_location = 'node', var_block_action = 'close')
E_IO = VTK_GEO_XML()
enddo

E_IO = VTK_END_XML()

if(p%glb%mpirank==0 .and. p%glb%mpisize>1 )then

    E_IO = PVTK_INI_XML(filename = './out/'//trim(p%glb%name)//'_'//name//'.pvts', mesh_topology = 'PStructuredGrid', &
                                 nx1=1,nx2=p%glb%node_x, ny1=1, ny2=p%glb%node_y, nz1=1, nz2=p%glb%node_z, tp='Float64')
    if( p%glb%split==1 )then
        do num = 0, p%glb%mpisize-1
            write(rank,'(I2.2)')num
            if(num>0)then
                E_IO = PVTK_GEO_XML(nx1=p%glb%mpis(num)-1, nx2=p%glb%mpie(num), ny1=1, ny2=p%glb%node_y, nz1=1, nz2=p%glb%node_z, source=rank//'_'//name//'.vts')
            else
                E_IO = PVTK_GEO_XML(nx1=p%glb%mpis(num)  , nx2=p%glb%mpie(num), ny1=1, ny2=p%glb%node_y, nz1=1, nz2=p%glb%node_z, source=rank//'_'//name//'.vts')
            endif
        enddo
    else if (p%glb%split==2)then
        do num = 0, p%glb%mpisize-1
            write(rank,'(I2.2)')num
            if(num>0)then
                E_IO = PVTK_GEO_XML(nx1=1,nx2=p%glb%node_x, ny1=p%glb%mpis(num)-1, ny2=p%glb%mpie(num), nz1=1, nz2=p%glb%node_z, source=rank//'_'//name//'.vts')
            else
                E_IO = PVTK_GEO_XML(nx1=1,nx2=p%glb%node_x, ny1=p%glb%mpie(num)  , ny2=p%glb%mpie(num), nz1=1, nz2=p%glb%node_z, source=rank//'_'//name//'.vts')
            endif
        enddo
    else 
        do num = 0, p%glb%mpisize-1
            write(rank,'(I2.2)')num
            if(num>0)then
                E_IO = PVTK_GEO_XML(nx1=1,nx2=p%glb%node_x, ny1=1, ny2=p%glb%node_y, nz1=p%glb%mpis(num)-1, nz2=p%glb%mpie(num), source=rank//'_'//name//'.vts')
            else
                E_IO = PVTK_GEO_XML(nx1=1,nx2=p%glb%node_x, ny1=1, ny2=p%glb%node_y, nz1=p%glb%mpis(num)  , nz2=p%glb%mpie(num), source=rank//'_'//name//'.vts')
            endif
        enddo
    endif

    E_IO = PVTK_DAT_XML(var_location = 'node', var_block_action = 'open')

    E_IO = PVTK_VAR_XML(varname='Concentration', tp='Float64', NC=1)
    E_IO = PVTK_VAR_XML(varname='Pressure', tp='Float64', NC=1)
    E_IO = PVTK_VAR_XML(varname='Density', tp='Float64', NC=1)
    E_IO = PVTK_VAR_XML(varname='Q', tp='Float64', NC=1)
    E_IO = PVTK_VAR_XML(varname='Omega', tp='Float64', NC=1)
    E_IO = PVTK_VAR_XML(varname='Lamb_div', tp='Float64', NC=1)

    E_IO = PVTK_VAR_XML(varname='Velocity', tp='Float64', NC=3)
    E_IO = PVTK_VAR_XML(varname='Vorticity', tp='Float64', NC=3)
    E_IO = PVTK_VAR_XML(varname='Lamb', tp='Float64', NC=3)
    E_IO = PVTK_VAR_XML(varname='Vort_adv', tp='Float64', NC=3)
    E_IO = PVTK_VAR_XML(varname='Vort_tws', tp='Float64', NC=3)
    E_IO = PVTK_VAR_XML(varname='Vort_baro', tp='Float64', NC=3)
    E_IO = PVTK_VAR_XML(varname='Vort_visc', tp='Float64', NC=3)

    E_IO = PVTK_DAT_XML(var_location = 'node', var_block_action = 'close')
    E_IO = PVTK_END_XML()
endif

p%glb%pid = p%glb%pid + 1
call p%sync

end subroutine
