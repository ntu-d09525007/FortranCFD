subroutine manager_ls_funs(p)
implicit none
class(manager) :: p
integer :: id, i, j
real(8) :: x, heavy, hp, pi, eps

eps = 1.0d-12
pi = dacos(-1.0_8)
    
!$omp parallel do private(i,j,x,heavy,hp)
do id = 0, p%glb%threads-1

    do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
    do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
    
        x = p%of(id)%loc%phi%now(i,j) / p%glb%ls_wid
        
        if( x > 1.0_8-eps )then
            heavy = 1.0_8
            hp = 0.0_8
        else if ( x < -1.0_8+eps )then
            heavy = 0.0_8
            hp = 0.0_8
        else
            heavy = 0.5_8 * (1.0_8 + x + dsin(pi*x) / pi )
            hp = 0.5_8 * ( 1.0_8 + dcos(pi*x) ) / p%glb%ls_wid 
        endif
        
        heavy = max(min(heavy,1.0d0),0.0d0) 
        
        p%of(id)%loc%heavy%now(i,j) = heavy
        p%of(id)%loc%delta%now(i,j) = hp
        p%of(id)%loc%sign%now(i,j) = 2.0_8*heavy-1.0_8
        
    end do
    end do
    
enddo
!$omp end parallel  do

end subroutine

subroutine manager_rho_mu(p)
implicit none
class(manager) :: p
integer :: id,i,j
real(8) :: heavy

call p%ls_funs

!$omp parallel do private(i,j,heavy)
do id = 0, p%glb%threads-1
          
    do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
    do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        
        heavy = p%of(id)%loc%vof%now(i,j)
        if( p%glb%method .ne. 3)heavy = p%of(id)%loc%heavy%now(i,j)
        
        p%of(id)%loc%rho%now(i,j) = heavy + p%glb%rho_12 * (1.0_8 - heavy )
        p%of(id)%loc%mu%now(i,j)  = heavy + p%glb%mu_12  * (1.0_8 - heavy )
        
    end do
    end do
 
enddo       
!$omp end parallel do
    
    
end subroutine

subroutine manager_ls_mv(p)
implicit none
class(manager) :: p
integer :: id, i, j
real(8) :: mass, vol, rho
real(8) :: dv

dv = p%glb%dx * p%glb%dy 

call p%rho_mu

!===========================  LS 

mass = 0.0_8; vol=0.0_8

!$omp parallel do private(i,j,rho), reduction(+:mass,vol)    
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
        rho = p%of(id)%loc%heavy%now(i,j) + p%glb%rho_12 * (1.0d0 - p%of(id)%loc%heavy%now(i,j))
        mass = mass + rho*p%of(id)%loc%heavy%now(i,j)*dv
        vol = vol + p%of(id)%loc%heavy%now(i,j)*dv
    
    enddo
    enddo

enddo
!$omp end parallel do

p%glb%mass = mass
p%glb%vol = vol

!===========================  VOF 

mass = 0.0_8; vol=0.0_8

!$omp parallel do private(i,j,rho), reduction(+:mass,vol)    
do id = 0, p%glb%threads-1
    
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        
        rho = p%of(id)%loc%vof%now(i,j) + p%glb%rho_12 * (1.0d0-p%of(id)%loc%vof%now(i,j))
        mass = mass + rho*p%of(id)%loc%vof%now(i,j)*dv
        vol = vol + p%of(id)%loc%vof%now(i,j)*dv
    
    enddo
    enddo

enddo
!$omp end parallel do

p%glb%massv = mass
p%glb%volv = vol

call p%sync

end subroutine

subroutine manager_surface_norms(p)
implicit none
class(manager) :: p
integer :: id,I,J
    
    !$omp parallel do private(i,j)
    do id = 0, p%glb%threads-1
        
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        
            call p%of(id)%loc%ccdsolvers%x%solve("ccd", p%of(id)%loc%phi%now(:,j),&
                        &p%of(id)%loc%normals%x%now(:,j),p%of(id)%loc%normals%xx%now(:,j) )
                                                                    
        enddo

        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        
            call p%of(id)%loc%ccdsolvers%y%solve("ccd", p%of(id)%loc%phi%now(i,:),&
                        &p%of(id)%loc%normals%y%now(i,:),p%of(id)%loc%normals%yy%now(i,:) )
                                                                            
        enddo
                
        !===========================================
        
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
        
            call p%of(id)%loc%ccdsolvers%y%solve("ccd", p%of(id)%loc%normals%x%now(i,:),&
                            &p%of(id)%loc%normals%xy%now(i,:),p%of(id)%loc%normals%curv%now(i,:) )
                                                                            
        enddo

        !===========================================
        
        call p%of(id)%bc(0,p%of(id)%loc%normals%x%now)
        call p%of(id)%bc(0,p%of(id)%loc%normals%xx%now)
        
        call p%of(id)%bc(0,p%of(id)%loc%normals%y%now)
        call p%of(id)%bc(0,p%of(id)%loc%normals%yy%now)
         
        call p%of(id)%bc(0,p%of(id)%loc%normals%xy%now)
     
    enddo       
    !$omp end parallel do
    
end subroutine

subroutine manager_surface_norms_sec(p)
implicit none
class(manager) :: p
integer :: id,I,J  

!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1

    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
    
        p%of(id)%loc%normals%x%now(i,j) = 0.5d0*( p%of(id)%loc%phi%now(i+1,j) - p%of(id)%loc%phi%now(i-1,j) )/p%glb%dx
        p%of(id)%loc%normals%y%now(i,j) = 0.5d0*( p%of(id)%loc%phi%now(i,j+1) - p%of(id)%loc%phi%now(i,j-1) )/p%glb%dy

        p%of(id)%loc%normals%xx%now(i,j) = ( p%of(id)%loc%phi%now(i+1,j) - 2.0d0*p%of(id)%loc%phi%now(i,j) - p%of(id)%loc%phi%now(i-1,j) )/p%glb%dx**2.0d0
        p%of(id)%loc%normals%yy%now(i,j) = ( p%of(id)%loc%phi%now(i,j+1) - 2.0d0*p%of(id)%loc%phi%now(i,j) - p%of(id)%loc%phi%now(i,j-1) )/p%glb%dy**2.0d0

        p%of(id)%loc%normals%xy%now(i,j) = ( p%of(id)%loc%phi%now(i+1,j+1)+p%of(id)%loc%phi%now(i-1,j-1) &
                                           & - p%of(id)%loc%phi%now(i-1,j+1)-p%of(id)%loc%phi%now(i+1,j-1) )/(4.0d0*p%glb%dx*p%glb%dy)
                                           
    end do
    end do

    !===========================================
    
    call p%of(id)%bc(0,p%of(id)%loc%normals%x%now)
    call p%of(id)%bc(0,p%of(id)%loc%normals%xx%now)
    
    call p%of(id)%bc(0,p%of(id)%loc%normals%y%now)
    call p%of(id)%bc(0,p%of(id)%loc%normals%yy%now)
     
    call p%of(id)%bc(0,p%of(id)%loc%normals%xy%now)
 
end do      
!$omp end parallel do
    
end subroutine

subroutine manager_curv(p)
implicit none
class(manager) :: p
integer :: id,i,j
real(8) :: fx,fxx,fy,fyy,fxy

call p%surface_norms

!$omp parallel do private(i,j,fx,fxx,fy,fyy,fxy)
do id = 0, p%glb%threads-1   

    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
    
        fx  = p%of(id)%loc%normals%x%now(i,j)
        fy  = p%of(id)%loc%normals%y%now(i,j)
        
        fxx = p%of(id)%loc%normals%xx%now(i,j)
        fyy = p%of(id)%loc%normals%yy%now(i,j)
        
        fxy = p%of(id)%loc%normals%xy%now(i,j)
        
        p%of(id)%loc%normals%curv%now(i,j) = ( fyy*fx**2.0d0+fxx*fy**2.0d0-2.0d0*fxy*fx*fy ) / (fx**2.0d0+fy**2.0d0+1.0d-12)**1.5d0   

    end do
    end do
    
    call p%of(id)%bc(0,p%of(id)%loc%normals%curv%now)
    
enddo
!$omp end parallel do

end subroutine