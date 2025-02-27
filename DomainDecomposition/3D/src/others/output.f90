subroutine output()
use all
implicit none
integer :: i,j,k,id,kk
real(8) :: eta,num,x,r,kh,exact,exact_vof,theta,pi
real(8) :: eta_vof,xx
real(8) :: nx,ny,nz,dz

pi=dacos(-1.0d0)

! level set method, loss of volume/mass in percentage
write(p%fil%ls_mv,*)p%glb%time,100.0d0*(p%glb%imass-p%glb%mass)/p%glb%imass,100.0d0*(p%glb%ivol-p%glb%vol)/p%glb%ivol

return

call p%surface_norms
call pt%normals%sync
    
eta=0.0d0
num=0.0d0
x=0.0d0
!$omp parallel do private(i,j,k,kk,nx,ny,nz,dz), reduction(+:eta,num,x)
do id = 0, p%glb%threads-1

    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        if( abs(p%of(id)%glb%x(i,j,k))<p%glb%dx*0.5d0 )then
            do k = p%of(id)%loc%ks, p%of(id)%loc%ke
            do j = p%of(id)%loc%js, p%of(id)%loc%je
                if(p%of(id)%loc%phi%now(i,j,k)*p%of(id)%loc%phi%now(i,j,k+1)<0.0d0)then
                    x = x + p%of(id)%glb%x(i,j,k)
                    num =  num + 1.0d0
                    do kk = k, k+1
                        nx = p%of(id)%loc%normals%x%now(i,j,kk)
                        ny = p%of(id)%loc%normals%y%now(i,j,kk)
                        nz = p%of(id)%loc%normals%z%now(i,j,kk)
                        dz = p%of(id)%loc%phi%now(i,j,kk) * dsqrt(nx*nx+ny*ny+nz*nz) / abs(nz) 
                        eta = eta + 0.5d0 * ( p%of(id)%glb%z(i,j,kk) + dz )
                    enddo
                endif
            enddo
            enddo
        endif
    enddo

enddo
!$omp end parallel do

x=x/num
eta=eta/num

eta_vof=0.0d0
num=0.0d0
xx=0.0d0
!$omp parallel do private(i,j,k,nx,ny,nz), reduction(+:eta_vof,num,xx)
do id = 0, p%glb%threads-1

    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        if( abs(p%of(id)%glb%x(i,j,k))<p%glb%dx*0.5d0 )then
            do k = p%of(id)%loc%ks, p%of(id)%loc%ke
            do j = p%of(id)%loc%js, p%of(id)%loc%je
                if( p%of(id)%loc%vof%now(i,j,k+1)<1.0d-8 .and. p%of(id)%loc%vof%now(i,j,k)>1.0d-8 )then
                    nx = p%of(id)%loc%normals%x%now(i,j,k)
                    ny = p%of(id)%loc%normals%y%now(i,j,k)
                    nz = p%of(id)%loc%normals%z%now(i,j,k)
                    xx = xx + p%of(id)%glb%x(i,j,k)
                    num =  num + 1.0d0
                    eta_vof = eta_vof + p%of(id)%glb%z(i,j,k-1) + p%glb%dz * find_pos(p%of(id)%loc%vof%now(i,j,k),Nx,Ny,Nz,0.5d0,0.5d0,0.0d0)
                endif
            enddo
            enddo
        endif
    enddo

enddo
!$omp end parallel do

xx=xx/num
eta_vof=eta_vof/num

kh = p%wa%k * abs(p%glb%zstart)

theta = p%wa%k * x - p%wa%w * p%glb%time
exact = p%wa%L * Stokes_wave_interface(theta,p%wa%steepness,kh) 

theta = p%wa%k * xx - p%wa%w * p%glb%time
exact_vof = p%wa%L * Stokes_wave_interface(theta,p%wa%steepness,kh) 

write(p%fil%wave,*)p%glb%time/(2.0d0*pi),eta,exact,eta_vof,exact_vof
write(*,'("LS:",2ES15.4)')eta,exact
write(*,'("VOF:",2ES15.4)')eta_vof,exact_vof

end subroutine

subroutine print_NS_info()
use all 
implicit none
        write(*,'("Divergence :",2ES15.4)')p%glb%vel_div,p%glb%vel_sdiv
        write(*,'("L2 norm    :",ES15.4)')p%glb%ns_l2f
        write(*,'("Linf norm  :",ES15.4)')p%glb%ns_linf
        write(*,*)''
        write(*,'("PPE iters  :",I15)')p%glb%piter
        write(*,'("PPE error  :",ES15.4)')p%glb%ppe_linf
        write(*,*)''        
end subroutine

subroutine print_LS_info()
use all
implicit none
        write(*,'("LS,  Loss of mass  (%) :",ES15.4)')100.0d0*(p%glb%imass-p%glb%mass)/p%glb%imass
        write(*,'("LS,  Loss of volume(%) :",ES15.4)')100.0d0*(p%glb%ivol-p%glb%vol)/p%glb%ivol
        write(*,'("LS,  redistance error  :",ES15.4)')p%glb%red_error
        write(*,*)''
        if(p%glb%method==3)then
            write(*,'("VOF, Loss of mass  (%) :",ES15.4)')100.0d0*(p%glb%imassv-p%glb%massv)/p%glb%imassv
            write(*,'("VOF, Loss of volume(%) :",ES15.4)')100.0d0*(p%glb%ivolv-p%glb%volv)/p%glb%ivolv
            write(*,*)''
        endif
end subroutine

subroutine print_CPU_info()
use all 
implicit none
real(8) :: total, totald

        call pt%cputime(totald)
        total = p%glb%ls_adv + p%glb%ls_red + p%glb%ns
        write(*,'("Total CPU time(s) :",F15.6)')total
        write(*,'(4A18)')"Inter. Adv.","Inter. Recon.","PPE","NS"
        write(*,'(F17.2,"%",F17.2,"%",F17.2,"%",F17.2,"%")')100.0d0*p%glb%ls_adv/total,100.0d0*p%glb%ls_red/total&
                                                &,100.0d0*p%glb%ppe/total,100.0d0*(p%glb%ns-p%glb%ppe)/total
        write(*,*)''
        write(*,'(A18,F17.2,"%")')"Data Sync:",totald/total*100.0d0
end subroutine
