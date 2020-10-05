subroutine output()
use all
implicit none
integer :: i,j,k,id
real(8) :: damfront, damh
    ! level set method, loss of volume/mass in percentage
    write(p%fil%ls_mv,*)p%glb%time,100.0d0*(p%glb%imass-p%glb%mass)/p%glb%imass,100.0d0*(p%glb%ivol-p%glb%vol)/p%glb%ivol

    damfront = 0.0d0; damh=0.0d0
    !$omp parallel do private(i,j,k), reduction(max:damfront,damh)
    do id = 0, p%glb%threads-1
        
        if( p%of(id)%loc%idz==0 )then
            k=1
            do j = p%of(id)%loc%js, p%of(id)%loc%je
            do i = p%of(id)%loc%is, p%of(id)%loc%ie
                if( p%of(id)%loc%phi%now(i,j,k)*p%of(id)%loc%phi%now(i+1,j,k) < 0.0d0 )then
                    damfront = max( damfront, p%glb%x(i) + &
                        p%glb%dx*abs(p%of(id)%loc%phi%now(i,j,k))/( abs(p%of(id)%loc%phi%now(i,j,k))+abs(p%of(id)%loc%phi%now(i+1,j,k))) )
                endif
            enddo
            enddo   
        endif

        if( p%of(id)%loc%idx==0 )then
            i=1
            do k = p%of(id)%loc%ks, p%of(id)%loc%ke
            do j = p%of(id)%loc%js, p%of(id)%loc%je
                if( p%of(id)%loc%phi%now(i,j,k)*p%of(id)%loc%phi%now(i,j,k+1) < 0.0d0 )then
                    damh = max( damh, p%glb%z(k) + &
                        p%glb%dz*abs(p%of(id)%loc%phi%now(i,j,k))/( abs(p%of(id)%loc%phi%now(i,j,k))+abs(p%of(id)%loc%phi%now(i,j,k+1))) )
                endif
            enddo
            enddo
        endif

    enddo
    !$omp end parallel do 

    write(p%fil%damdata, *)p%glb%time, damfront, damh

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
