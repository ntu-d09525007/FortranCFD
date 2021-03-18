subroutine output()
use all
implicit none
integer :: id,i,j
real(8) :: pi, eta, r, x

pi=dacos(-1.0d0)

! level set method, loss of volume/mass in percentage
write(p%fil%ls_mv,*)p%glb%time,100.0d0*(p%glb%imass-p%glb%mass)/p%glb%imass,100.0d0*(p%glb%ivol-p%glb%vol)/p%glb%ivol

eta=0.0d0
!$omp parallel do private(i,j)
do id = 0, p%glb%threads-1

    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        if( abs(p%glb%x(i,j)-2.5d0*2.0d0*pi)<p%glb%dx*0.5d0)then
            do j = p%of(id)%loc%js, p%of(id)%loc%je
                if(p%of(id)%loc%phi%now(i,j)*p%of(id)%loc%phi%now(i,j+1)<0.0d0)then
                    x = p%glb%x(i,j)
                    r = abs(p%of(id)%loc%phi%now(i,j))/(abs(p%of(id)%loc%phi%now(i,j))+abs(p%of(id)%loc%phi%now(i,j+1)))
                    eta = p%glb%y(i,j) + p%glb%dy * r
                endif
            enddo
        endif
    enddo

enddo
!$omp end parallel do

write(p%fil%wave,*)p%glb%time,eta,p%wa%L*dcos(p%wa%k*x-p%wa%w*p%glb%time)

write(*,'(2ES15.4)')eta,p%wa%L*dcos(p%wa%k*x-p%wa%w*p%glb%time)

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
