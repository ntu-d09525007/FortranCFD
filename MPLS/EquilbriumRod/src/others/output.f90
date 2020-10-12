subroutine output()
use all
implicit none
integer :: i,j,id
real(8) :: damfront, damh
    ! level set method, loss of volume/mass in percentage
    write(p%fil%ls_mv,*)p%glb%time,100.0d0*(p%glb%imass-p%glb%mass)/p%glb%imass,100.0d0*(p%glb%ivol-p%glb%vol)/p%glb%ivol

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

subroutine pressure_error
use all
implicit none
integer :: id,i,j
real(8) :: pout, pin, in, out, umax
real(8) :: L1, L2, curv

pout=0.0d0
out=0.0d0
umax=0.0d0
!$omp parallel do private(i,j), reduction(+:pout,out), reduction(max:umax)
do id = 0, p%glb%threads-1
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        umax = max( umax, p%of(id)%loc%nvel%x%now(i,j)**2.0d0+p%of(id)%loc%nvel%y%now(i,j)**2.0d0 )
        if (p%of(id)%loc%heavy%now(i,j) < 0.01 )then
            out = out + 1.0d0
            pout = pout + p%of(id)%loc%p%now(i,j)
        endif
    enddo
    enddo
enddo
!$omp end parallel do

pout = pout / out
umax = dsqrt( umax )

in=0.0d0;L1=0.0d0;L2=0.0d0
!$omp parallel do private(i,j), reduction(+:pin,in,L1,L2)
do id = 0, p%glb%threads-1
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        p%of(id)%loc%p%now(i,j) = p%of(id)%loc%p%now(i,j) - pout
        if( p%of(id)%loc%heavy%now(i,j)>0.99d0 )then
            pin = pin + p%of(id)%loc%p%now(i,j)
            in = in + 1.0d0
            L1 = L1 + p%glb%we*abs( p%of(id)%loc%p%now(i,j) - 1.0d0/p%glb%we )
            L2 = L2 + ( p%glb%we*abs( p%of(id)%loc%p%now(i,j) - 1.0d0/p%glb%we ) )**2.0d0
        endif
    enddo
    enddo
enddo
!$omp end parallel do

pin = pin / in

write(*,'(A,F8.5)')"Pressure inside the rod:",pin*p%glb%we
write(*,'(A,ES11.4)')"L1 error:",L1/in
write(*,'(A,ES11.4)')"L2 error:",dsqrt(L2/in)
write(*,'(A,ES11.4)')"Maximum Velocity:",umax

call p%curv

in=0.0d0;L1=0.0d0;L2=0.0d0;curv=0.0d0
!$omp parallel do private(i,j), reduction(+:in,L1,L2,curv)
do id = 0, p%glb%threads-1
    do j = p%of(id)%loc%js, p%of(id)%loc%je
    do i = p%of(id)%loc%is, p%of(id)%loc%ie
        if( abs(p%of(id)%loc%phi%now(i,j))<p%glb%ls_wid/2.0d0 )then
            in = in + 1.0d0
            curv = curv + abs(p%of(id)%loc%normals%curv%now(i,j))
            L1 = L1 + abs(abs(p%of(id)%loc%normals%curv%now(i,j))-1.0d0)
            L2 = L2 + abs(abs(p%of(id)%loc%normals%curv%now(i,j))-1.0d0)**2.0d0
        endif
    enddo
    enddo
enddo
!$omp end parallel do

write(*,'(A,F8.5)')"Curvature:",curv/in
write(*,'(A,ES11.4)')"L1 error:",L1/in
write(*,'(A,ES11.4)')"L2 error:",dsqrt(L2/in)

end subroutine
