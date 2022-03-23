subroutine output()
use all
implicit none
integer :: i,j,k,id,cnt,num,ks,ke,kk
real(8) :: rho, heavy, dv, x, eps, pi
real(8),dimension(:),allocatable :: mass
integer, dimension(1:p%glb%node_z) :: f
logical :: switch, finish

    ! level set method, loss of volume/mass in percentage
    write(p%fil%ls_mv,*)p%glb%time,100.0d0*(p%glb%imass-p%glb%mass)/p%glb%imass,100.0d0*(p%glb%ivol-p%glb%vol)/p%glb%ivol

    do k = 1, p%glb%node_z
        f(k) = -1
        do id = 0, p%glb%threads-1
            if( p%of(id)%loc%ks<=k .and. p%of(id)%loc%ke>=k)then
                do j = p%of(id)%loc%js, p%of(id)%loc%je
                do i = p%of(id)%loc%is, p%of(id)%Loc%ie
                    if( p%of(id)%loc%phi%now(i,j,k) > 1.0d-12 )then
                        f(k) = 1
                    endif
                enddo
                enddo
            endif
        enddo
    enddo

    cnt = 0
    do k = 1, p%glb%node_z-1
        if( f(k)*f(k+1) < 0 ) cnt = cnt + 1
    enddo
    cnt=cnt/2
    
    if(cnt<2)return

    allocate(mass(cnt))

    call p%rho_mu
    dv = p%glb%dx * p%glb%dy * p%glb%dz
    eps = 1.0d-12
    pi = dacos(-1.0_8)

    kk=1
    do num = 1, cnt

        switch=.false.
        finish = .false.

        do k = kk, p%glb%node_z
            if(finish)exit

            if( f(k)*f(k+1) < 0 )then
                if(.not.switch)then
                    ks=k+1
                    switch=.true.
                else
                    ke=k
                    finish=.true.
                    kk=k+1
                endif
            endif
        enddo

        mass(num)=0.0d0
        do k = ks, ke
            do id = 0, p%glb%threads-1
                if( p%of(id)%loc%ks<=k .and. p%of(id)%loc%ke>=k)then

                    do j = p%of(id)%loc%js, p%of(id)%loc%je
                    do i = p%of(id)%loc%is, p%of(id)%Loc%ie

                        mass(num)=mass(num)+p%of(id)%loc%rho%now(i,j,k)*p%of(id)%loc%heavy%now(i,j,k)*dv

                    enddo
                    enddo

                endif
            enddo
        enddo

    enddo

    if(p%glb%time<p%glb%dt)then
        p%glb%b1 = mass(1)
        p%glb%b2 = mass(2)

        write(*,*)p%glb%b1, p%glb%b2
    else
        write(*,*)cnt,(mass(1)/p%glb%b1-1.0)*100,(mass(2)/p%glb%b2-1.0)*100
        write(p%fil%bubble_mass,*)p%glb%time,(mass(1)/p%glb%b1-1.0)*100,(mass(2)/p%glb%b2-1.0)*100
    endif

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
