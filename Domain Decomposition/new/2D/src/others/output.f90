subroutine output()
use all
implicit none
integer :: i,j,id
real(8) :: damfront, damh
real(8) :: t, t1, t2

        ! level set method, loss of volume/mass in percentage
        call p%ls_mv; t=100.0d0*(p%glb%imass-p%glb%mass)/p%glb%imass
        !$omp parallel do private(i,j)
        do id = 0, p%glb%threads-1
        do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
        do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
            p%of(id)%loc%phi1%tmp(i,j) = p%of(id)%loc%phi%now(i,j)  ! store phi%now in phi1%tmp
            p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%phi1%now(i,j)
        enddo
        enddo
        enddo
        !$omp end parallel do

        call p%ls_mv; t1=100.0d0*(p%glb%imass/2-p%glb%mass)/p%glb%imass*2

        !$omp parallel do private(i,j)
        do id = 0, p%glb%threads-1
          do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
          do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
            p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%phi2%now(i,j)
          enddo
          enddo
        enddo
        !$omp end parallel do 

        call p%ls_mv; t2=100.0d0*(p%glb%imass/2-p%glb%mass)/p%glb%imass*2

        !$omp parallel do private(i,j)
        do id = 0, p%glb%threads-1
          do j = p%of(id)%loc%js-p%glb%ghc, p%of(id)%loc%je+p%glb%ghc
          do i = p%of(id)%loc%is-p%glb%ghc, p%of(id)%loc%ie+p%glb%ghc
            p%of(id)%loc%phi%now(i,j) = p%of(id)%loc%phi1%tmp(i,j)
          enddo
          enddo
        enddo
        !$omp end parallel do 

        write(p%fil%ls_mv,*)p%glb%time,t,t1,t2

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
