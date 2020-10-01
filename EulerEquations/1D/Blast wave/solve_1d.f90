SUBROUTINE SOLVE_1D(BTN)
use dat
IMPLICIT NONE
INTEGER :: BTN, i, cnt

!$omp parallel do
do i = -2, N+3
    u_old(i) = u(i)
    rho_old(i) = rho(i)
    ix_old(i) = ix(i)
    e_old(i) = e(i)
    p_old(i) = p(i)
    h_old(i) = h(i)
    a_old(i) = a(i)
end do
!$omp end parallel do

do cnt = 1, 4

    !$omp parallel do
    do i = -2, N+3
        u(i)   =   u_OLD(i)
        rho(i) = rho_OLD(i)
        ix(i)  =  ix_OLD(i)
        e(i)   =   e_OLD(i)
        p(i)   =   p_OLD(i)
        h(i)   =   h_OLD(i)
        a(i)   =   a_OLD(i)
    end do
    !$omp end parallel do
    
    CALL FIND_SPEED(btn)
    CALL RK3_RHO(BTN)
    CALL RK3_IX(BTN)
    CALL RK3_E(BTN)
    CALL RENEW

end do
   
END SUBROUTINE

subroutine find_speed(btn)
use dat
implicit none
integer :: i, btn
real(8) :: vc, va, vh

!if( btn==0 ) then
!    CALL weno_JS(rho,FR,FL,N)
!    CALL WENO_JS(U,GR,GL,N)
!    CALL WENO_JS(H,HR,HL,N)
!else if (btn==1)then
!    CALL CRWENO(rho,FR,FL,N)
!    CALL CRWENO(U,GR,GL,N)
!    CALL CRWENO(H,HR,HL,N)
!else IF(BTN==2)then
!    CALL CRWENO_LD(rho,FR,FL,N)
!    CALL CRWENO_LD(U,GR,GL,N)
!    CALL CRWENO_LD(H,HR,HL,N)
!ELSE
!    CALL oCRWENO_LD(rho,FR,FL,N)
!    CALL OCRWENO_LD(U,GR,GL,N)
!    CALL OCRWENO_LD(H,HR,HL,N)
!end if
!
!CALL limiter(rho,fr,fl,n)
!CALL LIMITER(U,GR,GL,N)
!CALL LIMITER(H,HR,HL,N)

 !$omp parallel private(vc,va,vh)
 do i = 0, N

    vc = (dsqrt(rho(i))*u(i)+dsqrt(rho(i+1))*u(i+1))/(dsqrt(rho(i))+dsqrt(rho(i+1)))
    vh = (dsqrt(rho(i))*h(i)+dsqrt(rho(i+1))*h(i+1))/(dsqrt(rho(i))+dsqrt(rho(i+1)))
    va = dsqrt((gamma-1.0d0)*(vh-0.5d0*vc**2.0d0))

    !vc = (dsqrt(FL(i))*GL(i)+dsqrt(FR(i))*GR(i))/(dsqrt(FL(i))+dsqrt(FR(i)))
    !vh = (dsqrt(FL(i))*HL(i)+dsqrt(GR(i))*HR(i))/(dsqrt(FR(i))+dsqrt(FR(i)))
    !va = dsqrt((gamma-1.0d0)*(vh-0.5d0*vc**2.0d0))


    sl(i) = vc-va
    sr(i) = vc+va 

 end do 
 !$omp end parallel 

end subroutine

subroutine limiter(f,fp,fm,n)
implicit none
real(8), dimension(-2:N+3) :: f,fp,fm,r
integer :: i,n
real(8) :: rd, rr

rd=0.75d0

!$omp parallel do
do i = 0, n
    R(i) = ( abs(2.0d0*(f(i+1)-f(i))*f(i)-f(i-1)) + 1.0d-13 )/ ( (f(i+1)-f(i))**2.0d0 + (f(i)-f(i-1))**2.0d0 + 1.0d-13 )
end do
!$omp end parallel do

!$omp parallel do private(rr)
do i = 0, n

    rr = min(r(i),r(i+1))
    rr = 0.5d0 * ( 1.0d0 + dtanh(3.0d0*(rr-rd)/max(rd,abs(rr-rd)))/dtanh(3.0d0) )
    rr = 1.0d0 - rr
    fp(i) = ( 1.0d0-rr ) * fp(i) + rr * (-f(i-1)+7.0d0*f(i)+7.0d0*f(i+1)-f(i+2))/12.0d0
    fm(i) = ( 1.0d0-rr ) * fm(i) + rr * (-f(i-1)+7.0d0*f(I)+7.0d0*f(i+1)-f(i+2))/12.0d0  
    !fp(i) = ( 1.0d0-rr ) * fp(i) + rr * 0.5d0*(f(i+1)+f(i))
    !fm(i) = ( 1.0d0-rr ) * fm(i) + rr * 0.5d0*(f(i+1)+f(i))  

end do
!$omp end parallel do



end subroutine
