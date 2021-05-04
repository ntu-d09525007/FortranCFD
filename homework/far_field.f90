subroutine task_run2(p,dx,name)
implicit none
class(task) :: p
real(8) :: dx
character(*) :: name
integer :: i

p%xstart=-20.0
p%xend=20.0

p%dx=dx
p%g=9.81
p%alpha=0.531
p%n=(p%xend-p%xstart)/p%dx + 1
p%name = trim(name)

allocate(p%A(p%n),p%B(p%n),p%C(p%n),p%S(p%n),p%x(p%n),p%diff(p%n))

call p%h%init(p%n)
call p%u%init(p%n)
call p%phi%init(p%n)
call p%eta%init(p%n)
call p%psi%init(p%n)
call p%hu%init(p%n)

!$omp parallel do 
do i =  1, p%n
    p%x(i) = p%xstart + (i-1)*dx
    p%h%now(i) = 0.3
    p%eta%now(i) = 0.1*dexp(-18.0*p%x(i)**2)
    p%u%now(i)=0.0
    p%phi%now(i) = 0.0
    p%diff(i) = 1.0d0
enddo
!$omp end parallel do 
call p%h%bc
call p%eta%bc
call p%phi%bc
call p%u%bc

p%A0 = -0.02088856896d0*derf(4.242640687d0*p%xstart) + 0.02088856896d0*derf(4.242640687d0*p%xend)
p%coef = (2.0d0/(dsqrt(p%g*0.3)* 0.3**2))**(1.0d0/3.0d0)

p%pltid=0
p%dt=0.9*p%dx/dsqrt(p%g*0.4)
p%t2s=10.0
p%t2p=2.5d0

call p%plot2

do 
    p%t=p%t+p%dt
    call p%solve
    call p%plot2
    if(p%t>p%t2s)exit

enddo

end subroutine

subroutine task_plot2(p)
! ==============================
! Plot the data with user-defined period
! ==============================
implicit none
class(task) :: p
integer :: i
character(2) :: name
real(8) :: x, xx, y
real(8) :: AI,BI,AD,BD

if(abs(p%t-p%pltid*p%t2p)>p%dt)return

write(name,'(i2.2)')p%pltid
open(unit=66,file='hw4_'//trim(p%name)//'_'//name//'.plt')


do i = 1, p%n
    if(p%t>0.0)then
        x=p%coef/(p%t)**(1.0d0/3.0d0)
        xx=x*( p%x(i) - dsqrt(p%g*0.3d0) * p%t )
        call AIRYA(xx,AI,BI,AD,BD)
        y=0.5d0*p%A0*p%coef*AI
    else
        y=0.0
    endif
    write(66,*)p%x(i),p%eta%now(i), y
enddo
close(unit=66)
p%pltid=p%pltid+1

end subroutine

SUBROUTINE AIRYA(X,AI,BI,AD,BD)

!       ======================================================
!       Purpose: Compute Airy functions and their derivatives
!       Input:   x  --- Argument of Airy function
!       Output:  AI --- Ai(x)
!                BI --- Bi(x)
!                AD --- Ai'(x)
!                BD --- Bi'(x)
!       Routine called:
!                AJYIK for computing Jv(x), Yv(x), Iv(x) and
!                Kv(x) with v=1/3 and 2/3
!       ======================================================

        IMPLICIT NONE
        REAL(8) :: X,AI,BI,AD,BD
        REAL(8) :: XA,PIR,C1,C2,SR3,Z,XQ
        REAL(8) :: VJ1,VJ2,VY1,VY2,VI1,VI2,VK1,VK2
        XA=DABS(X)
        PIR=0.318309886183891D0
        C1=0.355028053887817D0
        C2=0.258819403792807D0
        SR3=1.732050807568877D0
        Z=XA**1.5/1.5D0
        XQ=DSQRT(XA)
        CALL AJYIK(Z,VJ1,VJ2,VY1,VY2,VI1,VI2,VK1,VK2)
        IF (X.EQ.0.0D0) THEN
           AI=C1
           BI=SR3*C1
           AD=-C2
           BD=SR3*C2
        ELSE IF (X.GT.0.0D0) THEN
           AI=PIR*XQ/SR3*VK1
           BI=XQ*(PIR*VK1+2.0D0/SR3*VI1)
           AD=-XA/SR3*PIR*VK2
           BD=XA*(PIR*VK2+2.0D0/SR3*VI2)
        ELSE
           AI=0.5D0*XQ*(VJ1-VY1/SR3)
           BI=-0.5D0*XQ*(VJ1/SR3+VY1)
           AD=0.5D0*XA*(VJ2+VY2/SR3)
           BD=0.5D0*XA*(VJ2/SR3-VY2)
        ENDIF
        RETURN
        END SUBROUTINE


        SUBROUTINE AJYIK(X,VJ1,VJ2,VY1,VY2,VI1,VI2,VK1,VK2)

!       =======================================================
!       Purpose: Compute Bessel functions Jv(x) and Yv(x),
!                and modified Bessel functions Iv(x) and
!                Kv(x), and their derivatives with v=1/3,2/3
!       Input :  x --- Argument of Jv(x),Yv(x),Iv(x) and
!                      Kv(x) ( x ò 0 )
!       Output:  VJ1 --- J1/3(x)
!                VJ2 --- J2/3(x)
!                VY1 --- Y1/3(x)
!                VY2 --- Y2/3(x)
!                VI1 --- I1/3(x)
!                VI2 --- I2/3(x)
!                VK1 --- K1/3(x)
!                VK2 --- K2/3(x)
!       =======================================================

        IMPLICIT NONE
        REAL(8) :: X,VJ1,VJ2,VY1,VY2,VI1,VI2,VK1,VK2
        REAL(8) :: PI,RP2,GP1,GP2,GN,GN1,GN2,VV0,UU0
        INTEGER :: K0,K,L
        REAL(8) :: VL,VJL,R,A0,VV,PX,RP,QX,RQ,XK,CK,SK,B0,UJ1,UJ2,PV1,PV2
        REAL(8) :: C0,SUM,X2,VIL,VSL

        IF (X.EQ.0.0D0) THEN
           VJ1=0.0D0
           VJ2=0.0D0
           VY1=-1.0D+300
           VY2=1.0D+300
           VI1=0.0D0
           VI2=0.0D0
           VK1=-1.0D+300
           VK2=-1.0D+300
           RETURN
        ENDIF
        PI=3.141592653589793D0
        RP2=.63661977236758D0
        GP1=.892979511569249D0
        GP2=.902745292950934D0
        GN1=1.3541179394264D0
        GN2=2.678938534707747D0
        VV0=0.444444444444444D0
        UU0=1.1547005383793D0
        X2=X*X
        K0=12
        IF (X.GE.35.0) K0=10
        IF (X.GE.50.0) K0=8
        IF (X.LE.12.0) THEN
           DO 25 L=1,2
              VL=L/3.0D0
              VJL=1.0D0
              R=1.0D0
              DO 15 K=1,40
                 R=-0.25D0*R*X2/(K*(K+VL))
                 VJL=VJL+R
                 IF (DABS(R).LT.1.0D-15) GO TO 20
15            CONTINUE
20            A0=(0.5D0*X)**VL
              IF (L.EQ.1) VJ1=A0/GP1*VJL
              IF (L.EQ.2) VJ2=A0/GP2*VJL
25         CONTINUE
        ELSE
           DO 40 L=1,2
              VV=VV0*L*L
              PX=1.0D0
              RP=1.0D0
              DO 30 K=1,K0
                 RP=-0.78125D-2*RP*(VV-(4.0*K-3.0)**2.0)*(VV-  &
                    (4.0*K-1.0)**2.0)/(K*(2.0*K-1.0)*X2)
30               PX=PX+RP
              QX=1.0D0
              RQ=1.0D0
              DO 35 K=1,K0
                 RQ=-0.78125D-2*RQ*(VV-(4.0*K-1.0)**2.0)*(VV-  &
                    (4.0*K+1.0)**2.0)/(K*(2.0*K+1.0)*X2)
35               QX=QX+RQ
              QX=0.125D0*(VV-1.0)*QX/X
              XK=X-(0.5D0*L/3.0D0+0.25D0)*PI
              A0=DSQRT(RP2/X)
              CK=DCOS(XK)
              SK=DSIN(XK)
              IF (L.EQ.1) THEN
                 VJ1=A0*(PX*CK-QX*SK)
                 VY1=A0*(PX*SK+QX*CK)
              ELSE IF (L.EQ.2) THEN
                 VJ2=A0*(PX*CK-QX*SK)
                 VY2=A0*(PX*SK+QX*CK)
              ENDIF
40         CONTINUE
        ENDIF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF (X.LE.12.0D0) THEN
           DO 55 L=1,2
              VL=L/3.0D0
              VJL=1.0D0
              R=1.0D0
              DO 45 K=1,40
                 R=-0.25D0*R*X2/(K*(K-VL))
                 VJL=VJL+R
                 IF (DABS(R).LT.1.0D-15) GO TO 50
45            CONTINUE
50            B0=(2.0D0/X)**VL
              IF (L.EQ.1) UJ1=B0*VJL/GN1
              IF (L.EQ.2) UJ2=B0*VJL/GN2
55         CONTINUE
           PV1=PI/3.0D0
           PV2=PI/1.5D0
           VY1=UU0*(VJ1*DCOS(PV1)-UJ1)
           VY2=UU0*(VJ2*DCOS(PV2)-UJ2)
        ENDIF
        IF (X.LE.18.0) THEN
           DO 70 L=1,2
              VL=L/3.0D0
              VIL=1.0D0
              R=1.0D0
              DO 60 K=1,40
                 R=0.25D0*R*X2/(K*(K+VL))
                 VIL=VIL+R
                 IF (DABS(R).LT.1.0D-15) GO TO 65
60            CONTINUE
65            A0=(0.5D0*X)**VL
              IF (L.EQ.1) VI1=A0/GP1*VIL
              IF (L.EQ.2) VI2=A0/GP2*VIL
70         CONTINUE
        ELSE
           C0=DEXP(X)/DSQRT(2.0D0*PI*X)
           DO 80 L=1,2
              VV=VV0*L*L
              VSL=1.0D0
              R=1.0D0
              DO 75 K=1,K0
                 R=-0.125D0*R*(VV-(2.0D0*K-1.0D0)**2.0)/(K*X)
75               VSL=VSL+R
              IF (L.EQ.1) VI1=C0*VSL
              IF (L.EQ.2) VI2=C0*VSL
80         CONTINUE
        ENDIF
        IF (X.LE.9.0D0) THEN
           DO 95 L=1,2
              VL=L/3.0D0
               IF (L.EQ.1) GN=GN1
               IF (L.EQ.2) GN=GN2
               A0=(2.0D0/X)**VL/GN
               SUM=1.0D0
               R=1.0D0
               DO 85 K=1,60
                  R=0.25D0*R*X2/(K*(K-VL))
                  SUM=SUM+R
                  IF (DABS(R).LT.1.0D-15) GO TO 90
85             CONTINUE
90            IF (L.EQ.1) VK1=0.5D0*UU0*PI*(SUM*A0-VI1)
              IF (L.EQ.2) VK2=0.5D0*UU0*PI*(SUM*A0-VI2)
95         CONTINUE
        ELSE
           C0=DEXP(-X)*DSQRT(0.5D0*PI/X)
           DO 105 L=1,2
              VV=VV0*L*L
              SUM=1.0D0
              R=1.0D0
              DO 100 K=1,K0
                 R=0.125D0*R*(VV-(2.0*K-1.0)**2.0)/(K*X)
100              SUM=SUM+R
              IF (L.EQ.1) VK1=C0*SUM
              IF (L.EQ.2) VK2=C0*SUM
105        CONTINUE
        ENDIF
        RETURN
        END