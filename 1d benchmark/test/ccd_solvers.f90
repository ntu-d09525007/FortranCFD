module CCD_SOLVERS
USE MATRIX_SOLVER
implicit none

type working_matrix
integer :: is, ie
real(kind=ap) :: dx, dt
real(kind=ap),allocatable,dimension(:)   :: c3, SU, SSU, SD, SSD
real(kind=ap),ALLOCATABLE,DIMENSION(:,:) :: AU, AAU, BU, BBU
real(kind=ap),ALLOCATABLE,DIMENSION(:,:) :: AD, AAD, BD, BBD
contains
procedure alloc => working_matrix_alloc
procedure bc => working_matrix_bc
procedure decom => working_matrix_decom
procedure solve => working_matrix_solve
procedure init_srkccd_c3 => working_matrix_init_srkccd_c3
procedure init_srkccd => working_matrix_init_srkccd
procedure init_uccd => working_matrix_init_uccd
procedure init_ccd => working_matrix_init_ccd
procedure asgn_src_ccd => working_matrix_asgn_src_ccd
procedure asgn_src_uccd => working_matrix_asgn_src_uccd
procedure asgn_src_srkccd => working_matrix_asgn_src_srkccd
end type working_matrix

type ccd_root
integer :: is, ie
real(kind=ap) :: dx, dt
contains
procedure solve => root_solve
end type ccd_root

type ccd_manager
type(ccd_root) :: x,y,z
end type ccd_manager

contains

subroutine working_matrix_alloc(p,is,ie,dx,dt)
implicit none
class(working_matrix) :: p
integer :: is, ie
real(kind=ap) :: dx, dt

p%is = is
p%ie = ie
p%dx = dx
p%dt = dt

allocate(p%c3(is:ie))
allocate(p%AU(3,is:ie),p%AAU(3,is:ie),p%BU(3,is:ie),p%BBU(3,is:ie))
allocate(p%AD(3,is:ie),p%AAD(3,is:ie),p%BD(3,is:ie),p%BBD(3,is:ie))
allocate(p%SU(is:ie),p%SSU(is:ie),p%SD(is:ie),p%SSD(is:ie))

end subroutine

subroutine working_matrix_decom(p,twice)
implicit none
class(working_matrix) :: p
logical :: twice

call TWIN_DEC(p%AU, p%BU, p%AAU, p%BBU, p%is, p%ie)
if(twice)call TWIN_DEC(p%AD, p%BD, p%AAD, p%BBD, p%is, p%ie)

end subroutine

subroutine working_matrix_solve(p,twice)
implicit none
class(working_matrix) :: p
logical :: twice

call TWIN_BKS(p%au,p%bu,p%aau,p%bbu,p%su,p%ssu,p%is,p%ie)
if(twice)call TWIN_BKS(p%ad,p%bd,p%aad,p%bbd,p%sd,p%ssd,p%is,p%ie)

end subroutine

subroutine working_matrix_bc(p)
implicit none
class(working_matrix) :: p

p%AU(2,p%is) = 1.0_ap
p%BU(2,p%is) = 0.0_ap

p%AAU(2,p%is) = 0.0_ap
p%BBU(2,p%is) = 1.0_ap

p%AU(3,p%is) = 2.0_ap
p%BU(3,p%is) = - p%dx

p%AAU(3,p%is) = -2.5_ap / p%dx
p%BBU(3,p%is) = 8.5_ap

p%AU(2,p%ie) = 1.0_ap
p%BU(2,p%ie) = 0.0_ap

p%AAU(2,p%ie) = 0.0_ap
p%BBU(2,p%ie) = 1.0_ap

p%AU(1,p%ie) = 2.0_ap
p%BU(1,p%ie) = p%dx

p%AAU(1,p%ie) = 2.5_ap / p%dx
p%BBU(1,p%ie) = 8.5_ap

!---------------------------------------

p%AD(2,p%is) = 1.0_ap
p%BD(2,p%is) = 0.0_ap

p%AAD(2,p%is) = 0.0_ap
p%BBD(2,p%is) = 1.0_ap

p%AD(3,p%is) = 2.0_ap
p%BD(3,p%is) = - p%dx

p%AAD(3,p%is) = -2.5_ap / p%dx
p%BBD(3,p%is) = 8.5_ap

p%AD(2,p%ie) = 1.0_ap
p%BD(2,p%ie) = 0.0_ap

p%AAD(2,p%ie) = 0.0_ap
p%BBD(2,p%ie) = 1.0_ap

p%AD(1,p%ie) = 2.0_ap
p%BD(1,p%ie) = p%dx

p%AAD(1,p%ie) = 2.5_ap / p%dx
p%BBD(1,p%ie) = 8.5_ap

end subroutine

subroutine working_matrix_init_ccd(p)
implicit none
class(working_matrix) :: p
integer :: i

call p%bc

do i = p%is+1, p%ie-1

    p%AU(1,i) = 7.0_ap / 16.0_ap
    p%AU(2,i) = 1.0_ap
    p%AU(3,i) = 7.0_ap / 16.0_ap

    p%BU(1,i) = p%dx / 16.0_ap
    p%BU(2,i) = 0.0_ap
    p%BU(3,i) = - p%dx / 16.0_ap

    p%AAU(1,i) = - 9.0_ap / 8.0_ap / p%dx
    p%AAU(2,i) = 0.0_ap
    p%AAU(3,i) =   9.0_ap / 8.0_ap / p%dx

    p%BBU(1,i) = - 1.0_ap / 8.0_ap
    p%BBU(2,i) = 1.0_ap
    p%BBU(3,i) = - 1.0_ap / 8.0_ap

enddo

end subroutine

subroutine working_matrix_init_uccd(p)
implicit none
class(working_matrix) :: p
integer :: i
real(kind=ap) :: a1,b1,b2,b3

 A1 =  0.875_AP
 B1 =  0.1251282341599089_AP
 B2 = -0.2487176584009104_AP
 B3 =  0.0001282341599089_AP

call p%bc

do i = p%is+1, p%ie-1

    p%AU(1,i) = A1
    p%AU(2,i) = 1.0_ap
    p%AU(3,i) = 0.0_ap

    p%BU(1,i) = B1*p%dx
    p%BU(2,i) = B2*p%dx
    p%BU(3,i) = B3*p%dx

    p%AAU(1,i) = - 9.0_ap / 8.0_ap / p%dx
    p%AAU(2,i) = 0.0_ap
    p%AAU(3,i) =   9.0_ap / 8.0_ap / p%dx

    p%BBU(1,i) = - 1.0_ap / 8.0_ap
    p%BBU(2,i) = 1.0_ap
    p%BBU(3,i) = - 1.0_ap / 8.0_ap

    p%AD(1,i) = 0.0_ap
    p%AD(2,i) = 1.0_ap
    p%AD(3,i) = A1

    p%BD(1,i) = -B3*p%dx
    p%BD(2,i) = -B2*p%dx
    p%BD(3,i) = -B1*p%dx

    p%AAD(1,i) = - 9.0_ap / 8.0_ap / p%dx
    p%AAD(2,i) = 0.0_ap
    p%AAD(3,i) =   9.0_ap / 8.0_ap / p%dx

    p%BBD(1,i) = - 1.0_ap / 8.0_ap
    p%BBD(2,i) = 1.0_ap
    p%BBD(3,i) = - 1.0_ap / 8.0_ap

enddo

end subroutine

subroutine working_matrix_init_srkccd(p)
implicit none
class(working_matrix) :: p
integer :: i
real(kind=ap) :: a1,a3,b1,b3,c3

call p%bc

do i = p%is+1, p%ie-1

    c3 = p%c3(i)

    A1 =  131.0_ap / 128.0_ap - 5.0_ap * C3 / 8.0_ap
    A3 = - 19.0_ap / 128.0_ap + 5.0_ap * C3 / 8.0_ap
    
    B1 = 23.0_ap / 128.0_ap - C3 / 8.0_ap
    B3 =  7.0_ap / 128.0_ap - C3 / 8.0_ap

    p%AU(1,i) = A1
    p%AU(2,i) = 1.0_ap
    p%AU(3,i) = A3

    p%BU(1,i) = B1*p%dx
    p%BU(2,i) = 0.0_ap
    p%BU(3,i) = B3*p%dx

    p%AD(1,i) = A3
    p%AD(2,i) = 1.0_ap
    p%AD(3,i) = A1

    p%BD(1,i) = -B3*p%dx
    p%BD(2,i) = 0.0_ap
    p%BD(3,i) = -B1*p%dx

    p%AAU(1,i) = - 9.0_ap / 8.0_ap / p%dx
    p%AAU(2,i) = 0.0_ap
    p%AAU(3,i) =   9.0_ap / 8.0_ap / p%dx

    p%BBU(1,i) = - 1.0_ap / 8.0_ap
    p%BBU(2,i) = 1.0_ap
    p%BBU(3,i) = - 1.0_ap / 8.0_ap

    p%AAD(1,i) = - 9.0_ap / 8.0_ap / p%dx
    p%AAD(2,i) = 0.0_ap
    p%AAD(3,i) =   9.0_ap / 8.0_ap / p%dx

    p%BBD(1,i) = - 1.0_ap / 8.0_ap
    p%BBD(2,i) = 1.0_ap
    p%BBD(3,i) = - 1.0_ap / 8.0_ap


enddo

end subroutine

subroutine working_matrix_init_srkccd_c3(p,u)
implicit none
class(working_matrix) :: p
real(kind=ap), dimension(p%is:p%ie), intent(in) :: u
integer :: i
real(kind=ap) :: a,b,c,d,e,x

 A= 0.0001671_ap
 B= 7.43943_ap
 C= 0.840798_ap
 D= 0.00084915_ap
 E= -0.159194_ap

do i = p%is, p%ie
    x=abs(u(i))*p%dt/p%dx
    p%c3(i) = A*x**B+C*x**D+E
enddo

end subroutine

subroutine working_matrix_asgn_src_ccd(p,f)
implicit none
class(working_matrix) :: p
integer :: i
real(kind=ap), dimension(p%is:p%ie) :: f

p%SU(p%is)  = (-3.5_ap*f(p%is)+4.0_ap*f(p%is+1)-0.5_ap*f(p%is+2))/p%dx
P%SSU(p%is) = (34.0_ap/3.0_ap*f(p%is)-83.0_ap/4.0_ap*f(p%is+1)+10.0_ap*f(p%is+2)-7.0_ap/12.0_ap*f(p%is+3))/p%dx**2.0_ap

do i = p%is+1, p%ie-1

    p%su(i) = 15.0_ap/16.0_ap*(-f(i-1)+f(i+1))/p%dx
    p%ssu(i) = (3.0_ap*f(i-1)-6.0_ap*f(i)+3.0_ap*f(i+1))/p%dx**2.0_ap

enddo

p%SU(p%ie)  = -(-3.5_ap*f(p%ie)+4.0_ap*f(p%ie-1)-0.5_ap*f(p%ie-2))/p%dx
P%SSU(p%ie) = (34.0_ap/3.0_ap*f(p%ie)-83.0_ap/4.0_ap*f(p%ie-1)+10.0_ap*f(p%ie-2)-7.0_ap/12.0_ap*f(p%ie-3))/p%dx**2.0_ap

end subroutine

subroutine working_matrix_asgn_src_uccd(p,f)
implicit none
class(working_matrix) :: p
integer :: i
real(kind=ap), dimension(p%is:p%ie) :: f
real(kind=ap) :: c1,c2,c3

p%SU(p%is)  = (-3.5_ap*f(p%is)+4.0_ap*f(p%is+1)-0.5_ap*f(p%is+2))/p%dx
P%SSU(p%is) = (34.0_ap/3.0_ap*f(p%is)-83.0_ap/4.0_ap*f(p%is+1)+10.0_ap*f(p%is+2)-7.0_ap/12.0_ap*f(p%is+3))/p%dx**2.0_ap

p%SD(p%is)  = (-3.5_ap*f(p%is)+4.0_ap*f(p%is+1)-0.5_ap*f(p%is+2))/p%dx
P%SSD(p%is) = (34.0_ap/3.0_ap*f(p%is)-83.0_ap/4.0_ap*f(p%is+1)+10.0_ap*f(p%is+2)-7.0_ap/12.0_ap*f(p%is+3))/p%dx**2.0_ap

do i = p%is+1, p%ie-1

    c3 = -0.06096119008109_AP
    c2 = 1.99692238016218_AP
    c1 = -1.93596119008109_AP

    p%su(i) = (c1*f(i-1)+c2*f(i)+c3*f(i+1))/p%dx
    p%ssu(i) = (3.0_ap*f(i-1)-6.0_ap*f(i)+3.0_ap*f(i+1))/p%dx**2.0_ap

    p%sd(i) = -(c3*f(i-1)+c2*f(i)+c1*f(i+1))/p%dx
    p%ssd(i) = (3.0_ap*f(i-1)-6.0_ap*f(i)+3.0_ap*f(i+1))/p%dx**2.0_ap

enddo

p%SU(p%ie)  = -(-3.5_ap*f(p%ie)+4.0_ap*f(p%ie-1)-0.5_ap*f(p%ie-2))/p%dx
P%SSU(p%ie) = (34.0_ap/3.0_ap*f(p%ie)-83.0_ap/4.0_ap*f(p%ie-1)+10.0_ap*f(p%ie-2)-7.0_ap/12.0_ap*f(p%ie-3))/p%dx**2.0_ap

p%SD(p%ie)  = -(-3.5_ap*f(p%ie)+4.0_ap*f(p%ie-1)-0.5_ap*f(p%ie-2))/p%dx
P%SSD(p%ie) = (34.0_ap/3.0_ap*f(p%ie)-83.0_ap/4.0_ap*f(p%ie-1)+10.0_ap*f(p%ie-2)-7.0_ap/12.0_ap*f(p%ie-3))/p%dx**2.0_ap

end subroutine

subroutine working_matrix_asgn_src_srkccd(p,f)
implicit none
class(working_matrix) :: p
integer :: i
real(kind=ap), dimension(p%is:p%ie) :: f
real(kind=ap) :: c1,c2,c3

p%SU(p%is)  = (-3.5_ap*f(p%is)+4.0_ap*f(p%is+1)-0.5_ap*f(p%is+2))/p%dx
P%SSU(p%is) = (34.0_ap/3.0_ap*f(p%is)-83.0_ap/4.0_ap*f(p%is+1)+10.0_ap*f(p%is+2)-7.0_ap/12.0_ap*f(p%is+3))/p%dx**2.0_ap

p%SD(p%is)  = (-3.5_ap*f(p%is)+4.0_ap*f(p%is+1)-0.5_ap*f(p%is+2))/p%dx
P%SSD(p%is) = (34.0_ap/3.0_ap*f(p%is)-83.0_ap/4.0_ap*f(p%is+1)+10.0_ap*f(p%is+2)-7.0_ap/12.0_ap*f(p%is+3))/p%dx**2.0_ap

do i = p%is+1, p%ie-1

    c3 = p%c3(i)
    c2 = 15.0_ap/8.0_ap - 2.0_ap*c3
    c1 = c3 - 15.0_ap/8.0_ap

    p%su(i) = (c1*f(i-1)+c2*f(i)+c3*f(i+1))/p%dx
    p%ssu(i) = (3.0_ap*f(i-1)-6.0_ap*f(i)+3.0_ap*f(i+1))/p%dx**2.0_ap

    p%sd(i) = -(c3*f(i-1)+c2*f(i)+c1*f(i+1))/p%dx
    p%ssd(i) = (3.0_ap*f(i-1)-6.0_ap*f(i)+3.0_ap*f(i+1))/p%dx**2.0_ap

enddo

p%SU(p%ie)  = -(-3.5_ap*f(p%ie)+4.0_ap*f(p%ie-1)-0.5_ap*f(p%ie-2))/p%dx
P%SSU(p%ie) = (34.0_ap/3.0_ap*f(p%ie)-83.0_ap/4.0_ap*f(p%ie-1)+10.0_ap*f(p%ie-2)-7.0_ap/12.0_ap*f(p%ie-3))/p%dx**2.0_ap

p%SD(p%ie)  = -(-3.5_ap*f(p%ie)+4.0_ap*f(p%ie-1)-0.5_ap*f(p%ie-2))/p%dx
P%SSD(p%ie) = (34.0_ap/3.0_ap*f(p%ie)-83.0_ap/4.0_ap*f(p%ie-1)+10.0_ap*f(p%ie-2)-7.0_ap/12.0_ap*f(p%ie-3))/p%dx**2.0_ap

end subroutine


subroutine root_solve(p,type,f,fx,fxx,u)
implicit none
class(ccd_root) :: p
character(len=*) :: type
real(kind=ap), dimension(p%is:p%ie) :: f,fx,fxx
real(kind=ap), dimension(p%is:p%ie), optional :: u
type(working_matrix) :: work
integer :: i

call work%alloc(p%is,p%ie,p%dx,p%dt)

select case (trim(type))
    case ("ccd")
        call work%init_ccd
        call work%decom(.false.)
        call work%asgn_src_ccd(f)
        call work%solve(.false.)
        do i = p%is, p%ie
            fx(i) = work%su(i)
            fxx(i) = work%ssu(i)
        enddo
    case ("uccd")
        call work%init_uccd
        call work%decom(.true.)
        call work%asgn_src_uccd(f)
        call work%solve(.true.)
        do i = p%is, p%ie
            if( u(i)>0.0_ap )then
                fx(i) = work%su(i)
            else
                fx(i) = work%sd(i)
            endif
            fxx(i) = 0.5_ap*( work%ssu(i)+work%ssd(i) )
        enddo
    case ("srkccd")
        call work%init_srkccd_c3(u)
        call work%init_srkccd
        call work%decom(.true.)
        call work%asgn_src_srkccd(f)
        call work%solve(.true.)
        do i = p%is, p%ie
            if( u(i)>0.0_ap )then
                fx(i) = work%su(i)
            else
                fx(i) = work%sd(i)
            endif
            fxx(i) = 0.5_ap*( work%ssu(i)+work%ssd(i) )
        enddo
    case default
        stop "Something wrong when call CCD-type solvers"
end select

end subroutine

end module