module roots
implicit none

type dum_matrices
integer :: is, ie, js ,je
real(8),dimension(:,:),allocatable :: r,l,f,b,c,src
contains 
procedure alloc => dum_matrices_alloc
end type dum_matrices

type time_recorded
integer :: is, ie, js ,je
real(8),dimension(:,:),allocatable :: now, old, old2, tmp
contains
procedure alloc => time_recorded_alloc
procedure init => time_recorded_init
procedure switch => time_recorded_switch
end type time_recorded

type time_recorded_derivatives
integer :: is, ie, js ,je
type(time_recorded) :: x, y, xx, yy, xy, curv
contains
procedure alloc => time_recorded_derivatives_alloc
procedure switch => time_recorded_derivatives_switch
procedure init => time_recorded_derivatives_init
end type time_recorded_derivatives

type time_recorded_vec
integer :: is, ie, js ,je
type(time_recorded) :: x, y
contains
procedure alloc => time_recorded_vec_alloc
procedure init => time_recorded_vec_init
procedure switch => time_recorded_vec_switch
end type time_recorded_vec

!====================================================================================!
contains 
!====================================================================================!


subroutine dum_matrices_alloc(p,is,ie,js,je)
implicit none
class(dum_matrices) :: p
integer, intent(in) :: is,ie,js,jese

p%is = is; p%ie = ie
p%js = js; p%je = je

allocate( p%r(is:ie,js:je), p%l(is:ie,js:je), p%f(is:ie,js:je), &
        & p%b(is:ie,js:je), p%c(is:ie,js:je), p%src(is:ie,js:je) )

end subroutine 

subroutine time_recorded_alloc(p,is,ie,js,je)
implicit none
class(time_recorded) :: p
integer, intent(in) :: is,ie,js,je

p%is = is; p%ie = ie
p%js = js; p%je = je

allocate( p%now(is:ie,js:je), p%old(is:ie,js:je), p%old2(is:ie,js:je), p%tmp(is:ie,js:je) )

end subroutine

subroutine time_recorded_derivatives_alloc(p,is,ie,js,je)
implicit none
class(time_recorded_derivatives) :: p
integer, intent(in) :: is,ie,js,je

p%is = is; p%ie = ie
p%js = js; p%je = je

call p%x%alloc(is,ie,js,je)
call p%y%alloc(is,ie,js,je)
call p%z%alloc(is,ie,js,je)

call p%xx%alloc(is,ie,js,je)
call p%yy%alloc(is,ie,js,je)
call p%zz%alloc(is,ie,js,je)

call p%xy%alloc(is,ie,js,je)
call p%xz%alloc(is,ie,js,je)
call p%yz%alloc(is,ie,js,je)

call p%curv%alloc(is,ie,js,je)

end subroutine

subroutine time_recorded_vec_alloc(p,is,ie,js,je)
implicit none
class(time_recorded_vec) :: p
integer, intent(in) :: is,ie,js,je

p%is = is; p%ie = ie
p%js = js; p%je = je

call p%x%alloc(is,ie,js,je)
call p%y%alloc(is,ie,js,je)
call p%z%alloc(is,ie,js,je)

end subroutine

subroutine time_recorded_init(p)
implicit none
integer :: i,j
class(time_recorded) :: p

 do j = p%js, p%je
 do i = p%is, p%ie
    p%old(i,j) = p%now(i,j)
    p%old2(i,j) = p%now(i,j)
    p%tmp(i,j) = p%now(i,j)
 end do
 end do 
 
end subroutine

subroutine time_recorded_derivatives_init(p)
implicit none
class(time_recorded_derivatives) :: p

call p%x%init()
call p%y%init()

call p%xx%init()
call p%yy%init()

call p%xy%init()

call p%curv%init()

end subroutine

subroutine time_recorded_vec_init(p)
implicit none
class(time_recorded_vec) :: p

 call p%x%init
 call p%y%init
 
end subroutine

subroutine time_recorded_switch(p)
implicit none
integer :: i,j
class(time_recorded) :: p

do j = p%js, p%je
do i = p%is, p%ie
  p%old2(i,j) = p%old(i,j)
  p%old(i,j) = p%now(i,j)
end do
end do
 
end subroutine

subroutine time_recorded_derivatives_switch(p)
implicit none
class(time_recorded_derivatives) :: p

call p%x%switch
call p%y%switch

call p%xx%switch
call p%yy%switch

call p%xy%switch

call p%curv%switch
    
end subroutine

subroutine time_recorded_vec_switch(p)
implicit none
class(time_recorded_vec) :: p

call p%x%switch
call p%y%switch

end subroutine

end module roots

