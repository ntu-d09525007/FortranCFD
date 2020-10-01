module roots
implicit none

type dum_matrices
integer :: is, ie, js ,je, ks, ke
real(8),dimension(:,:,:),allocatable :: r,l,f,b,u,d,c,src
contains 
procedure alloc => dum_matrices_alloc
end type dum_matrices

type time_recorded
integer :: is, ie, js ,je, ks, ke
real(8),dimension(:,:,:),allocatable :: now, old, old2, tmp
contains
procedure alloc => time_recorded_alloc
procedure init => time_recorded_init
procedure switch => time_recorded_switch
end type time_recorded

type time_recorded_derivatives
integer :: is, ie, js ,je, ks, ke
type(time_recorded) :: x, y, z, xx, yy, zz, xy, xz, yz, curv
contains
procedure alloc => time_recorded_derivatives_alloc
procedure switch => time_recorded_derivatives_switch
procedure init => time_recorded_derivatives_init
end type time_recorded_derivatives

type time_recorded_vec
integer :: is, ie, js ,je, ks, ke
type(time_recorded) :: x, y, z
contains
procedure alloc => time_recorded_vec_alloc
procedure init => time_recorded_vec_init
procedure switch => time_recorded_vec_switch
end type time_recorded_vec

!====================================================================================!
contains 
!====================================================================================!

subroutine dum_matrices_alloc(p,is,ie,js,je,ks,ke)
implicit none
class(dum_matrices) :: p
integer, intent(in) :: is,ie,js,je,ks,ke

p%is = is; p%ie = ie
p%js = js; p%je = je
p%ks = ks; p%ke = ke

allocate( p%r(is:ie,js:je,ks:ke), p%l(is:ie,js:je,ks:ke), p%f(is:ie,js:je,ks:ke), &
		& p%b(is:ie,js:je,ks:ke), p%c(is:ie,js:je,ks:ke), p%src(is:ie,js:je,ks:ke), &
		& p%u(is:ie,js:je,ks:ke), p%d(is:ie,js:je,ks:ke) )

end subroutine 

subroutine time_recorded_alloc(p,is,ie,js,je,ks,ke)
implicit none
class(time_recorded) :: p
integer, intent(in) :: is,ie,js,je,ks,ke

p%is = is; p%ie = ie
p%js = js; p%je = je
p%ks = ks; p%ke = ke 

allocate( p%now(is:ie,js:je,ks:ke), p%old(is:ie,js:je,ks:ke), p%old2(is:ie,js:je,ks:ke), p%tmp(is:ie,js:je,ks:ke) )

end subroutine

subroutine time_recorded_derivatives_alloc(p,is,ie,js,je,ks,ke)
implicit none
class(time_recorded_derivatives) :: p
integer, intent(in) :: is,ie,js,je,ks,ke

p%is = is; p%ie = ie
p%js = js; p%je = je
p%ks = ks; p%ke = ke 

call p%x%alloc(is,ie,js,je,ks,ke)
call p%y%alloc(is,ie,js,je,ks,ke)
call p%z%alloc(is,ie,js,je,ks,ke)

call p%xx%alloc(is,ie,js,je,ks,ke)
call p%yy%alloc(is,ie,js,je,ks,ke)
call p%zz%alloc(is,ie,js,je,ks,ke)

call p%xy%alloc(is,ie,js,je,ks,ke)
call p%xz%alloc(is,ie,js,je,ks,ke)
call p%yz%alloc(is,ie,js,je,ks,ke)

call p%curv%alloc(is,ie,js,je,ks,ke)

end subroutine

subroutine time_recorded_vec_alloc(p,is,ie,js,je,ks,ke)
implicit none
class(time_recorded_vec) :: p
integer, intent(in) :: is,ie,js,je,ks,ke

p%is = is; p%ie = ie
p%js = js; p%je = je
p%ks = ks; p%ke = ke 

call p%x%alloc(is,ie,js,je,ks,ke)
call p%y%alloc(is,ie,js,je,ks,ke)
call p%z%alloc(is,ie,js,je,ks,ke)

end subroutine

subroutine time_recorded_init(p)
implicit none
integer :: i,j,k
class(time_recorded) :: p

 do k = p%ks, p%ke
 do j = p%js, p%je
 do i = p%is, p%ie
	p%old(i,j,k) = p%now(i,j,k)
	p%old2(i,j,k) = p%now(i,j,k)
	p%tmp(i,j,k) = p%now(i,j,k)
 end do
 end do 
 end do
 
end subroutine

subroutine time_recorded_derivatives_init(p)
implicit none
class(time_recorded_derivatives) :: p

	call p%x%init()
	call p%y%init()
	call p%z%init()
	
	call p%xx%init()
	call p%yy%init()
	call p%zz%init()
	
	call p%xy%init()
	call p%xz%init()
	call p%yz%init()
	
	call p%curv%init()

end subroutine

subroutine time_recorded_vec_init(p)
implicit none
class(time_recorded_vec) :: p

 call p%x%init
 call p%y%init
 call p%z%init
 
end subroutine

subroutine time_recorded_switch(p)
implicit none
integer :: i,j,k
class(time_recorded) :: p

 do k = p%ks, p%ke
 do j = p%js, p%je
 do i = p%is, p%ie
	p%old2(i,j,k) = p%old(i,j,k)
	p%old(i,j,k) = p%now(i,j,k)
 end do
 end do 
 end do
 
end subroutine

subroutine time_recorded_derivatives_switch(p)
implicit none
class(time_recorded_derivatives) :: p

	call p%x%switch
	call p%y%switch
	call p%z%switch
	
	call p%xx%switch
	call p%yy%switch
	call p%zz%switch
	
	call p%xy%switch
	call p%xz%switch
	call p%yz%switch
	
	call p%curv%switch
	
end subroutine

subroutine time_recorded_vec_switch(p)
implicit none
class(time_recorded_vec) :: p

  call p%x%switch
  call p%y%switch
  call p%z%switch

end subroutine

end module roots

