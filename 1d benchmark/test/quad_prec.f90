integer, parameter :: ap = selected_real_kind(p=24,r=1000)

interface sqrt
	module procedure my_sqrt
end interface

interface sin
	module procedure my_sin
end interface

interface cos
	module procedure my_cos
end interface

interface acos
	module procedure my_acos
end interface

interface asin
	module procedure my_asin
end interface

interface log
	module procedure my_log
end interface

contains

function my_sqrt(x) result(y)
real(kind=ap) :: x, y
 y = qsqrt(x)
end function 

function my_sin(x) result(y)
real(kind=ap) :: x, y
 y = qsin(x)
end function 

function my_cos(x) result(y)
real(kind=ap) :: x, y
 y = qcos(x)
end function 

function my_asin(x) result(y)
real(kind=ap) :: x, y
 y = qasin(x)
end function 

function my_acos(x) result(y)
real(kind=ap) :: x, y
 y = qacos(x)
end function 

function my_log(x) result(y)
real(kind=ap) :: x, y
 y = qlog(x)
end function 