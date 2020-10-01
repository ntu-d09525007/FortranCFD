include 'my_data.f90'

program main
use my_precision
use mydata
implicit none
integer,parameter :: num_of_test=2
integer :: i
type(pure_adv),allocatable :: dat(:)

allocate(dat(num_of_test))

w=0.5_ap
tol=1.0e-14_ap

IsPloting = .true.
IsShowing = .false.

do i = 1, num_of_test
  call dat(i)%init(i)
  call go(dat(i))
enddo


contains 

function wtime ( )
  use my_precision
  implicit none
  integer ( kind = 4 ) clock_max
  integer ( kind = 4 ) clock_rate
  integer ( kind = 4 ) clock_reading
  real ( kind = ap ) wtime
  call system_clock ( clock_reading, clock_rate, clock_max )
  wtime = real ( clock_reading, kind = ap ) &
        / real ( clock_rate, kind = ap )
  return
end function

subroutine go(p)
implicit none
type(pure_adv) :: p

 p%cpu_time = - wtime()
 
 do
	 
	 call p%find_dt
	 p%t = p%t + p%dt
	 call p%srk6
	 call p%plot
	 if(p%t>=p%t2s)exit
	 
 enddo
 
 p%cpu_time = p%cpu_time + wtime()
 
end subroutine

end program



