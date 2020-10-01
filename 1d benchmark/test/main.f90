include 'my_data.f90'

program main
use my_precision
use mydata
implicit none
integer,parameter :: num_of_test=4
integer :: i, id
type(pure_adv),allocatable :: dat(:)

allocate(dat(num_of_test))

w=0.5_ap
tol=1.0e-14_ap

IsPloting = .true.
IsShowing = .false.


do i = 1, num_of_test
  call dat(i)%init(1,i)
  call go(dat(i))

  if(i>1)then
          write(*,'(I5,2ES15.4,F8.3)')dat(i)%nx,dat(i)%dt,dat(i)%norm(2),log(dat(i)%norm(2)/dat(i-1)%norm(2))/log(dat(i)%dt/dat(i-1)%dt)
  else
          write(*,'(I5,2ES15.4)')dat(i)%nx,dat(i)%dt,dat(i)%norm(2)
  endif
  
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
	 if(p%t>p%t2s)exit
	 
 enddo
 
 p%cpu_time = p%cpu_time + wtime()
 
 call p%set(p%u_exact,p%c,p%t-p%t2s)
 call p%find_error(p%u,p%u_exact,p%norm(1),p%norm(2),p%norm(3))

end subroutine

end program



