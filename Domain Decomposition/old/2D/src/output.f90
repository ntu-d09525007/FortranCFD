subroutine output()
use all
implicit none

	! level set method, loss of volume/mass in percentage
	write(p%fil%ls_mv,*)p%glb%time,100.0d0*(p%glb%mass-p%glb%imass)/p%glb%imass,100.0d0*(p%glb%vol-p%glb%ivol)/p%glb%ivol

end subroutine