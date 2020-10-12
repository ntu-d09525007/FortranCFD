subroutine tsolver_roots_solve_srk4(p,err)
implicit none
class(tsolver_roots) :: p
integer :: i,j
real(8),intent(out) :: err
real(8) :: coef1_1, coef1_2
real(8) :: coef2_1, coef2_2
real(8) :: cc

cc = 1.0d0 / (2.0d0*dsqrt(3.0d0))

coef1_1 = 0.25d0; coef1_2 = 0.25d0 + cc
coef2_1 = 0.25d0 - cc; coef2_2 = 0.25d0

 err = 0.0_8

 do j = p%js, p%je
 do i = p%is, p%ie
 
    p%ss1(i,j) = p%s1(i,j)
    p%ss2(i,j) = p%s2(i,j)
    
    p%s1(i,j) = p%target(i,j) + p%dt*( coef1_1*p%l1(i,j) + coef1_2*p%l2(i,j)  )
    p%s2(i,j) = p%target(i,j) + p%dt*( coef2_1*p%l1(i,j) + coef2_2*p%l2(i,j)  )

    p%s1(i,j) = p%w * p%s1(i,j) + (1.0_8-p%w) * p%ss1(i,j) 
    p%s2(i,j) = p%w * p%s2(i,j) + (1.0_8-p%w) * p%ss2(i,j) 
    
 enddo
 enddo

 do j = p%js+p%ghc, p%je-p%ghc
 do i = p%is+p%ghc, p%ie-p%ghc
    
    err = max( err, abs( p%s1(i,j)-p%ss1(i,j) ), abs( p%s2(i,j)-p%ss2(i,j) ) )
    
 enddo
 enddo

end subroutine 

subroutine tsolver_data_solve_srk4(p,err)
implicit none
class(tsolver_data) :: p
real(8), intent(out) :: err
real(8) :: err1, err2

 call p%x%solve_srk4(err1)

 if( p%is_vector_solver )then
    call p%y%solve_srk4(err2)
    err = max( err1, err2 )
 else
    err = err1
 endif

end subroutine 

subroutine tsolver_roots_final_srk4(p)
implicit none
class(tsolver_roots) :: p
integer :: i,j

 do j = p%js, p%je
 do i = p%is, p%ie
     p%target(i,j) = p%target(i,j) + p%dt/2.0d0 * ( p%l1(i,j) + p%l2(i,j) )
 enddo
 enddo

end subroutine

subroutine tsolver_data_final_srk4(p)
implicit none
class(tsolver_data) :: p
integer :: i,j

 call p%x%final_srk4

 if( p%is_vector_solver )then
     call p%y%final_srk4
 endif
 
 p%is_vector_solver = .false.
 
 end subroutine