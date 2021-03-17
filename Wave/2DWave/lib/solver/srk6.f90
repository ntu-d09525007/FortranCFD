subroutine tsolver_roots_solve_srk6(p,err)
implicit none
class(tsolver_roots) :: p
integer :: i,j
real(8),intent(out) :: err
real(8) :: coef1_1, coef1_2, coef1_3
real(8) :: coef2_1, coef2_2, coef2_3
real(8) :: coef3_1, coef3_2, coef3_3
real(8) :: cc

cc = dsqrt(0.6_8)*0.5_8

coef1_1 = 5.0_8/36.0_8; coef1_2 = 2.0_8/9.0_8 + 2.0_8*cc/3.0_8; coef1_3 = 5.0_8/36.0_8 + cc/3.0_8
coef2_1 = 5.0_8/36.0_8 - 5.0_8*cc/12.0_8; coef2_2 = 2.0_8/9.0_8; coef2_3 = 5.0_8/36.0_8 + 5.0_8*cc/12.0_8 
coef3_1 = 5.0_8/36.0_8 - cc/3.0_8; coef3_2 = 2.0_8/9.0_8 - 2.0_8*cc/3.0_8; coef3_3 = 5.0_8/36.0_8

 err = 0.0_8

 do j = p%js, p%je
 do i = p%is, p%ie
 
    p%ss1(i,j) = p%s1(i,j)
    p%ss2(i,j) = p%s2(i,j)
    p%ss3(i,j) = p%s3(i,j)
    
    p%s1(i,j) = p%target(i,j) + p%dt*( coef1_1*p%l1(i,j) + coef1_2*p%l2(i,j) + coef1_3*p%l3(i,j) )
    p%s2(i,j) = p%target(i,j) + p%dt*( coef2_1*p%l1(i,j) + coef2_2*p%l2(i,j) + coef2_3*p%l3(i,j) )
    p%s3(i,j) = p%target(i,j) + p%dt*( coef3_1*p%l1(i,j) + coef3_2*p%l2(i,j) + coef3_3*p%l3(i,j) )
    
    p%s1(i,j) = p%w * p%s1(i,j) + (1.0_8-p%w) * p%ss1(i,j) 
    p%s2(i,j) = p%w * p%s2(i,j) + (1.0_8-p%w) * p%ss2(i,j) 
    p%s3(i,j) = p%w * p%s3(i,j) + (1.0_8-p%w) * p%ss3(i,j) 
    
 enddo
 enddo

 do j = p%js+p%ghc, p%je-p%ghc
 do i = p%is+p%ghc, p%ie-p%ghc
    
    err = max( err, abs( p%s1(i,j)-p%ss1(i,j) ), abs( p%s2(i,j)-p%ss2(i,j) ), abs( p%s3(i,j)-p%ss3(i,j) ) )
    
 enddo
 enddo
 
 
end subroutine

subroutine tsolver_data_solve_srk6(p,check)
implicit none
class(tsolver_data) :: p
logical, intent(out) :: check
real(8) :: err1, err2

p%error2 = p%error

call p%x%solve_srk6(err1)

if( p%is_vector_solver )then
    call p%y%solve_srk6(err2)
    p%error = max( err1, err2 )
else
    p%error = err1
endif

check = .false.

if( p%error / p%error2 < 0.9 )then
    if(p%error < 1.0d-14) check = .true.
else
    if(p%error < p%tol) check = .true.
endif

end subroutine

subroutine tsolver_roots_final_srk6(p)
implicit none
class(tsolver_roots) :: p
integer :: i,j

 do j = p%js, p%je
 do i = p%is, p%ie
     p%target(i,j) = p%target(i,j) + p%dt/18.0_8 * ( 5.0_8*p%l1(i,j) + 8.0_8*p%l2(i,j) + 5.0_8*p%l3(i,j) )
 enddo
 enddo

end subroutine

subroutine tsolver_data_final_srk6(p)
implicit none
class(tsolver_data) :: p

 call p%x%final_srk6

 if( p%is_vector_solver )then
     call p%y%final_srk6
 endif
 
 p%is_vector_solver = .false.
 
 end subroutine