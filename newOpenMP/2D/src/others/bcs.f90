subroutine bc(dat)
use all
implicit none
real(8),dimension(p%loc%is-p%glb%ghc:p%loc%ie+p%glb%ghc,&
                  p%loc%js-p%glb%ghc:p%loc%je+p%glb%ghc) :: dat
integer :: i,j

    !$omp parallel do 
    do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
    do i = 1, p%glb%ghc
        dat(p%loc%is-i,j) = dat(p%loc%is,j)
        dat(p%loc%ie+i,j) = dat(p%loc%ie,j)
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel do   
    do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
    do j = 1, p%glb%ghc
        dat(i,p%loc%js-j) = dat(i,p%loc%js)
        dat(i,p%loc%je+j) = dat(i,p%loc%je)
    enddo
    enddo
    !$omp end parallel do

end subroutine

subroutine velbc(u,v)
! doi.org/10.1063/1.1761178
use all
implicit none
real(8),dimension(p%loc%is-p%glb%ghc:p%loc%ie+p%glb%ghc,&
                  p%loc%js-p%glb%ghc:p%loc%je+p%glb%ghc) :: u,v
integer :: i,j,k
!==========================================
!  X-direction velocity boundary condition
!==========================================

if( p%glb%ubc(1) == 1 )then

    !$omp parallel do
    do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
    do i = 1, p%glb%ghc
        u(p%loc%is-i,j) = - u(p%loc%is-2+i,j)
        v(p%loc%is-i,j) = - v(p%loc%is-1+i,j)
    end do  
    u(p%loc%is-1,j) = 0.0d0
    end do
    !$omp end parallel do

else if ( p%glb%ubc(1) == 2)then

    !$omp parallel do 
    do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
    do i = 1, p%glb%ghc
        u(p%loc%is-i,j) = u(p%loc%is-2+i,j)
        v(p%loc%is-i,j) = v(p%loc%is-1+i,j)
    end do  
    u(p%loc%is-1,j) = 0.0d0
    end do
    !$omp end parallel do
            
endif

 if( p%glb%ubc(2) == 1 )then

    !$omp parallel do
    do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
    do i = 1, p%glb%ghc
        u(p%loc%ie+i,j) = - u(p%loc%ie-i,j)
        v(p%loc%ie+i,j) = - v(p%loc%ie+1-i,j)
        w(p%loc%ie+i,j) = - w(p%loc%ie+1-i,j)
    end do  
    u(p%loc%ie,j) = 0.0d0
    end do
    !$omp end parallel do

else if ( p%glb%ubc(2) == 2)then

    !$omp parallel do
    do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
    do i = 1, p%glb%ghc
        u(p%loc%ie+i,j) = u(p%loc%ie-i,j)
        v(p%loc%ie+i,j) = v(p%loc%ie+1-i,j)
        w(p%loc%ie+i,j) = w(p%loc%ie+1-i,j)
    end do  
    u(p%loc%ie,j) = 0.0d0
    end do
    !$omp end parallel do
            
endif

!==========================================
!  Y-direction velocity boundary condition
!==========================================

if( p%glb%vbc(1) == 1 )then

    !$omp parallel do 
    do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
    do j = 1, p%glb%ghc
        u(i,p%loc%js-j) = - u(i,p%loc%js-1+j)
        v(i,p%loc%js-j) = - v(i,p%loc%js-2+j)
        w(i,p%loc%js-j) = - w(i,p%loc%js-1+j)
    enddo
    v(i,p%loc%js-1) = 0.0d0
    enddo
    !$omp end parallel do
        
else if ( p%glb%vbc(1) == 2 )then

    !$omp parallel do 
    do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
    do j = 1, p%glb%ghc
        u(i,p%loc%js-j) = u(i,p%loc%js-1+j)
        v(i,p%loc%js-j) = v(i,p%loc%js-2+j)
        w(i,p%loc%js-j) = w(i,p%loc%js-1+j)
    enddo
    v(i,p%loc%js-1) = 0.0d0
    enddo
    !$omp end parallel do

endif

if( p%glb%vbc(2) == 1 )then

    !$omp parallel do 
    do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
    do j = 1, p%glb%ghc
        u(i,p%loc%je+j) = - u(i,p%loc%je+1-j)
        v(i,p%loc%je+j) = - v(i,p%loc%je-j)
        w(i,p%loc%je+j) = - w(i,p%loc%je+1-j)
    enddo
    v(i,p%loc%je,j) = 0.0d0
    enddo
    !$omp end parallel do
        
else if ( p%glb%vbc(2) == 2 )then

    !$omp parallel do
    do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
    do j = 1, p%glb%ghc
        u(i,p%loc%je+j) = u(i,p%loc%je+1-j)
        v(i,p%loc%je+j) = v(i,p%loc%je-j)
        w(i,p%loc%je+j) = w(i,p%loc%je+1-j)
    enddo
    v(i,p%loc%je,j) = 0.0d0
    enddo
    !$omp end parallel do
            
endif


end subroutine

subroutine nvelbc(u,v)
! doi.org/10.1063/1.1761178
use all
implicit none
real(8),dimension(p%loc%is-p%glb%ghc:p%loc%ie+p%glb%ghc,&
                  p%loc%js-p%glb%ghc:p%loc%je+p%glb%ghc) :: u,v
integer :: i,j

    !==========================================
    !  X-direction velocity boundary condition
    !==========================================

if( p%glb%ubc(1) == 1 )then

    !$omp parallel do 
    do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc       
    do i = 1, p%glb%ghc
        u(p%loc%is-i,j) = - u(p%loc%is-1+i,j)
        v(p%loc%is-i,j) = - v(p%loc%is-1+i,j)
        w(p%loc%is-i,j) = - w(p%loc%is-1+i,j)
    end do  
    end do
    !$omp end parallel do


else if ( p%glb%ubc(1) == 2)then

    !$omp parallel do 
    do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc       
    do i = 1, p%glb%ghc
        u(p%loc%is-i,j) = u(p%loc%is-1+i,j)
        v(p%loc%is-i,j) = v(p%loc%is-1+i,j)
        w(p%loc%is-i,j) = w(p%loc%is-1+i,j)
    end do  
    end do
    !$omp end parallel do
            
endif

if( p%glb%ubc(2) == 1 )then

    !$omp parallel do 
    do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc       
    do i = 1, p%glb%ghc
        u(p%loc%ie+i,j) = - u(p%loc%ie+1-i,j)
        v(p%loc%ie+i,j) = - v(p%loc%ie+1-i,j)
        w(p%loc%ie+i,j) = - w(p%loc%ie+1-i,j)
    end do  
    end do
    !$omp end parallel do

else if ( p%glb%ubc(2) == 2)then

    !$omp parallel do
    do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc       
    do i = 1, p%glb%ghc
        u(p%loc%ie+i,j) = u(p%loc%ie+1-i,j)
        v(p%loc%ie+i,j) = v(p%loc%ie+1-i,j)
        w(p%loc%ie+i,j) = w(p%loc%ie+1-i,j)
    end do  
    end do
    !$omp end parallel do
            
endif

!==========================================
!  Y-direction velocity boundary condition
!==========================================

if( p%glb%vbc(1) == 1 )then

    !$omp parallel do
    do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
    do j = 1, p%glb%ghc
        u(i,p%loc%js-j) = - u(i,p%loc%js-1+j)
        v(i,p%loc%js-j) = - v(i,p%loc%js-1+j)
        w(i,p%loc%js-j) = - w(i,p%loc%js-1+j)
    enddo
    enddo
    !$omp end parallel do
            
else if ( p%glb%vbc(1) == 2 )then

    !$omp parallel do
    do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
    do j = 1, p%glb%ghc
        u(i,p%loc%js-j) = u(i,p%loc%js-1+j)
        v(i,p%loc%js-j) = v(i,p%loc%js-1+j)
        w(i,p%loc%js-j) = w(i,p%loc%js-1+j)
    enddo
    enddo
    !$omp end parallel do
                
endif

if( p%glb%vbc(2) == 1 )then
    
    !$omp parallel do
    do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
    do j = 1, p%glb%ghc
        u(i,p%loc%je+j) = - u(i,p%loc%je+1-j)
        v(i,p%loc%je+j) = - v(i,p%loc%je+1-j)
        w(i,p%loc%je+j) = - w(i,p%loc%je+1-j)
    enddo
    enddo
    !$omp end parallel do
            
else if ( p%glb%vbc(2) == 2 )then

     !$omp parallel do
    do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
    do j = 1, p%glb%ghc
        u(i,p%loc%je+j) = u(i,p%loc%je+1-j)
        v(i,p%loc%je+j) = v(i,p%loc%je+1-j)
        w(i,p%loc%je+j) = w(i,p%loc%je+1-j)
    enddo
    enddo
    !$omp end parallel do
                
endif
    
end subroutine

subroutine find_stag_vel(u,v,us,vs)
use all
implicit none
real(8),dimension(p%loc%is-p%glb%ghc:p%loc%ie+p%glb%ghc,&
                  p%loc%js-p%glb%ghc:p%loc%je+p%glb%ghc) :: u,v,us,vs
integer :: i, j

!===============================================================
!  Us_t + Us*U_x +  V*U_y = 0
!  Vs_t +  U*V_x + Vs*V_y = 0
!===============================================================

!$omp parallel do collapse(2)
do j = p%loc%js, p%loc%je
do i = p%loc%is, p%loc%ie

    V(i,j) = 0.25d0*( Vs(i,j)+Vs(i,j-1)+Vs(i+1,j)+Vs(i+1,j-1) )
    U(i,j) = 0.25d0*( Us(i,j)+Us(i-1,j)+Us(i,j+1)+Us(i-1,j+1) )

end do
end do
!$omp end parallel do
 
!==========================================
!  X-direction velocity boundary condition
!==========================================

if( p%glb%ubc(1) == 1 )then

    !$omp parallel do 
    do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc       
    do i = 1, p%glb%ghc
        V(p%loc%is-i,j) = - V(p%loc%is-2+i,j)
        U(p%loc%is-i,j) = - U(p%loc%is-1+i,j)
    end do  
    V(p%loc%is-1,j) = 0.0d0
    end do
    !$omp end parallel do

else if ( p%glb%ubc(1) == 2)then

    !$omp parallel do
    do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc       
    do i = 1, p%glb%ghc
        V(p%loc%is-i,j) = V(p%loc%is-2+i,j)
        U(p%loc%is-i,j) = U(p%loc%is-1+i,j)
    end do  
    V(p%loc%is-1,j) = 0.0d0
    end do
    !$omp end parallel do
            
endif

if( p%glb%ubc(2) == 1 )then

    !$omp parallel do
    do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc       
    do i = 1, p%glb%ghc
        V(p%loc%ie+i,j) = - V(p%loc%ie-i,j)
        U(p%loc%ie+i,j) = - U(p%loc%ie+1-i,j)
    end do  
    V(p%loc%ie,j) = 0.0d0
    end do
    !$omp end parallel do

else if ( p%glb%ubc(2) == 2)then

    !$omp parallel do
    do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc       
    do i = 1, p%glb%ghc
        V(p%loc%ie+i,j) = V(p%loc%ie-i,j)
        U(p%loc%ie+i,j) = U(p%loc%ie+1-i,j)
    end do  
    V(p%loc%ie,j) = 0.0d0
    end do
    !$omp end parallel do
            
endif
            
!==========================================
!  Y-direction velocity boundary condition
!==========================================

if( p%glb%vbc(1) == 1 )then
    
    !$omp parallel do
    do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
    do j = 1, p%glb%ghc
        U(i,p%loc%js-j) = - U(i,p%loc%js-2+j)
        V(i,p%loc%js-j) = - V(i,p%loc%js-1+j)
    enddo
    U(i,p%loc%js-1) = 0.0d0
    enddo
    !$omp end parallel do
            
else if ( p%glb%vbc(1) == 2 )then

    !$omp parallel do 
    do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
    do j = 1, p%glb%ghc
        U(i,p%loc%js-j) = U(i,p%loc%js-2+j)
        V(i,p%loc%js-j) = V(i,p%loc%js-1+j)
    enddo
    U(i,p%loc%js-1) = 0.0d0
    enddo
    !$omp end parallel do
                
endif

if( p%glb%vbc(2) == 1 )then
    
    !$omp parallel do 
    do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc  
    do j = 1, p%glb%ghc     
        U(i,p%loc%je+j) = - U(i,p%loc%je-j)
        V(i,p%loc%je+j) = - V(i,p%loc%je+1-j)
    enddo
    U(i,p%loc%je) = 0.0d0
    enddo
    !$omp end parallel do
            
else if ( p%glb%vbc(2) == 2 )then

    !$omp parallel do 
    do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
    do j = 1, p%glb%ghc            
        U(i,p%loc%je+j) = U(i,p%loc%je-j)
        V(i,p%loc%je+j) = V(i,p%loc%je+1-j)
    enddo
    U(i,p%loc%je) = 0.0d0
    enddo
    !$omp end parallel do
                
endif
    
end subroutine
