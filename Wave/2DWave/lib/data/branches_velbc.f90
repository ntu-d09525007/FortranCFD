subroutine job_vel_bc(p,u,v)
! doi.org/10.1063/1.1761178
implicit none
class(job) :: p
integer :: i,j
real(8), dimension(p%loc%is-p%glb%ghc:p%loc%ie+p%glb%ghc,&
                  &p%loc%js-p%glb%ghc:p%loc%je+p%glb%ghc) :: u,v
real(8) :: src
!==========================================
!  X-direction velocity boundary condition
!==========================================

if( p%loc%idx == 0 .and. .not. p%glb%xper )then

    if( p%glb%ubc(1) == 1 )then
    
        do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
        do i = 1, p%glb%ghc
            u(p%loc%is-i,j) = - u(p%loc%is-2+i,j)
            v(p%loc%is-i,j) = - v(p%loc%is-1+i,j)
        end do  
        u(p%loc%is-1,j) = 0.0d0
        end do
    
    else if ( p%glb%ubc(1) == 2)then

        do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
        do i = 1, p%glb%ghc
            u(p%loc%is-i,j) = u(p%loc%is-2+i,j)
            v(p%loc%is-i,j) = v(p%loc%is-1+i,j)
        end do  
        u(p%loc%is-1,j) = 0.0d0
        end do

    else if ( p%glb%ubc(1) == 3)then

        do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc

            src = p%loc%vel%x%old(p%loc%is-1,j)*(u(p%loc%is+1,j)-u(p%loc%is,j))/p%glb%dx*p%glb%dt

            u(p%loc%is-1,j) = p%loc%vel%x%old(p%loc%is-1,j) - src

            do i = 2, p%glb%ghc
                u(p%loc%is-i,j) = 2.0d0*u(p%loc%is-i+1,j)-u(p%loc%is-i+2,j)
            end do

            do i = 1, p%glb%ghc
                v(p%loc%is-i,j) = - v(p%loc%is-1+i,j)
            end do

        end do
                
    endif

endif

if ( p%loc%idx == p%glb%grid_x-1 .and. .not. p%glb%xper )then

     if( p%glb%ubc(2) == 1 )then
    
        do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
        do i = 1, p%glb%ghc
            u(p%loc%ie+i,j) = - u(p%loc%ie-i,j)
            v(p%loc%ie+i,j) = - v(p%loc%ie+1-i,j)
        end do  
        u(p%loc%ie,j) = 0.0d0
        end do
    
    else if ( p%glb%ubc(2) == 2)then

        do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
        do i = 1, p%glb%ghc
            u(p%loc%ie+i,j) = u(p%loc%ie-i,j)
            v(p%loc%ie+i,j) = v(p%loc%ie+1-i,j)
        end do  
        u(p%loc%ie,j) = 0.0d0
        end do

    else if ( p%glb%ubc(2) == 3)then

        do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc
            
            src = p%loc%vel%x%old(p%loc%ie,j)*(u(p%loc%ie-1,j)-u(p%loc%ie-2,j))/p%glb%dx*p%glb%dt

            u(p%loc%ie,j) = p%loc%vel%x%old(p%loc%ie,j) - src

            do i = 1, p%glb%ghc
                u(p%loc%ie+i,j) = 2.0d0*u(p%loc%ie+i-1,j)-u(p%loc%ie+i-2,j)
                v(p%loc%ie+i,j) = - v(p%loc%ie+1-i,j)
            end do

        end do 
                
    endif
        
endif

!==========================================
!  Y-direction velocity boundary condition
!==========================================

if( p%loc%idy==0 .and. .not. p%glb%yper )then

    if( p%glb%vbc(1) == 1 )then

        do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
        do j = 1, p%glb%ghc
            u(i,p%loc%js-j) = - u(i,p%loc%js-1+j)
            v(i,p%loc%js-j) = - v(i,p%loc%js-2+j)
        enddo
        v(i,p%loc%js-1) = 0.0d0
        enddo
            
    else if ( p%glb%vbc(1) == 2 )then

        do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
        do j = 1, p%glb%ghc
            u(i,p%loc%js-j) = u(i,p%loc%js-1+j)
            v(i,p%loc%js-j) = v(i,p%loc%js-2+j)
        enddo
        v(i,p%loc%js-1) = 0.0d0
        enddo

    else if ( p%glb%vbc(1) == 3 )then

        do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc

            src = p%loc%vel%y%old(i,p%loc%js-1)*(v(i,p%loc%js+1)-v(i,p%loc%js))/p%glb%dy*p%glb%dt

            v(i,p%loc%js-1) = p%loc%vel%y%old(i,p%loc%js-1) - src

            do j = 2, p%glb%ghc
                v(i,p%loc%js-j) = 2.0d0*v(i,p%loc%js-j+1) - v(i,p%loc%js-j+2)
            enddo

            do j = 1, p%glb%ghc
                u(i,p%loc%js-j) = - u(i,p%loc%js-1+j)
            enddo

        enddo

    endif

endif
    
if (p%loc%idy==p%glb%grid_y-1 .and. .not. p%glb%yper )then

    if( p%glb%vbc(2) == 1 )then
    
        do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
        do j = 1, p%glb%ghc
            u(i,p%loc%je+j) = - u(i,p%loc%je+1-j)
            v(i,p%loc%je+j) = - v(i,p%loc%je-j)
        enddo
        v(i,p%loc%je) = 0.0d0
        enddo
            
    else if ( p%glb%vbc(2) == 2 )then

        do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
        do j = 1, p%glb%ghc
            u(i,p%loc%je+j) = u(i,p%loc%je+1-j)
            v(i,p%loc%je+j) = v(i,p%loc%je-j)
        enddo
        v(i,p%loc%je) = 0.0d0
        enddo

    else if ( p%glb%vbc(2) == 3 )then

        do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc

            src = p%loc%vel%y%old(i,p%loc%je)*(v(i,p%loc%je-1)-v(i,p%loc%je-2))/p%glb%dy*p%glb%dt

            v(i,p%loc%je) = p%loc%vel%y%old(i,p%loc%je) - src

            do j = 1, p%glb%ghc
                u(i,p%loc%je+j) = - u(i,p%loc%je+1-j)
                v(i,p%loc%je+j) = 2.0d0*v(i,p%loc%je+j-1) - v(i,p%loc%je+j-2)
            enddo

        enddo

    endif
    
endif

end subroutine
                  
subroutine job_nvel_bc(p,u,v)
! doi.org/10.1063/1.1761178
implicit none
class(job) :: p
real(8), dimension(p%loc%is-p%glb%ghc:p%loc%ie+p%glb%ghc,&
                  &p%loc%js-p%glb%ghc:p%loc%je+p%glb%ghc) :: u,v
integer :: i,j

!==========================================
!  X-direction velocity boundary condition
!==========================================

if( p%loc%idx == 0 .and. .not. p%glb%xper )then

    if( p%glb%ubc(1) == 1 )then
    
        do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc       
        do i = 1, p%glb%ghc
            u(p%loc%is-i,j) = - u(p%loc%is-1+i,j)
            v(p%loc%is-i,j) = - v(p%loc%is-1+i,j)
        end do
        end do
    
    else if ( p%glb%ubc(1) == 2)then

        do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc       
        do i = 1, p%glb%ghc
            u(p%loc%is-i,j) = u(p%loc%is-1+i,j)
            v(p%loc%is-i,j) = v(p%loc%is-1+i,j)
        end do
        end do

    else if ( p%glb%ubc(1) == 3)then

        do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc       
        do i = 1, p%glb%ghc
            u(p%loc%is-i,j) = u(p%loc%is,j)
            v(p%loc%is-i,j) = v(p%loc%is,j)
        end do
        end do

    endif

endif 

if ( p%loc%idx == p%glb%grid_x-1 .and. .not. p%glb%xper )then

    if( p%glb%ubc(2) == 1 )then
    
        do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc       
        do i = 1, p%glb%ghc
            u(p%loc%ie+i,j) = - u(p%loc%ie+1-i,j)
            v(p%loc%ie+i,j) = - v(p%loc%ie+1-i,j)
        end do
        end do
    
    else if ( p%glb%ubc(2) == 2)then

        do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc       
        do i = 1, p%glb%ghc
            u(p%loc%ie+i,j) = u(p%loc%ie+1-i,j)
            v(p%loc%ie+i,j) = v(p%loc%ie+1-i,j)
        end do
        end do

    else if ( p%glb%ubc(2) == 3)then

        do j = p%loc%js-p%glb%ghc, p%loc%je+p%glb%ghc       
        do i = 1, p%glb%ghc
            u(p%loc%ie+i,j) = u(p%loc%ie,j)
            v(p%loc%ie+i,j) = v(p%loc%ie,j)
        end do
        end do

    endif
        
endif

!==========================================
!  Y-direction velocity boundary condition
!==========================================

if( p%loc%idy==0 .and. .not. p%glb%yper )then

    if( p%glb%vbc(1) == 1 )then
        
        do j = 1, p%glb%ghc
        do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
            u(i,p%loc%js-j) = - u(i,p%loc%js-1+j)
            v(i,p%loc%js-j) = - v(i,p%loc%js-1+j)
        enddo
        enddo
                
    else if ( p%glb%vbc(1) == 2 )then

        do j = 1, p%glb%ghc
        do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
            u(i,p%loc%js-j) = u(i,p%loc%js-1+j)
            v(i,p%loc%js-j) = v(i,p%loc%js-1+j)
        enddo
        enddo

    else if ( p%glb%vbc(1) == 3 )then

        do j = 1, p%glb%ghc
        do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
            u(i,p%loc%js-j) = u(i,p%loc%js)
            v(i,p%loc%js-j) = v(i,p%loc%js)
        enddo
        enddo

    endif

endif

if (p%loc%idy==p%glb%grid_y-1 .and. .not. p%glb%yper )then

    if( p%glb%vbc(2) == 1 )then
        
        do j = 1, p%glb%ghc
        do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
            u(i,p%loc%je+j) = - u(i,p%loc%je+1-j)
            v(i,p%loc%je+j) = - v(i,p%loc%je+1-j)
        enddo
        enddo
                
    else if ( p%glb%vbc(2) == 2 )then

        do j = 1, p%glb%ghc
        do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
            u(i,p%loc%je+j) = u(i,p%loc%je+1-j)
            v(i,p%loc%je+j) = v(i,p%loc%je+1-j)
        enddo
        enddo

    else if ( p%glb%vbc(2) == 3 )then

        do j = 1, p%glb%ghc
        do i = p%loc%is-p%glb%ghc, p%loc%ie+p%glb%ghc
            u(i,p%loc%je+j) = u(i,p%loc%je)
            v(i,p%loc%je+j) = v(i,p%loc%je)
        enddo
        enddo

    endif

endif
    
end subroutine