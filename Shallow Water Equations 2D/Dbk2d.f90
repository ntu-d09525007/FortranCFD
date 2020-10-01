
       parameter(m=1001,m1=1000,n=101,n1=100,nt=1810000)
	   common     reposedry,reposewet

       dimension  dxc(m),dyc(n),x(m,n),y(m,n),an(m,n),zb(m,n)
       dimension  bednet(m,n),dzb(m,n),bedslopex(m,n),bedslopey(m,n)
       dimension  qx1(m,n),qy1(m,n),h1(m,n),zs1(m,n), &
                  u1(m,n),v1(m,n),uv1(m,n),c1(m,n)
       dimension  qx2(m,n),qy2(m,n),h2(m,n),zs2(m,n), &
                  u2(m,n),v2(m,n),c2(m,n),idix(m,n),kdry(m,n),zbmin(m,n)
       dimension  fmass(m1,n1),fxmom(m1,n1),fymom(m1,n1),fsed(m1,n1)
       dimension  gmass(m1,n1),gxmom(m1,n1),gymom(m1,n1),gsed(m1,n1)
       dimension  alenii(m,n,8),xxrep(m,n,8),yyrep(m,n,8),   &
                  volrep(m,n,8),idixrep(m,n,8)   
       dimension  dzdxmin(m,n), dzdymin(m,n)

	   CHARACTER(30) Fname1	!..............................by He
	   
	   
	   OPEN (10,file="2DbkInput.dat",status='old')
	   OPEN (20,file="2DbkOutput.dat")

       epsilon=0.00000001
       hmin=0.001

       nwsl=1

!!!!       dt0=0.02
       time=0.0

	   read (10,*)
	   read (10,*) chanl,chanwid,xdam
	   read (10,*)
	   read (10,*) damcrest,damheight,depthcrest,depthBreach,slopeup,slopedwn
	   read (10,*)
	   read (10,*) dt0,aMann,diased,poro
	   read (10,*)
	   read (10,*) reposedry,reposewet

!	   write (*,*) chanl,chanwid,xdam
!	   write (*,*) damcrest,damheight,depthcrest,slopeup,slopedwn
!	   write (*,*) dt0,aMann,diased,poro
!	   write (*,*) reposedry,reposewet	   
!	   pause

!!!       chanl=48.63
!!!       chanwid=10.1
!!!       xdam=36.13			!742,744,746
!!!       damcrest=0.2
!       damcrest=0.3
!!!       damheight=0.5
!!!       depthcrest=0.01
!!!       slopeup=0.588235
!!!       slopedwn=0.588235
!       slopeup=0.5
!       slopedwn=0.5

       do i=1,m
         dxc(i)=chanl/real(m-1)
       end do
      
       write(*,*) " dxc(3) , dxc(7)", dxc(3) , dxc(7)
       pause
       
       do j=1,n
         dyc(j)=chanwid/real(n-1)
       end do

       x(1,1)=0.0
       y(1,1)=0.0
       do i=2,m
         x(i,1)=x(i-1,1)+dxc(i)
         y(i,1)=y(i-1,1)
       end do
       do i=1,m
       do j=2,n
         x(i,j)=x(i,j-1)
         y(i,j)=y(i,j-1)+dyc(j)
       end do
       end do

      !$omp parallel do
       do i=1,m
       do j=1,n
         zb(i,j)=0.0
!         an(i,j)=0.015
!         an(i,j)=0.018	!0.02		!0.018
		 an(i,j)=aMann
       end do
       end do
       !$omp end parallel do

       xuptoe=xdam-damcrest/2.0-damheight/slopeup
       xupcorner=xdam-damcrest/2.0
       xdwntoe=xdam+damcrest/2.0+damheight/slopedwn
       xdwncorner=xdam+damcrest/2.0
	   
	   xtoeLength=damcrest+damheight/slopeup+damheight/slopedwn											!........by He 2009
	   ixtoeUp=xuptoe/dxc(1)
	   ixtoeDwn=1+(xuptoe+xtoeLength)/dxc(1)
	   idammid=1+xdam/dxc(1)
	   idamup=(xdam-damcrest/2.0)/dxc(1)
	   idamdwn=1+(xdam+damcrest/2.0)/dxc(1)
	   write(*,*) ixtoeUp,ixtoeDwn,idammid,idamup,idamdwn,n/2
!!	   xdwntoeL1=xdwntoe
       
	   xdam=xdam-xuptoe
       xupcorner=xupcorner-xuptoe
       xdwntoe=xdwntoe-xuptoe
       xdwncorner=xdwncorner-xuptoe
       
       !$omp parallel do
       do i=1,m
       do j=1,n
         x(i,j)=x(i,j)-xuptoe
       end do
       end do
       !$omp end parallel do
       
       xuptoe=0.0

       !$omp parallel do
       do i=1,m
       do j=1,n
	   
         if(x(i,j) .lt. xuptoe) then
            zb(i,j)=0.0
         else if(x(i,j) .ge. xuptoe .and. x(i,j) .le. xupcorner) then
            zb(i,j)=(x(i,j)-xuptoe)*slopeup
         else if(x(i,j) .ge. xupcorner .and. x(i,j) .le. xdwncorner) then
            zb(i,j)=damheight
!			write(*,*) i,j,x(i,j),xupcorner,xdwncorner,zb(i,j)
         else if(x(i,j) .ge. xdwncorner .and. x(i,j) .le. xdwntoe) then
            zb(i,j)=(xdwntoe-x(i,j))*slopedwn
         else if(x(i,j) .gt. xdwntoe) then
            zb(i,j)=0.0
         endif
         
         zbmin(i,j)=0.0      !nonerodible
         
       end do
       end do
       !$omp end parallel do

       !$omp parallel do
       do i=1,m
       do j=1,n
!!	   if(j.ge.n/2-2.and.j.le.n/2+3) then
!!           zb(i,j)=min(zb(i,j),damheight-0.02)
!!         endif
!!         if(j.eq.n/2-3.or.j.eq.n/2+4) then
!!           zb(i,j)=min(zb(i,j),damheight-0.01)
!!         endif
!!	     if(j.ge.n/2-0.and.j.le.n/2+2) then
!!           zb(i,j)=min(zb(i,j),damheight-0.01)
!!         endif
!!!!!         if(j.eq.n/2+1) then
		 if(j > n/2-3 .and. j < n/2+3) then
!!!!!           zb(i,j)=min(zb(i,j),damheight-0.02)		!........by He 200905
           zb(i,j)=min(zb(i,j),damheight-depthBreach)
!			 write(*,*) i,j,x(i,j),zb(i,j)
         end if

!         if(j.eq.n/2-4.or.j.eq.n/2+5) then
!           zb(i,j)=min(zb(i,j),damheight-0.01)
!         endif
!         if(j.eq.n/2-5.or.j.eq.n/2+6) then
!           zb(i,j)=min(zb(i,j),damheight-0.01)
!         endif
       end do
       end do
       !$omp end parallel do

!...... 
       !$omp parallel do
       do i=2,m-1
   	   do j=2,n-1
         xxrep(i,j,1)   =   x(i-1,j-1)
         yyrep(i,j,1)   =   y(i-1,j-1)
         volrep(i,j,1)  =  dxc(i-1)*dyc(j-1)
         xxrep(i,j,2)   =   x(i,j-1)
         yyrep(i,j,2)   =   y(i,j-1)
         volrep(i,j,2)  =  dxc(i)*dyc(j-1)
         xxrep(i,j,3)   =   x(i+1,j-1)
         yyrep(i,j,3)   =   y(i+1,j-1)
         volrep(i,j,3)  =  dxc(i+1)*dyc(j-1)
         xxrep(i,j,4)   =   x(i-1,j)
         yyrep(i,j,4)   =   y(i-1,j)
         volrep(i,j,4)  =  dxc(i-1)*dyc(j)
         xxrep(i,j,5)   =   x(i+1,j)
         yyrep(i,j,5)   =   y(i+1,j)
         volrep(i,j,5)  =  dxc(i+1)*dyc(j)
         xxrep(i,j,6)   =   x(i-1,j+1)
         yyrep(i,j,6)   =   y(i-1,j+1)
         volrep(i,j,6)  =  dxc(i-1)*dyc(j+1)
         xxrep(i,j,7)   =   x(i,j+1)
         yyrep(i,j,7)   =   y(i,j+1)
         volrep(i,j,7)  =  dxc(i)*dyc(j+1)
         xxrep(i,j,8)   =   x(i+1,j+1)
         yyrep(i,j,8)   =   y(i+1,j+1)
         volrep(i,j,8)  =  dxc(i+1)*dyc(j+1)
      end do
      end do
      !$omp end parallel do
      
      !$omp parallel do private(kl)
       do i=2,m-1
   	   do j=2,n-1
   	    
         do kl=1,8
            alenii(i,j,kl)=sqrt( (x(i,j)-xxrep(i,j,kl))**2  &
                                +(y(i,j)-yyrep(i,j,kl))**2)         
         end do
!!!!         if(i.gt.700.and.i.le.800) then			!........by He 2009  	   ixtoeUp,ixtoeDwn
         if(i .ge. ixtoeUp .and. i .le. ixtoeDwn+20) then
            do kl=1,8
               idixrep(i,j,kl)=1
            end do
            idix(i,j)=1
         else
            do kl=1,8
               idixrep(i,j,kl)=0
            end do
            idix(i,j)=0
         endif
         
       end do
       end do
       !$omp end parallel do

!.......................................................output initial inform....................
	call Output0(m,n,time,x,y,zb,zs2,h2,u2,v2,qx2,qy2,c2)
!................................................................................................
       it=0
       call bedrepos0(m,n,x,y,zb,idix,dxc,dyc,it,alenii,   &
                      xxrep,yyrep,volrep,idixrep)
!       call bedrepos(m,n,x,y,zb,dxc,dyc,it)



       !$omp parallel do
       do i=1,m
       do j=1,n
!         qx1(i,j)=0.00122
         qx1(i,j)=0.0
         qy1(i,j)=0.0
         
         if(x(i,j) .gt. xdam) then
            zs1(i,j)=max(0.010,epsilon+zb(i,j))
!           zs1(i,j)=max(1.0,epsilon+zb(i,j))
            zs1(i,j)=min(zs1(i,j),damheight-depthcrest)
         else
            zs1(i,j)=damheight-depthcrest
         endif
         
         h1(i,j)=max(epsilon,zs1(i,j)-zb(i,j))
         u1(i,j)=qx1(i,j)/h1(i,j)
         v1(i,j)=qy1(i,j)/h1(i,j)
         uv1(i,j)=sqrt(u1(i,j)**2+v1(i,j)**2)
         c1(i,j)=0.0
         
         if(h1(i,j) .lt. hmin) then                            
  	        kdry(i,j)=0
         else
            kdry(i,j)=1
         endif
         
       end do
       end do
       !$omp end parallel do

!!!!       diased=0.00025		!...by He 200905
       spec=2.65
!       poro=0.29
!!!!       poro=0.36	!0.35	!...by He 200905
       rhos=2650.0
       rhow=1000.0

       rhob=rhos*(1.0-poro)+rhow*poro 

        amiu=0.000001
        amiud=amiu*1000.0/(diased*1000.0)
        wset=sqrt((13.95*amiud)**2+1.09*(spec-1.0)*9.81*diased)-13.95*amiud

        write(*,*) wset



!....................................................Start Time Loop......................................
       do it=1,nt

!          do i=1,m
!            u1(i)=q1(i)/a1(i)
!            cstar=2.0/2650.0*u1(i)**3/9.8/h1(i)/wset
!            bednet(i)=1.0*wset*(1.0-min(0.5,c1(i)))**4*(c1(i)-cstar)
!          enddo

!          do i=1,m
!            u1(i)=q1(i)/a1(i)
!            alpha=min(2.0,(1.0-poro)/(c1(i)+0.000001))
!            drate=alpha*c1(i)*wset*(1.0-min(0.8,alpha*c1(i)))**2
!            erate=0.015*(an(i)**2*u1(i)**2/h1(i)**0.333/1.65/diased
!     &                   -0.045) *u1(i)/h1(i)/diased**0.2
!            erate=max(0.0,erate)
!            bednet(i)=drate-erate
!            if(x(i).gt.xdwntoe) bednet(i)=max(0.0,bednet(i))
!c        bednet(i)=0.0
!          enddo

         cqb=sqrt(1.65*9.81*diased**3)
         ustarcr2=0.03*1.65*9.81*diased
         
         !$omp parallel do private(uv,ustar2,alamda,ustar2x,ustar2y,exm,wset1,qsstar,Codune,dune,qbstar,qtstar,&
         !$omp &slopx,slopy,aCoqtlimit,qtstarlimit)
         do i=1,m
         do j=1,n
          if(kdry(i,j) .eq. 1) then
            uv=uv1(i,j)
!            if((i.gt.1.and.i.lt.m).and.(j.gt.1.and.j.lt.n)) then
!              if(h1(i,j).gt.hmin) then
!                uv=sqrt( (0.5*(fmass(i,j)+fmass(i-1,j))/h1(i,j))**2
!     &                 +(0.5*(gmass(i,j)+gmass(i,j-1))/h1(i,j))**2 )
!               write(*,*) i,j, uv
!              endif
!            endif
            ustar2=an(i,j)**2*9.81*uv**2/h1(i,j)**0.333333
            
            if(i .eq. 1) then
               bedslopex(i,j)=(zb(i+1,j)-zb(i,j))/dxc(i)
            else if(i .gt. 1 .and. i .lt. m) then
               bedslopex(i,j)=(zb(i+1,j)-zb(i-1,j))/2.0/dxc(i)
            else if(i .eq. m) then
               bedslopex(i,j)=(zb(i,j)-zb(i-1,j))/dxc(i)
            endif
            
            if(j .eq. 1) then
               bedslopey(i,j)=(zb(i,j+1)-zb(i,j))/dyc(j)
            else if(j .gt. 1 .and. j .lt. n) then
               bedslopey(i,j)=(zb(i,j+1)-zb(i,j-1))/2.0/dyc(j)
            else if(j .eq. n) then
               bedslopey(i,j)=(zb(i,j)-zb(i,j-1))/dyc(j)
            endif
!            if(bedslope(i).lt.0.0) then
!               alamda=1.0+0.22*(ustar2/ustarcr2)**0.15
!     &                        *exp(2.0*abs(bedslope(i)/0.625))
!            else
                alamda=1.0
!            endif
            ustar2x=ustar2*u1(i,j)/(uv+epsilon)   &
                           -alamda*ustarcr2*bedslopex(i,j)/0.65
            ustar2y=ustar2*v1(i,j)/(uv+epsilon)   &
                           -alamda*ustarcr2*bedslopey(i,j)/0.65
            ustar2=sqrt(ustar2x**2+ustar2y**2)
!            wset1=wset
!            wset1=wset*(1.0-min(0.65,c1(i,j)))**4
            exm=4.7-1.324*alog10(wset*diased/amiu/5.0)
            exm=min( 4.7, max(2.3,exm) )
            wset1=wset*(1.0-min(0.65,c1(i,j)))**exm

            qsstar=0.0000262*cqb*(max(0.0,ustar2/ustarcr2-1.0)    &
                            *uv/wset1)**1.74
!.....................................................................by He
			Codune=20.0		!20.0
            dune=(diased**0.166667/Codune/an(i,j))**1.5

            qbstar=0.0053*cqb*(max(0.0,dune*ustar2/ustarcr2-1.0))**2.2
            qtstar=qbstar+qsstar

            if((i .gt. 1 .and. i .lt. m).and.(j .gt. 1 .and. j .lt. n)) then
               slopx=(zb(i+1,j)-zb(i-1,j))/dxc(i)/2.0
               slopy=(zb(i,j+1)-zb(i,j-1))/dyc(j)/2.0
               qtstar=qtstar*sqrt(1.0+slopx**2+slopy**2)
            endif
!.....................................................................by He
			aCoqtlimit=0.5			!0.5
            qtstarlimit=aCoqtlimit*uv*h1(i,j)*(1.0-poro)
            
            if(qtstar .gt. qtstarlimit) qtstar=qtstarlimit
      
            bednet(i,j)=(c1(i,j)*uv*h1(i,j)-qtstar)/0.05

            if((i.gt.1.and.i.lt.m).and.(j.gt.1.and.j.lt.n)) then
              if(kdry(i+1,j)+kdry(i-1,j)+kdry(i,j+1)+kdry(i,j-1).eq.4) then
                if(zb(i,j).lt.zbmin(i,j).and.bednet(i,j).lt.0.0)bednet(i,j)=0.0
              endif
            endif
            
            if(i.eq.m) then
                if(zb(i,j) .lt. zbmin(i,j) .and. bednet(i,j) .lt. 0.0) bednet(i,j)=0.0
            endif

            if(h1(i,j).lt.4.0*hmin) bednet(i,j)=0.0
            
          else
          
            bednet(i,j)=0.0
          
          endif
          
         end do
         end do
         !$omp end parallel do

         dt=dt0
         !$omp parallel do reduction(min:dt), private(dtij)
         do i=1,m
         do j=1,n
           if(kdry(i,j) .eq. 1) then
             dtij=0.1*h1(i,j)*(1.0-poro)/(abs(bednet(i,j))+0.01*epsilon)  ! Bed Change <0.1h
             !if(dtij .lt. dt)    write(*,*) i,j, h1(i,j) 
             if(dtij .lt. dt) dt=dtij
             dtij=0.5/( (abs(u1(i,j))+sqrt(9.81*h1(i,j)))/dxc(i)   & 
                        +(abs(v1(i,j))+sqrt(9.81*h1(i,j)))/dyc(j) )    !CFL
             if(dtij.lt.dt) dt=dtij
           endif
         end do
         end do
         !$omp end parallel do
         time=time+dt

         if(it/10*10.eq.it) write(*,*) it, dt, time
         
         !$omp parallel do
         do i=1,m
         do j=1,n
            dzb(i,j)=dt*bednet(i,j)/(1.0-poro)
            zb(i,j)=zb(i,j)+dzb(i,j)
         end do
         end do
         !$omp end parallel do

         call bedrepos(m,n,x,y,zb,zs2,idix,dxc,dyc,it,alenii,   &
                       xxrep,yyrep,volrep,idixrep,hmin)

!         call cflux_Ying(m,n,hmin,h1,zb,zs1,fmass,gmass,fxmom,gxmom,
!     &       fymom,gymom,fsed,gsed,epsilon,qx1,qy1,u1,v1,c1,dxc,dyc,dt)
        
        
      !    call cflux_HLL(m,n,hmin,h1,zb,zs1,fmass,gmass,fxmom,gxmom,   &
      !       fymom,gymom,fsed,gsed,epsilon,qx1,qy1,u1,v1,c1,dxc,dyc,dt)
             
        
          call cflux_HLL_3(m,n,hmin,h1,zb,zs1,fmass,gmass,fxmom,gxmom,   &
             fymom,gymom,fsed,gsed,epsilon,qx1,qy1,u1,v1,c1,dxc,dyc,dt)      
             
             
       !  call cflux_Weno(m,n,hmin,h1,zb,zs1,fmass,gmass,fxmom,gxmom,   &
       !      fymom,gymom,fsed,gsed,epsilon,qx1,qy1,u1,v1,c1,dxc,dyc,dt)    
             
         !$omp parallel do
         do i=2,m-1
         do j=2,n-1
            h2(i,j)=h1(i,j)-dt/dxc(i)*(fmass(i,j)-fmass(i-1,j))   &
                           -dt/dyc(j)*(gmass(i,j)-gmass(i,j-1))   &
                           -0*dt*bednet(i,j)/(1.0-poro)					!.......by He 200905
            h2(i,j)=max(hmin,h2(i,j))
            zs2(i,j)=zb(i,j)+h2(i,j)
         enddo
         enddo
         !$omp end parallel do
         
         !$omp parallel do
         do j=1,n
!           zs2(1,j)=zs2(2,j)
            zs2(1,j)=1.5*zs2(2,j)-0.5*zs2(3,j)
            h2(1,j) =max(hmin, zs2(1,j)-zb(1,j))
!           zs2(m,j)=zs1(m,j)
!           zs2(m,j)=zs2(m-1,j)
            zs2(m,j)=1.5*zs2(m-1,j)-0.5*zs2(m-2,j)
            h2(m,j) =max(hmin, zs2(m,j)-zb(m,j))
!!!!!!		h2(m,j)=uv1(m,j)**2/9.81
         end do
         !$omp end parallel do
         
         !$omp parallel do
         do i=1,m
            zs2(i,1)=zs2(i,2)
!!            zs2(i,1)=2.0*zs2(i,2)-zs2(i,3)
!            if(h1(i,1).gt.2.0*hmin.and.h1(i,2).gt.2.0*hmin) then
!              zs2(i,1)=max(2.0*zs2(i,2)-zs2(i,3),zb(i,1)+epsilon)
!            else
!              zs2(i,1)=max(zs2(i,2),zb(i,1)+epsilon)
!            endif
            h2(i,1) =max(hmin, zs2(i,1)-zb(i,1))

            zs2(i,n)=zs2(i,n-1)
!!            zs2(i,n)=2.0*zs2(i,n-1)-zs2(i,n-2)
!            if(h1(i,n).gt.2.0*hmin.and.h1(i,n-1).gt.2.0*hmin) then
!              zs2(i,n)=max(2.0*zs2(i,n-1)-zs2(i,n-2),zb(i,n)+epsilon)
!            else
!              zs2(i,n)=max(zs2(i,n-1),zb(i,n)+epsilon)
!            endif
            h2(i,n) =max(hmin, zs2(i,n)-zb(i,n))
         end do
         !$omp end parallel do

         if(nwsl.eq.2) call dzdxyminmod(m,n,dzdxmin,dzdymin,x,y,zs2)

!..............X-Momentum
        !$omp parallel do
         do j=1,n
!            qx2(1,j)=0.01
!!            qx2(1,j)=0.0025
           qx2(1,j)=0.0025				!.............by He
!			qx2(1,j)=0.07/10.0
!            qx2(1,j)=0.0
         end do
         !$omp end parallel do

        !$omp parallel do private(zsup,zsdown,dzdx,rho)
         do i=2,m-1
         do j=2,n-1
         
           if(h2(i,j).le.hmin) then
             qx2(i,j)=0.0
!           else if(h2(i+1,j).le.hmin.and.zb(i+1,j).gt.zs2(i,j)) then
!             qx2(i,j)=0.0
!           else if(h2(i-1,j).le.hmin.and.zb(i-1,j).gt.zs2(i,j)) then
!             qx2(i,j)=0.0
           else
           
             if(nwsl.eq.1) then
               zsup  =zs2(i-1,j)
               zsdown=zs2(i+1,j)
!               if(zs2(i-1,j).lt.zb(i,j)) zsup  =zb(i,j)
!               if(zs2(i+1,j).lt.zb(i,j)) zsdown=zb(i,j)
               if(h2(i-1,j).le.hmin.and.zs2(i-1,j).gt.zs2(i,j)) zsup=zs2(i,j)           
               if(h2(i+1,j).le.hmin.and.zs2(i+1,j).gt.zs2(i,j)) zsdown=zs2(i,j)  
                        
	           if( qx1(i,j) .gt. epsilon .and. qx1(i+1,j) .gt. -epsilon ) then             
                 dzdx=zsdown-zs2(i,j)
               else if( qx1(i,j) .lt. -epsilon .and. qx1(i-1,j) .lt. epsilon ) then   
                 dzdx=zs2(i,j)-zsup
	           else                                              
                 dzdx=(zsdown-zsup)/2.0
               endif
               
             endif
             
             if(nwsl.eq.2) then
               zsup  =0.5*( zs2(i-1,j)+0.5*dxc(i-1)*dzdxmin(i-1,j)   &
                         +zs2(i  ,j)-0.5*dxc(i  )*dzdxmin(i  ,j) ) 
               zsdown=0.5*( zs2(i  ,j)+0.5*dxc(i  )*dzdxmin(i  ,j)   &
                         +zs2(i+1,j)-0.5*dxc(i+1)*dzdxmin(i+1,j) ) 
!              if(h2(i-1,j).le.hmin.and.zs2(i-1,j).gt.zs2(i,j))   &
!                                    zsup  =zs2(i,j)           
!              if(h2(i+1,j).le.hmin.and.zs2(i+1,j).gt.zs2(i,j))   &
!                                    zsdown=zs2(i,j)           
               if(i.eq.2  ) zsup  =zs2(1,j)
               if(i.eq.m-1) zsdown=zs2(m,j)
               dzdx=zsdown-zsup
               if(h2(i-1,j).le.hmin.and.zs2(i-1,j).gt.zs2(i,j))dzdx=0.0           
               if(h2(i+1,j).le.hmin.and.zs2(i+1,j).gt.zs2(i,j))dzdx=0.0     
                     
             endif

!             if(h2(i-1,j).le.hmin.and.zb(i-1,j).gt.zs2(i,j)) 
!     &                               dzdx=0.0           
!             if(h2(i+1,j).le.hmin.and.zb(i+1,j).gt.zs2(i,j)) 
!     &                               dzdx=0.0           

             rho=rhos*c1(i,j)+rhow*(1.0-c1(i,j)) 
             qx2(i,j)=qx1(i,j)-dt/dxc(i)*(fxmom(i,j)-fxmom(i-1,j))  &
                              -dt/dyc(j)*(gxmom(i,j)-gxmom(i,j-1))  &
                             -9.81*h2(i,j)*dt/dxc(i)*dzdx           &
!c     &             -9.81*dt*an(i,j)**2*uv1(i,j)*u1(i,j)/h1(i,j)**0.3333
                             -dt*0.5*9.81*h1(i,j)*h1(i,j)*(rhos-rhow)&
                           *(c1(i+1,j)-c1(i-1,j))/2.0/dxc(i)/rho     &
                     +0*dt*u1(i,j)*(rhob-rho)/rho*bednet(i,j)/(1.0-poro)			!!.......by He 200905
             qx2(i,j)=qx2(i,j)/(1.0+9.81*dt*an(i,j)**2*uv1(i,j)/h1(i,j)**1.3333)
!             qx2(i,j)=qx1(i,j)-dt/dxc(i)*(fxmom(i,j)-fxmom(i-1,j))
!     &                         -dt/dyc(j)*(gxmom(i,j)-gxmom(i,j-1))
!     &             -9.81*h2(i,j)*dt/dxc(i)*dzdx
!     &             -9.81*dt*an(i,j)**2*uv1(i,j)*u1(i,j)/h1(i,j)**0.3333
!     &             -dt*0.5*9.81*h1(i,j)*h1(i,j)*(rhos-rhow)
!     &                    *(c1(i+1,j)-c1(i-1,j))/2.0/dxc(i)/rho
!     &             +dt*u1(i,j)*(rhob-rho)/rho*bednet(i,j)/(1.0-poro)
           endif
         end do
         end do
         !$omp end parallel do
         
         !$omp parallel do
         do j=1,n
!           qx2(m,j)=qx2(m-1,j)
			qx2(m,j)=SQRT(9.81*h2(m,j)**3.0)
         end do
         !$omp end parallel do
         
         !$omp parallel do
         do i=1,m
!            qx2(i,1)=0.0
!            qx2(i,n)=0.0
            qx2(i,1)=qx2(i,2)
            qx2(i,n)=qx2(i,n-1)
         end do
         !$omp end parallel do

!..............Y-Momentum
         !$omp parallel do
         do j=1,n
            qy2(1,j)=0.0
         end do
         !$omp end parallel do

        !$omp parallel do private(zsup,zsdown,dzdy,rho)
         do i=2,m-1
         do j=2,n-1
         
           if(h2(i,j) .le. hmin) then
             qy2(i,j)=0.0
!           else if(h2(i,j+1).le.hmin.and.zb(i,j+1).gt.zs2(i,j)) then
!             qy2(i,j)=0.0
!           else if(h2(i,j-1).le.hmin.and.zb(i,j-1).gt.zs2(i,j)) then
!             qy2(i,j)=0.0
           else

             if(nwsl.eq.1) then
               zsup  =zs2(i,j-1)
               zsdown=zs2(i,j+1)
!               if(zs2(i,j-1).lt.zb(i,j)) zsup  =zb(i,j)
!               if(zs2(i,j+1).lt.zb(i,j)) zsdown=zb(i,j)
               if(h2(i,j-1).le.hmin.and.zs2(i,j-1).gt.zs2(i,j))zsup  =zs2(i,j)           
               if(h2(i,j+1).le.hmin.and.zs2(i,j+1).gt.zs2(i,j))zsdown=zs2(i,j)     
                     
	           if(qy1(i,j) .gt. epsilon .and. qy1(i,j+1) .gt. -epsilon) then             
                 dzdy=zsdown-zs2(i,j)
               else if(qy1(i,j) .lt. -epsilon .and. qy1(i,j-1) .lt. epsilon) then   
                 dzdy=zs2(i,j)-zsup
  	           else                                              
                 dzdy=(zsdown-zsup)/2.0
               endif
             endif
             
             if(nwsl .eq. 2) then
               zsup  =0.5*( zs2(i,j-1)+0.5*dyc(j-1)*dzdymin(i,j-1)   &
                         +zs2(i,j  )-0.5*dyc(j  )*dzdymin(i,j  ) ) 
               zsdown=0.5*( zs2(i,j  )+0.5*dyc(j  )*dzdymin(i,j  )   &
                         +zs2(i,j+1)-0.5*dyc(j+1)*dzdymin(i,j+1) ) 
!               if(h2(i,j-1).le.hmin.and.zs2(i,j-1).gt.zs2(i,j))    &
!                                    zsup  =zs2(i,j)           
!               if(h2(i,j+1).le.hmin.and.zs2(i,j+1).gt.zs2(i,j))    &
!                                    zsdown=zs2(i,j)           
               if(j .eq. 2  ) zsup  =zs2(i,1)
               if(j .eq. n-1) zsdown=zs2(i,n)
               dzdy=zsdown-zsup
               
               if(h2(i,j-1).le.hmin.and.zs2(i,j-1).gt.zs2(i,j))dzdy=0.0           
               if(h2(i,j+1).le.hmin.and.zs2(i,j+1).gt.zs2(i,j))dzdy=0.0 
                         
             endif

!             if(h2(i,j-1).le.hmin.and.zb(i,j-1).gt.zs2(i,j)) 
!     &                              dzdy=0.0           
!             if(h2(i,j+1).le.hmin.and.zb(i,j+1).gt.zs2(i,j)) 
!     &                              dzdy=0.0          

             rho=rhos*c1(i,j)+rhow*(1.0-c1(i,j)) 
             qy2(i,j)=qy1(i,j)-dt/dxc(i)*(fymom(i,j)-fymom(i-1,j))   &
                              -dt/dyc(j)*(gymom(i,j)-gymom(i,j-1))   &
                              -9.81*h2(i,j)*dt/dyc(j)*dzdy           &
!     &             -9.81*dt*an(i,j)**2*uv1(i,j)*v1(i,j)/h1(i,j)**0.3333
                              -dt*0.5*9.81*h1(i,j)*h1(i,j)*(rhos-rhow)&
                               *(c1(i,j+1)-c1(i,j-1))/2.0/dyc(j)/rho  &
                    +0*dt*v1(i,j)*(rhob-rho)/rho*bednet(i,j)/(1.0-poro)			!.......by He 200905
             qy2(i,j)=qy2(i,j) /(1.0+9.81*dt*an(i,j)**2*uv1(i,j)/h1(i,j)**1.3333)
!             qy2(i,j)=qy1(i,j)-dt/dxc(i)*(fymom(i,j)-fymom(i-1,j))
!     &                         -dt/dyc(j)*(gymom(i,j)-gymom(i,j-1))
!     &             -9.81*h2(i,j)*dt/dyc(j)*dzdy
!     &             -9.81*dt*an(i,j)**2*uv1(i,j)*v1(i,j)/h1(i,j)**0.3333
!     &             -dt*0.5*9.81*h1(i,j)*h1(i,j)*(rhos-rhow)
!     &                    *(c1(i,j+1)-c1(i,j-1))/2.0/dyc(j)/rho
!     &             +dt*v1(i,j)*(rhob-rho)/rho*bednet(i,j)/(1.0-poro)

           endif
           
         end do
         end do
         !$omp end parallel do
        
         !$omp parallel do
         do j=1,n
            qy2(m,j)=qy2(m-1,j)
         end do
         !$omp end parallel do
         
         !$omp parallel do
         do i=1,m
            qy2(i,1)=0.0
            qy2(i,n)=0.0
         end do
         !$omp end parallel do

!.........Sediment
         
         !$omp parallel do
         do j=1,n
            c2(1,j)=0.0
         end do
         !$omp end parallel do
         
         !$omp parallel do
         do i=2,m-1
         do j=2,n-1
            c2(i,j)=(c1(i,j)*h1(i,j)-dt/dxc(i)*(fsed(i,j)-fsed(i-1,j))   &
                                    -dt/dyc(j)*(gsed(i,j)-gsed(i,j-1))   &
                                    -dt*bednet(i,j) )/h2(i,j)
            c2(i,j)=max(0.0,c2(i,j))
            c2(i,j)=min(1-poro,c2(i,j))
         end do
         end do
         !$omp end parallel do
         
         !$omp parallel do
         do j=1,n
            c2(m,j)=c2(m-1,j)
         end do
         !$omp end parallel do
         
         !$omp parallel do
         do i=1,m
            c2(i,1)=c2(i,2)
            c2(i,n)=c2(i,n-1)
         end do
         !$omp end parallel do

         !$omp parallel do
         do i=1,m
         do j=1,n
            if(h2(i,j).le.hmin) then
               u2(i,j)=0.0
               v2(i,j)=0.0
            else
               u2(i,j)=qx2(i,j)/h2(i,j)
               v2(i,j)=qy2(i,j)/h2(i,j)
            endif
            u1(i,j)=u2(i,j)        
            v1(i,j)=v2(i,j)        
            uv1(i,j)=sqrt(u1(i,j)**2+v1(i,j)**2)
            qx1(i,j)=qx2(i,j)        
            qy1(i,j)=qy2(i,j)        
            zs1(i,j)=zs2(i,j)        
            h1(i,j)=h2(i,j)        
            c1(i,j)=c2(i,j)
            if(h1(i,j).lt.hmin) then 
  	         kdry(i,j)=0
            else
               kdry(i,j)=1
            endif
         end do
         end do
         !$omp end parallel do

         if(it/240*240-it.eq.0) then
            write(99,*) it,time
!            do i=m/2-100,m/2+50
!............................................................by He
			qsum1=0.0		
			qsum2=0.0		!................................by He
			qsum3=0.0		!................................by He
			WidthBreach742=0.
			WidthBreach744=0.
			WidthBreach2=0.
			WidthBreach1=0.
!.................................................................
!!	   ixtoeUp,	   ixtoeDwn,	   idammid,	   idamup,	   idamdwn
!!!!!!!!!            do i=700,850			!.............by He 200905
			do i=ixtoeUp,ixtoeDwn
!            do i=150,400
              qsum=0.0
			  
			  !$omp parallel do reduction(+:qsum), private(WidthBreach1,WidthBreach2,WidthBreach742,WidthBreach744)
              do j=2,n-1
                qsum=qsum+qx2(i,j)*dyc(j)
!.............................................................by He
!!!				If (i==742) Then
				If (i==idammid) Then
					IF (zb(i,j) <= damheight .and. abs(zb(i,j-1)-damheight) < 1.E-10 .and.  zb(i,j+1) < zb(i,j) ) WidthBreach1=y(i,j)
					IF (zb(i,j) <= damheight .and. abs(zb(i,j+1)-damheight) < 1.E-10 .and. zb(i,j)>zb(i,j-1) ) WidthBreach2=y(i,j)
					WidthBreach742=WidthBreach2-WidthBreach1
				Endif
				If (i==idamup+1) Then
					IF (zb(i,j) <= damheight .and. abs(zb(i,j-1)-damheight) < 1.E-10 .and.  zb(i,j+1) < zb(i,j) ) WidthBreach1=y(i,j)
					IF (zb(i,j) <= damheight .and. abs(zb(i,j+1)-damheight) < 1.E-10 .and. zb(i,j)>zb(i,j-1) ) WidthBreach2=y(i,j)
					WidthBreach744=WidthBreach2-WidthBreach1
				Endif
!..................................................................
              enddo
			  !$omp end parallel do
			  
              write(99,*) x(i,n/2),qsum,zb(i,n/2),zs2(i,n/2)
!              write(99,*) x(i,n/2),qx2(i,n/2),zb(i,n/2),zs2(i,n/2)
!.............................................................by He
!			  IF (i==742) qsum1=qsum
!			  IF (i==735) qsum2=qsum
!			  IF (i==730) qsum3=qsum
			  IF (i==idammid) qsum1=qsum
			  IF (i==idamup+1) qsum2=qsum
			  IF (i==idamdwn) qsum3=qsum
!..................................................................
            enddo
!...................................................................by He
			  WRITE(20,980) time,x(744,n/2+1),qsum1,qsum2,qsum3,zb(744,n/2+1),zb(744,n/2+1),zs2(740,n/2+1)*1000,zs2(742,n/2+1)*1000,zs2(744,n/2+1)*1000,WidthBreach742,WidthBreach744
980			  FORMAT(12F18.12)
!........................................................................
         endif

         numio=it/6000
         if(it/6000*6000-it.eq.0) then
!           numio=it/24000
!          if(it/24000*24000-it.eq.0) then
           write(20+numio,*) 'time=', time
           do i=1,m
             write(*,*) x(i,n/2),zs2(i,n/2),zb(i,n/2)
!              write(20+numio,*) x(i,n/2),zb(i,n/2),zs2(i,n/2),qx1(i,n/2)
             do j=1,n
               write(20+numio,*) x(i,j),y(i,j),zb(i,j),zs2(i,j),h2(i,j),   &
                              u2(i,j),v2(i,j),qx2(i,j),qy2(i,j)
             enddo
           enddo
         endif

!.................................................................by He			
         Ioutput=80*10
         if(MOD(it,Ioutput).eq.1) then
!           numio=it/24000
!          if(it/24000*24000-it.eq.0) then
!           write(20+numio,*) 'time=', time
		   WRITE(Fname1,5000) it/Ioutput
5000	   FORMAT('Output/Tec2D',I5.5,'.dat')
		   OPEN(300,FILE=Fname1)
		   write(300,*)	'title= " Time= ', time, ' s " ' 
		   write(300,*) 'VARIABLES ="X","Y","Z","Zs","h","U","V","W","Qx","Qy","Concen"'
		   write(300,*) 'ZONE T="contour"',' I=',m,' J=',n,' F=point'
           do j=1,n
		     do i=1,m
             
               write(300,9980) x(i,j),y(i,j),zb(i,j),zs2(i,j),h2(i,j),   &
                              u2(i,j),v2(i,j),0.0,qx2(i,j),qy2(i,j),c2(i,j)
             enddo
           enddo
9980	   format(2F12.4,9E18.8)
		   CLOSE(300)	
		 endif
!................................................................................

       enddo

	   close(10)
	   close(20)
       end

!#################################################################################################################

       subroutine cflux_HLL(m,n,hmin,h1,zb,zs1,fmass,gmass,fxmom,gxmom,   &
              fymom,gymom,fsed,gsed,epsilon,qx1,qy1,u1,v1,c1,dxc,dyc,dt)
       dimension  dxc(m),dyc(n)
       dimension  qx1(m,n),qy1(m,n),h1(m,n),zs1(m,n),      &
                  u1(m,n),v1(m,n),c1(m,n),zb(m,n)
       dimension  fmass(m-1,n-1),fxmom(m-1,n-1),fymom(m-1,n-1)
       dimension  gmass(m-1,n-1),gxmom(m-1,n-1),gymom(m-1,n-1)
       dimension  fsed(m-1,n-1),gsed(m-1,n-1)

      !$omp parallel do
       do j=2,n-1
         fmass(1,j)=qx1(1,j)
         fxmom(1,j)=qx1(1,j)*u1(1,j)
         fymom(1,j)=qx1(1,j)*v1(1,j)
         fsed(1,j) =qx1(1,j)*c1(1,j)
         fmass(m-1,j)=qx1(m,j)
         fxmom(m-1,j)=qx1(m,j)*u1(m,j)
         fymom(m-1,j)=qx1(m,j)*v1(m,j)
         fsed(m-1,j) =qx1(m,j)*c1(m,j)
       end do
       !$omp end parallel do
       
       !$omp parallel do private(ux1,uxr,hx1,hxr,ax1,axr,hxstar,uxstar,aqxr,aqxl,sxl,sxr,&
       !$omp &fmassl,fmassr,fxmoml,fxmomr,fymoml,fymomr)
       do j=2,n-1
       do i=2,m-2
       
            uxl=u1(i  ,j)
            uxr=u1(i+1,j)
            hxl=h1(i,j)
            hxr=h1(i+1,j)
            axl=sqrt(9.81*hxl)
            axr=sqrt(9.81*hxr)
            hxstar=0.5*(hxl+hxr)-0.25*(uxr-uxl)*(hxl+hxr)/(axl+axr)
            uxstar=0.5*(uxl+uxr)-(hxr-hxl)*(axl+axr)/(hxl+hxr)
              
	        if(hxstar .gt. hxl) then
               aqxl=sqrt(0.5*(hxstar+hxl)*hxstar/hxl**2)
            else
               aqxl=1.0
            endif
            
	        if(hxstar .gt. hxr) then
               aqxr=sqrt(0.5*(hxstar+hxr)*hxstar/hxr**2)
            else
               aqxr=1.0
            endif
            
            sxl=uxl-axl*aqxl
            sxr=uxr+axr*aqxr

            fmassl=qx1(i,  j)              
            fmassr=qx1(i+1,j)              
            fxmoml=fmassl*u1(i,  j)
            fxmomr=fmassr*u1(i+1,j)
            fymoml=fmassl*v1(i,  j)
            fymomr=fmassr*v1(i+1,j)
!            fxmoml=qx1(i,  j)*u1(i,  j)
!            fxmomr=qx1(i+1,j)*u1(i+1,j)
!            fymoml=qx1(i,  j)*v1(i,  j)
!            fymomr=qx1(i+1,j)*v1(i+1,j)

            if(sxl.gt.0.0) then
               fmass(i,j)=fmassl             
               fxmom(i,j)=fxmoml
               fymom(i,j)=fymoml
            else if(sxl.le.0.0.and.sxr.ge.0.0) then
               fmass(i,j)=( sxr*fmassl - sxl*fmassr    &
                            + sxl*sxr*(zs1(i+1,j)-zs1(i,j)) )/(sxr-sxl)             
               fxmom(i,j)=( sxr*fxmoml - sxl*fxmomr    &
                            + sxl*sxr*(qx1(i+1,j)-qx1(i,j)) )/(sxr-sxl)
               fymom(i,j)=( sxr*fymoml - sxl*fymomr    &
                            + sxl*sxr*(qy1(i+1,j)-qy1(i,j)) )/(sxr-sxl)
            else if(sxr.lt.0.0) then
               fmass(i,j)=fmassr             
               fxmom(i,j)=fxmomr
               fymom(i,j)=fymomr
            endif

            if(fmass(i,j).gt.0.0) then
               fsed(i,j)=fmass(i,j)*c1(i,j) 
            else if(fmass(i,j).lt.0.0) then
               fsed(i,j)=fmass(i,j)*c1(i+1,j) 
            else
               fsed(i,j)=fmass(i,j)*0.5*(c1(i,j)+c1(i+1,j)) 
            endif

            if(h1(i,j).le.hmin.and.zb(i,j).gt.zs1(i+1,j)) then
               fmass(i,j)=0.0
               fxmom(i,j)=0.0
               fymom(i,j)=0.0
               fsed(i,j) =0.0
	        endif 
	   	    if(h1(i+1,j).le.hmin.and.zb(i+1,j).gt.zs1(i,j)) then
               fmass(i,j)=0.0
               fxmom(i,j)=0.0
               fymom(i,j)=0.0
               fsed(i,j) =0.0
            endif
	  	    if(h1(i,j).le.hmin.and.h1(i+1,j).le.hmin) then
               fmass(i,j)=0.0
               fxmom(i,j)=0.0
               fymom(i,j)=0.0
               fsed(i,j) =0.0
            endif
            
         end do
         end do
         !$omp end parallel do

		!$omp parallel do
        do i=2,m-1
         gmass(i,1)=qy1(i,1)
         gxmom(i,1)=qy1(i,1)*u1(i,1)
         gymom(i,1)=qy1(i,1)*v1(i,1)
         gsed(i,1) =qy1(i,1)*c1(i,1)
         gmass(i,n-1)=qy1(i,n)
         gxmom(i,n-1)=qy1(i,n)*u1(i,n)
         gymom(i,n-1)=qy1(i,n)*v1(i,n)
         gsed(i,n-1) =qy1(i,n)*c1(i,n)
		enddo
		!$omp end parallel do
		
		!$omp parallel do private(uyl,uyr,hyl,hyr,ayl,ayr,hystar,uystar,aqyl,aqyr,syl,syr,&
		!$omp &gmassl,gmassr,gxmoml,gxmomr,gymoml,gymomr)
		do i=2,m-1
         do j=2,n-2
            uyl=v1(i,j  )
            uyr=v1(i,j+1)
            hyl=h1(i,j)
            hyr=h1(i+1,j)
            ayl=sqrt(9.81*hyl)
            ayr=sqrt(9.81*hyr)
            hystar=0.5*(hyl+hyr)-0.25*(uyr-uyl)*(hyl+hyr)/(ayl+ayr)
            uystar=0.5*(uyl+uyr)-(hyr-hyl)*(ayl+ayr)/(hyl+hyr)
              
	        if(hystar.gt.hyl) then
               aqyl=sqrt(0.5*(hystar+hyl)*hystar/hyl**2)
            else
               aqyl=1.0
            endif
			
	      if(hystar.gt.hyr) then
               aqyr=sqrt(0.5*(hystar+hyr)*hystar/hyr**2)
            else
               aqyr=1.0
            endif
            syl=uyl-ayl*aqyl
            syr=uyr+ayr*aqyr

            gmassl=qy1(i,j  )              
            gmassr=qy1(i,j+1)              
            gxmoml=gmassl*u1(i,j  )
            gxmomr=gmassr*u1(i,j+1)
            gymoml=gmassl*v1(i,j  )
            gymomr=gmassr*v1(i,j+1)
!            gxmoml=qy1(i,j  )*u1(i,j  )
!            gxmomr=qy1(i,j+1)*u1(i,j+1)
!            gymoml=qy1(i,j  )*v1(i,j  )
!            gymomr=qy1(i,j+1)*v1(i,j+1)

            if(syl.gt.0.0) then
               gmass(i,j)=gmassl             
               gxmom(i,j)=gxmoml
               gymom(i,j)=gymoml
            elseif(syl.le.0.0.and.syr.ge.0.0) then
               gmass(i,j)=( syr*gmassl - syl*gmassr      &
                            + syl*syr*(zs1(i,j+1)-zs1(i,j)) )/(syr-syl)             
               gxmom(i,j)=( syr*gxmoml - syl*gxmomr      &
                            + syl*syr*(qx1(i,j+1)-qx1(i,j)) )/(syr-syl)
               gymom(i,j)=( syr*gymoml - syl*gymomr      &
                            + syl*syr*(qy1(i,j+1)-qy1(i,j)) )/(syr-syl)
            else if(syr.lt.0.0) then
               gmass(i,j)=gmassr             
               gxmom(i,j)=gxmomr
               gymom(i,j)=gymomr
            endif

            if(gmass(i,j).gt.0.0) then
               gsed(i,j)=gmass(i,j)*c1(i,j) 
            else if(gmass(i,j).lt.0.0) then
               gsed(i,j)=gmass(i,j)*c1(i,j+1) 
            else
               gsed(i,j)=gmass(i,j)*0.5*(c1(i,j)+c1(i,j+1)) 
            endif

            if(h1(i,j).le.hmin.and.zb(i,j).gt.zs1(i,j+1)) then
               gmass(i,j)=0.0
               gxmom(i,j)=0.0
               gymom(i,j)=0.0
               gsed(i,j) =0.0
	        endif
		    if(h1(i,j+1).le.hmin.and.zb(i,j+1).gt.zs1(i,j)) then
               gmass(i,j)=0.0
               gxmom(i,j)=0.0
               gymom(i,j)=0.0
               gsed(i,j) =0.0
	        endif
            if(h1(i,j).le.hmin.and.h1(i,j+1).le.hmin) then
               gmass(i,j)=0.0
               gxmom(i,j)=0.0
               gymom(i,j)=0.0
               gsed(i,j) =0.0
	        endif
	     enddo
         enddo
		!$omp end parallel do 
		
		
      return
      end
      
      
      subroutine cflux_HLL_1(m,n,hmin,h1,zb,zs1,fmass,gmass,fxmom,gxmom,   &
              fymom,gymom,fsed,gsed,epsilon,qx1,qy1,u1,v1,c1,dxc,dyc,dt)
       dimension  dxc(m),dyc(n)
       dimension  qx1(m,n),qy1(m,n),h1(m,n),zs1(m,n),      &
                  u1(m,n),v1(m,n),c1(m,n),zb(m,n)
       dimension  fmass(m-1,n-1),fxmom(m-1,n-1),fymom(m-1,n-1)
       dimension  gmass(m-1,n-1),gxmom(m-1,n-1),gymom(m-1,n-1)
       dimension  fsed(m-1,n-1),gsed(m-1,n-1)
       
       
       dimension  up(-2:m+3,-2:n+3),um(-2:m+3,-2:n+3)
       dimension  fi(-2:m+3,-2:n+3), ff(-2:m+3,-2:n+3)
       
       dimension  fin(-2:m+3,-2:n+3) , uh(-2:m+3,-2:n+3)

       dimension  dfp(-2:m+3,-2:n+3),dfm(-2:m+3,-2:n+3)
       
       
        dimension  c_x(-2:m+3,-2:n+3)

		!$omp parallel do 
        do j=2,n-1
         fmass(1,j)=qx1(1,j)
         fxmom(1,j)=qx1(1,j)*u1(1,j)
         fymom(1,j)=qx1(1,j)*v1(1,j)
         fsed(1,j) =qx1(1,j)*c1(1,j)
         fmass(m-1,j)=qx1(m,j)
         fxmom(m-1,j)=qx1(m,j)*u1(m,j)
         fymom(m-1,j)=qx1(m,j)*v1(m,j)
         fsed(m-1,j) =qx1(m,j)*c1(m,j)
        enddo 
		!$omp end parallel do 
         
		 !$omp parallel do private(uxl,uxr,hxl,hxr,axl,axr,hxstar,uxstar,aqxl,aqxr,sxl,sxr,&
		 !$omp &fmassl,fmassr,fxmoml,fxmomr,fymoml,fymomr)
         do j=2,n-1
         do i=2,m-2
            uxl=u1(i  ,j)
            uxr=u1(i+1,j)
            hxl=h1(i,j)
            hxr=h1(i+1,j)
            axl=sqrt(9.81*hxl)
            axr=sqrt(9.81*hxr)
            hxstar=0.5*(hxl+hxr)-0.25*(uxr-uxl)*(hxl+hxr)/(axl+axr)
            uxstar=0.5*(uxl+uxr)-(hxr-hxl)*(axl+axr)/(hxl+hxr)
              
	        if(hxstar.gt.hxl) then
               aqxl=sqrt(0.5*(hxstar+hxl)*hxstar/hxl**2)
            else
               aqxl=1.0
            endif
	        if(hxstar.gt.hxr) then
               aqxr=sqrt(0.5*(hxstar+hxr)*hxstar/hxr**2)
            else
               aqxr=1.0
            endif
            sxl=uxl-axl*aqxl
            sxr=uxr+axr*aqxr

            fmassl=qx1(i,  j)              
            fmassr=qx1(i+1,j)              
            fxmoml=fmassl*u1(i,  j)
            fxmomr=fmassr*u1(i+1,j)
            fymoml=fmassl*v1(i,  j)
            fymomr=fmassr*v1(i+1,j)


            if(sxl.gt.0.0) then
               fmass(i,j)=fmassl             
               fxmom(i,j)=fxmoml
               fymom(i,j)=fymoml
            elseif(sxl.le.0.0.and.sxr.ge.0.0) then
               fmass(i,j)=( sxr*fmassl - sxl*fmassr    &
                            + sxl*sxr*(zs1(i+1,j)-zs1(i,j)) )/(sxr-sxl)             
               fxmom(i,j)=( sxr*fxmoml - sxl*fxmomr    &
                            + sxl*sxr*(qx1(i+1,j)-qx1(i,j)) )/(sxr-sxl)
               fymom(i,j)=( sxr*fymoml - sxl*fymomr    &
                            + sxl*sxr*(qy1(i+1,j)-qy1(i,j)) )/(sxr-sxl)
            else if(sxr.lt.0.0) then
               fmass(i,j)=fmassr             
               fxmom(i,j)=fxmomr
               fymom(i,j)=fymomr
            endif
            
        enddo
        enddo 
        !$omp end parallel do 
        
      
      
	   !   do i= 2,m-2
	   !   do j= 2,n-2
	   !     c_x(i,j)=0.5*( uh(i+1,j)+ uh(i,j) )     
	   !   enddo
	   !   enddo
            
        
      
	  ! 
	  !   do j=2,n-1
    !   do i=2,m-2
    !  
    !   	
    !   !	if( c_x(i,j) >=0.0 	) then
	  !   !   fmass(i,j)= c_x(i,j)*(1.0/12.0*(-fin(i-1,j)+7.0*fin(i,j)+7.0*fin(i+1,j)-fin(i+2,j))  -phyn(dfp(i-2,j),dfp(i-1,j),dfp(i,j),dfp(i+1,j)) )
	  !   !	       
	  !   !	else      
	  !   !	       
	  !   !  fmass(i,j)=c_x(i,j)*(1.0/12.0*(-fin(i-1,j)+7.0*fin(i,j)+7.0*fin(i+1,j)-fin(i+2,j)) 	 + phyn(dfm(i+2,j),dfm(i+1,j),dfm(i,j),dfm(i-1,j))    )
	  !   !	       
	  !   !	endif
	  !   
	  !        fmass(i,j) =  fmass(i,j)
	  !         
	  !  	       
	  !    enddo
	  !    enddo
            
            
            
    !     do i=-2,m+3
	  !     do j=-2,n+3
	  !       fin(i,j)= 0.0
	  !       uh(i,j)=  0.0
	  !     enddo
	  !     enddo 
    !     
    !     
    !     do i=1,m
	  !     do j=1,n
	  !       fin(i,j)= h1(i,j)
	  !       uh(i,j)=  u1(i,j)
	  !     
	  !     enddo
	  !     enddo
    !     
    !   
    !     CALL BC2D(uh,M,N)
    !     CALL BC2D(fin,M,N)
    !   
    !   
	  !     do i=-2,m+3
	  !     do j=-2,n+3
	  !       up(i,j)=0.5*(uh(i,j)+abs(uh(i,j)))
	  !       um(i,j)=0.5*(uh(i,j)-abs(uh(i,j)))
	  !     !  vp(i,j)=0.5*(vh(i,j)+abs(vh(i,j)))
	  !     !  vm(i,j)=0.5*(vh(i,j)-abs(vh(i,j)))
	  !     enddo
	  !     enddo
    !         
    !     
    !    do i=0,m+1 
    !    do j=1,n
    !   
	  !       ff(i,j)=uh(i,j)*fin(i,j)
	  !     enddo
	  !     enddo
	  !     
	  !     ! do j=1,n
	  !     do i=-1,m+1
	  !   	do j=1,n
	  !       dfp(i,j)=up(i+1,j)*fi(i+1,j)-up(i,j)*fi(i,j)
	  !       dfm(i,j)=um(i+1,j)*fi(i+1,j)-um(i,j)*fi(i,j)
	  !     enddo
	  !     enddo
	  !     
	  !   
    !    do i=2,m-2
    !    do j=2,n-1 	
	  !       fmass(i,j)=1.0/12.0*(-ff(i-1,j)+7.0*ff(i,j)+7.0*ff(i+1,j)-ff(i+2,j)) &
	  !   	       -phyn(dfp(i-2,j),dfp(i-1,j),dfp(i,j),dfp(i+1,j)) &
	  !   	       +phyn(dfm(i+2,j),dfm(i+1,j),dfm(i,j),dfm(i-1,j))
	  !     enddo
	  !     enddo    
    !         
    !         
    !         
    !      do j=2,n-1
    !      do i=2,m-2
    !         if(fmass(i,j).gt.0.0) then
    !            fsed(i,j)=fmass(i,j)*c1(i,j) 
    !         else if(fmass(i,j).lt.0.0) then
    !            fsed(i,j)=fmass(i,j)*c1(i+1,j) 
    !         else
    !            fsed(i,j)=fmass(i,j)*0.5*(c1(i,j)+c1(i+1,j)) 
    !         endif
    !
    !         if(h1(i,j).le.hmin.and.zb(i,j).gt.zs1(i+1,j)) then
    !            fmass(i,j)=0.0
    !            fxmom(i,j)=0.0
    !            fymom(i,j)=0.0
    !            fsed(i,j) =0.0
	  !       endif 
	  !  	    if(h1(i+1,j).le.hmin.and.zb(i+1,j).gt.zs1(i,j)) then
    !            fmass(i,j)=0.0
    !            fxmom(i,j)=0.0
    !            fymom(i,j)=0.0
    !            fsed(i,j) =0.0
    !         endif
	  ! 	    if(h1(i,j).le.hmin.and.h1(i+1,j).le.hmin) then
    !            fmass(i,j)=0.0
    !            fxmom(i,j)=0.0
    !            fymom(i,j)=0.0
    !            fsed(i,j) =0.0
    !         endif
    !    enddo
    !    enddo

		!$omp parallel do 
        do i=2,m-1
         gmass(i,1)=qy1(i,1)
         gxmom(i,1)=qy1(i,1)*u1(i,1)
         gymom(i,1)=qy1(i,1)*v1(i,1)
         gsed(i,1) =qy1(i,1)*c1(i,1)
         gmass(i,n-1)=qy1(i,n)
         gxmom(i,n-1)=qy1(i,n)*u1(i,n)
         gymom(i,n-1)=qy1(i,n)*v1(i,n)
         gsed(i,n-1) =qy1(i,n)*c1(i,n)
		end do
		!$omp end parallel do 
		
		!$omp parallel do(uyl,uyr,hyl,hyr,ayl,ayr,hystar,uystar,aqyl,aqyr,syl,syr,&
		!$omp &gmassl,gmassr,gxmoml,gxmomr,gymoml,gymomr)
		do i=2,m-1
         do j=2,n-2
            uyl=v1(i,j  )
            uyr=v1(i,j+1)
            hyl=h1(i,j)
            hyr=h1(i+1,j)
            ayl=sqrt(9.81*hyl)
            ayr=sqrt(9.81*hyr)
            hystar=0.5*(hyl+hyr)-0.25*(uyr-uyl)*(hyl+hyr)/(ayl+ayr)
            uystar=0.5*(uyl+uyr)-(hyr-hyl)*(ayl+ayr)/(hyl+hyr)
              
	      if(hystar.gt.hyl) then
               aqyl=sqrt(0.5*(hystar+hyl)*hystar/hyl**2)
            else
               aqyl=1.0
            endif
	      if(hystar.gt.hyr) then
               aqyr=sqrt(0.5*(hystar+hyr)*hystar/hyr**2)
            else
               aqyr=1.0
            endif
            syl=uyl-ayl*aqyl
            syr=uyr+ayr*aqyr

            gmassl=qy1(i,j  )              
            gmassr=qy1(i,j+1)              
            gxmoml=gmassl*u1(i,j  )
            gxmomr=gmassr*u1(i,j+1)
            gymoml=gmassl*v1(i,j  )
            gymomr=gmassr*v1(i,j+1)
!            gxmoml=qy1(i,j  )*u1(i,j  )
!            gxmomr=qy1(i,j+1)*u1(i,j+1)
!            gymoml=qy1(i,j  )*v1(i,j  )
!            gymomr=qy1(i,j+1)*v1(i,j+1)

            if(syl.gt.0.0) then
               gmass(i,j)=gmassl             
               gxmom(i,j)=gxmoml
               gymom(i,j)=gymoml
            elseif(syl.le.0.0.and.syr.ge.0.0) then
               gmass(i,j)=( syr*gmassl - syl*gmassr      &
                            + syl*syr*(zs1(i,j+1)-zs1(i,j)) )/(syr-syl)             
               gxmom(i,j)=( syr*gxmoml - syl*gxmomr      &
                            + syl*syr*(qx1(i,j+1)-qx1(i,j)) )/(syr-syl)
               gymom(i,j)=( syr*gymoml - syl*gymomr      &
                            + syl*syr*(qy1(i,j+1)-qy1(i,j)) )/(syr-syl)
            else if(syr.lt.0.0) then
               gmass(i,j)=gmassr             
               gxmom(i,j)=gxmomr
               gymom(i,j)=gymomr
            endif

            if(gmass(i,j).gt.0.0) then
               gsed(i,j)=gmass(i,j)*c1(i,j) 
            else if(gmass(i,j).lt.0.0) then
               gsed(i,j)=gmass(i,j)*c1(i,j+1) 
            else
               gsed(i,j)=gmass(i,j)*0.5*(c1(i,j)+c1(i,j+1)) 
            endif

            if(h1(i,j).le.hmin.and.zb(i,j).gt.zs1(i,j+1)) then
               gmass(i,j)=0.0
               gxmom(i,j)=0.0
               gymom(i,j)=0.0
               gsed(i,j) =0.0
	        endif
		    if(h1(i,j+1).le.hmin.and.zb(i,j+1).gt.zs1(i,j)) then
               gmass(i,j)=0.0
               gxmom(i,j)=0.0
               gymom(i,j)=0.0
               gsed(i,j) =0.0
	        endif
            if(h1(i,j).le.hmin.and.h1(i,j+1).le.hmin) then
               gmass(i,j)=0.0
               gxmom(i,j)=0.0
               gymom(i,j)=0.0
               gsed(i,j) =0.0
	        endif
	     enddo
         enddo
		 !$omp end parallel do 

      return
      end
      
      
      
      !#################################################################################################################

     subroutine cflux_HLL_2(m,n,hmin,h1,zb,zs1,fmass,gmass,fxmom,gxmom,   &
              fymom,gymom,fsed,gsed,epsilon,qx1,qy1,u1,v1,c1,dxc,dyc,dt)
       dimension  dxc(m),dyc(n)
       dimension  qx1(m,n),qy1(m,n),h1(m,n),zs1(m,n),      &
                  u1(m,n),v1(m,n),c1(m,n),zb(m,n)
       dimension  fmass(m-1,n-1),fxmom(m-1,n-1),fymom(m-1,n-1)
       dimension  gmass(m-1,n-1),gxmom(m-1,n-1),gymom(m-1,n-1)
       dimension  fsed(m-1,n-1),gsed(m-1,n-1)
       
       
       dimension  c_x(m-1,n-1)
       dimension  dfp(0:m,0:n) , dfm(0:m,0:n)
       dimension  hh1(0:m,0:n) 


       REAL CFL ,dx , a1 , a2 , a3

		!$omp parallel do 
       do j=2,n-1
         fmass(1,j)=qx1(1,j)
         fxmom(1,j)=qx1(1,j)*u1(1,j)
         fymom(1,j)=qx1(1,j)*v1(1,j)
         fsed(1,j) =qx1(1,j)*c1(1,j)
         fmass(m-1,j)=qx1(m,j)
         fxmom(m-1,j)=qx1(m,j)*u1(m,j)
         fymom(m-1,j)=qx1(m,j)*v1(m,j)
         fsed(m-1,j) =qx1(m,j)*c1(m,j)
       enddo   
	   !$omp end parallel do 
         
		 !$omp parallel do private(uxl,uxr,hxl,hxr,axl,axr,hxstar,uxstar,aqxl,aqxr,sxl,sxr,
		 !$omp &fmassl,fmassr,fxmoml,fxmomr,fymoml,fymomr)
         do j=2,n-1
         do i=2,m-2
            uxl=u1(i  ,j)
            uxr=u1(i+1,j)
            hxl=h1(i,j)
            hxr=h1(i+1,j)
            axl=sqrt(9.81*hxl)
            axr=sqrt(9.81*hxr)
            hxstar=0.5*(hxl+hxr)-0.25*(uxr-uxl)*(hxl+hxr)/(axl+axr)
            uxstar=0.5*(uxl+uxr)-(hxr-hxl)*(axl+axr)/(hxl+hxr)
              
	        if(hxstar.gt.hxl) then
               aqxl=sqrt(0.5*(hxstar+hxl)*hxstar/hxl**2)
            else
               aqxl=1.0
            endif
	        if(hxstar.gt.hxr) then
               aqxr=sqrt(0.5*(hxstar+hxr)*hxstar/hxr**2)
            else
               aqxr=1.0
            endif
            sxl=uxl-axl*aqxl
            sxr=uxr+axr*aqxr

            fmassl=qx1(i,  j)              
            fmassr=qx1(i+1,j)              
            fxmoml=fmassl*u1(i,  j)
            fxmomr=fmassr*u1(i+1,j)
            fymoml=fmassl*v1(i,  j)
            fymomr=fmassr*v1(i+1,j)
!            fxmoml=qx1(i,  j)*u1(i,  j)
!            fxmomr=qx1(i+1,j)*u1(i+1,j)
!            fymoml=qx1(i,  j)*v1(i,  j)
!            fymomr=qx1(i+1,j)*v1(i+1,j)

            if(sxl.gt.0.0) then
               fmass(i,j)=fmassl             
               fxmom(i,j)=fxmoml
               fymom(i,j)=fymoml
            elseif(sxl.le.0.0.and.sxr.ge.0.0) then
               fmass(i,j)=( sxr*fmassl - sxl*fmassr    &
                            + sxl*sxr*(zs1(i+1,j)-zs1(i,j)) )/(sxr-sxl)             
               fxmom(i,j)=( sxr*fxmoml - sxl*fxmomr    &
                            + sxl*sxr*(qx1(i+1,j)-qx1(i,j)) )/(sxr-sxl)
               fymom(i,j)=( sxr*fymoml - sxl*fymomr    &
                            + sxl*sxr*(qy1(i+1,j)-qy1(i,j)) )/(sxr-sxl)
            else if(sxr.lt.0.0) then
               fmass(i,j)=fmassr             
               fxmom(i,j)=fxmomr
               fymom(i,j)=fymomr
            endif

      enddo
      enddo
	  !$omp end parallel do 
      
      
        
	     !   do j=2,n-1
       !   do i=2,m-2
	     !     c_x(i,j)=0.5*( u1(i+1,j)+ u1(i,j) )     
	     !   enddo
	     !   enddo
       !  
       !  
       !  
       !     do j=2,n-1
       !     do i=1,m-1
	     !        hh1(i,j)= h1(i,j)
	     !        hh1(i,j)= h1(i,j)
	     !      enddo
	     !      enddo
	     !      
	     !      
	     !     do j=2,n-1
	     !     	  hh1(1,j)= 2.0*hh1(1,j) - hh1(2,j) 
	     !        hh1(1,j)= 2.0*hh1(m-1,j) - hh1(m-2,j)
	     !     	
       !     	  hh1(0,j)= hh1(1,j)
	     !        hh1(m,j)= hh1(m-1,j)
	     !      enddo
	     !      
       !  
       !  
       !     do j=2,n-1
       !     do i=0,m
	     !        dfp(i,j)= hh1(i+1,j)- hh1(i,j)
	     !        dfm(i,j)= hh1(i+1,j)- hh1(i,j)
	     !      enddo
	     !      enddo
       !  
       !      A1 =    -1.0
       !      A2 =    6.0
       !      A3 =    3.0
	     !     
       !  
       !  
	     !    do j=2,n-1
       !    do i=4,m-4
       !   
       !    	
       !    	if( c_x(i,j) >=0.0 	) then
	     !      ! fmass(i,j)= c_x(i,j)*(1.0/12.0*(-h1(i-1,j)+7.0*h1(i,j)+7.0*h1(i+1,j)-h1(i+2,j))  -phyn(dfp(i-2,j),dfp(i-1,j),dfp(i,j),dfp(i+1,j)) )
	     !      !
	     !      !   fmass(i,j)= c_x(i,j)*(   1.0/3.0*hh1(i-2,j)- 7.0/6.0*hh1(i-1,j) + 11.0/6.0*hh1(i,j)  ) 
	     !       
	     !      !  fmass(i,j)= c_x(i,j)*(   -1.0/2.0*hh1(i-1,j) + 3.0/2.0*hh1(i,j)  ) 
	     !      !   fmass(i,j)= c_x(i,j)*(   hh1(i,j)  ) 
	     !      
	     !        !  fmass(i,j)= (   1.0/3.0*hh1(i-2,j)- 7.0/6.0*hh1(i-1,j) + 11.0/6.0*hh1(i,j)  )  
	     !        ! fmass(i,j)=   -1.0/2.0*hh1(i-1,j) + 3.0/2.0*hh1(i,j)    
	     !      
	     !        fmass(i,j) = 1.0/8.0*( A1*hh1(I-1,j) + A2*hh1(I,j) + A3*hh1(I+1,j)  )
	     !      
	     !    	else      
	     !    	       
	     !       ! fmass(i,j)=c_x(i,j)*(1.0/12.0*(-h1(i-1,j)+7.0*h1(i,j)+7.0*h1(i+1,j)-h1(i+2,j)) 	 + phyn(dfm(i+2,j),dfm(i+1,j),dfm(i,j),dfm(i-1,j))    )
	     !       ! fmass(i,j)=  c_x(i,j)*(   -1.0/3.0*hh1(i+2,j)+ 7.0/6.0*hh1(i+1,j) - 11.0/6.0*hh1(i,j)  ) 	
	     !       ! fmass(i,j)=  c_x(i,j)*(   hh1(i+1,j) )
	     !       
	     !         fmass(i,j)=   hh1(i+1,j) 
	     !    	       
	     !    	endif
	     !    
	     !         
	     !          
	     !   	       
	     !     enddo
	     !     enddo
       !      
       !      
       !      dx = dxc(4) - dxc(3)
       !      
       !      do j=2,n-1
       !      do i=4,m-4 
       !    
       ! 
       !       IF( c_x(I,j)>= 0.0 ) THEN
       !                                !                              P         U           D
       !         CFL = ABS(c_x(I,j) * Dt / Dx)
       !         fmass(I,j) = FV_Phi_COMPUTE_UL(CFL , fmass(I,j) , hh1(I,j) , hh1(I-1,j) , hh1(I+1,j))
       !       
       !        !fmass(i,j)=   hh1(i,j)   
       !        
       !       ELSE 
       !       	
       !        !	CFL = ABS(Uc(I) * Dt / Dx)
       !        !	SsSS(I) = FV_Phi_COMPUTE_UL(CFL , SSsS(I) , U_old(I+1) , U_old(I+2) , U_old(I))
       !       	
       !       	  fmass(i,j)=   hh1(i+1,j) 
       !       	  
       !       END IF
       !       
       !     enddo
       !     enddo
       !     
       !     
       !      do j=2,n-1
       !      do i=4,m-4 
       !
       !         fmass(I,j) =  fmass(I,j)*c_x(I,j)
       !
       !     enddo
       !     enddo
       
       
       
       

       
			!$omp parallel do
           do j=2,n-1
           do i=2,m-2
	           c_x(i,j)=0.5*( u1(i+1,j)+ u1(i,j) )     
	         enddo
	         enddo
			 !$omp end parallel do 
          
          
            
			!$omp parallel do 
	        do j=2,n-1
            do i=4,m-4
           
            	
            	if( c_x(i,j) >=0.0 	) then
	            ! fmass(i,j)= c_x(i,j)*(1.0/12.0*(-h1(i-1,j)+7.0*h1(i,j)+7.0*h1(i+1,j)-h1(i+2,j))  -phyn(dfp(i-2,j),dfp(i-1,j),dfp(i,j),dfp(i+1,j)) )
	            !
	            !   fmass(i,j)= c_x(i,j)*(   1.0/3.0*hh1(i-2,j)- 7.0/6.0*hh1(i-1,j) + 11.0/6.0*hh1(i,j)  ) 
	             
	            !  fmass(i,j)= c_x(i,j)*(   -1.0/2.0*hh1(i-1,j) + 3.0/2.0*hh1(i,j)  ) 
	            !   fmass(i,j)= c_x(i,j)*(   hh1(i,j)  ) 
	            
	              !  fmass(i,j)= (   1.0/3.0*hh1(i-2,j)- 7.0/6.0*hh1(i-1,j) + 11.0/6.0*hh1(i,j)  )  
	              ! fmass(i,j)=   -1.0/2.0*hh1(i-1,j) + 3.0/2.0*hh1(i,j)    
	            
	             fxmom(i,j)=   qx1(i,j)  !-1.0/2.0*qx1(i-1,j) + 3.0/2.0*qx1(i,j) 
	            
	          	else      
	          	       
	             ! fmass(i,j)=c_x(i,j)*(1.0/12.0*(-h1(i-1,j)+7.0*h1(i,j)+7.0*h1(i+1,j)-h1(i+2,j)) 	 + phyn(dfm(i+2,j),dfm(i+1,j),dfm(i,j),dfm(i-1,j))    )
	             ! fmass(i,j)=  c_x(i,j)*(   -1.0/3.0*hh1(i+2,j)+ 7.0/6.0*hh1(i+1,j) - 11.0/6.0*hh1(i,j)  ) 	
	             ! fmass(i,j)=  c_x(i,j)*(   hh1(i+1,j) )
	             
	              fxmom(i,j)=  qx1(i+1,j)
	          	       
	          	endif
	          
	               
	                
	         	       
	           enddo
	           enddo
              !$omp end parallel do
              
        !      dx = dxc(4) - dxc(3)
        !      
        !      do j=2,n-1
        !      do i=4,m-4 
        !    
        ! 
        !       IF( c_x(I,j)>= 0.0 ) THEN
        !                                !                              P         U           D
        !         CFL = ABS(c_x(I,j) * Dt / Dx)
        !         fmass(I,j) = FV_Phi_COMPUTE_UL(CFL , fmass(I,j) , hh1(I,j) , hh1(I-1,j) , hh1(I+1,j))
        !       
        !        !fmass(i,j)=   hh1(i,j)   
        !        
        !       ELSE 
        !       	
        !        !	CFL = ABS(Uc(I) * Dt / Dx)
        !        !	SsSS(I) = FV_Phi_COMPUTE_UL(CFL , SSsS(I) , U_old(I+1) , U_old(I+2) , U_old(I))
        !       	
        !       	  fmass(i,j)=   hh1(i+1,j) 
        !       	  
        !       END IF
        !       
        !     enddo
        !     enddo
             
             !$omp parallel do 
              do j=2,n-1
              do i=4,m-4 
        
                 fxmom(I,j) =  fxmom(I,j)*c_x(I,j)
        
             enddo
             enddo
			 !$omp end parallel do 
     
     
     
     
       	 

		!$omp parallel do 
       do j=2,n-1
       do i=2,m-2
            if(fmass(i,j).gt.0.0) then
               fsed(i,j)=fmass(i,j)*c1(i,j) 
            else if(fmass(i,j).lt.0.0) then
               fsed(i,j)=fmass(i,j)*c1(i+1,j) 
            else
               fsed(i,j)=fmass(i,j)*0.5*(c1(i,j)+c1(i+1,j)) 
            endif

            if(h1(i,j).le.hmin.and.zb(i,j).gt.zs1(i+1,j)) then
               fmass(i,j)=0.0
               fxmom(i,j)=0.0
               fymom(i,j)=0.0
               fsed(i,j) =0.0
	        endif 
	   	    if(h1(i+1,j).le.hmin.and.zb(i+1,j).gt.zs1(i,j)) then
               fmass(i,j)=0.0
               fxmom(i,j)=0.0
               fymom(i,j)=0.0
               fsed(i,j) =0.0
            endif
	  	    if(h1(i,j).le.hmin.and.h1(i+1,j).le.hmin) then
               fmass(i,j)=0.0
               fxmom(i,j)=0.0
               fymom(i,j)=0.0
               fsed(i,j) =0.0
            endif
       
       enddo
       enddo
	   !$omp end parallel do 
       
       
		!$omp parallle do
        do i=2,m-1
         gmass(i,1)=qy1(i,1)
         gxmom(i,1)=qy1(i,1)*u1(i,1)
         gymom(i,1)=qy1(i,1)*v1(i,1)
         gsed(i,1) =qy1(i,1)*c1(i,1)
         gmass(i,n-1)=qy1(i,n)
         gxmom(i,n-1)=qy1(i,n)*u1(i,n)
         gymom(i,n-1)=qy1(i,n)*v1(i,n)
         gsed(i,n-1) =qy1(i,n)*c1(i,n)
		end do 
		!$omp end parallle do 
		
		!$omp parallel do private(uyl,uyr,hyl,hyr,ayl,ayr,hystar,uystar,aqyl,aqyr,syl,syr,&
		!$omp &gmassl,gmassr,gxmoml,gxmomr,gymoml,gymomr)
		do i=2,m-1
         do j=2,n-2
            uyl=v1(i,j  )
            uyr=v1(i,j+1)
            hyl=h1(i,j)
            hyr=h1(i+1,j)
            ayl=sqrt(9.81*hyl)
            ayr=sqrt(9.81*hyr)
            hystar=0.5*(hyl+hyr)-0.25*(uyr-uyl)*(hyl+hyr)/(ayl+ayr)
            uystar=0.5*(uyl+uyr)-(hyr-hyl)*(ayl+ayr)/(hyl+hyr)
              
	      if(hystar.gt.hyl) then
               aqyl=sqrt(0.5*(hystar+hyl)*hystar/hyl**2)
            else
               aqyl=1.0
            endif
	      if(hystar.gt.hyr) then
               aqyr=sqrt(0.5*(hystar+hyr)*hystar/hyr**2)
            else
               aqyr=1.0
            endif
            syl=uyl-ayl*aqyl
            syr=uyr+ayr*aqyr

            gmassl=qy1(i,j  )              
            gmassr=qy1(i,j+1)              
            gxmoml=gmassl*u1(i,j  )
            gxmomr=gmassr*u1(i,j+1)
            gymoml=gmassl*v1(i,j  )
            gymomr=gmassr*v1(i,j+1)
!            gxmoml=qy1(i,j  )*u1(i,j  )
!            gxmomr=qy1(i,j+1)*u1(i,j+1)
!            gymoml=qy1(i,j  )*v1(i,j  )
!            gymomr=qy1(i,j+1)*v1(i,j+1)

            if(syl.gt.0.0) then
               gmass(i,j)=gmassl             
               gxmom(i,j)=gxmoml
               gymom(i,j)=gymoml
            elseif(syl.le.0.0.and.syr.ge.0.0) then
               gmass(i,j)=( syr*gmassl - syl*gmassr      &
                            + syl*syr*(zs1(i,j+1)-zs1(i,j)) )/(syr-syl)             
               gxmom(i,j)=( syr*gxmoml - syl*gxmomr      &
                            + syl*syr*(qx1(i,j+1)-qx1(i,j)) )/(syr-syl)
               gymom(i,j)=( syr*gymoml - syl*gymomr      &
                            + syl*syr*(qy1(i,j+1)-qy1(i,j)) )/(syr-syl)
            else if(syr.lt.0.0) then
               gmass(i,j)=gmassr             
               gxmom(i,j)=gxmomr
               gymom(i,j)=gymomr
            endif

            if(gmass(i,j).gt.0.0) then
               gsed(i,j)=gmass(i,j)*c1(i,j) 
            else if(gmass(i,j).lt.0.0) then
               gsed(i,j)=gmass(i,j)*c1(i,j+1) 
            else
               gsed(i,j)=gmass(i,j)*0.5*(c1(i,j)+c1(i,j+1)) 
            endif

            if(h1(i,j).le.hmin.and.zb(i,j).gt.zs1(i,j+1)) then
               gmass(i,j)=0.0
               gxmom(i,j)=0.0
               gymom(i,j)=0.0
               gsed(i,j) =0.0
	        endif
		    if(h1(i,j+1).le.hmin.and.zb(i,j+1).gt.zs1(i,j)) then
               gmass(i,j)=0.0
               gxmom(i,j)=0.0
               gymom(i,j)=0.0
               gsed(i,j) =0.0
	        endif
            if(h1(i,j).le.hmin.and.h1(i,j+1).le.hmin) then
               gmass(i,j)=0.0
               gxmom(i,j)=0.0
               gymom(i,j)=0.0
               gsed(i,j) =0.0
	        endif
	     enddo
         enddo
		 !$omp end parallel do 

      return
      end
      
      
      subroutine cflux_HLL_3(m,n,hmin,h1,zb,zs1,fmass,gmass,fxmom,gxmom,   &
              fymom,gymom,fsed,gsed,epsilon,qx1,qy1,u1,v1,c1,dxc,dyc,dt)
       dimension  dxc(m),dyc(n)
       dimension  qx1(m,n),qy1(m,n),h1(m,n),zs1(m,n),      &
                  u1(m,n),v1(m,n),c1(m,n),zb(m,n)
       dimension  fmass(m-1,n-1),fxmom(m-1,n-1),fymom(m-1,n-1)
       dimension  gmass(m-1,n-1),gxmom(m-1,n-1),gymom(m-1,n-1)
       dimension  fsed(m-1,n-1),gsed(m-1,n-1)
       
       
       dimension  c_x(m-1,n-1)
       dimension  dfp(0:m,0:n) , dfm(0:m,0:n)
       dimension  hh1(0:m,0:n) 


       REAL CFL ,dx , a1 , a2 , a3
       
       
      
       
		!$omp parallel do 
       do j=2,n-1
         fmass(1,j)=qx1(1,j)
         fxmom(1,j)=qx1(1,j)*u1(1,j)
         fymom(1,j)=qx1(1,j)*v1(1,j)
         fsed(1,j) =qx1(1,j)*c1(1,j)
         fmass(m-1,j)=qx1(m,j)
         fxmom(m-1,j)=qx1(m,j)*u1(m,j)
         fymom(m-1,j)=qx1(m,j)*v1(m,j)
         fsed(m-1,j) =qx1(m,j)*c1(m,j)
       enddo   
       !$omp end parallel do 
       
       
       
         !$omp parallel do private(uxl,uxr,hxl,hxr,axl,axr,hxstar,uxstar,aqxl,aqxr,sxl,sxr,
		 !$omp &fmassl,fmassr,fxmoml,fxmomr,fymoml,fymomr)
         do j=2,n-1
         do i=2,m-2
            uxl=u1(i  ,j)
            uxr=u1(i+1,j)
            hxl=h1(i,j)
            hxr=h1(i+1,j)
            axl=sqrt(9.81*hxl)
            axr=sqrt(9.81*hxr)
            hxstar=0.5*(hxl+hxr)-0.25*(uxr-uxl)*(hxl+hxr)/(axl+axr)
            uxstar=0.5*(uxl+uxr)-(hxr-hxl)*(axl+axr)/(hxl+hxr)
              
	        if(hxstar.gt.hxl) then
               aqxl=sqrt(0.5*(hxstar+hxl)*hxstar/hxl**2)
            else
               aqxl=1.0
            endif
	        if(hxstar.gt.hxr) then
               aqxr=sqrt(0.5*(hxstar+hxr)*hxstar/hxr**2)
            else
               aqxr=1.0
            endif
            sxl=uxl-axl*aqxl
            sxr=uxr+axr*aqxr

            fmassl=qx1(i,  j)              
            fmassr=qx1(i+1,j)              
            fxmoml=fmassl*u1(i,  j)
            fxmomr=fmassr*u1(i+1,j)
            fymoml=fmassl*v1(i,  j)
            fymomr=fmassr*v1(i+1,j)
!            fxmoml=qx1(i,  j)*u1(i,  j)
!            fxmomr=qx1(i+1,j)*u1(i+1,j)
!            fymoml=qx1(i,  j)*v1(i,  j)
!            fymomr=qx1(i+1,j)*v1(i+1,j)

            if(sxl.gt.0.0) then
               fmass(i,j)=fmassl             
               fxmom(i,j)=fxmoml
               fymom(i,j)=fymoml
            elseif(sxl.le.0.0.and.sxr.ge.0.0) then
               fmass(i,j)=( sxr*fmassl - sxl*fmassr    &
                            + sxl*sxr*(zs1(i+1,j)-zs1(i,j)) )/(sxr-sxl)             
               fxmom(i,j)=( sxr*fxmoml - sxl*fxmomr    &
                            + sxl*sxr*(qx1(i+1,j)-qx1(i,j)) )/(sxr-sxl)
               fymom(i,j)=( sxr*fymoml - sxl*fymomr    &
                            + sxl*sxr*(qy1(i+1,j)-qy1(i,j)) )/(sxr-sxl)
            else if(sxr.lt.0.0) then
               fmass(i,j)=fmassr             
               fxmom(i,j)=fxmomr
               fymom(i,j)=fymomr
            endif

      enddo
      enddo
      !$omp end parallel do 
      
    !    write(*,*) " fmass" ,  fmass
    !     pause
      
        !$omp parallel do 
	    do j=2,n-1
        do i=2,m-2
	       c_x(i,j)=0.5*( u1(i+1,j)+ u1(i,j) )     
	    enddo
	    enddo
		!$omp end parallel do 
      
       
         
      
          A1 =    -1.0
          A2 =    6.0
          A3 =    3.0
	       
      
      !$omp parallel do
	    do j=2,n-1
        do i=4,m-4
       
        	
        	if( c_x(i,j) >0.0 	) then
	        ! fmass(i,j)= c_x(i,j)*(1.0/12.0*(-h1(i-1,j)+7.0*h1(i,j)+7.0*h1(i+1,j)-h1(i+2,j))  -phyn(dfp(i-2,j),dfp(i-1,j),dfp(i,j),dfp(i+1,j)) )
	        !
	        !   fmass(i,j)= c_x(i,j)*(   1.0/3.0*hh1(i-2,j)- 7.0/6.0*hh1(i-1,j) + 11.0/6.0*hh1(i,j)  ) 
	         
	        !  fmass(i,j)= c_x(i,j)*(   -1.0/2.0*hh1(i-1,j) + 3.0/2.0*hh1(i,j)  ) 
	        !   fmass(i,j)= c_x(i,j)*(   hh1(i,j)  ) 
	        
	         ! fxmom(i,j)= (   1.0/3.0*qx1(i-2,j)- 7.0/6.0*qx1(i-1,j) + 11.0/6.0*qx1(i,j)  )  
	         
	         
	        !   fxmom(i,j)=   qx1(i,j)  
	            fxmom(i,j)=   -1.0/2.0*qx1(i-1,j) + 3.0/2.0*qx1(i,j)    
	        
	        !  fxmom(i,j) = 1.0/8.0*( A1*qx1(I-1,j) + A2*qx1(I,j) + A3*qx1(I+1,j)  )
	        
	      	else      
	      	       
	         ! fmass(i,j)=c_x(i,j)*(1.0/12.0*(-h1(i-1,j)+7.0*h1(i,j)+7.0*h1(i+1,j)-h1(i+2,j)) 	 + phyn(dfm(i+2,j),dfm(i+1,j),dfm(i,j),dfm(i-1,j))    )
	         ! fmass(i,j)=  c_x(i,j)*(   -1.0/3.0*hh1(i+2,j)+ 7.0/6.0*hh1(i+1,j) - 11.0/6.0*hh1(i,j)  ) 	
	         ! fmass(i,j)=  c_x(i,j)*(   hh1(i+1,j) )
	         
	          
	         !  fxmom(i,j)= (   11.0/6.0*qx1(i+1,j)- 7.0/6.0*qx1(i+2,j) + 1.0/3.0*qx1(i+3,j)  ) 
	         
	          !  fxmom(i,j)=   qx1(i+1,j) 
	      	     
	      	    fxmom(i,j)=  3.0/2.0*qx1(i+1,j) - 1.0/2.0*qx1(i+2,j)
	      	       
	      	endif
	      
	           
	            
	     	       
	       enddo
	       enddo
           !$omp end parallel do 
          
           dx = dxc(9) !- dxc(9)
           
        !   write(*,*) "dx" ,dx
        !   pause
          
		  !$omp parallel do private(CFL)
          do j=2,n-1
          do i=4,m-4 
        
       
           IF( c_x(I,j)> 0.0 ) THEN
                                    !                              P         U           D
             CFL = ABS(c_x(I,j) * Dt / Dx)
             
           !  write(*,*) "DT" , DT , DX
           !  pause
             fxmom(I,j) = FV_Phi_COMPUTE_UL(CFL ,  fxmom(I,j) , qx1(I,j) , qx1(I-1,j) , qx1(I+1,j))
           
         !     write(*,*) "fxmom" , fxmom
         !     pause
           
           
           !  fxmom(i,j)=   qx1(i,j)   
            
           ELSE
           	
            	CFL = ABS(c_x(I,j) * Dt / Dx) 
             	fxmom(I,j) = FV_Phi_COMPUTE_UL(CFL , fxmom(I,j) , qx1(I+1,j) , qx1(I+2,j) , qx1(I,j))
            
           !	
           !	   fxmom(i,j)=   qx1(i+1,j) 
          
       
           	  
           END IF
           
         enddo
         enddo
		 !$omp end parallel do 
         
        !$omp parallel do
          do j=2,n-1
          do i=4,m-4 
     
             fxmom(I,j) =  fxmom(I,j)*c_x(I,j)
     
         enddo
         enddo
		!$omp end parallel do
     
     
        
     
     
       	 

		!$omp parallel do
       do j=2,n-1
       do i=2,m-2
            if(fmass(i,j).gt.0.0) then
               fsed(i,j)=fmass(i,j)*c1(i,j) 
            else if(fmass(i,j).lt.0.0) then
               fsed(i,j)=fmass(i,j)*c1(i+1,j) 
            else
               fsed(i,j)=fmass(i,j)*0.5*(c1(i,j)+c1(i+1,j)) 
            endif

            if(h1(i,j).le.hmin.and.zb(i,j).gt.zs1(i+1,j)) then
               fmass(i,j)=0.0
               fxmom(i,j)=0.0
               fymom(i,j)=0.0
               fsed(i,j) =0.0
	        endif 
	   	    if(h1(i+1,j).le.hmin.and.zb(i+1,j).gt.zs1(i,j)) then
               fmass(i,j)=0.0
               fxmom(i,j)=0.0
               fymom(i,j)=0.0
               fsed(i,j) =0.0
            endif
	  	    if(h1(i,j).le.hmin.and.h1(i+1,j).le.hmin) then
               fmass(i,j)=0.0
               fxmom(i,j)=0.0
               fymom(i,j)=0.0
               fsed(i,j) =0.0
            endif
       
       enddo
       enddo
	   !$omp end parallel do
       
       
		!$omp parallel do
        do i=2,m-1
         gmass(i,1)=qy1(i,1)
         gxmom(i,1)=qy1(i,1)*u1(i,1)
         gymom(i,1)=qy1(i,1)*v1(i,1)
         gsed(i,1) =qy1(i,1)*c1(i,1)
         gmass(i,n-1)=qy1(i,n)
         gxmom(i,n-1)=qy1(i,n)*u1(i,n)
         gymom(i,n-1)=qy1(i,n)*v1(i,n)
         gsed(i,n-1) =qy1(i,n)*c1(i,n)\
		end do
		!$omp end parallel do 
		
		!$omp parallel do private(uyl,uyr,hyl,hyr,ayl,ayr,hystar,uystar,aqyl,aqyr,syl,syr,&
		!$omp &gmassl,gmassr,gxmoml,gxmomr,gymoml,gymomr)
		do i=2,m-1
        do j=2,n-2
            uyl=v1(i,j  )
            uyr=v1(i,j+1)
            hyl=h1(i,j)
            hyr=h1(i+1,j)
            ayl=sqrt(9.81*hyl)
            ayr=sqrt(9.81*hyr)
            hystar=0.5*(hyl+hyr)-0.25*(uyr-uyl)*(hyl+hyr)/(ayl+ayr)
            uystar=0.5*(uyl+uyr)-(hyr-hyl)*(ayl+ayr)/(hyl+hyr)
              
	      if(hystar.gt.hyl) then
               aqyl=sqrt(0.5*(hystar+hyl)*hystar/hyl**2)
            else
               aqyl=1.0
            endif
	      if(hystar.gt.hyr) then
               aqyr=sqrt(0.5*(hystar+hyr)*hystar/hyr**2)
            else
               aqyr=1.0
            endif
            syl=uyl-ayl*aqyl
            syr=uyr+ayr*aqyr

            gmassl=qy1(i,j  )              
            gmassr=qy1(i,j+1)              
            gxmoml=gmassl*u1(i,j  )
            gxmomr=gmassr*u1(i,j+1)
            gymoml=gmassl*v1(i,j  )
            gymomr=gmassr*v1(i,j+1)
!            gxmoml=qy1(i,j  )*u1(i,j  )
!            gxmomr=qy1(i,j+1)*u1(i,j+1)
!            gymoml=qy1(i,j  )*v1(i,j  )
!            gymomr=qy1(i,j+1)*v1(i,j+1)

            if(syl.gt.0.0) then
               gmass(i,j)=gmassl             
               gxmom(i,j)=gxmoml
               gymom(i,j)=gymoml
            elseif(syl.le.0.0.and.syr.ge.0.0) then
               gmass(i,j)=( syr*gmassl - syl*gmassr      &
                            + syl*syr*(zs1(i,j+1)-zs1(i,j)) )/(syr-syl)             
               gxmom(i,j)=( syr*gxmoml - syl*gxmomr      &
                            + syl*syr*(qx1(i,j+1)-qx1(i,j)) )/(syr-syl)
               gymom(i,j)=( syr*gymoml - syl*gymomr      &
                            + syl*syr*(qy1(i,j+1)-qy1(i,j)) )/(syr-syl)
            else if(syr.lt.0.0) then
               gmass(i,j)=gmassr             
               gxmom(i,j)=gxmomr
               gymom(i,j)=gymomr
            endif

            if(gmass(i,j).gt.0.0) then
               gsed(i,j)=gmass(i,j)*c1(i,j) 
            else if(gmass(i,j).lt.0.0) then
               gsed(i,j)=gmass(i,j)*c1(i,j+1) 
            else
               gsed(i,j)=gmass(i,j)*0.5*(c1(i,j)+c1(i,j+1)) 
            endif

            if(h1(i,j).le.hmin.and.zb(i,j).gt.zs1(i,j+1)) then
               gmass(i,j)=0.0
               gxmom(i,j)=0.0
               gymom(i,j)=0.0
               gsed(i,j) =0.0
	        endif
		    if(h1(i,j+1).le.hmin.and.zb(i,j+1).gt.zs1(i,j)) then
               gmass(i,j)=0.0
               gxmom(i,j)=0.0
               gymom(i,j)=0.0
               gsed(i,j) =0.0
	        endif
            if(h1(i,j).le.hmin.and.h1(i,j+1).le.hmin) then
               gmass(i,j)=0.0
               gxmom(i,j)=0.0
               gymom(i,j)=0.0
               gsed(i,j) =0.0
	        endif
	     enddo
         enddo
		 !$omp end parallel do 

      return
      end
      
     subroutine cflux_weno(m,n,hmin,h1,zb,zs1,fmass,gmass,fxmom,gxmom,   &
              fymom,gymom,fsed,gsed,epsilon,qx1,qy1,u1,v1,c1,dxc,dyc,dt)
      
      
      
       dimension  dxc(m),dyc(n)
       dimension  qx1(m,n),qy1(m,n),h1(m,n),zs1(m,n),      &
                  u1(m,n),v1(m,n),c1(m,n),zb(m,n)
       dimension  fmass(m-1,n-1),fxmom(m-1,n-1),fymom(m-1,n-1)
       dimension  gmass(m-1,n-1),gxmom(m-1,n-1),gymom(m-1,n-1)
       dimension  fsed(m-1,n-1),gsed(m-1,n-1)


         
      
       !  call fflux(fh,fi,u,m,n)
          



	  !$omp parallel do 
      do j=2,n-1
         fmass(1,j)=qx1(1,j)
         fxmom(1,j)=qx1(1,j)*u1(1,j)
         fymom(1,j)=qx1(1,j)*v1(1,j)
         fsed(1,j) =qx1(1,j)*c1(1,j)
         fmass(m-1,j)=qx1(m,j)
         fxmom(m-1,j)=qx1(m,j)*u1(m,j)
         fymom(m-1,j)=qx1(m,j)*v1(m,j)
         fsed(m-1,j) =qx1(m,j)*c1(m,j)
       enddo   
       !$omp end parallel do 
	   
       call fflux(fmass,h1,u1,m,n)
       call fflux(fxmom,qx1,u1,m,n)
       call fflux(fymom,qx1,v1,m,n)   
         
		!$omp parallel do 
        do j=2,n-1  
        do i=2,m-2   
         if(fmass(i,j).gt.0.0) then
               fsed(i,j)=fmass(i,j)*c1(i,j) 
            else if(fmass(i,j).lt.0.0) then
               fsed(i,j)=fmass(i,j)*c1(i+1,j) 
            else
               fsed(i,j)=fmass(i,j)*0.5*(c1(i,j)+c1(i+1,j)) 
         endif
       
       enddo
       enddo  
        !$omp end parallel do  
	
		!$omp parallel do 
       do j=2,n-1  
       do i=2,m-2
           

            if(h1(i,j).le.hmin.and.zb(i,j).gt.zs1(i+1,j)) then
               fmass(i,j)=0.0
               fxmom(i,j)=0.0
               fymom(i,j)=0.0
               fsed(i,j) =0.0
	        endif 
	   	    if(h1(i+1,j).le.hmin.and.zb(i+1,j).gt.zs1(i,j)) then
               fmass(i,j)=0.0
               fxmom(i,j)=0.0
               fymom(i,j)=0.0
               fsed(i,j) =0.0
            endif
	  	    if(h1(i,j).le.hmin.and.h1(i+1,j).le.hmin) then
               fmass(i,j)=0.0
               fxmom(i,j)=0.0
               fymom(i,j)=0.0
               fsed(i,j) =0.0
            endif
       enddo
       enddo
       !$omp end parallel do 
       
       
		!$omp parallel do 
        do i=2,m-1
         gmass(i,1)=qy1(i,1)
         gxmom(i,1)=qy1(i,1)*u1(i,1)
         gymom(i,1)=qy1(i,1)*v1(i,1)
         gsed(i,1) =qy1(i,1)*c1(i,1)
         gmass(i,n-1)=qy1(i,n)
         gxmom(i,n-1)=qy1(i,n)*u1(i,n)
         gymom(i,n-1)=qy1(i,n)*v1(i,n)
         gsed(i,n-1) =qy1(i,n)*c1(i,n)
		end do 
		!$omp end parallel do 
		
		!$omp parallel do private(uyl,uyr,hyl,hyr,ayl,ayr,hystar,uystar,aqyl,aqyr,syl,syr,&
		!$omp &gmassl,gmassr,gxmoml,gxmomr,gymoml,gymomr)
		do i=2,m-1
        do j=2,n-2
            uyl=v1(i,j  )
            uyr=v1(i,j+1)
            hyl=h1(i,j)
            hyr=h1(i+1,j)
            ayl=sqrt(9.81*hyl)
            ayr=sqrt(9.81*hyr)
            hystar=0.5*(hyl+hyr)-0.25*(uyr-uyl)*(hyl+hyr)/(ayl+ayr)
            uystar=0.5*(uyl+uyr)-(hyr-hyl)*(ayl+ayr)/(hyl+hyr)
              
	      if(hystar.gt.hyl) then
               aqyl=sqrt(0.5*(hystar+hyl)*hystar/hyl**2)
            else
               aqyl=1.0
            endif
	      if(hystar.gt.hyr) then
               aqyr=sqrt(0.5*(hystar+hyr)*hystar/hyr**2)
            else
               aqyr=1.0
            endif
            syl=uyl-ayl*aqyl
            syr=uyr+ayr*aqyr

            gmassl=qy1(i,j  )              
            gmassr=qy1(i,j+1)              
            gxmoml=gmassl*u1(i,j  )
            gxmomr=gmassr*u1(i,j+1)
            gymoml=gmassl*v1(i,j  )
            gymomr=gmassr*v1(i,j+1)
!            gxmoml=qy1(i,j  )*u1(i,j  )
!            gxmomr=qy1(i,j+1)*u1(i,j+1)
!            gymoml=qy1(i,j  )*v1(i,j  )
!            gymomr=qy1(i,j+1)*v1(i,j+1)

            if(syl.gt.0.0) then
               gmass(i,j)=gmassl             
               gxmom(i,j)=gxmoml
               gymom(i,j)=gymoml
            elseif(syl.le.0.0.and.syr.ge.0.0) then
               gmass(i,j)=( syr*gmassl - syl*gmassr      &
                            + syl*syr*(zs1(i,j+1)-zs1(i,j)) )/(syr-syl)             
               gxmom(i,j)=( syr*gxmoml - syl*gxmomr      &
                            + syl*syr*(qx1(i,j+1)-qx1(i,j)) )/(syr-syl)
               gymom(i,j)=( syr*gymoml - syl*gymomr      &
                            + syl*syr*(qy1(i,j+1)-qy1(i,j)) )/(syr-syl)
            else if(syr.lt.0.0) then
               gmass(i,j)=gmassr             
               gxmom(i,j)=gxmomr
               gymom(i,j)=gymomr
            endif

            if(gmass(i,j).gt.0.0) then
               gsed(i,j)=gmass(i,j)*c1(i,j) 
            else if(gmass(i,j).lt.0.0) then
               gsed(i,j)=gmass(i,j)*c1(i,j+1) 
            else
               gsed(i,j)=gmass(i,j)*0.5*(c1(i,j)+c1(i,j+1)) 
            endif

            if(h1(i,j).le.hmin.and.zb(i,j).gt.zs1(i,j+1)) then
               gmass(i,j)=0.0
               gxmom(i,j)=0.0
               gymom(i,j)=0.0
               gsed(i,j) =0.0
	        endif
		    if(h1(i,j+1).le.hmin.and.zb(i,j+1).gt.zs1(i,j)) then
               gmass(i,j)=0.0
               gxmom(i,j)=0.0
               gymom(i,j)=0.0
               gsed(i,j) =0.0
	        endif
            if(h1(i,j).le.hmin.and.h1(i,j+1).le.hmin) then
               gmass(i,j)=0.0
               gxmom(i,j)=0.0
               gymom(i,j)=0.0
               gsed(i,j) =0.0
	        endif
	     enddo
         enddo
		 !$omp end parallel do

      return
      end



      subroutine cflux_Ying(m,n,hmin,h1,zb,zs1,fmass,gmass,fxmom,gxmom,   &
              fymom,gymom,fsed,gsed,epsilon,qx1,qy1,u1,v1,c1,dxc,dyc,dt)
       dimension  dxc(m),dyc(n)
       dimension  qx1(m,n),qy1(m,n),h1(m,n),zs1(m,n),    &
                  u1(m,n),v1(m,n),c1(m,n),zb(m,n)
       dimension  fmass(m-1,n-1),fxmom(m-1,n-1),fymom(m-1,n-1)
       dimension  gmass(m-1,n-1),gxmom(m-1,n-1),gymom(m-1,n-1)
       dimension  fsed(m-1,n-1),gsed(m-1,n-1)

		!$omp parallel do 
        do j=2,n-1
         fmass(1,j)=qx1(1,j)
         fxmom(1,j)=qx1(1,j)*u1(1,j)
         fymom(1,j)=qx1(1,j)*v1(1,j)
         fsed(1,j) =qx1(1,j)*c1(1,j)
         fmass(m-1,j)=qx1(m,j)
         fxmom(m-1,j)=qx1(m,j)*u1(m,j)
         fymom(m-1,j)=qx1(m,j)*v1(m,j)
         fsed(m-1,j) =qx1(m,j)*c1(m,j)
		end do 
		!$omp end parallel do
		
		!$omp parallel do private(vif,hif)
        do j=2,n-1
        do i=2,m-2
		
	    if(qx1(i+1,j).gt.epsilon.and.qx1(i,j).gt.epsilon) then 
           fmass(i,j)=qx1(i,j)
           fxmom(i,j)=qx1(i,j)*u1(i,j)
           fymom(i,j)=qx1(i,j)*v1(i,j)
           fsed(i,j) =qx1(i,j)*c1(i,j)
	    else if(qx1(i+1,j).lt.-epsilon.and.qx1(i,j).lt.-epsilon) then
           fmass(i,j)=qx1(i+1,j)
           fxmom(i,j)=qx1(i+1,j)*u1(i+1,j)
           fymom(i,j)=qx1(i+1,j)*v1(i+1,j)
           fsed(i,j) =qx1(i+1,j)*c1(i+1,j)
        else
           vif=-9.8*dt/dxc(i)*(zs1(i+1,j)-zs1(i,j))
	       hif=(h1(i+1,j)+h1(i,j))*0.5
           fmass(i,j)=(qx1(i+1,j)+qx1(i,j))/2.0+vif*hif         
           fxmom(i,j)=(qx1(i+1,j)*u1(i+1,j)+qx1(i,j)*u1(i,j))/3.0
           fymom(i,j)=(qx1(i+1,j)*v1(i+1,j)+qx1(i,j)*v1(i,j))/3.0
           fsed(i,j) =(qx1(i+1,j)*c1(i+1,j)+qx1(i,j)*c1(i,j))/2.0
        endif
		
        if(h1(i,j).le.hmin.and.zb(i,j).gt.zs1(i+1,j)) then
           fmass(i,j)=0.0
           fxmom(i,j)=0.0
           fymom(i,j)=0.0
           fsed(i,j) =0.0
	    endif 
		if(h1(i+1,j).le.hmin.and.zb(i+1,j).gt.zs1(i,j)) then
           fmass(i,j)=0.0
           fxmom(i,j)=0.0
           fymom(i,j)=0.0
           fsed(i,j) =0.0
        endif
		if(h1(i,j).le.hmin.and.h1(i+1,j).le.hmin) then
           fmass(i,j)=0.0
           fxmom(i,j)=0.0
           fymom(i,j)=0.0
           fsed(i,j) =0.0
        endif
		
        enddo
        enddo
		!$omp end parallel do 
			
		!$omp parallel do 
        do i=2,m-1
         gmass(i,1)=qy1(i,1)
         gxmom(i,1)=qy1(i,1)*u1(i,1)
         gymom(i,1)=qy1(i,1)*v1(i,1)
         gsed(i,1) =qy1(i,1)*c1(i,1)
         gmass(i,n-1)=qy1(i,n)
         gxmom(i,n-1)=qy1(i,n)*u1(i,n)
         gymom(i,n-1)=qy1(i,n)*v1(i,n)
         gsed(i,n-1) =qy1(i,n)*c1(i,n)
		enddo
		!$omp end parallel do 
		
		!$omp parallel do private(vif,hif)
        do i=2,m-1
        do j=2,n-2
		
        if(qy1(i,j+1).gt.epsilon.and.qy1(i,j).gt.epsilon) then 
           gmass(i,j)=qy1(i,j)
           gxmom(i,j)=qy1(i,j)*u1(i,j)
           gymom(i,j)=qy1(i,j)*v1(i,j)
           gsed(i,j) =qy1(i,j)*c1(i,j)
	    else if(qy1(i,j+1).lt.-epsilon.and.qy1(i,j).lt.-epsilon) then
           gmass(i,j)=qy1(i,j+1)
           gxmom(i,j)=qy1(i,j+1)*u1(i,j+1)
           gymom(i,j)=qy1(i,j+1)*v1(i,j+1)
           gsed(i,j) =qy1(i,j+1)*c1(i,j+1)
        else
           vif=-9.8*dt/dyc(j)*(zs1(i,j+1)-zs1(i,j))
	       hif=(h1(i,j+1)+h1(i,j))*0.5
           gmass(i,j)=(qy1(i,j+1)+qy1(i,j))/2.0+vif*hif         
           gxmom(i,j)=(qy1(i,j+1)*u1(i,j+1)+qy1(i,j)*u1(i,j))/3.0
           gymom(i,j)=(qy1(i,j+1)*v1(i,j+1)+qy1(i,j)*v1(i,j))/3.0
           gsed(i,j) =(qy1(i,j+1)*c1(i,j+1)+qy1(i,j)*c1(i,j))/2.0
        endif
        if(h1(i,j).le.hmin.and.zb(i,j).gt.zs1(i,j+1)) then
           gmass(i,j)=0.0
           gxmom(i,j)=0.0
           gymom(i,j)=0.0
           gsed(i,j) =0.0
	    endif
		if(h1(i,j+1).le.hmin.and.zb(i,j+1).gt.zs1(i,j)) then
           gmass(i,j)=0.0
           gxmom(i,j)=0.0
           gymom(i,j)=0.0
           gsed(i,j) =0.0
	    endif
          if(h1(i,j).le.hmin.and.h1(i,j+1).le.hmin) then
           gmass(i,j)=0.0
           gxmom(i,j)=0.0
           gymom(i,j)=0.0
           gsed(i,j) =0.0
	    endif
	    enddo
        enddo
		!$omp end parallel do 

      return
	end


!     =================================================================
!     This subroutine handles the steep slope case
!           Made by Dr. Weiming Wu, ..............
      subroutine bedrepos(m,n,x,y,zb,zs,idix,dxc,dyc,it,alenii,   &
                          xxrep,yyrep,volrep,idixrep,hmin)
!    repose angle is considered here to check the bed profile
!    if the bed slope is larger than the repose angle, the bed will collapse 
!    due to the gravity  (Horizontal 2-D or 3D)
!     =================================================================
      common	reposedry,reposewet
	  dimension x(m,n),y(m,n),zb(m,n),zs(m,n),idix(m,n),dxc(m),dyc(n)
      dimension alenii(m,n,8),xxrep(m,n,8),yyrep(m,n,8),volrep(m,n,8),idixrep(m,n,8)   
      dimension slpbed(8),zbrep(8),dzbrep(8),zsrep(8),repose(8)  

!      reposedry=5.67
!!!!!!      reposedry=5.145	!6.5			!2.5			!.......by He 200905
!!!!!!      reposewet=0.65		!0.65

      iter=0
 101  nodesum=0
      iter=iter+1

      ntimesed=it+iter

      if(mod(ntimesed,4).eq.0) then
         istart=2
         ifinal=m-1
         iskip =1
         jstart=2
         jfinal=n-1
         jskip =1
      elseif(mod(ntimesed,4).eq.1) then
         istart=m-1
         ifinal=2
         iskip =-1
         jstart=2
         jfinal=n-1
         jskip =1
      elseif(mod(ntimesed,4).eq.2) then
         istart=2
         ifinal=m-1
         iskip =1
         jstart=n-1
         jfinal=2
         jskip =-1
      elseif(mod(ntimesed,4).eq.3) then
         istart=m-1
         ifinal=2
         iskip =-1
         jstart=n-1
         jfinal=2
         jskip =-1
      endif

	    !$omp parallel do private(zbrep,noderep1,noderep2,arearep,alenrep,zbavrep,kl,dzbrepii), reduction(+:nodesum)
        do i=istart,ifinal,iskip
        do j=jstart,jfinal,jskip

         if(idix(i,j).eq.1) then
            zbrep(1)   =  zb(i-1,j-1)
            zbrep(2)   =  zb(i,j-1)
            zbrep(3)   =  zb(i+1,j-1)
            zbrep(4)   =  zb(i-1,j)
            zbrep(5)   =  zb(i+1,j)
            zbrep(6)   =  zb(i-1,j+1)
            zbrep(7)   =  zb(i,j+1)
            zbrep(8)   =  zb(i+1,j+1)
            zsrep(1)   =  zs(i-1,j-1)
            zsrep(2)   =  zs(i,j-1)
            zsrep(3)   =  zs(i+1,j-1)
            zsrep(4)   =  zs(i-1,j)
            zsrep(5)   =  zs(i+1,j)
            zsrep(6)   =  zs(i-1,j+1)
            zsrep(7)   =  zs(i,j+1)
            zsrep(8)   =  zs(i+1,j+1)

 105        noderep1=0
            arearep=0.0
            alenrep=0.0
            zbavrep=0.0
			
            do kl=1,8
              slpbed(kl)=(zb(i,j)-zbrep(kl))/alenii(i,j,kl)
              if(zs(i,j).gt.zb(i,j)+hmin) then
                repose(kl)=reposewet
              else
                repose(kl)=reposedry
              endif
!              if((zs(i,j).lt.zb(i,j)+hmin).
!     &                   and.(zsrep(kl).lt.zbrep(kl)+hmin)) then
!                repose(kl)=reposedry
!              else
!                repose(kl)=reposewet
!              endif
              if(slpbed(kl).gt.repose(kl)+0.0001    &
                         .and.idixrep(i,j,kl).eq.1) then 
                noderep1=noderep1+1
                arearep=arearep+volrep(i,j,kl)
                alenrep=alenrep+alenii(i,j,kl)*volrep(i,j,kl)*repose(kl)
                zbavrep=zbavrep+zbrep(kl)*volrep(i,j,kl)
              endif
            enddo

            if(noderep1.gt.0) then
             dzbrepii=(-zb(i,j)*arearep+zbavrep+alenrep)/(arearep+dxc(i)*dyc(j))
!              write(*,*) 'L,',i,j, dzbrepii
             do kl=1,8
               if(slpbed(kl).gt.repose(kl)+0.0001      &
                         .and.idixrep(i,j,kl).eq.1) then
                 dzbrep(kl)=zb(i,j)+dzbrepii-zbrep(kl)   &
                                   -repose(kl)*alenii(i,j,kl)
               else
                 dzbrep(kl)=0.0
               endif
             enddo

             do kl=1,8
               zbrep(kl)=zbrep(kl)+dzbrep(kl)*0.5
             enddo
             zb(i,j)=zb(i,j)+dzbrepii*0.5
            endif

            nodesum=nodesum+noderep1

!  ....................................................................
 106        noderep2=0
            arearep=0.0
            alenrep=0.0
            zbavrep=0.0
            do kl=1,8
              slpbed(kl)=(zb(i,j)-zbrep(kl))/alenii(i,j,kl)
              if(zsrep(kl).gt.zbrep(kl)+hmin) then
                repose(kl)=reposewet
              else
                repose(kl)=reposedry
              endif
!              if((zs(i,j).lt.zb(i,j)+hmin).
!     &                   and.(zsrep(kl).lt.zbrep(kl)+hmin)) then
!                repose(kl)=reposedry
!              else
!                repose(kl)=reposewet
!              endif
              if(slpbed(kl).lt.-repose(kl)-0.0001          &
                               .and.idixrep(i,j,kl).eq.1) then 
                noderep2=noderep2+1
                arearep=arearep+volrep(i,j,kl)
                alenrep=alenrep+alenii(i,j,kl)*volrep(i,j,kl)*repose(kl)
                zbavrep=zbavrep+zbrep(kl)*volrep(i,j,kl)
              endif
            enddo

            if(noderep2.gt.0) then
             dzbrepii=(-zb(i,j)*arearep+zbavrep-alenrep)    &
                         /(arearep+dxc(i)*dyc(j))
             do kl=1,8
              if(slpbed(kl).lt.-repose(kl)-0.0001        &
                               .and.idixrep(i,j,kl).eq.1) then
                dzbrep(kl)=zb(i,j)+dzbrepii-zbrep(kl)     &
                                  +repose(kl)*alenii(i,j,kl)
              else
                dzbrep(kl)=0.0
              endif
             enddo

             do kl=1,8
               zbrep(kl)=zbrep(kl)+dzbrep(kl)*0.5
             enddo
               zb(i,j)=zb(i,j)+dzbrepii*0.5
            endif

            nodesum=nodesum+noderep2

            zb(i-1,j-1)=zbrep(1)
            zb(i  ,j-1)=zbrep(2)
            zb(i+1,j-1)=zbrep(3)
            zb(i-1,j  )=zbrep(4)
            zb(i+1,j  )=zbrep(5)
            zb(i-1,j+1)=zbrep(6)
            zb(i  ,j+1)=zbrep(7)
            zb(i+1,j+1)=zbrep(8)

         endif

      enddo
      enddo
	  !$omp end parallel do 
	  
!         write(*,*) 'nodesum=', nodesum

      if(iter.gt.10) goto 102
      if(nodesum.ge.1) goto 101 
 102  continue

      return
      end


!     =================================================================
!     This subroutine handles the steep slope case
!           Made by Dr. Weiming Wu, ..............
      subroutine bedrepos0(m,n,x,y,zb,idix,dxc,dyc,it,alenii,   &
                          xxrep,yyrep,volrep,idixrep)
!    repose angle is considered here to check the bed profile
!    if the bed slope is larger than the repose angle, the bed will collapse 
!    due to the gravity  (Horizontal 2-D or 3D)
!     =================================================================
      dimension x(m,n),y(m,n),zb(m,n),idix(m,n),dxc(m),dyc(n)
      dimension alenii(m,n,8),xxrep(m,n,8),yyrep(m,n,8),volrep(m,n,8),idixrep(m,n,8)   
      dimension slpbed(8),zbrep(8),dzbrep(8)  

!!      repose=0.65
	  repose=1.0

      iter=0
 101  nodesum=0
      iter=iter+1

      ntimesed=it+iter

      if(mod(ntimesed,4).eq.0) then
         istart=2
         ifinal=m-1
         iskip =1
         jstart=2
         jfinal=n-1
         jskip =1
      elseif(mod(ntimesed,4).eq.1) then
         istart=m-1
         ifinal=2
         iskip =-1
         jstart=2
         jfinal=n-1
         jskip =1
      elseif(mod(ntimesed,4).eq.2) then
         istart=2
         ifinal=m-1
         iskip =1
         jstart=n-1
         jfinal=2
         jskip =-1
      elseif(mod(ntimesed,4).eq.3) then
         istart=m-1
         ifinal=2
         iskip =-1
         jstart=n-1
         jfinal=2
         jskip =-1
      endif

	  !$omp parallel do private(zbrep,noderep1,noderep2,arearep,alenrep,zbavrep,kl,dzbrepii), reduction(+:nodesum)
      do i=istart,ifinal,iskip
      do j=jstart,jfinal,jskip

         if(idix(i,j).eq.1) then
            zbrep(1)   =  zb(i-1,j-1)
            zbrep(2)   =  zb(i,j-1)
            zbrep(3)   =  zb(i+1,j-1)
            zbrep(4)   =  zb(i-1,j)
            zbrep(5)   =  zb(i+1,j)
            zbrep(6)   =  zb(i-1,j+1)
            zbrep(7)   =  zb(i,j+1)
            zbrep(8)   =  zb(i+1,j+1)

 105        noderep1=0
            arearep=0.0
            alenrep=0.0
            zbavrep=0.0
            do kl=1,8
              slpbed(kl)=(zb(i,j)-zbrep(kl))/alenii(i,j,kl)
              if(slpbed(kl).gt.repose+0.0001       &
                         .and.idixrep(i,j,kl).eq.1) then 
                noderep1=noderep1+1
                arearep=arearep+volrep(i,j,kl)
                alenrep=alenrep+alenii(i,j,kl)*volrep(i,j,kl)
                zbavrep=zbavrep+zbrep(kl)*volrep(i,j,kl)
              endif
            enddo

            if(noderep1.gt.0) then
             dzbrepii=(-zb(i,j)*arearep+zbavrep+repose*alenrep)     &
                      /(arearep+dxc(i)*dyc(j))
!              write(*,*) 'L,',i,j, dzbrepii
             do kl=1,8
               if(slpbed(kl).gt.repose+0.0001              &
                          .and.idixrep(i,j,kl).eq.1) then
                 dzbrep(kl)=zb(i,j)+dzbrepii-zbrep(kl)     &
                                    -repose*alenii(i,j,kl)
               else
                 dzbrep(kl)=0.0
               endif
             enddo

             do kl=1,8
               zbrep(kl)=zbrep(kl)+dzbrep(kl)
             enddo
             zb(i,j)=zb(i,j)+dzbrepii
            endif

            nodesum=nodesum+noderep1

!  ....................................................................
 106        noderep2=0
            arearep=0.0
            alenrep=0.0
            zbavrep=0.0
            do kl=1,8
              slpbed(kl)=(zb(i,j)-zbrep(kl))/alenii(i,j,kl)
              if(slpbed(kl).lt.-repose-0.0001      &
                              .and.idixrep(i,j,kl).eq.1) then 
                noderep2=noderep2+1
                arearep=arearep+volrep(i,j,kl)
                alenrep=alenrep+alenii(i,j,kl)*volrep(i,j,kl)
                zbavrep=zbavrep+zbrep(kl)*volrep(i,j,kl)
              endif
            enddo

            if(noderep2.gt.0) then
             dzbrepii=(-zb(i,j)*arearep+zbavrep-repose*alenrep)   &
                         /(arearep+dxc(i)*dyc(j))
             do kl=1,8
              if(slpbed(kl).lt.-repose-0.0001      &
                              .and.idixrep(i,j,kl).eq.1) then
                dzbrep(kl)=zb(i,j)+dzbrepii-zbrep(kl)   &
                                  +repose*alenii(i,j,kl)
              else
                dzbrep(kl)=0.0
              endif
             enddo

             do kl=1,8
               zbrep(kl)=zbrep(kl)+dzbrep(kl)
             enddo
             zb(i,j)=zb(i,j)+dzbrepii
            endif

            nodesum=nodesum+noderep2

            zb(i-1,j-1)=zbrep(1)
            zb(i  ,j-1)=zbrep(2)
            zb(i+1,j-1)=zbrep(3)
            zb(i-1,j  )=zbrep(4)
            zb(i+1,j  )=zbrep(5)
            zb(i-1,j+1)=zbrep(6)
            zb(i  ,j+1)=zbrep(7)
            zb(i+1,j+1)=zbrep(8)

         endif

      enddo
      enddo
	  !$omp end parallel do 
!         write(*,*) 'nodesum=', nodesum

      if(iter.gt.50) goto 102
      if(nodesum.ge.1) goto 101 
 102  continue

      return
      end
!.##################################################################################
      subroutine dzdxyminmod(m,n,dzdxmin,dzdymin,x,y,zs)
       dimension  x(m,n),y(m,n),zs(m,n)
       dimension  dzdxmin(m,n), dzdymin(m,n)
        
		!$omp parallel do 
        do j=1,n
         dzdxmin(1,j)=(zs(2,j)-zs(1,j))/(x(2,j)-x(1,j))          
         dzdxmin(m,j)=(zs(m,j)-zs(m-1,j))/(x(m,j)-x(m-1,j))          
		enddo
		!$omp end parallel do 
		
		!$omp parallel do private(dzdxup,dzdxlw)
        do j=1,n
         do i=2,m-1
            dzdxup=(zs(i  ,j)-zs(i-1,j))/(x(i  ,j)-x(i-1,j))          
			dzdxlw=(zs(i+1,j)-zs(i  ,j))/(x(i+1,j)-x(i  ,j))          
			if(dzdxup*dzdxlw.gt.0.0) then
			   if(dzdxup.gt.0.0) then
			      dzdxmin(i,j)=min(dzdxup, dzdxlw)
			   else
			      dzdxmin(i,j)=max(dzdxup, dzdxlw)
               endif
            else
			   dzdxmin(i,j)=0.0
            endif
         enddo
	   enddo
	   !$omp end parallel do 

		!$omp parallel do 
       do i=1,m
         dzdymin(i,1)=(zs(i,2)-zs(i,1))/(y(i,2)-y(i,1))          
         dzdymin(i,n)=(zs(i,n)-zs(i,n-1))/(y(i,n)-y(i,n-1)) 
		enddo
		!$omp end parallel do 
		
		!$omp parallel do private(dzdyup,dzdylw)
         do j=2,n-1
            dzdyup=(zs(i,j  )-zs(i,j-1))/(y(i,j  )-y(i,j-1))          
			dzdylw=(zs(i,j+1)-zs(i  ,j))/(y(i,j+1)-y(i  ,j))          
			if(dzdyup*dzdylw.gt.0.0) then
			   if(dzdyup.gt.0.0) then
			      dzdymin(i,j)=min(dzdyup, dzdylw)
			   else
			      dzdymin(i,j)=max(dzdyup, dzdylw)
               endif
            else
			   dzdymin(i,j)=0.0
            endif
         enddo
	   enddo
	   !$omp end parallel do 

       return
      end

!#############################################################################
	Subroutine Output0(m,n,time,x,y,zb,zs2,h2,u2,v2,qx2,qy2,c2)
       
	   dimension  x(m,n),y(m,n),zb(m,n),zs2(m,n),h2(m,n),u2(m,n),v2(m,n),qx2(m,n),qy2(m,n),c2(m,n)

!..................................................Output Initial Information............................
		   OPEN(301,FILE="Z0.dat")
		   write(301,*)	'title= " Time= ', time, ' s " ' 
		   write(301,*) 'VARIABLES ="X","Y","Z","Zs","h","U","V","W","Qx","Qy","Concen"'
		   write(301,*) 'ZONE T="contour"',' I=',m,' J=',n,' F=point'
           do j=1,n
		     do i=1,m
             
               write(301,9970) x(i,j),y(i,j),zb(i,j),zs2(i,j),h2(i,j),   &
                              u2(i,j),v2(i,j),0.0,qx2(i,j),qy2(i,j),c2(i,j)
             enddo
           enddo
9970	   format(2F12.4,9E18.8)
		   CLOSE(301)
!...........................................................................................................

	end subroutine
	
	
     !	 subroutine fflux(fh,ffi,uu,m,n)
     !  
     !!  subroutine fflux(fmass,h1,u1,m,n)
     !  
     !  
     !!
     !!     fflux --- give fi, u, u+, u- to calculate fh at i+1/2
     !!     f, dfp, dfm are dummy array
     !!
     !		                 
     !    INTEGER :: m , n
     !
     !    
     !	  REAL :: fh(1:m-1,1:n-1), ffi(1:m,1:n) ,  uu(1:m,1:n) 
     !	  
     !	  REAL :: up(0:m+1,0:n+1), um(0:m+1,0:n+1) , dfp(0:m+1,0:n+1)  ,  dfm(0:m+1,0:n+1)
     !	  
     !	  REAL :: f(1:m,1:n)
     !
     !		REAL :: fi(0:m+1,0:n+1),u(0:m+1,0:n+1)
     !
     !    
     !    INTEGER :: I , J
     !    
     !    	
     !    do i=1,m
     !	  do j=1,n
     !	    u(i,j)= uu(i,j)  
     !	  enddo
     !	  enddo
     !	  
     !	  
     !	  do j=1,n-1
     !	    u(1,j)= 2.0*uu(1,j) -  uu(2,j) 
     !	    u(m,j)= 2.0*uu(m,j) -  uu(m-1,j)
     !	    
     !	    
     !	    u(0,j)=  u(1,j) 
     !	    u(m+1,j)= u(m,j) 
     !	    
     !	  enddo
     !	  
     !	  
     !	  do i=1,m
     !	  do j=1,n-1
     !	    fi(i,j)= ffi(i,j)  
     !	  enddo
     !	  enddo
     !	  
     !	  
     !	  do j=1,n-1    
     !	    fi(1,j)= 2.0*ffi(1,j) -  ffi(2,j) 
     !	    fi(m,j)= 2.0*ffi(m,j) -  ffi(m-1,j)
     !	    
     !	    
     !	    fi(0,j)=  fi(1,j) 
     !	    fi(m+1,j)= fi(m,j) 
     !	  enddo
     !	
     !    	
     !    	
     !    	
     !    do i=1,m 		                 
     !	  do j=1,n-1  
     !	    f(i,j)=u(i,j)*fi(i,j)
     !	  enddo
     !	  enddo
     !	  
     !	  
     ! 
     !	  
     !	  
     !	  
     !	  
     !	  
     !	  do i=0,m+1
     !	  do j=1,n-1
     !	    up(i,j)=0.5*(u(i,j)+abs(u(i,j)))
     !	    um(i,j)=0.5*(u(i,j)-abs(u(i,j)))
     !	  enddo
     !	  enddo
     !	  
     !	  
     !	  do i=0 , m
     !		do j=1 , n-1
     !	    dfp(i,j)=up(i+1,j)*fi(i+1,j)-up(i,j)*fi(i,j)
     !	    dfm(i,j)=um(i+1,j)*fi(i+1,j)-um(i,j)*fi(i,j)
     !	  enddo
     !	  enddo
     !	  
     !	
     !    do i=2,m-2
     !    do j=1,n-1 	
     !     	
     !	    fh(i,j)=1.0/12.0*(-f(i-1,j)+7.0*f(i,j)+7.0*f(i+1,j)-f(i+2,j)) &
     !		       -phyn(dfp(i-2,j),dfp(i-1,j),dfp(i,j),dfp(i+1,j)) &
     !		       +phyn(dfm(i+2,j),dfm(i+1,j),dfm(i,j),dfm(i-1,j))
     !	  enddo
     !	  enddo
     !
     !	  end subroutine fflux
	  
	 subroutine fflux(fh,fi,u,up,um,dfp,dfm,m,n)
!
!     fflux --- give fi, u, u+, u- to calculate fh at i+1/2
!     f, dfp, dfm are dummy array
!
		                 
    INTEGER,INTENT(IN) :: m , n

    INTEGER,PARAMETER  :: DP = SELECTED_REAL_KIND(P=15)
	  REAL(KIND=DP) :: fh(-2:m+3,-2:n+3),fi(-2:m+3,-2:n+3),u(-2:m+3,-2:n+3), &
	                   up(-2:m+3,-2:n+3),um(-2:m+3,-2:n+3),f(-2:m+3,-2:n+3), &
		                 dfp(-2:m+3,-2:n+3),dfm(-2:m+3,-2:n+3)

    REAL(KIND=DP) ::  uu(1:m,1:n) , ffi(1:m,1:n)
    
    
    INTEGER :: I , J
    
  !   do j=1,n
  !   do i=1,m
	!    u(i,j)=uu(i,j)
	!    fi(i,j)=ffi(i,j)
	!   enddo
	!   enddo
    
    
      !$omp parallel do 	                 
	  do j=1,n
      do i=0,m+1
	    f(i,j)=u(i,j)*fi(i,j)
	  enddo
	  enddo
	  !$omp end parallel do 
	  
	  !$omp parallel do 
	  ! do j=1,n
	  do i=-1,m+1
		do j=1,n
	    dfp(i,j)=up(i+1,j)*fi(i+1,j)-up(i,j)*fi(i,j)
	    dfm(i,j)=um(i+1,j)*fi(i+1,j)-um(i,j)*fi(i,j)
	  enddo
	  enddo
	  !$omp end parallel do 
	  
	!  do j=1,n
	!$omp parallel do 
     do i=2,m-2
     do j=1,n 	
	    fh(i,j)=1.0_DP/12.0_DP*(-f(i-1,j)+7.0_DP*f(i,j)+7.0_DP*f(i+1,j)-f(i+2,j)) &
		       -phyn(dfp(i-2,j),dfp(i-1,j),dfp(i,j),dfp(i+1,j)) &
		       +phyn(dfm(i+2,j),dfm(i+1,j),dfm(i,j),dfm(i-1,j))
	  enddo
	  enddo
	  !$omp end parallel do 
	  
	  end subroutine fflux 
	  
	  
	  FUNCTION PHYN(A,B,C,D)
    INTEGER,PARAMETER  :: DP = SELECTED_REAL_KIND(P=15)      
    REAL(KIND=DP) :: PHYN
    REAL(KIND=DP) :: A,B,C,D,EPS,ALP0,ALP1,ALP2,W0,W2,IS0,IS1,IS2
        EPS=1.0E-09_DP
        IS0=13.0_DP*(A-B)**2+3.0_DP*(A-3.0_DP*B)**2
        IS1=13.0_DP*(B-C)**2+3.0_DP*(B+C)**2
        IS2=13.0_DP*(C-D)**2+3.0_DP*(3.0_DP*C-D)**2
		ALP0=1.0_DP/(EPS+IS0)**2
		ALP1=6.0_DP/(EPS+IS1)**2
		ALP2=3.0_DP/(EPS+IS2)**2
		W0=ALP0/(ALP0+ALP1+ALP2)
		W2=ALP2/(ALP0+ALP1+ALP2)
		PHYN=W0/3.0_DP*(A-2.0_DP*B+C)+(W2-0.5_DP)/6.0_DP*(B-2.0_DP*C+D)

	  END FUNCTION PHYN
	  
	  
	SUBROUTINE BC2D(U,M,N)
	implicit none 
    INTEGER,PARAMETER  :: DP = SELECTED_REAL_KIND(P=15)
	INTEGER,INTENT(IN) :: m , n
    REAL(KIND=DP) :: U(-2:m+3,-2:n+3)
	INTEGER :: I , J
   
    !$omp parallel do 
    DO I=2,M-1
        U(I,1)=2.0_DP*U(I,2)-U(I,3)
        U(I,0)=2.0_DP*U(I,1)-U(I,2)
        U(I,-1)=2.0_DP*U(I,0)-U(I,1)
        U(I,-2)=2.0*U(I,-1)-U(I,0)
        U(I,N)=2.0_DP*U(I,N-1)-U(I,N-2)
        U(I,N+1)=2.0_DP*U(I,N)-U(I,N-1)
        U(I,N+2)=2.0_DP*U(I,N+1)-U(I,N)
        U(I,N+3)=2.0*U(I,N+2)-U(I,N+1)
    ENDDO
	!$omp end parallel do 
	
	!$omp parallel do 
    DO J=2,N-1
        U(1,J)=2.0_DP*U(2,J)-U(3,J)
        U(0,J)=2.0_DP*U(1,J)-U(2,J)
        U(-1,J)=2.0_DP*U(0,J)-U(1,J)
        U(-2,J)=2.0*U(-1,J)-U(0,J)
        U(M,J)=2.0_DP*U(M-1,J)-U(M-2,J)
        U(M+1,J)=2.0_DP*U(M,J)-U(M-1,J)
        U(M+2,J)=2.0_DP*U(M+1,J)-U(M,J)
        U(M+3,J)=2.0*U(M+2,J)-U(M+1,J)
    ENDDO
	!$omp parallel do 

!           DO J=1,N
!     
!           FI1(1,J)=FI1(2,J)    
!           FI1(M,J)=FI1(M-1,J)     
!           fi1(0,J)=fi1(1,J)
!           fi1(-1,J)=fi1(1,J)
!           fi1(-2,J)=fi1(1,J)
!           fi1(M+1,J)=fi1(M,J)
!           fi1(M+2,J)=fi1(M,J)
!           fi1(M+3,J)=fi1(M,J)
!         END DO
!   !    J - DIR
!        DO I=-2,M+3
!     
!           
!           fi1(I,1)=fi1(I,2)  
!           fi1(I,N)=fi1(I,N-1)   
!           fi1(I,0)=fi1(I,1)
!           fi1(I,-1)=fi1(I,1)
!           fi1(I,-2)=fi1(I,1)
!           fi1(I,N+1)=fi1(I,N)
!           fi1(I,N+2)=fi1(I,N)
!           fi1(I,N+3)=fi1(I,N)
!         END DO
    
END SUBROUTINE BC2D


     !  FUNCTION FV_Phi_COMPUTE_UL(CFL , Phi_f , Phi_P , Phi_U , Phi_D)
     !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !  !Function to compute the bounded phi_e by SMART-like method.
     !  !
     !  !
     !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !    INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(p=15)
     !    REAL(KIND=DP),INTENT(IN ) :: CFL , Phi_f , Phi_p , Phi_U , phi_D
     !    REAL(KIND=DP) :: FV_Phi_COMPUTE_UL , EPS
     !    REAL(KIND=DP) :: phi_p_hat , phi_f_hat
     !    REAL(KIND=DP) :: DEL , C , DEL_1 , phi_max , phi_min
     ! 
     !  !  EPS = 1.0E-20_DP
     !    EPS = 1.0E-3_DP
     !    DEL = Phi_D - Phi_U
     !    C   = CFL
     ! 
     ! 
     !   ! to avoid zero value
     !    IF(ABS(DEL) < EPS) DEL = EPS
     ! 
     !    IF(C < EPS) C = EPS
     !    
     !    
     ! 
     !    DEL_1 = 1.0_DP / DEL
     ! 
     !    phi_p_hat = (Phi_p - Phi_U) * DEL_1
     ! 
     ! 
     !   ! write(*,*) " phi_p_hat" ,  phi_p_hat
     !   ! pause
     ! 
     ! 
     !    IF( phi_p_hat < 0.0_DP .OR. phi_p_hat > 1.0_DP) THEN
     !        FV_Phi_COMPUTE_UL = Phi_P
     !        
     !      ! write(*,*) " FV_Phi_COMPUTE_UL = Phi_P" ,  FV_Phi_COMPUTE_UL
     !      ! pause
     !      
     !        
     !    ELSE
     !      FV_Phi_COMPUTE_UL = Phi_f
     ! 
     !     !  write(*,*) " FV_Phi_COMPUTE_UL = Phi_f" ,  FV_Phi_COMPUTE_UL
     !     !  pause
     ! 
     ! 
     !      phi_f_hat = (Phi_f - Phi_U) * DEL_1
     ! 
     !      phi_max = phi_p_hat / C
     !      phi_max = MIN(phi_max , 1.0_DP)   
     !      phi_min = phi_p_hat               !! the lower limit
     ! 
     !      IF(phi_f_hat > phi_max) FV_Phi_COMPUTE_UL = phi_max*DEL + Phi_U
     !      IF(phi_f_hat < phi_min) FV_Phi_COMPUTE_UL = phi_P
     ! 
     !    END IF
     ! 
     !  END FUNCTION FV_Phi_COMPUTE_UL
 
 
   !   FUNCTION FV_Phi_COMPUTE_UL(CFL , Phi_f , Phi_P , Phi_U , Phi_D)
   !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !   !Function to compute the bounded phi_e by SMART-like method.
   !   !
   !   !
   !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !    ! INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(p=15)
   !     REAL CFL , Phi_f , Phi_p , Phi_U , phi_D
   !     REAL FV_Phi_COMPUTE_UL , EPS
   !     REAL phi_p_hat , phi_f_hat
   !     REAL DEL , C , DEL_1 , phi_max , phi_min
   !  
   !   ! EPS = 1.0E-20_DP
   !    ! EPS = 1.0E-1
   !    
   !     eps = 0.05
   !     DEL = Phi_D - Phi_U
   !     C   = CFL
   !  
   !   
   !     ! to avoid zero value
   !      IF(ABS(DEL) < EPS) DEL = EPS
   !   
   !      IF(C < EPS) C = EPS
   !      
   !       
   !   
   !       DEL_1 = 1.0 / DEL
   !   
   !       phi_p_hat = (Phi_p - Phi_U) * DEL_1
   !   
   !   
   !      ! write(*,*) " phi_p_hat" ,  phi_p_hat
   !      ! pause
   !   
   !    
   !       IF( phi_p_hat < 0.0 .OR. phi_p_hat > 1.0) THEN
   !           FV_Phi_COMPUTE_UL = Phi_P
   !           
   !         ! write(*,*) " FV_Phi_COMPUTE_UL = Phi_P" ,  FV_Phi_COMPUTE_UL
   !         ! pause
   !         
   !           
   !       ELSE
   !         FV_Phi_COMPUTE_UL = Phi_f
   !   
   !        !  write(*,*) " FV_Phi_COMPUTE_UL = Phi_f" ,  FV_Phi_COMPUTE_UL
   !        !  pause
   !   
   !   
   !         phi_f_hat = (Phi_f - Phi_U) * DEL_1
   !   
   !         phi_max = phi_p_hat / C
   !         phi_max = MIN(phi_max , 1.0)   
   !         phi_min = phi_p_hat               !! the lower limit
   !   
   !         IF(phi_f_hat > phi_max) FV_Phi_COMPUTE_UL = phi_max*DEL + Phi_U
   !         IF(phi_f_hat < phi_min) FV_Phi_COMPUTE_UL = phi_P
   !   
   !       END IF
   !   
   !     END FUNCTION FV_Phi_COMPUTE_UL
        
        
        
        
      FUNCTION FV_Phi_COMPUTE_UL(CFL , Phi_f , Phi_P , Phi_U , Phi_D)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Function to compute the bounded phi_e by SMART-like method.
      !
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(p=15)
        REAL CFL , Phi_f , Phi_p , Phi_U , phi_D
        REAL FV_Phi_COMPUTE_UL , EPS
        REAL phi_p_hat , phi_f_hat
        REAL DEL , C , DEL_1 , phi_max , phi_min
     
      ! EPS = 1.0E-20_DP
       ! EPS = 1.0E-1
       
        eps = 0.05
        DEL = Phi_D - Phi_U
        C   = CFL
     
      
        ! to avoid zero value
         IF(ABS(DEL) < EPS) DEL = EPS
      
         IF(C < EPS) C = EPS
         
          
      
          DEL_1 = 1.0 / DEL
      
          phi_p_hat = (Phi_p - Phi_U) * DEL_1
      
      
         ! write(*,*) " phi_p_hat" ,  phi_p_hat
         ! pause
      
       
          IF( phi_p_hat < ( 0.0 + 1.0E-10) .OR. phi_p_hat > (1.0 - 1.0E-10)) THEN
              FV_Phi_COMPUTE_UL = Phi_P
              
              
              
            ! write(*,*) " FV_Phi_COMPUTE_UL = Phi_P" ,  FV_Phi_COMPUTE_UL
            ! pause
            
              
          ELSE
            FV_Phi_COMPUTE_UL = Phi_f
      
           !  write(*,*) " FV_Phi_COMPUTE_UL = Phi_f" ,  FV_Phi_COMPUTE_UL
           !  pause
      
      
            phi_f_hat = (Phi_f - Phi_U) * DEL_1
      
            phi_max = phi_p_hat / C
            phi_max = MIN(phi_max , 1.0)   
            phi_min = phi_p_hat               !! the lower limit
      
            IF(phi_f_hat > phi_max) FV_Phi_COMPUTE_UL = phi_max*DEL + Phi_U
            IF(phi_f_hat < phi_min) FV_Phi_COMPUTE_UL = phi_P
      
          END IF
      
        END FUNCTION FV_Phi_COMPUTE_UL
      
      
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 