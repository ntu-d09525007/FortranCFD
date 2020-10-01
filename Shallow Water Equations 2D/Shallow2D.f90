program main 
!$ use omp_lib
implicit none 
integer :: i,j,it,kl,ioutput
integer :: m,m1,n,n1,nt,threads,cnt

real(8) :: chanl,chanwid,xdam
real(8) :: damcrest,damheight,depthcrest,depthBreach,widthBreach,slopeup,slopedwn
real(8) :: dt,dt0,dtij,aMann,diased,poro
real(8) :: reposedry,reposewet
real(8) :: epsilon,hmin,nwsl,time,time2stop,time2plot
real(8) :: xuptoe,xupcorner,xdwntoe,xdwncorner,xtoeLength
real(8) :: ixtoeUp,ixtoeDwn,idammid,idamup,idamdwn, jbreach
real(8) :: spec,rhos,rhow,rhob,amiu,amiud,wset 

real(8) :: cqb,uv,ustar2,ustarcr2,alamda,ustar2x,ustar2y,exm,wset1,qsstar,Codune,dune,qbstar,qtstar
real(8) :: slopx,slopy,aCoqtlimit,qtstarlimit
real(8) :: zsup,zsdown,dzdx,dzdy,rho
real(8) :: qsum,qsum1,qsum2,qsum3
real(8) :: WidthBreach742,WidthBreach744,WidthBreach1,WidthBreach2
real(8) :: breach_ystart, breach_yend
real(8) :: cputime

real(8),dimension(:),allocatable :: dxc,dyc
real(8),dimension(:,:),allocatable :: x,y,an,zb
real(8),dimension(:,:),allocatable :: bednet,dzb,bedslopex,bedslopey
real(8),dimension(:,:),allocatable :: qx1,qy1,h1,zs1,u1,v1,uv1,c1
real(8),dimension(:,:),allocatable :: qx2,qy2,h2,zs2,u2,v2,c2,idix,kdry,zbmin
real(8),dimension(:,:),allocatable :: fmass,fxmom,fymom,fsed
real(8),dimension(:,:),allocatable :: gmass,gxmom,gymom,gsed
real(8),dimension(:,:,:),allocatable :: alenii,xxrep,yyrep,volrep,idixrep
real(8),dimension(:,:),allocatable :: dzdxmin,dzdymin


	m1=1000
	m=m1+1
	   
	n1=400
	n=n1+1
	   
	time2stop=60.0d0*25.0d0
	time2plot=10.0d0
	
	allocate( dxc(m),dyc(n),x(m,n),y(m,n),an(m,n),zb(m,n) )
	allocate( bednet(m,n),dzb(m,n),bedslopex(m,n),bedslopey(m,n) )
	allocate( qx1(m,n),qy1(m,n),h1(m,n),zs1(m,n),u1(m,n),v1(m,n),uv1(m,n),c1(m,n) )
	allocate( qx2(m,n),qy2(m,n),h2(m,n),zs2(m,n),u2(m,n),v2(m,n),c2(m,n),idix(m,n),kdry(m,n),zbmin(m,n) )
	allocate( fmass(m1,n1),fxmom(m1,n1),fymom(m1,n1),fsed(m1,n1) )
	allocate( gmass(m1,n1),gxmom(m1,n1),gymom(m1,n1),gsed(m1,n1) )
	allocate( alenii(m,n,8),xxrep(m,n,8),yyrep(m,n,8),volrep(m,n,8),idixrep(m,n,8) )
	allocate( dzdxmin(m,n), dzdymin(m,n) )

open(unit=111,file='cpu time.txt')

!$ do cnt = 5, 5

	threads=1
	!$ threads = 64!2**(cnt)


	!$ call omp_set_dynamic(0)
	!$ call omp_set_num_threads(threads)
       
 
	OPEN (10,file="2DbkInput.txt",status='old')
	OPEN (99,FILE="Breach Width.plt")
	!OPEN (20,file="2DbkOutput.dat")

    epsilon=0.00000001
    hmin=0.0001
    nwsl=1
    time=0.0

	read (10,*)
	read (10,*) chanl,chanwid,xdam
	read (10,*)
	read (10,*) damcrest,damheight,depthcrest,depthBreach,widthBreach,slopeup,slopedwn
	read (10,*)
	read (10,*) dt0,aMann,diased,poro
	read (10,*)
	read (10,*) reposedry,reposewet

    do i=1,m
        dxc(i)=chanl/real(m-1)
    end do
           
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
		an(i,j)=aMann
    end do
    end do
    !$omp end parallel do

    xuptoe=xdam-damcrest/2.0-damheight/slopeup
    xupcorner=xdam-damcrest/2.0
    xdwntoe=xdam+damcrest/2.0+damheight/slopedwn
    xdwncorner=xdam+damcrest/2.0
   
    xtoeLength=damcrest+damheight/slopeup+damheight/slopedwn
    ixtoeUp=xuptoe/dxc(1)
    ixtoeDwn=1+(xuptoe+xtoeLength)/dxc(1)
    idammid=1+xdam/dxc(1)
    idamup=(xdam-damcrest/2.0)/dxc(1)
    idamdwn=1+(xdam+damcrest/2.0)/dxc(1)

    jbreach=WidthBreach/2.0/dyc(1)
       
    xdam=xdam-xuptoe
    xupcorner=xupcorner-xuptoe
    xdwntoe=xdwntoe-xuptoe
    xdwncorner=xdwncorner-xuptoe


    !$omp parallel do num_threads(threads)
    do j=1,n
    do i=1,m
        x(i,j)=x(i,j)-xuptoe
    end do
    end do
    !$omp end parallel do
       
    xuptoe=0.0

    !$omp parallel do num_threads(threads)
    do j=1,n
    do i=1,m
   
        if(x(i,j) .lt. xuptoe) then
           zb(i,j)=0.0
        else if(x(i,j) .ge. xuptoe .and. x(i,j) .le. xupcorner) then
           zb(i,j)=(x(i,j)-xuptoe)*slopeup
        else if(x(i,j) .ge. xupcorner .and. x(i,j) .le. xdwncorner) then
           zb(i,j)=damheight
        else if(x(i,j) .ge. xdwncorner .and. x(i,j) .le. xdwntoe) then
           zb(i,j)=(xdwntoe-x(i,j))*slopedwn
        else if(x(i,j) .gt. xdwntoe) then
           zb(i,j)=0.0
        endif
         
        zbmin(i,j)=0.0      !nonerodible
         
    end do
    end do
    !$omp end parallel do

    !$omp parallel do num_threads(threads)
    do j=1,n
    do i=1,m
        if(j > n/2-jbreach .and. j < n/2+jbreach) then
           zb(i,j)=min(zb(i,j),damheight-depthBreach)
        end if
    end do
    end do
    !$omp end parallel do

    !$omp parallel do num_threads(threads)
    do j=2,n-1
    do i=2,m-1
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
      
    !$omp parallel do num_threads(threads)
    do j=2,n-1
    do i=2,m-1
       
        do kl=1,8
            alenii(i,j,kl)=sqrt( (x(i,j)-xxrep(i,j,kl))**2  &
                                +(y(i,j)-yyrep(i,j,kl))**2)         
        end do

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

    it=0


    call Output(m,n,time,x,y,zb,zs2,h2,u2,v2,qx2,qy2,c2,it,ioutput)
    ioutput=1
    call bedrepos0(threads,m,n,x,y,zb,idix,dxc,dyc,it,alenii,xxrep,yyrep,volrep,idixrep)

    !$omp parallel do num_threads(threads)
    do j=1,n
    do i=1,m
        qx1(i,j)=0.0
        qy1(i,j)=0.0
         
        if(x(i,j) .gt. xdam) then
           zs1(i,j)=max(epsilon,hmin+zb(i,j))
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

    spec=2.65
    rhos=2650.0
    rhow=1000.0

    rhob=rhos*(1.0-poro)+rhow*poro 

    amiu=0.000001
    amiud=amiu*1000.0/(diased*1000.0)
    wset=sqrt((13.95*amiud)**2+1.09*(spec-1.0)*9.81*diased)-13.95*amiud


!$ cputime = -omp_get_wtime()
!....................................................Start Time Loop......................................
    do
      
		if(time>time2stop)exit
      
        cqb=sqrt(1.65*9.81*diased**3)
        ustarcr2=0.03*1.65*9.81*diased
         
        !$omp parallel do private(uv,ustar2,alamda,ustar2x,ustar2y,exm,wset1,qsstar,Codune,dune,qbstar,qtstar,&
        !$omp &slopx,slopy,aCoqtlimit,qtstarlimit), num_threads(threads)
		do j=1,n
		do i=1,m
			if(kdry(i,j) .eq. 1) then
				uv=uv1(i,j)
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
			
				alamda=1.0
				ustar2x=ustar2*u1(i,j)/(uv+epsilon)-alamda*ustarcr2*bedslopex(i,j)/0.65
				ustar2y=ustar2*v1(i,j)/(uv+epsilon)-alamda*ustarcr2*bedslopey(i,j)/0.65
				ustar2=sqrt(ustar2x**2+ustar2y**2)

                !wset1=wset
                wset1=wset*(1.0-min(0.65,c1(i,j)))**4.0

				!exm=4.7-1.324*log10(wset*diased/amiu/5.0)
				!exm=min( 4.7, max(2.3,exm) )
				!wset1=wset*(1.0-min(0.65,c1(i,j)))**exm
				qsstar=0.0000262*cqb*(max(0.0,ustar2/ustarcr2-1.0)*uv/wset1)**1.74
				Codune=20.0
				dune=(diased**0.166667/Codune/an(i,j))**1.5
				qbstar=0.0053*cqb*(max(0.0,dune*ustar2/ustarcr2-1.0))**2.2
				qtstar=qbstar+qsstar

				if((i .gt. 1 .and. i .lt. m).and.(j .gt. 1 .and. j .lt. n)) then
					slopx=(zb(i+1,j)-zb(i-1,j))/dxc(i)/2.0
					slopy=(zb(i,j+1)-zb(i,j-1))/dyc(j)/2.0
					qtstar=qtstar*sqrt(1.0+slopx**2+slopy**2)
				endif

				aCoqtlimit=0.5
				qtstarlimit=aCoqtlimit*uv*h1(i,j)*(1.0-poro)
            
				if(qtstar .gt. qtstarlimit) qtstar=qtstarlimit
      
				bednet(i,j)=(c1(i,j)*uv*h1(i,j)-qtstar)/0.05 ! adaption length

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
        !$omp parallel do reduction(min:dt), private(dtij), num_threads(threads)
		do j=1,n
		do i=1,m
			if(kdry(i,j) .eq. 1) then
				dtij=0.1*h1(i,j)*(1.0-poro)/(abs(bednet(i,j))+0.01*epsilon)
				if(dtij .lt. dt) dt=dtij
				dtij=0.5/( (abs(u1(i,j))+sqrt(9.81*h1(i,j)))/dxc(i)+(abs(v1(i,j))+sqrt(9.81*h1(i,j)))/dyc(j) )   
				if(dtij.lt.dt) dt=dtij
			endif
        end do
        end do
        !$omp end parallel do

        time=time+dt
     
        !$omp parallel do num_threads(threads)
		do j=1,n
		do i=1,m
			dzb(i,j)=dt*bednet(i,j)/(1.0-poro)
            zb(i,j)=zb(i,j)+dzb(i,j)
        end do
        end do
        !$omp end parallel do

        call bedrepos(threads,m,n,x,y,zb,zs2,idix,dxc,dyc,it,alenii,xxrep,yyrep,volrep,idixrep,hmin,reposedry,reposewet)
   
        call cflux_HLL(threads,m,n,hmin,h1,zb,zs1,fmass,gmass,fxmom,gxmom,fymom,gymom,fsed,gsed,epsilon,qx1,qy1,u1,v1,c1,dxc,dyc,dt)
   
        !$omp parallel do num_threads(threads)
		do j=2,n-1
		do i=2,m-1
            h2(i,j)=h1(i,j)-dt/dxc(i)*(fmass(i,j)-fmass(i-1,j))   &
                           -dt/dyc(j)*(gmass(i,j)-gmass(i,j-1))   &
                           -0*dt*bednet(i,j)/(1.0-poro)					
            h2(i,j)=max(hmin,h2(i,j))
            zs2(i,j)=zb(i,j)+h2(i,j)
        enddo
        enddo
        !$omp end parallel do
         
        !$omp parallel do num_threads(threads)
        do j=1,n
            zs2(1,j)=zs2(2,j)
            !zs2(1,j)=1.5*zs2(2,j)-0.5*zs2(3,j)
            h2(1,j) =max(hmin, zs2(1,j)-zb(1,j))
            zs2(m,j)=zs2(m-1,j)
            !zs2(m,j)=1.5*zs2(m-1,j)-0.5*zs2(m-2,j)
            h2(m,j) =max(hmin, zs2(m,j)-zb(m,j))
        end do
        !$omp end parallel do
         
        !$omp parallel do num_threads(threads)
        do i=1,m
            zs2(i,1)=zs2(i,2)
            h2(i,1) =max(hmin, zs2(i,1)-zb(i,1))
            zs2(i,n)=zs2(i,n-1)
            h2(i,n) =max(hmin, zs2(i,n)-zb(i,n))
        end do
        !$omp end parallel do

        if(nwsl.eq.2) call dzdxyminmod(threads,m,n,dzdxmin,dzdymin,x,y,zs2)

!..............X-Momentum
        !$omp parallel do num_threads(threads)
        do j=1,n
           qx2(1,j)=0.07/chanwid
         end do
        !$omp end parallel do

        !$omp parallel do private(zsup,zsdown,dzdx,rho), num_threads(threads)
		do j=2,n-1
		do i=2,m-1
         
			if(h2(i,j).le.hmin) then
				qx2(i,j)=0.0
			else
           
				if(nwsl.eq.1) then
					zsup  =zs2(i-1,j)
					zsdown=zs2(i+1,j)
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
					zsup  =0.5*( zs2(i-1,j)+0.5*dxc(i-1)*dzdxmin(i-1,j)+zs2(i  ,j)-0.5*dxc(i  )*dzdxmin(i  ,j) ) 
					zsdown=0.5*( zs2(i  ,j)+0.5*dxc(i  )*dzdxmin(i  ,j)+zs2(i+1,j)-0.5*dxc(i+1)*dzdxmin(i+1,j) )         
					if(i.eq.2  ) zsup  =zs2(1,j)
					if(i.eq.m-1) zsdown=zs2(m,j)
					dzdx=zsdown-zsup
					if(h2(i-1,j).le.hmin.and.zs2(i-1,j).gt.zs2(i,j))dzdx=0.0           
					if(h2(i+1,j).le.hmin.and.zs2(i+1,j).gt.zs2(i,j))dzdx=0.0     
				endif

				rho=rhos*c1(i,j)+rhow*(1.0-c1(i,j)) 
				qx2(i,j)=qx1(i,j)-dt/dxc(i)*(fxmom(i,j)-fxmom(i-1,j))  &
								-dt/dyc(j)*(gxmom(i,j)-gxmom(i,j-1))  &
								-9.81*h2(i,j)*dt/dxc(i)*dzdx           &
								-dt*0.5*9.81*h1(i,j)*h1(i,j)*(rhos-rhow)*(c1(i+1,j)-c1(i-1,j))/2.0/dxc(i)/rho     &
								+0*dt*u1(i,j)*(rhob-rho)/rho*bednet(i,j)/(1.0-poro)			
				qx2(i,j)=qx2(i,j)/(1.0+9.81*dt*an(i,j)**2*uv1(i,j)/h1(i,j)**1.3333)
			endif
        end do
        end do
        !$omp end parallel do
         
        !$omp parallel do num_threads(threads)
        do j=1,n
			qx2(m,j)=qx2(m-1,j)
			!qx2(m,j)=SQRT(9.81*h2(m,j)**3.0)
        end do
        !$omp end parallel do
         
        !$omp parallel do num_threads(threads)
        do i=2,m
            !qx2(i,1)=0.0
            !qx2(i,n)=0.0
            qx2(i,1)=qx2(i,2)
            qx2(i,n)=qx2(i,n-1)
        end do
        !$omp end parallel do

!..............Y-Momentum
        !$omp parallel do num_threads(threads)
        do j=1,n
            qy2(1,j)=0.0
        end do
        !$omp end parallel do

        !$omp parallel do private(zsup,zsdown,dzdy,rho), num_threads(threads)
		do j=2,n-1
		do i=2,m-1
         
			if(h2(i,j) .le. hmin) then
				qy2(i,j)=0.0
			else
				if(nwsl.eq.1) then
					zsup  =zs2(i,j-1)
					zsdown=zs2(i,j+1)
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
					zsup  =0.5*( zs2(i,j-1)+0.5*dyc(j-1)*dzdymin(i,j-1)+zs2(i,j  )-0.5*dyc(j  )*dzdymin(i,j  ) ) 
					zsdown=0.5*( zs2(i,j  )+0.5*dyc(j  )*dzdymin(i,j  )+zs2(i,j+1)-0.5*dyc(j+1)*dzdymin(i,j+1) )          
					if(j .eq. 2  ) zsup  =zs2(i,1)
					if(j .eq. n-1) zsdown=zs2(i,n)
					dzdy=zsdown-zsup           
					if(h2(i,j-1).le.hmin.and.zs2(i,j-1).gt.zs2(i,j))dzdy=0.0           
					if(h2(i,j+1).le.hmin.and.zs2(i,j+1).gt.zs2(i,j))dzdy=0.0         
				endif
			
				rho=rhos*c1(i,j)+rhow*(1.0-c1(i,j)) 
				qy2(i,j)=qy1(i,j)-dt/dxc(i)*(fymom(i,j)-fymom(i-1,j))   &
								-dt/dyc(j)*(gymom(i,j)-gymom(i,j-1))   &
								-9.81*h2(i,j)*dt/dyc(j)*dzdy           &
								-dt*0.5*9.81*h1(i,j)*h1(i,j)*(rhos-rhow)*(c1(i,j+1)-c1(i,j-1))/2.0/dyc(j)/rho  &
								+0*dt*v1(i,j)*(rhob-rho)/rho*bednet(i,j)/(1.0-poro)			
				qy2(i,j)=qy2(i,j) /(1.0+9.81*dt*an(i,j)**2*uv1(i,j)/h1(i,j)**1.3333)
            endif
           
        end do
        end do
        !$omp end parallel do
        
        !$omp parallel do num_threads(threads)
        do j=1,n
			qy2(m,j)=qy2(m-1,j)
        end do
        !$omp end parallel do
         
        !$omp parallel do num_threads(threads)
        do i=1,m
            qy2(i,1)=0.0
            qy2(i,n)=0.0
        end do
        !$omp end parallel do

!.........Sediment
         
        !$omp parallel do num_threads(threads)
        do j=1,n
			c2(1,j)=0.0
        end do
        !$omp end parallel do
         
        !$omp parallel do num_threads(threads)
		do j=2,n-1
		do i=2,m-1
            c2(i,j)=(c1(i,j)*h1(i,j)-dt/dxc(i)*(fsed(i,j)-fsed(i-1,j))   &
                                    -dt/dyc(j)*(gsed(i,j)-gsed(i,j-1))   &
                                    -dt*bednet(i,j) )/h2(i,j)
            c2(i,j)=max(0.0,c2(i,j))
            c2(i,j)=min(1-poro,c2(i,j))
        end do
        end do
        !$omp end parallel do
         
        !$omp parallel do num_threads(threads)
        do j=1,n
            c2(m,j)=c2(m-1,j)
        end do
        !$omp end parallel do
         
        !$omp parallel do num_threads(threads)
        do i=1,m
            c2(i,1)=c2(i,2)
            c2(i,n)=c2(i,n-1)
        end do
        !$omp end parallel do

        !$omp parallel do num_threads(threads)
		do j=1,n
		do i=1,m
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
		
		qsum1=0.0		
		qsum2=0.0		
		qsum3=0.0		
		WidthBreach742=0.
		WidthBreach744=0.
		WidthBreach2=0.
		WidthBreach1=0.
			
		i=idammid
		qsum=0.0
		do j=2,n-1
			qsum=qsum+qx2(i,j)*dyc(j)
			IF (zb(i,j) <= damheight .and. abs(zb(i,j-1)-damheight) < 1.E-10 .and. zb(i,j+1) < zb(i,j) ) WidthBreach1=y(i,j)
			IF (zb(i,j) <= damheight .and. abs(zb(i,j+1)-damheight) < 1.E-10 .and. zb(i,j)>zb(i,j-1) ) WidthBreach2=y(i,j)
			WidthBreach742=WidthBreach2-WidthBreach1
		enddo
	 
        if( abs(time-time2plot*Ioutput)<dt ) then
			call Output(m,n,time,x,y,zb,zs2,h2,u2,v2,qx2,qy2,c2,it,Ioutput)
			Ioutput=Ioutput+1
		endif
!................................................................................

		
		write(99,*)time,WidthBreach742,qsum
		write(*,'(ES11.4," ",F10.2," ",F8.4," ",F8.4)') dt, time,WidthBreach742,qsum
		
		
		
	enddo

	
	!$ cputime = cputime + omp_get_wtime()
	
	close(10)
	close(20)
	
	write(111,*)threads,cputime

!$ end do 

contains 

!#################################################################################################################

subroutine cflux_HLL(threads,m,n,hmin,h1,zb,zs1,fmass,gmass,fxmom,gxmom,fymom,gymom,fsed,gsed,epsilon,qx1,qy1,u1,v1,c1,dxc,dyc,dt)
implicit none 
integer :: i,j,m,n,threads
real(8) :: dxc(m),dyc(n)
real(8) :: qx1(m,n),qy1(m,n),h1(m,n),zs1(m,n)
real(8) :: u1(m,n),v1(m,n),c1(m,n),zb(m,n)
real(8) :: fmass(m-1,n-1),fxmom(m-1,n-1),fymom(m-1,n-1)
real(8) :: gmass(m-1,n-1),gxmom(m-1,n-1),gymom(m-1,n-1)
real(8) :: fsed(m-1,n-1),gsed(m-1,n-1)
real(8) :: hmin,dt,epsilon
real(8) :: uxl,uxr,hxl,hxr,axl,axr,hxstar,uxstar,aqxr,aqxl,sxl,sxr
real(8) :: uyl,uyr,hyl,hyr,ayl,ayr,hystar,uystar,aqyl,aqyr,syl,syr
real(8) :: fmassl,fmassr,fxmoml,fxmomr,fymoml,fymomr
real(8) :: gmassl,gmassr,gxmoml,gxmomr,gymoml,gymomr
real(8),dimension(-2:m+3,-2:n+3) :: f,fp,fm,hp,hm

    !$omp parallel do num_threads(threads)
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
	
	!call weno_flux(h1,f,fp,fm,hp,hm,m,n)
       
    !$omp parallel do private(uxl,uxr,hxl,hxr,axl,axr,hxstar,uxstar,aqxr,aqxl,sxl,sxr,&
    !$omp &fmassl,fmassr,fxmoml,fxmomr,fymoml,fymomr), num_threads(threads)
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

        if(sxl.gt.0.0) then
            fmass(i,j)=fmassl             
            fxmom(i,j)=fxmoml
            fymom(i,j)=fymoml
        else if(sxl.le.0.0.and.sxr.ge.0.0) then
            fmass(i,j)=( sxr*fmassl-sxl*fmassr+sxl*sxr*(zs1(i+1,j)-zs1(i,j)) )/(sxr-sxl)             
            fxmom(i,j)=( sxr*fxmoml-sxl*fxmomr+sxl*sxr*(qx1(i+1,j)-qx1(i,j)) )/(sxr-sxl)
            fymom(i,j)=( sxr*fymoml-sxl*fymomr+sxl*sxr*(qy1(i+1,j)-qy1(i,j)) )/(sxr-sxl)
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

		!$omp parallel do num_threads(threads)
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
		!$omp &gmassl,gmassr,gxmoml,gxmomr,gymoml,gymomr), num_threads(threads)
		do j=2,n-1
		do i=2,m-1
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

            if(syl.gt.0.0) then
				gmass(i,j)=gmassl             
				gxmom(i,j)=gxmoml
				gymom(i,j)=gymoml
            elseif(syl.le.0.0.and.syr.ge.0.0) then
				gmass(i,j)=( syr*gmassl - syl*gmassr + syl*syr*(zs1(i,j+1)-zs1(i,j)) )/(syr-syl)             
				gxmom(i,j)=( syr*gxmoml - syl*gxmomr + syl*syr*(qx1(i,j+1)-qx1(i,j)) )/(syr-syl)
				gymom(i,j)=( syr*gymoml - syl*gymomr + syl*syr*(qy1(i,j+1)-qy1(i,j)) )/(syr-syl)
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

end subroutine

subroutine weno_flux(ff,f,fp,fm,hp,hm,m,n)
implicit none
integer :: i,j,m,n
real(8),dimension(-2:m+3,-2:n+3) :: f,fp,fm,hp,hm
real(8),dimension(m,n) :: ff
real(8) :: a1,a2,a3,b1,b2,b3,w1,w2,w3,eps

	!$omp parallel do
	do j = 1, n
	do i = 1, m
		f(i,j) = ff(i,j)
	end do
	end do 
	!$omp end parallel do 
	
	!$omp parallel do
	do i = 1, m
	do j = 1, 3
		f(i,1-j) = f(i,1)
		f(i,n+j) = f(i,n)
	end do
	end do
	!$omp end parallel do 

	!$omp parallel do
	do i = 1, 3
	do j = 1, n
		f(1-i,j) = f(1,j)
		f(n+i,j) = f(n,j)
	end do
	end do
	!$omp end parallel do
	
	do j = 1, n
		!call weno_js(f(:,j),fp(:,j),fm(:,j),m)
		call ocrweno(f(:,j),fp(:,j),fm(:,j),m)
	end do

	do i = 1, m
		!call weno_js(f(i,:),hp(i,:),hm(i,:),n)
		call ocrweno(f(i,:),hp(i,:),hm(i,:),n)
	end do	
	

end subroutine
	
subroutine weno_js(f,fp,fm,Nx)
implicit none
integer :: Nx, i
real(8),dimension(-2:Nx+3) :: f, fp, fm
real(8) :: a1,a2,a3,b1,b2,b3,w1,w2,w3,eps

eps = 1.0d-12

!$OMP PARALLEL DO PRIVATE(B1,B2,B3,A1,A2,A3,W1,W2,W3)
do i = 0, Nx
	
	b1 = 13.0d0*(f(i-2)-2.0*f(i-1)+f(i))**2 + 3.0d0*(f(i-2)-4*f(i-1)+3.0*f(i))**2
	b2 = 13.0d0*(f(i-1)-2.0*f(i)+f(i+1))**2 + 3.0d0*(f(i-1)-f(i+1))**2
	b3 = 13.0d0*(f(i)-2.0*f(i+1)+f(i+2))**2 + 3.0d0*(3.0*f(i)-4.0*f(i+1)+f(i+2))**2
	
	a1 = 1.0d0/(EPS+b1)**2
	a2 = 6.0d0/(EPS+b2)**2
	a3 = 3.0d0/(EPS+b3)**2
	
	w1 = a1/(a1+a2+a3)
	w2 = a2/(a1+a2+a3)
	w3 = a3/(a1+a2+a3)
	
	fm(i) = w1/3.0d0*f(i-2) - (7.0d0*w1+w2)/6.0d0*f(i-1) + (11.0d0*w1+5.0d0*w2+2.0d0*w3)/6.0d0*f(i) &
			+ (2.0d0*w2+5.0d0*w3)/6.0d0*f(i+1) - w3/6.0d0*f(i+2)
			
end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(B1,B2,B3,A1,A2,A3,W1,W2,W3)	
do i = 0, Nx
	
	b3 = 13.0*(f(i-1)-2.0*f(i)  +f(i+1))**2 + 3.0*(f(i-1)-4.0*f(i)+3.0*f(i+1))**2
	b2 = 13.0*(f(i)  -2.0*f(i+1)+f(i+2))**2 + 3.0*(f(i)-f(i+2))**2
	b1 = 13.0*(f(i+1)-2.0*f(i+2)+f(i+3))**2 + 3.0*(3.0*f(i+1)-4.0*f(i+2)+f(i+3))**2
	
  a1 = 1.0/(EPS+b1)**2
  a2 = 6.0/(EPS+b2)**2
  a3 = 3.0/(EPS+b3)**2	
	
  w1 = a1 / (a1+a2+a3)
  w2 = a2 / (a1+a2+a3)
  w3 = a3 / (a1+a2+a3)
 	
 	fp(i) = w3*(-f(i-1)+5.0*f(i)+2.0*f(i+1))/6.0 &
 				 +w2*(2.0*f(i)+5.0*f(i+1)-f(i+2))/6.0 &
 				 +w1*(11.0*f(i+1)-7.0*f(i+2)+2.0*f(i+3))/6.0	
	
	
end do
!$OMP END PARALLEL DO
end subroutine

subroutine weno_bc(f,fp,fm,Nx)
implicit none
integer :: Nx, i
real(8),dimension(-2:Nx+3) :: f, fp, fm
real(8) :: a1,a2,a3,b1,b2,b3,w1,w2,w3,eps

eps = 1.0d-12

!$OMP PARALLEL DO PRIVATE(B1,B2,B3,A1,A2,A3,W1,W2,W3)
do i = 0, Nx, Nx
	
	b1 = 13.0*(f(i-2)-2.0*f(i-1)+f(i))**2 + 3.0*(f(i-2)-4*f(i-1)+3.0*f(i))**2
	b2 = 13.0*(f(i-1)-2.0*f(i)+f(i+1))**2 + 3.0*(f(i-1)-f(i+1))**2
	b3 = 13.0*(f(i)-2.0*f(i+1)+f(i+2))**2 + 3.0*(3.0*f(i)-4.0*f(i+1)+f(i+2))**2
	
	a1 = 1.0*(1.0+abs(b3-b1)/(EPS+b1))
	a2 = 6.0*(1.0+abs(b3-b1)/(EPS+b2))
	a3 = 3.0*(1.0+abs(b3-b1)/(EPS+b3))
	
	w1 = a1/(a1+a2+a3)
	w2 = a2/(a1+a2+a3)
	w3 = a3/(a1+a2+a3)
	
	fm(i) = w1/3.0*f(i-2) - (7.0*w1+w2)/6.0*f(i-1) + (11.0*w1+5.0*w2+2.0*w3)/6.0*f(i) &
			+ (2.0*w2+5.0*w3)/6.0*f(i+1) - w3/6.0*f(i+2)
			
end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(B1,B2,B3,A1,A2,A3,W1,W2,W3)	
do i = 0, Nx, Nx
	
	b3 = 13.0*(f(i-1)-2.0*f(i)  +f(i+1))**2 + 3.0*(f(i-1)-4.0*f(i)+3.0*f(i+1))**2
	b2 = 13.0*(f(i)  -2.0*f(i+1)+f(i+2))**2 + 3.0*(f(i)-f(i+2))**2
	b1 = 13.0*(f(i+1)-2.0*f(i+2)+f(i+3))**2 + 3.0*(3.0*f(i+1)-4.0*f(i+2)+f(i+3))**2
	
  a1 = 1.0 * ( 1.0 + abs(b1-b3)/(EPS+b1) )
  a2 = 6.0 * ( 1.0 + abs(b1-b3)/(EPS+b2) )
  a3 = 3.0 * ( 1.0 + abs(b1-b3)/(EPS+b3) )	
	
  w1 = a1 / ( a1+a2+a3)
  w2 = a2 / ( a1+a2+a3)
  w3 = a3 / ( a1+a2+a3)
 	
 	fp(i) = w3*(-f(i-1)+5.0*f(i)+2.0*f(i+1))/6.0 &
 				 +w2*(2.0*f(i)+5.0*f(i+1)-f(i+2))/6.0 &
 				 +w1*(11.0*f(i+1)-7.0*f(i+2)+2.0*f(i+3))/6.0	
	
	
end do
!$OMP END PARALLEL DO
end subroutine

subroutine ocrweno(f,fp,fm,Nx)
implicit none
integer :: Nx, i
real(8),dimension(-2:Nx+3) :: f, fp, fm
real(8) :: a1,a2,a3,b1,b2,b3,w1,w2,w3,c1,c2,c3,eps
real(8),dimension(1:Nx-1) :: A,B,C,S

eps = 1.0d-12

	!c1 = 2.0
	!c2 = 5.0
	!c3 = 3.0
	
	c1 = 0.2089141306d0
	c2 = 0.4999999998d0
	c3 = 0.2910858692d0

 call weno_bc(f,fp,fm,Nx)

!$OMP PARALLEL DO PRIVATE(B1,B2,B3,A1,A2,A3,W1,W2,W3)
do i = 1, Nx-1
	
	b1 = 13.0*(f(i-2)-2.0*f(i-1)+f(i))**2   + 3.0*(    f(i-2)-4*f(i-1)+3.0*f(i))**2
	b2 = 13.0*(f(i-1)-2.0*f(i)  +f(i+1))**2 + 3.0*(    f(i-1)  -f(i+1))**2
	b3 = 13.0*(f(i)  -2.0*f(i+1)+f(i+2))**2 + 3.0*(3.0*f(i)-4.0*f(i+1)+f(i+2))**2
	
	a1 = c1*(  1.0 +  abs(b3-b1) / (EPS+b1)  )
	a2 = c2*(  1.0 +  abs(b3-b1) / (EPS+b2)  )
	a3 = c3*(  1.0 +  abs(b3-b1) / (EPS+b3)  )
	
	w1 = a1/(a1+a2+a3)
	w2 = a2/(a1+a2+a3)
	w3 = a3/(a1+a2+a3)	
	
	A(i) = (2.0*w1+w2)/3.0
	B(i) = (w1+2.0*(w2+w3))/3.0
	C(i) =  w3/3.0
	S(i) = w1/6.0*f(i-1) + (5.0*(w1+w2)+w3)/6.0*f(i) + (w2+5.0*w3)/6.0*f(i+1)
	
end do
!$OMP END PARALLEL DO

  S(1) = S(1) - A(1)*fm(0)
  S(Nx-1) = S(Nx-1) - C(Nx-1)*fm(Nx)
  
  
  call solve_tridiagonal(A,B,C,S,fm(1:Nx-1),1,Nx-1)
  
!$OMP PARALLEL DO PRIVATE(B1,B2,B3,A1,A2,A3,W1,W2,W3)  
do i = 1, Nx-1
	
	b3 = 13.0*(f(i-1)-2.0*f(i)  +f(i+1))**2 + 3.0*(f(i-1)-4.0*f(i)+3.0*f(i+1))**2
	b2 = 13.0*(f(i)  -2.0*f(i+1)+f(i+2))**2 + 3.0*(f(i)-f(i+2))**2
	b1 = 13.0*(f(i+1)-2.0*f(i+2)+f(i+3))**2 + 3.0*(3.0*f(i+1)-4.0*f(i+2)+f(i+3))**2
	
	a1 = c1*(1.0+abs(b3-b1)/(EPS+b1))
	a2 = c2*(1.0+abs(b3-b1)/(EPS+b2))
	a3 = c3*(1.0+abs(b3-b1)/(EPS+b3))
	
	w1 = a1/(a1+a2+a3)
	w2 = a2/(a1+a2+a3)
	w3 = a3/(a1+a2+a3)	
	
	A(i) = (w3)/3.0
	B(i) = (w1+2.0*(w2+w3))/3.0
	C(i) = (w2+2.0*w1)/3.0
	S(i) = (5.0*w3+w2)/6.0*f(i) + (w3+5.0*(w2+w1))/6.0*f(i+1) + w1/6.0*f(i+2)
	
end do
!$OMP END PARALLEL DO

  S(1) = S(1) - A(1)*fp(0)
  S(Nx-1) = S(Nx-1) - C(Nx-1)*fp(Nx)
  
  call solve_tridiagonal(A,B,C,S,fp(1:Nx-1),1,Nx-1)
  
  
end subroutine

 subroutine solve_tridiagonal(A,B,C,S,X,M,N)
!cccccccccccccccccccccccccccccccccc
!
! A, first  coefficient matrix
! B, second coefficient matrix
! C, third  coefficient matrix
! S, sorce matrix
! X, solution
! [M,N], index domain
!
!cccccccccccccccccccccccccccccccccc
implicit none
INTEGER,PARAMETER :: DP = SELECTED_REAL_KIND(p=15)
integer :: M, N, i
real(kind=dp),dimension(m:n) :: A,B,C,S,X

 C(M) = C(M)/ B(M)
 S(M) = S(M)/ B(M)
 
 do i = M+1, N, 1
	 
	 C(i) = C(i) / ( B(I)-A(I)*C(I-1) )
	 S(I) = (S(I) - A(i)*S(I-1))/(B(i)-A(I)*C(i-1))

 end do
 
 X(N) = S(N)
 
 do i = N-1, M, -1
	 
	 X(i) = S(i) - C(i)*X(i+1)	 
	 
 end do


end subroutine


!     =================================================================
!     This subroutine handles the steep slope case
!           Made by Dr. Weiming Wu, ..............
subroutine bedrepos(threads,m,n,x,y,zb,zs,idix,dxc,dyc,it,alenii,   &
                    xxrep,yyrep,volrep,idixrep,hmin,reposedry,reposewet)
!    repose angle is considered here to check the bed profile
!    if the bed slope is larger than the repose angle, the bed will collapse 
!    due to the gravity  (Horizontal 2-D or 3D)
!     =================================================================
implicit none
integer :: i,j,m,n,iter,nodesum,it,ntimesed,kl,threads
integer :: istart,ifinal,iskip,jstart,jfinal,jskip
real(8) :: reposedry,reposewet,hmin
real(8) :: x(m,n),y(m,n),zb(m,n),zs(m,n),idix(m,n),dxc(m),dyc(n)
real(8) :: alenii(m,n,8),xxrep(m,n,8),yyrep(m,n,8),volrep(m,n,8),idixrep(m,n,8)   
real(8) :: slpbed(8),zbrep(8),dzbrep(8),zsrep(8),repose(8)  
real(8) :: noderep1,noderep2,arearep,alenrep,zbavrep,dzbrepii

    iter=0
101	nodesum=0
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

	!$omp parallel do private(zbrep,noderep1,noderep2,arearep,alenrep,zbavrep,kl,dzbrepii), reduction(+:nodesum), num_threads(threads)
		do j=jstart,jfinal,jskip       
		do i=istart,ifinal,iskip

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

			noderep1=0
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

				if(slpbed(kl).gt.repose(kl)+0.0001 .and. idixrep(i,j,kl).eq.1) then 
					noderep1=noderep1+1
					arearep=arearep+volrep(i,j,kl)
					alenrep=alenrep+alenii(i,j,kl)*volrep(i,j,kl)*repose(kl)
					zbavrep=zbavrep+zbrep(kl)*volrep(i,j,kl)
				endif
            enddo

            if(noderep1.gt.0) then
			
				dzbrepii=(-zb(i,j)*arearep+zbavrep+alenrep)/(arearep+dxc(i)*dyc(j))
				
				do kl=1,8
					if(slpbed(kl).gt.repose(kl)+0.0001 .and. idixrep(i,j,kl).eq.1) then
						dzbrep(kl)=zb(i,j)+dzbrepii-zbrep(kl)-repose(kl)*alenii(i,j,kl)
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
			noderep2=0
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

				if(slpbed(kl).lt.-repose(kl)-0.0001 .and. idixrep(i,j,kl).eq.1) then 
					noderep2=noderep2+1
					arearep=arearep+volrep(i,j,kl)
					alenrep=alenrep+alenii(i,j,kl)*volrep(i,j,kl)*repose(kl)
					zbavrep=zbavrep+zbrep(kl)*volrep(i,j,kl)
				endif
				
            enddo

            if(noderep2.gt.0) then
			
				dzbrepii=(-zb(i,j)*arearep+zbavrep-alenrep)/(arearep+dxc(i)*dyc(j))
				
				do kl=1,8
				
					if(slpbed(kl).lt.-repose(kl)-0.0001 .and. idixrep(i,j,kl).eq.1) then
						dzbrep(kl)=zb(i,j)+dzbrepii-zbrep(kl)+repose(kl)*alenii(i,j,kl)
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

		if(iter.gt.10) goto 102
		if(nodesum.ge.1) goto 101 
 102  continue

end subroutine

!     =================================================================
!     This subroutine handles the steep slope case
!           Made by Dr. Weiming Wu, ..............
subroutine bedrepos0(threads,m,n,x,y,zb,idix,dxc,dyc,it,alenii,   &
                     xxrep,yyrep,volrep,idixrep)
!    repose angle is considered here to check the bed profile
!    if the bed slope is larger than the repose angle, the bed will collapse 
!    due to the gravity  (Horizontal 2-D or 3D)
!     =================================================================
implicit none
integer :: i,j,m,n,iter,nodesum,it,ntimesed,kl,threads
integer :: istart,ifinal,iskip,jstart,jfinal,jskip
real(8) :: x(m,n),y(m,n),zb(m,n),idix(m,n),dxc(m),dyc(n)
real(8) :: alenii(m,n,8),xxrep(m,n,8),yyrep(m,n,8),volrep(m,n,8),idixrep(m,n,8)   
real(8) :: slpbed(8),zbrep(8),dzbrep(8)
real(8) :: repose,noderep1,noderep2,arearep,alenrep,zbavrep,dzbrepii

	repose=1.0

    iter=0
101	nodesum=0
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

	!$omp parallel do private(zbrep,noderep1,noderep2,arearep,alenrep,zbavrep,kl,dzbrepii), reduction(+:nodesum), num_threads(threads)
	do j=jstart,jfinal,jskip
    do i=istart,ifinal,iskip
    

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
				if(slpbed(kl).gt.repose+0.0001 .and. idixrep(i,j,kl).eq.1) then 
					noderep1=noderep1+1
					arearep=arearep+volrep(i,j,kl)
					alenrep=alenrep+alenii(i,j,kl)*volrep(i,j,kl)
					zbavrep=zbavrep+zbrep(kl)*volrep(i,j,kl)
				endif
            enddo

            if(noderep1.gt.0) then
				dzbrepii=(-zb(i,j)*arearep+zbavrep+repose*alenrep)/(arearep+dxc(i)*dyc(j))

				do kl=1,8
					if(slpbed(kl).gt.repose+0.0001 .and. idixrep(i,j,kl).eq.1) then
						dzbrep(kl)=zb(i,j)+dzbrepii-zbrep(kl)-repose*alenii(i,j,kl)
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

 106        noderep2=0
            arearep=0.0
            alenrep=0.0
            zbavrep=0.0
			
            do kl=1,8
				slpbed(kl)=(zb(i,j)-zbrep(kl))/alenii(i,j,kl)
				if(slpbed(kl).lt.-repose-0.0001 .and. idixrep(i,j,kl).eq.1) then 
					noderep2=noderep2+1
					arearep=arearep+volrep(i,j,kl)
					alenrep=alenrep+alenii(i,j,kl)*volrep(i,j,kl)
					zbavrep=zbavrep+zbrep(kl)*volrep(i,j,kl)
				endif
            enddo

            if(noderep2.gt.0) then
			
				dzbrepii=(-zb(i,j)*arearep+zbavrep-repose*alenrep)/(arearep+dxc(i)*dyc(j))
			 
				do kl=1,8
					if(slpbed(kl).lt.-repose-0.0001 .and. idixrep(i,j,kl).eq.1) then
						dzbrep(kl)=zb(i,j)+dzbrepii-zbrep(kl)+repose*alenii(i,j,kl)
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

      if(iter.gt.50) goto 102
      if(nodesum.ge.1) goto 101 
 102  continue

end subroutine

!.##################################################################################
subroutine dzdxyminmod(threads,m,n,dzdxmin,dzdymin,x,y,zs)
implicit none
integer :: i,j,m,n,threads
real(8) :: x(m,n),y(m,n),zs(m,n)
real(8) :: dzdxmin(m,n), dzdymin(m,n)
real(8) :: dzdxup,dzdxlw,dzdyup,dzdylw
	  
	!$omp parallel do num_threads(threads)
    do j=1,n
        dzdxmin(1,j)=(zs(2,j)-zs(1,j))/(x(2,j)-x(1,j))          
        dzdxmin(m,j)=(zs(m,j)-zs(m-1,j))/(x(m,j)-x(m-1,j))          
	enddo
	!$omp end parallel do 
		
	!$omp parallel do private(dzdxup,dzdxlw), num_threads(threads)
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

	!$omp parallel do num_threads(threads)
    do i=1,m
        dzdymin(i,1)=(zs(i,2)-zs(i,1))/(y(i,2)-y(i,1))          
        dzdymin(i,n)=(zs(i,n)-zs(i,n-1))/(y(i,n)-y(i,n-1)) 
	enddo
	!$omp end parallel do 
		
	!$omp parallel do private(dzdyup,dzdylw), num_threads(threads)
	do j=2,n-1
	do i=1,m
	
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


end subroutine

!#############################################################################
Subroutine Output(m,n,time,x,y,zb,zs2,h2,u2,v2,qx2,qy2,c2,it,Ioutput)
implicit none
real(8) :: time
integer :: i,j,m,n,it,ioutput
real(8) :: x(m,n),y(m,n),zb(m,n),zs2(m,n),h2(m,n),u2(m,n),v2(m,n),qx2(m,n),qy2(m,n),c2(m,n)
character(3) :: name

	WRITE(name,'(I3.3)') Ioutput
	OPEN(300,FILE='Output_'//name//'.tec')
			
	write(300,*) 'title= " Time= ', time, ' s " ' 
	write(300,*) 'VARIABLES ="X","Y","Z","Zs","h","U","V","W","Qx","Qy","Concen"'
	write(300,*) 'ZONE T="contour"',' I=',m,' J=',n,' F=point'
			
	do j=1,n
	do i=1,m
		write(300,9980) x(i,j),y(i,j),zb(i,j),zs2(i,j),h2(i,j),u2(i,j),v2(i,j),0.0d0,qx2(i,j),qy2(i,j),c2(i,j)
	enddo
	enddo
			
9980  	format(2F12.4,9E18.8)

	CLOSE(300)	
			
end subroutine

function wtime ( )
  implicit none
  integer ( kind = 4 ) clock_max
  integer ( kind = 4 ) clock_rate
  integer ( kind = 4 ) clock_reading
  real ( kind = 8 ) wtime
  call system_clock ( clock_reading, clock_rate, clock_max )
  wtime = real ( clock_reading, kind = 8 ) &
        / real ( clock_rate, kind = 8 )
  return
end function

end program 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
