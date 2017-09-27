        program lammpsdump2xyz
        implicit none
        character a*3
        integer nframes,nmaxatms,NPt,NWAT,K,L,m,n,p,t
		REAL*8 XXCM,YYCM,ZZCM
		REAL*8 XXCP,YYCP,ZZCP
		REAL*8 XX,YY,RSQ,SUM1
		REAL*8 SUMX, SUMY, SUMZ, SUMHIST,TOTSUM
		REAL*8 SUM1X, SUM1Y, SUM1Z, SUM2Z
		REAL*8 XSET,YSET,ZSET ,DENS
		REAL*8 XPSET,YPSET,ZPSET
		REAL*8 CMX, CMY, CMZ,AREA, CM
		REAL*8 CPX, CPY, CPZ, CP1Z
	    REAL*8 PI,BOX,DZ
		INTEGER ZBIN, RBIN,Nzb
		INTEGER MAXZBIN, MAXRBIN
		INTEGER ZFLAG, RFLAG

		PARAMETER (PI=3.14159265359)
		PARAMETER (ZPSET=57.359738495410234)
		PARAMETER (BOX=98.0D0)

!        parameter (nframes=10,nmaxatms=2e6)
        !parameter (nframes=8000,nmaxatms=11000)
		PARAMETER (MAXZBIN=5000)
		PARAMETER (MAXRBIN=5000)
		INTEGER HIST(MAXZBIN,MAXRBIN),HIST2(MAXRBIN),HIST1(MAXZBIN)
        !parameter (nframes=2,nmaxatms=11000)
        parameter (nframes=700,nmaxatms=11000)
        integer i,j,atm_id(nmaxatms),typ(nmaxatms),time(nframes)
        integer natms,ind,atmn_id(nmaxatms),typn(nmaxatms)
		integer NWATS(nframes),ZZ(nmaxatms)
		REAL*8 R(MAXRBIN), H(MAXZBIN)
		REAL*8 EDS(MAXZBIN)
        double precision x(nmaxatms),y(nmaxatms),z(nmaxatms)
        double precision ax(nmaxatms),ay(nmaxatms),az(nmaxatms)
        double precision ox(nmaxatms),oy(nmaxatms),oz(nmaxatms)
        double precision sox(nmaxatms),soy(nmaxatms),soz(nmaxatms)
        double precision hx(nmaxatms),hy(nmaxatms),hz(nmaxatms)
        double precision px(nmaxatms),py(nmaxatms),pz(nmaxatms)
        double precision shx(nmaxatms),shy(nmaxatms),shz(nmaxatms)
        double precision sh1x(nmaxatms),sh1y(nmaxatms),sh1z(nmaxatms)
        double precision chrg(nmaxatms),Z1CP(nmaxatms),Z1CM(nmaxatms)
        double precision XCM(nmaxatms),YCM(nmaxatms),ZCM(nmaxatms)
        double precision XCP(nmaxatms),YCP(nmaxatms),ZCP(nmaxatms)

        open(unit=11,file='positions.out')
!        open(unit=20,file='traj.xyz')
!
		AREA=7.0
		DZ=0.05
		RBIN=900
		ZBIN=1500
		natms = 10341

!
		DO I = 1, RBIN
		 R(I) = DSQRT(I*AREA/PI)
		ENDDO
!
		DO I = 1, ZBIN
		   H(I) = I*DZ
		ENDDO

		DO I = 1, ZBIN
		 DO J = 1, RBIN
		  HIST(I,J) = 0
		  ENDDO
		ENDDO

		DO I = 1, ZBIN
		 HIST1(I) = 0
		ENDDO

		DO I = 1, RBIN
		 HIST2(I) = 0
		ENDDO


		SUM1 = 0
        ind = 0
        do t = 1,nframes
!        ind=ind+1
!        write(*,*)ind
!        read(11,*)
!        read(11,*)time(k)
!        read(11,*)
!        read(11,*)natms
!        read(11,*)
!        read(11,*)
!        read(11,*)
!        read(11,*)
!        read(11,*)
		do j = 1,natms
		read(11,*)ax(j),ay(j),az(j)
		enddo
!
		do j=1,natms
		do i = 1,natms
!		if(atm_id(i).eq.j) then
!		atmn_id(j)=atm_id(i)
!		typn(j)=typ(i)
!		ax(j)=x(i)*BOX
!		ay(j)=y(i)*BOX
!		az(j)=z(i)*BOX
!		endif
		enddo
		enddo
!
		n=0
		m=0
		p=0
		do i = 1,natms
        if(i.le.2197)then
		n=n+1
        ox(n)=ax(i)
		oy(n)=ay(i)
		oz(n)=az(i)
        endif
!
        if(i.gt.2197 .and. i.le.6591)then
		m= m+1
		hx(m)=ax(i)
		hy(m)=ay(i)
		hz(m)=az(i)
		!write(22,*) m,atmn_id(i),typn(i),hx(m),hy(m),hz(m)
        endif
		if(i.gt.6591)then
		p=p+1
		px(p)=ax(i)
		py(p)=ay(i)
		pz(p)=az(i)
		endif
        enddo
		NWAT=n
		NPt=p
!
		do i = 1,NWAT
        sox(i)=ox(i)
		soy(i)=oy(i)
		!soz(i)=oz(i)-4.92
		soz(i)=oz(i)
!
		shx(i)=hx(2*i-1)
		shy(i)=hy(2*i-1)
		!shz(i)=hz(2*i-1)-4.92
		shz(i)=hz(2*i-1)
!
		sh1x(i)=hx(2*i)
		sh1y(i)=hy(2*i)
		!shz(i)=hz(2*i)-4.92
		sh1z(i)=hz(2*i)
		enddo
		!write(45,*) sox(1),shx(2),shx(3)
!
		do i = 1,NWAT
		XXCM=(sox(i)*16.0+shx(i)*1.008+sh1x(i)*1.008)
		YYCM=(soy(i)*16.0+shy(i)*1.008+sh1y(i)*1.008)
		ZZCM=(soz(i)*16.0+shz(i)*1.008+sh1z(i)*1.008)

		XCM(i) = XXCM / 18.016
		YCM(i) = YYCM / 18.016
		ZCM(i) = ZZCM / 18.016
		enddo
!
		do i = 1,NPt
		XCP(i) = px(i)
		YCP(i) = py(i)
		ZCP(i) = pz(i)
		enddo

		SUMX = 0.0D0
		SUM1X = 0.0D0
		SUMY = 0.0D0
		SUM1Y = 0.0D0
		SUMZ = 0.0D0
		SUM1Z = 0.0D0
		SUM2Z = 0.0D0
		do j = 1,NWAT
		  SUMX = SUMX + XCM(j)
		  SUMY = SUMY + YCM(j)
		  SUMZ = SUMZ + ZCM(j)
		enddo

		do j = 1,NPt
		  SUM1X = SUM1X + XCP(j)
		  SUM1Y = SUM1Y + YCP(j)
		  SUM1Z = SUM1Z + ZCP(j)
		enddo
		    CMX = SUMX/NWAT
		    CMY = SUMY/NWAT
		    CMZ = SUMZ/NWAT

		!    CPX = SUM1X/NPt
		!    CPY = SUM1Y/NPt
		!    CPZ = SUM1Z/NPt

		WRITE(42,*) (t*400), CMX, CMY, CMZ
		WRITE(43,*) (t*400), CMX*0.5292, CMY*0.5292, CMZ*0.5292

		enddo

        stop
        end
