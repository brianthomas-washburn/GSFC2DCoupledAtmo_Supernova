c
c    /misc/kah04/fleming/2dclinux> em radlib_9ga_nosad.f &
c    [3] 22717
c
      subroutine radstart
      include 'PARAM.INC'
      include 'COMMONC.INC'
      include 'COMMOND.INC'
      include 'COMMONR.INC'
      include 'aero.inc'
      logical cldfl,rinit
      common/clouds/ctop(19),cfrc(19,12),cldlat(19),ctop1(n$),
     * ncld,cldfl

      dimension jlbedo(37)
      data ndi/24/
      data rinit/.false./
      mz=45      !mzrad
      mz1=mz+1
c
      data jlbedo/84,79,72,59,43,28,19,17,12,7*8,3*7,2*8,2*9,5*10,11,
     * 15,19,26,35,45,53,58,61/
c
         ple(1)=0.180000e-02
         ple(2)=0.245944e-02
         ple(3)=0.336047e-02
         ple(4)=0.459160e-02
         ple(5)=0.627376e-02
         ple(6)=0.857219e-02
         ple(7)=0.117126e-01
         ple(8)=0.160036e-01
         ple(9)=0.218667e-01
         ple(10)=0.298776e-01
         ple(11)=0.408234e-01
         ple(12)=0.557793e-01
         ple(13)=0.762144e-01
         ple(14)=0.104136e+00
         ple(15)=0.142287e+00
         ple(16)=0.194414e+00
         ple(17)=0.265639e+00
         ple(18)=0.362957e+00
         ple(19)=0.495928e+00
         ple(20)=0.677613e+00
         ple(21)=0.925861e+00
         ple(22)=0.126506e+01
         ple(23)=0.172851e+01
         ple(24)=0.236176e+01
         ple(25)=0.322701e+01
         ple(26)=0.440924e+01
         ple(27)=0.602459e+01
         ple(28)=0.823173e+01
         ple(29)=0.112475e+02
         ple(30)=0.153680e+02
         ple(31)=0.209982e+02
         ple(32)=0.286910e+02
         ple(33)=0.392021e+02
         ple(34)=0.535639e+02
         ple(35)=0.731874e+02
         ple(36)=0.100000e+03
         ple(37)=0.135000e+03
         ple(38)=0.180000e+03
         ple(39)=0.250000e+03
         ple(40)=0.340000e+03
         ple(41)=0.435800e+03
         ple(42)=0.548700e+03
         ple(43)=0.661600e+03
         ple(44)=0.774500e+03
         ple(45)=0.887400e+03
         ple(46)=0.100000e+04
 
      do 13 j=1,37
      lat(j)=-90.+ (j-1)*5.
   13 d3(j)=real(jlbedo(j))
 
      call inter(d3x,yp,n$,d3,lat,37)
 
      do 14 j=1,n$
  14  ilbedo(j)=nint(d3x(j))
c
c
c...initialize model pressures
c
      ps=1000.
      his=0.0
      do 5435 n=1,mz
      dpl(n)=ple(n+1)-ple(n)
      pl(n)=sqrt(ple(n)*ple(n+1))
 5435 continue
       
c
c...full grid of model and edge pressures
c
      mztot=mz+mz1
      do 5436 n=1,mz
      nn=2*n
 5436 ptot(nn)=pl(n)
      do 5437 n=1,mz1
      nn=2*n-1
 5437 ptot(nn)=ple(n)
      do 5438 n=1,91
      zjr(n)=h*alog(1000./ptot(n))/1.e5 ! jr z's in km
 5438 zjrrev(92-n)=zjr(n)               ! jr z's reversed in km
      
      do 5439 n=1,m$
 5439 zprev(mp$-n)=zp(n) ! mrs zp's reversed
      do 5440 n=2,90,2
 5440 zjrmid(n/2)=zjrrev(n)

      zptop=0.
      do 5441 k=1,m$
         zptop=amax1( zptop, zp(k) )
 5441 continue

      do 5442 k=1,91 
         if( zjrrev(k) .le. zptop) mzrad2=k
 5442 continue

      mzrad=mzrad2/2 
      if(mzrad .lt. mz) then
        print *,'mzrad,mz = ',mzrad,mz
cjer        stop2
      endif
      mz=mzrad  
c
c...initialize radiation quantities
c
 8000 format( 'in radstart zptop, mzrad ,mzrad2, zjrtop ',
     >  f10.4,1x,i4,i4,f10.4)
      write( 6,8000) zptop,mzrad,mzrad2,zjrrev(mzrad2)


      aerfl=.false.
      volcano=.false.
      cldfl=.false.
      if(.not.rinit) then
        call radini(mz,yp,yy)
        rinit=.true.
      endif
 
      return
      end


      subroutine radcal(decmrs, YIN, ZIN)

c     called from getrad  

c     set up intervals and heights for rad model

      include 'aero.inc'
      include 'PARAM.INC'
      include 'COMMONC.INC'
      include 'COMMONR.INC'
      include 'COMMOND.INC'

      integer N$22
      PARAMETER (N$22=N$-22)
      
ccelf      include 'test_value.inc'
ccelf      logical nantest
      
      logical flag1
      
c  COMMOND.INC has common/time/ with day, month, etc.
cjer quantities for radiative flux at tropopause
c    for radiative forcing studies
      logical flag 
      real ptrop(n$),lwflux(n$),swflux(n$),netflx(n$)
      real irflx(46),solflx(46), yin(ns$),zin(ms$)
      common /radflx/ ntrop, solflx, irflx

      logical lprint,rad,plot,cogli,solar,ir,sphere,cldfl
      common/cool/coolr(48),sflux
      common/sun/decang,cosz,fday,cogli
      common/hght/fcloud,fclear,ntopf,
     1 ntopt,rhostd,dlat(n$),long,rlat(n$)
      common/mix/o3mix(48),co2mix(48),o2mix(48)
      dimension hetnet(n$,45)
cjer stuff for psc and subvisible cirrus heating
      common / pscld/ipsc(45),pscabs(45),pscext(45),psch(45),psc96(45),
     * psc15(45)
      common / pscpass/psctotal(N$,45),pscdyn(N$,M$) 

      common/cyrzero/yrzero

      common/solid /chemlat(ns$),sh2o(ns$,ms$),shno3(ns$,ms$)
      common/solid1/ssh2o(N$,ms$),sshno3(N$,ms$)
      dimension soltmp1(ns$),soltmp2(N$)
    

      common/clouds/ctop(19),cfrc(19,12),cldlat(19),ctop1(n$),
     * ncld,cldfl
      dimension cfrc1(n$),cldmon(19)
      dimension cfrc2(n$)
      data cfrc2/0.7,0.7,0.5,0.3,0.2,9*0.0,0.1,0.2,0.3,0.35,0.45,
     * 0.45,0.3,0.1,n$22*0.0/
ccelf originally, above cfrc2 assumes n$=36 (No. latitudes in dynamics), now pad w/ 0's using N$22 - EF 1/08
 
      real los
      data los/2.69e19/
c
      data idy/1/
c
      data solar/.true./, ir/.true./, cogli/.true./, sphere/.false./
      data ndi/24/

      logical phinit
      data phinit /.false./
      !integer yrzero, NS$, MS$, N$22
      integer yrzero

      mz=mzrad
      mz1=mz+1
 
      fnh=ndi
      decang=decmrs
 
c      write(6,20) decang
   20 format(/,' decang = ',f8.2,//)

C     Set Initial Year Number
      if (.not. phinit) then
         yrzero = nyear
         phinit = .true.
         endif

C     Get Year Number Relative to Start Year (Same as in multjr.f).
      ideltyr = nyear - yrzero

c      write(6,6010) ideltyr, month,nday
c 6010 format(' ideltyr, month, day = ',3i4)

c...select cloud data for month
c      if(cldfl) then
c        do 27 jj=1,ncld
c        cldmon(jj)=cfrc(jj,month)
c   27   continue
c        call inter(cfrc1,llat,n$,cldmon,cldlat,ncld)        
c      end if 

c...interpolate solid h2o and hno3 from chemistry lats to dynamics lats

      do jj=1,N$
      do kk=1,ms$
        ssh2o(jj,kk)=0.0
        sshno3(jj,kk)=0.0
      enddo
      enddo
      kflg=0
      do jj=1,ns$
      do kk=1,ms$
        if(sh2o(jj,kk).ne.0.0) kflg=1
      enddo
      enddo

      if(kflg.eq.1) then
            
c  loop over chemistry altitudes
        do 30 kk=1,ms$

           do 28 jj=1,ns$
              soltmp1(jj)=sh2o(jj,kk)
   28      continue
           call inter(soltmp2,yp,N$,soltmp1,chemlat,ns$)
           do 31 jj=1,N$
              if(soltmp2(jj) .lt. 0.0) soltmp2(jj)=0.0
              ssh2o(jj,kk)=soltmp2(jj)
   31      continue

           do 32 jj=1,ns$
              soltmp1(jj)=shno3(jj,kk)
   32      continue
           call inter(soltmp2,yp,N$,soltmp1,chemlat,ns$)
           do 33 jj=1,N$
              if(soltmp2(jj) .lt. 0.0) soltmp2(jj)=0.0
              sshno3(jj,kk)=soltmp2(jj)
   33      continue


   30   continue
      end if

c
c...beginning of lat loop
 
      do 90 j=1,n$
c      flag=.false.
c      if(dlat(j) .ge. 60. .and. dlat(j) .le. 61.) flag=.true.
c      do 90 j=19,19
c       write(6,6000) j,dlat(j)
6000  format(i6,' latitude = ',f8.2)

      do 650 n=1,mz
      cooljr(j,n)=0.0
  650 heatav(j,n)=0.0
 
c...get temps and density at midpoints and edges
 
      do 62 n=1,mztot
      patm=ptot(n)/1013.25
      rhotot(n)=patm*.3531/ttot(j,n)
   62 continue
      do 82 n=1,mz1
      nn=2*n-1
   82 tle(n)=ttot(j,nn)
      do 84 n=1,mz
      nn=2*n
      tl(n)=ttot(j,nn)
      rhojr(n)=rhotot(nn)
   84 continue
 
c...get height at midpoints
 
      t1=tle(mz1)
      t2=tl(mz)
      pr1=ple(mz1)
      pr2=pl(mz)
      hi(mz)=his-(t1+t2)*.5/34.17*alog(pr2/pr1)
      do 88 n=2,mz
      nn=mz-n+1
      pr1=pl(nn+1)
      pr2=pl(nn)
      t1=tl(nn+1)
      t2=tl(nn)
      hi(nn)=hi(nn+1)-(t1+t2)*.5/34.17*alog(pr2/pr1)
   88 continue
c
c...get heights at edges
c
      hie(mz1)=0.0
      do 8085 n=1,mz
      nn=mz-n+1
      t1=tl(nn)
      t2=tle(nn)
      pr1=pl(nn)
      pr2=ple(nn)
      hie(nn)=hi(nn)-(t1+t2)*.5/34.17*alog(pr2/pr1)
 8085 continue
c
c...surface quantities
c
      tgjr=tle(mz1)
      idum=ilbedo(j)
      rsurf=real(idum)/100.

c  get aerosol optical depths
      lprint=.false.
c      if(j.eq.19.and.nyear.ge.2.and.month.ge.4.and.idayx.eq.1) 
c     * lprint=.true.      
      if(aerfl) call aeros1(mz,ideltyr,month,nday,j,lprint)
c  get clouds
      do 652 n=1,mz
  652 cloud(n)=0.0
      if(cldfl) then
cjer if use cloud data set
c         ictop = nearidx(ctop1(j),ple,mz)
c         cloud(ictop) = abs(cfrc1(j))
c         write(6,656) ple(ictop)
c  656 format(' cloud top pressure = ',f10.1)
c      end if
cjer
c if don't use cloud data set
c
        ictop = nearidx(700.,ple,mz)
        cloud(ictop) = 0.5
        ictop = nearidx(250.,ple,mz)
        cloud(ictop) = cfrc2(j)
cjer
      end if


c
c...get mixing ratios
c
c      lprint=.true.
      lprint=.false.
c      print *,'before radred'
      call radred(mz,j,lprint)
c      print *,'after radred'

c...get psc and sub_visible optical depths for heating calculation
      lprint=.false.
c      lprint=.true.
        call psc(j,dlat(j),lprint, YIN, ZIN)
        do 800 kk=1,45
        psctotal(j,kk)=pscext(kk)
  800   continue
c      endif

c find tropopause level
c minimum temp definition
c start looking at 50 mb
      do n=1,mz
        nstart=n
        if(ple(n+1) .gt. 50.) go to 5008
      enddo
 5008 continue
      
      do n=nstart,mz1
        ntrop=n
        if(tle(n+1).gt.tle(n)) go to 5009
      enddo
 5009 continue
      ptrop(j)=ple(ntrop)

c
c      print *,'before htrdly'
      if(.not.solar) go to 6300

c      write(6,8) dlat(j),cogli,sphere
c    8 format('latitude = ',f8.2,' cogli = ',l4,'   sphere = ',l4)

      if(cogli) go to 8600
 
c...compute diurnal average solar heating   (not cogli)
 
      jalb=j
      fday=1.0
      sinf=sin(rlat(j))
      cosf=cos(rlat(j))
      sind=sin(decang*pi/180.)
      cosd=cos(decang*pi/180.)
      do 654 id=1,ndi
      hour=id*(24./fnh)
      hcos=cos(hour*pi/12.)
      cosz=sinf*sind+cosf*cosd*hcos
      angle=acos(cosz)
c      write(6,657) id,hour,cosz
  657 format(' id,hour,cosz,angle:',i4,f6.1,f8.4)
      angle=angle*180./pi
      lprint=.false.
      if(j.eq.19) lprint=.true.
      call htrdly(mz,j,lprint,sphere)
      do 653 n=1,mz
  653 heatav(j,n)=heatav(j,n)+as(n)
  654 continue
      do 655 n=1,mz
  655 heatav(j,n)=heatav(j,n)/fnh
      goto 6300
 
c  cogli and borucki diurnal average
 
 8600 continue
      lprint=.false.
c      if(j.eq.19) lprint=.true.
c
C      PRINT*, "Doing Cogli ..."
c      IF (fday .NE. fday) THEN
c       nantest=test_nan(fday)
c       if(nantest) then
c         PRINT*, "FDAY NAN at COGLI:"
c         PRINT*, "   FDAY = ", FDAY
c         PRINT*, "   COSZ = ", COSZ
c         STOP
c         ENDIF
      call htrdly(mz,j,lprint,sphere)
      do 6530 n=1,mz
      heatav(j,n)=as(n)
 6530 continue
 6300 continue
c
c      print *,'before irdriv'
      if(.not.ir) go to 658
c
c...compute longwave cooling
c
      lprint=.false.
c      lprint=.true.
c      if(j.eq.18) lprint=.true.
      call irdriv(mz,lprint)
      do 651 n=1,mz
      cooljr(j,n)=coolr(n+3)
c      nantest=test_nan(cooljr(j,n))
c      if(nantest) then
c        print *,'cooljr nan in radcal: j,n,coolr = ',j,n,coolr(n+3)
c      endif	
  651 continue
c      print *,'after irdriv'
      olr(j) = tir
      
      do 110 n=1,mz
      hetnet(j,n)=heatav(j,n)-cooljr(j,n)
  110 continue
  658 continue

      lprint=.false.
c      if(j.eq.19) lprint=.true.
      if(.not. lprint) go to 705
      write(6,5000) j,dlat(j)
 5000 format(/,' heating and cooling for j = ',i3,' lat = ',f8.2)
      write(6,5001) 
 5001 format(2x,'n',6x,'pl',5x,'tl',9x,'solar',2x,'ir',5x,'net',
     * 4x,'cloud',4x,'height')
      do 5002 n=1,mz
      write(6,5003) n,pl(n),tl(n),heatav(j,n),cooljr(j,n),
     * hetnet(j,n),cloud(n),hi(n)
 5002 continue
 5003 format(' ',i3,f10.4,1x,f7.2,2x,3(f7.3,1x),f6.2,2x,f6.2)
  705 continue
  

c      lprint=.true.
      lprint=.false.
      if(.not.lprint) go to 706
      write(6,5004)
 5004 format(//,' solar (down), ir (up), net flux (up)(W/M**2)')
      do 5005 n=1,mz1
      write(6,5006) n,ple(n),tle(n),solflx(n),irflx(n),netflx(n)
 5005 continue
 5006 format(i4,f10.4,3x,f10.2,2x,f10.2,2x,f10.2,2x,f10.2)
  706 continue

      swflux(j)=solflx(ntrop)
      lwflux(j)=irflx(ntrop)
      netflx(j)=irflx(ntrop)-solflx(ntrop)

      if(.not.lprint) go to 5012
      write(6,5011)
 5011 format(//,' tropopause p, t, and fluxes',/)
      write(6,5010) ple(ntrop),tle(ntrop),solflx(ntrop),irflx(ntrop),
     *  netflx(j)
 5010 format(5f10.2)
 5012 continue


   90 continue  !...end of lat loop
   
   
c      stop1
      
c write our tropopause fluxes for radiative forcing
c      if(idayx .eq. 15) then
c       fiyr=float(iyr)
c       fmon=float(month)
c       write(69) fiyr,fmon
c       write(69) ptrop,swflux,lwflux,netflx
c      endif


      return
      end
 
      subroutine radini(mz,lat,xlat)

c     this routine initializes all the values needed for the
c     radiation parameterization
 
      include 'PARAM.INC'
      include 'COMMONR1.INC'
      include 'aero.inc'
 
      logical cldfl
      real los
      real lat(n$),xlat(n$)
      common/clouds/ctop(19),cfrc(19,12),cldlat(19),ctop1(n$),
     * ncld,cldfl
      common/hght/fcloud,fclear,ntopf,
     1 ntopt,rhostd,dlat(n$),long,rlat(n$)
      common/h2o/gl(30,20),b250(41,2),coeff(30,20,2)
      common/init/nup1,nsb,nc,npt,nw,nt,im
      dimension taustd(49,49),tauco2(48,48),fsn(19600)
      common/co2/tauco2,taustd,fsn
      logical co2var
      common/co2mr/co2vmr(45),co2lw(48),co2var

      data ctop/19*350./
      data pi/3.1415926/
 
      do 5 jj=1,n$
      dlat(jj)=lat(jj)
      rlat(jj)=xlat(jj)
    5 continue
      rhostd=1.2254e-3
      nup1=20
      npt=mz + 3
      nsb=1
      nw=22
      nt=25


c   ...b250 is integrated planck function for h2o bands
c   ...i=1,2 refers to band center and band wing regions
c   ...it=1,nt refers to nt temperatures
      read(1,*)
c      write(6,201)
  201 format(' reading b250 values')
      read(1,33)((b250(it,i),it=1,nt),i=1,2)
c      write(6,*)((b250(it,i),it=1,nt),i=1,2)
c   ...gl is h2o transmittance at 250 deg.
      read(1,*)
c      write(6,202)
  202 format(' reading gl values')
      read(1,34)((gl(iw,iu),iw=1,nw),iu=1,nup1)
c     write(6,34)((gl(iw,iu),iw=1,nw),iu=1,nup1)
c   ...coeff are coefficients for temperature expansion of transm.
      read(1,*)
c     write(6,203)
  203 format(' reading coeff values')
      read(1,34)(((coeff(iw,iu,i),iw=1,nw),iu=1,nup1),i=1,2)
   33 format(11f10.2)
   34 format(11f7.4)
      close(unit=1)
c     write(6,34)(((coeff(iw,iu,i),iw=1,nw),iu=1,nup1),i=1,2)
 
c      read data for fixed co2 transmission model
c
      mpt=49
      mptm1=mpt-1
      read(15,369) taustd,fsn
369   format(8e14.8)
      close(unit=15)

c  read aerosol and cloud flags - NO VOLCANO here - EF
      open (unit=200, file='radiation_back.dat', 
     >                        status='old',form='formatted')
      read(200,*) aerfl
      read(200,*) volcano
      read(200,*) cldfl
      read(200,*) co2var
      if(co2var) read(200,*) co2vmr
      close (unit=200)

      write(6,6005) aerfl,volcano,cldfl
 6005 format(//,' read radiation_back.dat:',/,' aerfl = ',l4,
     >          ' volcano = ', l4,' cldfl = ',l4,//)
      if(.not.co2var) print *, 'FIXED CO2'
      if(co2var) print *, 'VARIABLE CO2'
      if(co2var) print *,'chemistry CO2 interactive with radiation'
      print *, '  '

cjer not using this cloud data set now
c      if(cldfl) then
c        open (unit=200,file='cloud_fr.dat',status='old',
c     c        form='formatted') 
c        read(200,*) ncld
c        read(200,*) cfrc
c        close (unit=200)
c        do 14 j=1,ncld
c        cldlat(j)=-90. + (j-1)*10.
c   14 continue
c        call inter(ctop1,llat,n$,ctop,cldlat,ncld)
c      end if

      return
      end
 
 
      block data
c
      include 'COMMONR1.INC'
c
      data tcond/37*0.0,
     *               1.0,2.0,4.0,6.0,6.0,8.0,8.0,8.0/
      data tpene/38*0.0,7*8.0/
      data tlowl/16./
      data tmidl/8.0/
      data fk/0.107,0.104,0.073,0.044,0.025/
      data xk/0.005,0.041,0.416,4.752,72.459/
      data nfk/5/
      end
 
      subroutine radred(mzdum,jlat,lprint)
 
      include 'aero.inc'
      include 'PARAM.INC'
      include 'COMMONC.INC'
      
ccelf      include 'test_value.inc'
ccelf      logical nantest1,nantest2
      logical co2var
 
c...this routine prepares vertical arrays for radiation calculation
 
      dimension rhot(48),tmp(48),pmb(48),z(48),wat(48),
     1 pr1(m$),oz1(m$),watcol(m$),co21(m$),cloud1(48)
      dimension vmrh2o(48)
 
      common/ozot/ vmro3(48)
      common/hght/fcloud,fclear,ntopf,
     1 ntopt,rhostd,dlat(n$),long,rlat(n$)
      logical lprint,cogli
      common/sun/decang,cosz,fday,cogli
      common/mix/o3mix(48),co2mix(48),o2mix(48)
      common/co2mr/co2vmr(45),co2lw(48),co2var
 
      common/o3/co2m(n$,m$),o3m(n$,m$),h2om(n$,m$),wlat(18),zwat(23),
     1watz(18,23,12),xlsbuv(18),zsbuv(17),sbuv(18,17,12)
 
      include 'COMMONR1.INC'
 
      common/cirrus/icir(45),hcir(45),sigcir(45),taucir(45),epsa(45),
     1 eps11(45),tauvis(45),prmax
 
      real los,mubar,mu2bar
      data los/2.69e19/
      data twopi/6.28319/
      data co2frc/3.3e-4/
      data bcirr/1.6e-4/, t0cir/82.5/
      data rhostd/1.2254e-3/
      mzd1=mzdum+1
      mz=mzdum+2
      mzp1=mz+1
 
c...reverse top down arrays from gcm
 
      do 1012 n=1,mzdum
      nn=mzdum+1-n
      pmb(n+1)=pl(nn)
      tmp(n+1)=tl(nn)
      z(n+1)=hi(nn)
      rhot(n+1)=rhojr(nn)
      cloud1(n+1)=cloud(nn)
 1012 continue
 
c...surface quantities from gcm
 
      tmp(1)=tgjr
      z(1)=his
      pmb(1)=ple(mzd1)
      rhot(1)=.3531*ple(mzd1)/(tgjr*1013.25)
      wat(1)=shle(mzd1)
c
c...add 2 extra levels at top
c
      delog=alog(pl(1)) - alog(pl(2))
      pmb(mz)=exp(alog(pl(1)) + delog)
      pmb(mzp1)=exp(alog(pl(1)) + 2.*delog)
      pl1=alog(pmb(mzdum))
      pl2=alog(pmb(mzd1))
      tl1=tmp(mzdum)
      tl2=tmp(mzd1)
      do 5 ip=mz,mzp1
      prl=alog(pmb(ip))
      tmp(ip)=tl2+(tl1-tl2)*(prl-pl2)/(pl1-pl2)
      patm=pmb(ip)/1013.25
      rhot(ip)=patm*.3531/tmp(ip)
      z(ip)=z(ip-1)-(tmp(ip)+tmp(ip-1))*.5/34.17*alog(pmb(ip)/pmb(ip-1))
    5 continue
 
c      get o3, h2o and co2 from the model
c      m$ is number of vertical levels in dynamics
 
       do 31 k=1,m$
       watcol(mp$-k)=h2om(jlat,k)
       oz1(mp$-k)=o3m(jlat,k)
       if(co2var) co21(mp$-k)=co2m(jlat,k)
31     pr1(mp$-k)=1000.*rho(k)
       nno=m$
 
c      get o3 mixing ratios
 
       call conmix(vmro3,pmb,mzp1,oz1,pr1,nno,lprint)
 
c      water vapor  mixing ratios
 
       call conmix(vmrh2o,pmb,mzp1,watcol,pr1,nno,lprint)
 
      do 150 n=1,mzp1
      o3mix(n)=vmro3(n)*rhot(n)/rhostd*1.e-6
  150 continue
      do 25 n=1,mzrad
      shl(mzrad+1-n)=vmrh2o(n+1)*0.622e-6
25    continue
c      get co2 mixing ratios
 
c      if(co2var) call conmix(co2vmr,pmb,mzp1,co21,pr1,nno,lprint)
      if(co2var) call conmix(co2vmr,pl,mzdum,co21,pr1,nno,lprint)
 
ccccc arrays for htrdly (relative to loscmidt's number)  cccccccccccccc
      
      do 10 n=1,mzdum
      if(co2var) co2mix(n)=co2vmr(n)*rhot(n)/rhostd
      if(.not. co2var) co2mix(n)=co2frc*rhot(n)/rhostd
      o2mix(n)=0.20947*rhot(n)/rhostd
   10 continue
   
      do 11 n=mzd1,mzp1
      if(co2var) co2mix(n)=co2mix(mzdum)
      if(.not. co2var) co2mix(n)=co2frc*rhot(n)/rhostd
      o2mix(n)=0.20947*rhot(n)/rhostd
   11 continue

      if(.not.cogli) go to 190
 
cccccccccccccccccccc   diurnal average quantities  ccccccccccccccccccc
c...to: half day length (cogley & borucki, 1976)
c...mubar: average direction cosine of sun zenith angle
c...mu2bar: average of squared direction cosines of sun zenith angle
 
       phi=rlat(jlat)
       delta=decang/57.296

C       PRINT*, "In RADRED:"
C       PRINT*, "   DECANG = ", DECANG
C       PRINT*, "    DELTA = ", DELTA

       sind=sin(delta)
       acb=sin(phi)*sind
       bcb=cos(phi)*cos(delta)
       if(bcb.eq.0..and.phi.gt.0..and.sind.gt.0.) go to 120
       if(bcb.eq.0..and.phi.gt.0..and.sind.le.0.) go to 200
       if(bcb.eq.0..and.phi.lt.0..and.sind.lt.0.) go to 120
       if(bcb.eq.0..and.phi.lt.0..and.sind.ge.0.) go to 200
       ratio=acb/bcb
       if(ratio.ge.1.0) go to 120
       if(ratio.le.-1.0) go to 200
       t0=acos(-acb/bcb)/twopi
       go to 210
c...sun does not rise
  200  cosz=0.0
       fday=0.0
       go to 190
c...sun does not set
  120  t0=0.5
  210  continue
       mubar=2.*(acb*t0+bcb*sin(twopi*t0)/twopi)
       mu2bar=2.*(t0*acb**2+2.*acb*bcb*sin(twopi*t0)/
     *           twopi+bcb**2*(0.5*t0+sin(2.*twopi*t0)/
     *           (4.*twopi)))

C     Added 9/23/96 by PEM.  Above calculation allows MU2BAR to be negative
C     under some conditions (through roundoff-type errors, I believe).  This
C     is non-physical and will cause B4 to be NAN.  MU2BAR is restricted to
C     be no less than 1.0E-06.

      IF (MU2BAR .LE. 1.0E-06) MU2BAR = 1.0E-06

            b4=sqrt(2.*t0/mu2bar)
            a4=0.5*mubar*b4
            cosz=1./b4
            fday=2.*a4

c      IF ((cosz .NE. cosz) .OR. (fday .NE. fday)) THEN
c        nantest1=test_nan( cosz )
c        nantest2=test_nan( fday )
c        if(nantest1 .or. nantet2 ) then
        
c         PRINT*, "#1 NAN in RADRED:"
c         PRINT*, "    COSZ = ", COSZ
c         PRINT*, "    FDAY = ", FDAY
c         PRINT*, "      B4 = ", B4
c         PRINT*, "      A4 = ", A4
c         PRINT*, "      T0 = ", T0
c         PRINT*, "   MUBAR = ", MUBAR
c         PRINT*, "  MU2BAR = ", MU2BAR
c         PRINT*, "     ACB = ", ACB
c         PRINT*, "     BCB = ", BCB
c         PRINT*, "   TWOPI = ", TWOPI
c         STOP
c         ENDIF

C     End Added 9/23/95 by PEM.

  190 continue
 
c  cloud top, total and fractional cloudiness parameters
c  this code is taken from solar1 of the glas gcm
c    cloud fraction is the maximum of the fractions of all cloud types
 
      fcloud = 0.0
      fclear = 1.0
      ntopt = mzdum + 1
      ntopf = mzdum + 1
      icld=1
         do 240 n=icld,mzdum
         xx = cloud(n)
         if (xx.lt.0.01)  go to 240
         if (xx.gt.0.99)  go to 230
      fc = amax1(xx,fcloud)
      fcloud = fc
      fclear = 1.0 - fcloud
         if (ntopf.lt.mzdum)  go to 240
      ntopf = n
         go to 240
  230    continue
         if (ntopt.lt.mzdum)  go to 240
      ntopt = n
      fclear = 0.0
  240          continue
               if (fclear.gt.0.99)  go to 300
      if (fcloud.lt.0.01)  fcloud=1.00
  300 continue
 
c cirrus quantities
      prmax=400.
      do 500 n=1,mzdum
      icir(n)=0
      if(ple(n+1).lt.prmax .and. ple(n) .ge.100.0 
     * .and.cloud(n).gt.0.0) icir(n)=1
  500 continue
      do 520 n=1,mzdum
      if(icir(n).ne.1) go to 520
      hcir(n)=hie(n)-hie(n+1)
      taucir(n)=0.0
      tauvis(n)=0.0
      if(tl(n).lt.190.5) go to 520
      tt= tl(n)-273.+t0cir
      sigcir(n)= bcirr*tt*tt
      taucir(n)=sigcir(n)*hcir(n)
      tauvis(n)=2.*taucir(n)
c   beam emittance
      epsa(n)=1.0 - exp(-taucir(n))
c   diffusivity factor of 2.2 gives eps agreeing with platt and harsh
      eps11(n)=1.0 - exp(-2.2*taucir(n))
  520 continue
 
      if(.not.lprint) go to 425
      write(6,394)
  394 format(/,' printout from radred')
      write(6,395)
  395 format(' ','i ','n',3x,'z(n)',6x,'press',8x,'temp',5x,
     * 'vmrh2o(ppm)',2x,
     * 2x,'vmro3(ppm)',5x,'cloud',/)
      do 400 n=1,mzp1
      write(6,410)n,z(n),pmb(n),tmp(n),vmrh2o(n),vmro3(n),
     *  cloud1(n)
  410 format(' ',i3,f8.2,2x,f10.4,f10.3,2x,e11.4,2x,f8.2,2x,f6.2)
  400 continue
      write(6,*) 'co2vmr'
      do n=1,mzdum
        write(6,*) n,co2vmr(n)
      enddo	
      write(6,415) ntopt,ntopf
  415 format(//,' ntopt,ntopf:',2i5)
      write(6,416) fclear,fcloud
  416 format(' fclear,fcloud:',2f10.2)
      do 419 n=1,mzdum
      if(icir(n).eq.1) write(6,418) tl(n),hcir(n),sigcir(n),taucir(n),
     *  epsa(n),eps11(n),tauvis(n)
  418 format(/,6x,'tl',4x,3x,'hcir',3x,2x,'sigcir',2x,2x,'taucir',
     *  2x,3x,'epsa',3x,3x,'eps11',2x,'tauvis:',/,
     *  2f10.2,2f10.4,f10.4,2f10.4)
  419 continue
      write(6,420) cogli,decang,cosz,fday
  420 format(1x,'cogli,decang,cosz,fday:',l4,f8.3,f10.6,f8.2)
  425 continue
       return
       end
 
      subroutine conmix(cmix,pmb,mzp1,oz,pr,npo3,lprint)
 
c this routine calculates the constituent mixing ratios
c by logarithmic interpolation onto the radiation model grids
 
c       input:
c       pmb(mzp1)  are the model pressure levels
c       oz(npo3) is the constituent volume mixing ratio at npo3 pressure
c       pr(npo3) are the pressure levels
c       lprint = diagnostic switch
 
c       output:
c       cmix(mzp1) the interpolated mixing ratio values

      include 'PARAM.INC'
 
      dimension pr(m$),oz(m$),plog(m$),olog(m$),
     * ppt(m$),ozpt(m$)
      dimension cmix(48),pmb(48),pres(48),o3log(48)
      logical lprint
 
      do 11 n=1,npo3
      olog(n)=alog(oz(n))
      plog(n)=alog(pr(n))
11    continue
 
 
      do 101 n=1,mzp1
101   pres(n)=alog(pmb(n))
 
      call inter(o3log,pres,mzp1,olog,plog,npo3)
 
      do 102 n=1,mzp1
102   cmix(n)=exp(o3log(n))
 
      return
      end
 
      subroutine htrdly(nzdum,jlat,lprint,sphere)
c...main solar heating routine
      include 'PARAM.INC'
      include 'aero.inc'
      
      real solflx(46),irflx(46)
      common/radflx/ntrop, solflx, irflx
      
      real cnf(48),heat(48,8),n0,m,i1,i2,ls,ll,o3(48),o2(48),
     *no2(48),no3(48),delz(48),ptoth(48),
     * n2(48),k1(48),k2(48),k3(48),k6(48),haeff(48),hueff(48),
     * cph(48)
      real tmp(48),pmb(48),z(49),rhot(48),totht(48),wat(48)
      common/hght/fcloud,fclear,ntopf,
     1 ntopt,rhostd,dlat(n$),long,rlat(n$)
      common/mix/o3mix(48),co2mix(48),o2mix(48)
      logical lprint,cogli
      logical sphere,twilit
      common/sun/decang,cosz,fday,cogli
      dimension heato3(48),heato2(48),co2tmp(48)
      dimension flux(45,8),flxsum(8)
      dimension b4(50)
      
      include 'COMMONR1.INC'
 
      data rearth/6371.e5/
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      data n0/2.69e19/,sigc/2.85e-21/,fc/3.7e5/,sigha/8.8e-18/,fha/4580.
     */
      data ae1/.08/,ae2/2.6e-4/
      data sigsrc/1.e-17/,fsrc/1.1/
      data flux/360*0.0/
c    data for new huggins model
	data em1/0.14516/,shu1/3.3184/,fa1/646.20/,con11/3.91342e-19/
        data em2/0.10058/,shu2/6.9003e-6/,fa2/373.64/
        data con21/6.03748e-13/,con22/2.68669e-13/
c this calculates the solar heating rates for o3 and o2 see strobel
 
      pi=3.1415926
      nz=nzdum+2
      nzd1=nzdum+1
      nzpp1=nz+1
      nz1=nzpp1+1
      nb=6
      nb1=7
      nbt=8
      if(lprint) write(6,394) fclear,fcloud,ntopt,ntopf
  394 format(/,' printout from htrdly',/,'    fclear,fcloud,ntopt,
     *ntopf:',2f4.1,2i4)
      sg=0.
      do 396 n=1,nzdum
  396 as(n)=0.0
      if(cosz.lt.-0.13) return
      if(.not.sphere .and. cosz.lt.0.001) return
      secz=1./cosz
      sinz=sqrt(1.-cosz*cosz)
      ncl=nzdum +2 - 0.5*(ntopt+ntopf)
      nc1=ncl-1
      if(cosz.lt.0.001) go to 392
c...get heating due to h2o
      if(.not.aerfl) call solar5(nzdum,jlat,fday,cosz,lprint)
      if(aerfl) call solaer(nzdum,jlat,fday,cosz,lprint)
  392 continue
      twilit=.false.
      if( (sphere) .and. (cosz.lt.0.1)  ) twilit=.true.
c     if(twilit) write(6,6005) twilit
 6005 format(' twilit = ',l4)
      if(lprint) write(6,89) fclear,fcloud,ntopt,ntopf
   89 format(' after return from solar5',/,'  fclear,fcloud,ntopt,
     *ntopf:',2f4.1,2i4)

ccelf      ntime=ntime+ntime2-ntime1
ccelf      time=real(ntime)/100.
c     write(6,88) time
   88 format(' total time in solar5 is ',f7.2,' seconds'/)

c...loop over longitudes
c
      if(lprint) write(6,398) jlat,rsurf,ncl,cosz,fday
  398 format(' jlat,rsurf, ncl,cosz,fday:',i5,f10.3,i5,f10.6,f10.3)
c...invert arrays passed from gcm
      do 1012 n=1,nzdum
       nn=nzdum+1-n
       pmb(n+1)=pl(nn)
       tmp(n+1)=tl(nn)
       z(n+1)=hi(nn)
       rhot(n+1)=rhojr(nn)
 1012 continue
c
c...add 2 levels at top and extrapolate temperature
c
      delog=alog(pl(1))-alog(pl(2))
      pmb(nz)=exp(alog(pl(1))+delog)
      pmb(nzpp1)=exp(alog(pl(1))+2.*delog)
      pl1=alog(pmb(nzdum))
      pl2=alog(pmb(nzd1))
      tl1=tmp(nzdum)
      tl2=tmp(nzd1)
      do 393 ip=nzd1,nzpp1
      prl=alog(pmb(ip))
      tmp(ip)=tl2+(tl1-tl2)*(prl-pl2)/(pl1-pl2)
      patm=pmb(ip)/1013.25
      rhot(ip)=patm*.3531/tmp(ip)
      z(ip)=z(ip-1)-(tmp(ip)+tmp(ip-1))*.5/34.17*alog(pmb(ip)/pmb(ip-1))
  393 continue
      tmp(1)=tgjr
      z(1)=his
      pmb(1)=ple(nzd1)
      rhot(1)=.3531*ple(nzd1)/(tgjr*1013.25)
      do  1  n=1,nzpp1
      z(n)=z(n)*1.e5
      k1(n)=2.e-11*exp(107./tmp(n))
      k2(n)=2.9e-11*exp(67./tmp(n))
      k3(n)=2.e-15
      k6(n)=2.2e-18*(tmp(n)/300.)**0.8
      ptoth(n)=pmb(n)/1013.25
    1 o2mix(n)=0.20947*rhot(n)/rhostd
      do  2  n=1,nzpp1
      o2(n)=n0*o2mix(n)
      o3(n)=n0*o3mix(n)
    2 n2(n)=.78084*n0*rhot(n)/rhostd
c...number densities assumed to have exponential form for column
c...density integration. 5.2, 3.1 km are approximate scale heights
      do  3  n=2,nzpp1
      i=nz1-n
    3 delz(i)=z(i+1)-z(i)
      delz(nzpp1)=delz(nz)
c...section for spherical earth
                if (.not.(sphere.and.twilit)) go to 601
                z(nz1)=z(nzpp1) + delz(nzpp1)
                if(cosz.lt.0.) go to 607
c... the sun is above the horizon
c     write(6,6006)
 6006 format(' the sun is above the horizon')
                do 305 i=1,nzpp1
                ip1=i+1
                b4(i)=0.
                do 304 j=ip1,nz1
                path=(rearth+z(j))**2-
     1          (rearth+z(i))**2*sinz**2
                path=sqrt(path)-(rearth+z(i))*cosz
c...total path from j to i
 304            b4(j)=path
                do 306 j=ip1,nz1
                jrev=ip1+nz1-j
c...incremental path through j (a function of i)
  306           b4(jrev)=b4(jrev)-b4(jrev-1)
c               if (i.eq.1) write (6,309) cosz,b4
c309c          format(' cosz=',e11.3/' b4=',(9e11.3))
                no2(nzpp1)=o2(nzpp1)*5.2e5*b4(nz1)/delz(nzpp1)
                no3(nzpp1)=o3(nzpp1)*3.1e5*b4(nz1)/delz(nzpp1)
                if(i.eq.nzpp1) goto 305
                no2(i)=no2(nzpp1)
                no3(i)=no3(nzpp1)
                do 214 n=i,nz
                r1=o2(n+1)/o2(n)
                r2=o3(n+1)/o3(n)
                if(r1.eq.1.0) r1=1.001
                if(r2.eq.1.0) r2=1.001
c...total column density above i
                no2(i)=no2(i)+o2(n)*(r1-1.)/alog(r1)*b4(n+1)
                no3(i)=no3(i)+o3(n)*(r2-1.)/alog(r2)*b4(n+1)
  214           continue
  305           continue
                goto 602
c...sun is below the horizon
  607           continue
c     write(6,6007)
 6007 format(' the sun is below the horizon')
                no2(1)=1.e30
                no3(1)=1.e25
                do 380 i=2,nzpp1
                rtem=(rearth+z(i))*sinz
                no2(i)=1.e30
                no3(i)=1.e25
c...if path intersects earth
                if(rtem .le.(rearth+z(1)))goto 380
                do 371 j=1,nzpp1
                jtem=j
c...determine level above which path doesn't intersect earth
                if((rearth+z(j)) .gt. rtem) goto 372
  371           continue
  372           j=jtem
                jp1=j+1
                b4(j)=0.
                do 374 k=jp1,nz1
                path=(rearth+z(k))**2-(rearth+z(j))**2
                path=sqrt(path)
  374           b4(k)=path
                do 376 k=jp1,nz1
                krev=jp1+nz1-k
                b4(krev)=b4(krev)-b4(krev-1)
  376           continue
c                write (6,391) cosz,b4
c  391           format(' 391 cosz=',e11.3/' b4=',(9e11.3))
                no2(nzpp1)=o2(nzpp1)*5.2e5*b4(nz1)/delz(nzpp1)
                no3(nzpp1)=o3(nzpp1)*3.1e5*b4(nz1)/delz(nzpp1)
                if(i.eq.nzpp1) goto 380
                no2(i)=no2(nzpp1)
                no3(i)=no3(nzpp1)
                do 378 n=j,nz
                ifac=1.
                if(n.le.i) ifac=2
                r1=o2(n+1)/o2(n)
                r2=o3(n+1)/o3(n)
                if (r1.eq.1.) r1=0.001
                if (r2.eq.1.) r2=0.001
                no2(i)=no2(i)+o2(n)*(r1-1.)/alog(r1)*b4(n+1)*ifac
                no3(i)=no3(i)+o3(n)*(r2-1.)/alog(r2)*b4(n+1)*ifac
 378            continue
 380            continue
                goto 602
 601            continue
      no2(nzpp1)=o2(nzpp1)*5.2e5*secz
      no3(nzpp1)=o3(nzpp1)*3.1e5*secz
      do  4  n=2,nzpp1
      i=nz1-n
      r1=o2(i+1)/o2(i)
      r2=o3(i+1)/o3(i)
      if (r1 .eq. 1.0) r1=1.001
      if (r2 .eq. 1.0) r2=1.001
      no2(i)=no2(i+1)+delz(i)*o2(i)*(r1-1.)/alog(r1)*secz
    4 no3(i)=no3(i+1)+delz(i)*o3(i)*(r2-1.)/alog(r2)*secz
 602  continue

c...ozone absorption
c  albt is the albedo of the reflecting region for the chappuis
c  band absorption.  it includes the effective albedo of the
c  lower atmosphere, albc, and the ground albedo, rsurf. (lacis &
c  hansen)
 
c     clear sky case
 
      albc = 0.219/(1.0+0.816*cosz)
 
c     cloudy sky case
 
      albc = albc*fclear + rcloud*(1.-fclear)
      albt=albc+(1.0-albc)*0.856*rsurf/(1.0-0.144*rsurf)
      if(lprint) write(6,3000) fclear,rcloud,albc,albt
 3000 format(' fclear,rcloud,albc,albt:',4f10.3)

      do  7  n=2,nzpp1
      pq1=1./(1.+ae1/(k3(n)*n2(n)))
      pq2=1.+ae2/(k6(n)*o2(n))
      haeff(n)=0.2*(2.39+1.63*(k1(n)*n2(n)+k2(n)*o2(n)*pq1)/
     *(k1(n)*n2(n)+k2(n)*o2(n))+0.98/pq2)
      hueff(n)=0.25*(3.02+0.98/pq2)
      a2=sigc*no3(ncl)
      a1=sigc*no3(n)
      a3=sigc*(no3(ncl)-no3(n))
      a3=a3*cosz
      attc=exp(-amin1(50.,a2))
      att=exp(-amin1(50.,a1))
c    chappuis bands
                heat(n,1)=o3(n)*fc*sigc* att *fday
               if(twilit.and.(cosz .gt.0.)) heat(n,1)=heat(n,1)+o3(n)
     1       *fc*sigc *(2.*albt*e2(a3,0)*attc*cosz)
               if(.not.twilit) heat(n,1)=heat(n,1)+o3(n)*fc*sigc
     1                *(2.*albt*e2(a3,0)*attc*cosz)*fday
      
c     heat(n,1)=o3(n)*fc*sigc*(att+2.*albt*e2(a3)*attc/secz)*fday
c   hartley region
      att=exp(-amin1(50.,sigha*no3(n)))
      heat(n,2)=fday*o3(n)*fha*sigha*att
      heat(n,2)=heat(n,2)*haeff(n)
c   huggins bands
      heat(n,3)=(1.-exp(-shu1*con11*no3(n)))*fa1/em1
     1 + (exp(-shu2*con22*no3(n))-exp(-shu2*con21*no3(n)))*fa2/em2
      heat(n,3)=heat(n,3)*o3(n)/no3(n)*fday*hueff(n)
c...oxygen absorption
c    herzberg continuum (ozone & oxygen)
      a1=amin1(50.,(6.6e-24*no2(n)+4.9e-18*no3(n)))
      heat(n,4)=(7.9e-21*o2(n)+5.88e-15*o3(n))*fday*exp(-a1)
c   schumann-runge continuum (net heating)
      a1=amin1(50.,2.9e-19*no2(n))
      a2=amin1(50.,1.7e-18*no2(n))
      a3=amin1(50.,1.15e-17*no2(n))
      a5=amin1(50.,sigsrc*no2(n))
      heat(n,5)=(o2(n)/no2(n))*(0.98*exp(-a1)-0.55*exp(-a2)-0.43*exp(-a3
     *))*fday
      heat(n,5)=heat(n,5)+o2(n)*fsrc*sigsrc*fday*exp(-a5)
c   schumann-runge bands (total heating)
      if (no2(n) .lt. 1.0e18) go to 21
      heat(n,6)=fday*o2(n)/(.143*no2(n)+9.64e8*sqrt(no2(n)))
      go to 7
   21 heat(n,6)=9.03e-19*fday*o2(n)
    7 continue
c...compute heating rates for o3 and o2
      do  5  n=2,nzpp1
c...8.604e-3 contains 1/cp and conversion to deg/day
    5 cnf(n)=8.604e-3/rhot(n)
      do  9  n=2,nzpp1
      do  9  i=1,nb
    9 heat(n,i)=heat(n,i)*cnf(n)
c...get heating due to co2
      do 11 n=1,nzpp1
      co2tmp(n)=co2mix(n)
   11 continue
      if(cosz.gt.0.0)
     *call co2wdy(tmp,ptoth,rhot,co2tmp,rhostd,z,nzpp1,fday,cosz,heat)
c
c...reduce heating due to o3,o2,and co2 below cloud top
      if(fcloud.lt.0.01) go to 85
      do 80 n=2,nc1
      do 80 i=1,nb1
   80 heat(n,i)=heat(n,i)*(1.-fcloud)
   85 continue
c
c...sum over all bands except h2o
      do  93  n=2,nzd1
      totht(n)=0.
      do  93  i=1,nb1
  93  totht(n)=totht(n)+heat(n,i)
c   section for printout
      if(.not.lprint) go to 2000
c      write(6,201)
c  201 format(10x,'press',10x,'haeff',10x,'hueff',/)
c      do 203 n=2,nzd1
c  203 write(6,202) pmb(n),haeff(n),hueff(n)
c  202 format(e11.5,2e15.4)
      do 200 n=2,nzd1
      heato3(n)=heat(n,1)+heat(n,2)+heat(n,3)
      heato2(n)=heat(n,4)+heat(n,5)+heat(n,6)
  200 continue
      write(6,100)
  100 format('0','daily solar heating rates,band-by-band and total for'
     x ,' ozone and oxygen, in k/day  ')
      write(6,101)
  101 format ('0',9x,'z',13x,'press',8x,'cha',8x,'ha',8x,'hu',4x,
     *'o3tot',7x,'hzcont',7x,'src',7x,'srb',7x,'o2tot')
      do 96 n=2,nzd1
   96 write(6,97) n,z(n),pmb(n),(heat(n,i),i=1,3),heato3(n),
     * (heat(n,i),i=4,6),heato2(n)
   97 format(' ',i3,e11.4,2x,e11.5,2x,9e11.4)
      write(6,102)
  102 format('1',/,' ','daily solar heating rates, by species and total'
     x ,'in k/day ',/)
      write(6,103)
  103 format(' ',9x,'z',11x,'press',8x,'o3tot',7x,'o2tot',7x,'co2',7x,
     * '  sum ')
      do 104 n=2,nzd1
      write(6,97)n,z(n),pmb(n),heato3(n),heato2(n),heat(n,7),
     * totht(n)
  104 continue
 2000 continue
c..  end of printout section
c...reverse heating array to pass back to gcm
cjer solar flux stuff 
      do 1013 n=1,nzdum
      flux(n,8)=as(n)
       nn=nzdum-n+2
       as(n)=totht(nn)  +  as(n)
      do 1010 i=1,nb1
      flux(n,i)=heat(nn,i)
 1010 continue
 1013 continue
      if(twilit .and. cosz.le.0.0) go to 6000
      fac=3.6*24./10030.*.98
      do 1016 i=1,nbt
      flxsum(i)=0.0
      do 1015 n=1,nzdum
      flux(n,i)=flux(n,i)/fac*(ple(n+1)-ple(n))
 1015 flxsum(i)=flxsum(i) + flux(n,i)
 1016 continue
      solflx(1)=1365.*fday
      do 1017 n=2,nzd1
      solflx(n)=solflx(n-1)
      do 1017 i=1,nbt
      solflx(n)=solflx(n) - flux(n-1,i)*1.e-3
 1017 continue
c      if(.not.lprint) go to 1056
c      write(6,1020)
c 1020 format(//,' flux differences',/,
c     * 22x,'cha',9x,'har',9x,'hug',9x,'hzc',9x,'src',9x,'srb',9x,'co2',
c     * 9x,'h2o',/)
c      do 1050 n=1,nzdum
c      write(6,1025) n,pl(n),(flux(n,i),i=1,nbt)
c 1025 format(i4,f12.4,8f12.2)
c 1050 continue
c      write(6,1055)  flxsum
c 1055 format(/,' sum',12x,8f12.2)
c      write(6,1100) 
c 1100 format(//,' net solar flux down (W/M**2)',/)
c      do 1105 n=1,nzd1
c      write(6,1025) n,ple(n),solflx(n)
c 1105 continue


 1056 continue
c      s0=1365.
c      sgcosz=s0*fday*cosz
c      if(.not.lprint) go to 1060
c      write(6,1057) s0
c 1057 format(/,' incident solar flux, normal to beam  (w/m**2):',f10.3)
c      write(6,1058) sgcosz
c 1058 format(/,' diurnal average solar flux, normal to surface ',f10.3)
c 1060 continue
c...convert surface flux in h2o regions to w/m**2
c      sg1=sg*.48449074
c      if(lprint) write(6,1070) sg1
c 1070 format(/,' surface flux in h2o regions:',f10.3)
c      flxco2=flxsum(7)*1.e-3
c      if(lprint) write(6,1071) flxco2
c 1071 format(/,' flux absorbed by co2:',f10.3)
c      sg1=sg1 - flxco2*(1.-rsurf)
c      if(lprint) write(6,1072) sg1
c 1072 format(/,' surface flux in h2o regions with co2 abs:',f10.3)
c..surface flux absorbed in visible and uv
c      flxo3 = 0.0
c      do 1080 n=1,6
c 1080 flxo3 = flxo3 + flxsum(n)*1.e-3
c      if(lprint) write(6,1085) flxo3
c 1085 format(/,' flux abs.by atmosphere in o3 and o2 regions:',f10.3)
c..clear sky
c...lacis & hansen albedo
c      rclr=0.28/(1.+6.43*cosz)
c      sg2 = (0.647*sgcosz - flxo3 - sgcosz*rclr)*(1.-rsurf)/
c     * (1.-0.0685*rsurf)
c...cloudy sky
c      sg2=fclear*sg2 + (1.-fclear)*(0.647*sgcosz - flxo3)*
c     * (1.-rcloud)*(1.-rsurf)/(1.-0.144*rsurf)
c      if(lprint) write(6,1088) rclr,rcloud
c 1088 format(/,' rayleigh scattering albedo:',f10.3,/,
c     * ' cloud albedo:',f10.3)
c      if(lprint) write(6,1090) sg2
c 1090 format(/,' surface flux in uv & vis:',f10.3)
c      sgtot=sg1 + sg2
c      if(lprint) write(6,1095) sgtot
c 1095 format(/,' total surface flux:',f10.3)
c      sg=sgtot
 6000 continue
      return
      end
c
      subroutine co2wdy(t,ptot,prs,co2mix,rhostd,z,imax,fday,cosz,
     * heat)
c     calculates daily diurnally averaged
c     solar heating by co2
      real t(48),ptot(48),prs(48),co2mix(48),z(48),
     *heat(48,8),fco2(48)
      data cmix/1.0/
      do  1  n=1,imax
    1 fco2(n)=co2mix(n)*rhostd/prs(n)
      ia=imax-1
      do 1000 i=1,ia
      ic=imax-i
      dely=z(ic+1)-z(ic)
      p=prs(ic)*2.87*t(ic)
      sw=sqrt(293.0/t(ic))
 1000 continue
      do  1001  ia=2,imax
      a1=prs(ia)*2.87*t(ia)
      a1=(a1**0.84)
      azen=35.0/sqrt((1224.0*(cosz**2.0))+1.)
      zen=1./azen
      a1=a1*sqrt(cmix)
      sco2=zen*sqrt(azen)*((.396/(1.0+ (84.5*a1*sqrt(azen))))+(.167/(1.0
     1+(9.0*a1*sqrt(azen)))))*fday
      sco2=sco2+(zen*sqrt(azen)*fday*0.04/(1.0+(3.6*a1*sqrt(azen))))
      sco2=sco2*(t(ia)/220.0)
      sco2=sco2*sqrt(cmix)
      sco2=sco2*fco2(ia)/3.3e-4
      heat(ia,7)=sco2
 1001 continue
      return
      end
 
      function e2(q,kkk)
      real t,a(5),a1,a2,b1,b2,e1,at
      data a/.99999193,-.24991055,.05519968,-.00976004,.00107857
     *  /
      data a1/2.334733  /,a2/.250621  /,b1/3.330657  /,b2/1.681534  /

c add kkk to change storage somewhat 
      
      t=q
      if (t .ge. 20.  ) go to 6
      if (t .le. 5.0e-5) go to 5
      if (t .ge. 1.0  ) go to 1
      e1=0.0
      do  2  i=1,5
  2   e1=e1+a(i)*t**i
      e1=e1-alog(t)-.57721566
      go to 4
  1   at=(t**2+a1*t+a2)/(t**2+b1*t+b2)
      e1=at/(t*exp(t))
  4   e2=exp(-t)-t*e1
      if(kkk.eq.1) print *,'e2 = ',e2
      return
  5   e2=1.0
      return
  6   e2=0.
      return
      end
 
      function e3(q)
      t=q
      if (t .le. 5.0d-5) go to 5
      e3=0.5*(exp(-amin1(50.,t))-t*e2(q,0))
      return
  5   e3=0.5
      return
      end
c
      subroutine solar5 (nlay,jlat,fday,cosz,lprint)
 
c  compute solar heating due to water vapor
c  lacis and hansen parameterization
      include 'aero.inc'
      include 'PARAM.INC'
      include 'COMMONR1.INC'
      
ccelf      include 'test_value.inc'
ccelf      logical nantest
      
      common/hght/fcloud,fclear,ntopf,
     1 ntopt,rhostd,dlat(n$),long,rlat(n$)
 
      common/cirrus/icir(45),hcir(45),sigcir(45),taucir(45),epsa(45),
     * eps11(45),tauvis(45),prmax

      common/pscld/ipsc(45),pscabs(45),pscext(45),psch(45),psc96(45),
     * psc15(45)

c  radiation and source term fields
 
      logical lprint
      dimension heat(45)
      data gt0p0/120.1612/, delta/.0001/
c fraction of solar flux absorbed by an amount of h2o (cm)
      awater(x) = 2.9*x/((1. + 141.51*x)**.635 + 5.925*x)
c * * *
c
      if(lprint) write(6,40)
   40 format(' in solar5')
c...don't do h2o solar computation above 1 mb
      nlay1 = nlay+1
      ilay=2
      do 41 n=1,nlay
      al(n)=0.0
      swil(n)=0.0
      swale(n)=0.0
   41 continue
      al(nlay1)=0.0
      swale(nlay1)=0.0
      do 50 n=1,nlay
      ilay=n
      if(pl(n).gt.1.) go to 55
   50 continue
   55 continue
      ilay1=ilay+1
      rcloud=0.0
      rmean=0.0
c l&h magnification factor
      cosmag = 35.0/sqrt(1224.0*cosz**2 + 1.0)
c * * *
c  partition of incident flux subject to scattering
      s0=1365./.48449074
      s0=s0*fday
      sgcosz= s0*cosz
      scosz=s0/cosmag

c      IF (SCOSZ .NE. SCOSZ) THEN
c       nantest=test_nan( scosz )
c       if(nantest) then
c         PRINT*, "#4 SCOSZ NAN in SOLAR5:"
c         PRINT*, "     FDAY = ", FDAY
c         PRINT*, "       S0 = ", S0
c         PRINT*, "   COSMAG = ", COSMAG
c         PRINT*, "     COSZ = ", COSZ
c         STOP
c         ENDIF

c
c  scaled water vapor content above each layer edge
      db = ple(ilay)**2
      swale(ilay) = db*shl(ilay-1)/(sqrt(tl(ilay-1))*gt0p0)
c     write(6,105) swale(ilay)
  105 format(' swale(ilay):',e15.3)
               do 120 n=ilay,nlay
      m = n + 1
      da = ple(m)**2
      w = (da - db)*shl(n)/(sqrt(tl(n))*gt0p0)
      swil(n) = w
      swale(m) = swale(n) + w
      db = da
  120 continue
c      write(6,122)
c  122 format(' shl, scaled water vapor in layer and column')
c      do 125 n=1,nlay
c      write(6,124) n,pl(n),shl(n),swil(n),swale(n+1)
c  124 format(i4,f10.2,3e15.3)
c  125 continue
c 
      do 130 n=ilay,nlay
  130 taul(n)=0.0

      if(fclear.lt.0.01) go to 250
c
c  absorption by water vapor in clear atmosphere
c
c  l & h eqs. 25-27
c     write(6,202)
  202 format(' direct beam - pl,w,da,db,al')
c
      w = swale(ilay)*cosmag
      db = awater(w)
               do 210 n=ilay,nlay
      w = swale(n+1)*cosmag
      da = awater(w)
c  al is the fraction of solar flux absorbed in a layer
      al(n) = da - db
c     write(6,205) n,pl(n),w,da,db,al(n)
  205 format(i4,f10.2,4e15.3)
      db = da
  210          continue
c fraction of solar flux transmitted
      trans = 1.0 - db
c  fraction of solar flux reflected
      rf = trans*rsurf
c     write(6,212) trans,rsurf,rf
  212 format(' trans, rsurf, rf:',3e15.3)
      al(nlay1) = (trans - 0.647)*(1.0 - rsurf)
c            =(.353 - db)*(1-rsurf)
c            = fraction of total solar flux absorbed by ground in
c                    h2o regions
               if (rf.lt.0.001)  go to 230
c
c  now do absorption of reflected radiation
c     write(6,213)
  213 format(/,' diffuse part - pl,w,da,db,diff,al')
c
      ww = w*(1.0 + 1.66/cosmag)
               do 220 n=1,nlay
      m = nlay1 - n
      w = ww - 1.66*swale(m)
      da = awater(w)
      diff=(da-db)*rf
      al(m) = al(m) + (da-db)*rf
c     write(6,215) m,pl(m),w,da,db,diff,al(m)
  215 format(i4,f10.2,5e15.3)
      db=da
  220          continue
  230          continue
      aclear = fclear*scosz
               do 240 n=ilay,nlay
      as(n) = aclear*al(n)
  240          continue
      sg = fclear*sgcosz*al(nlay1)
c 
          if (fclear.gt.0.99)  go to 300
c 
c  absorption by water vapor in cloudy atmosphere
 250  continue
c
c  get scattering optical depth of layers
c
      do 150 n=ilay,nlay
      taul(n) = 0.0
      if( cloud(n).gt.0.01 )  then 
        taul(n)=tcond(n)
        if(icir(n).eq.1) taul(n)=tauvis(n)
        if(ipsc(n).eq.1) taul(n)=pscext(n)
      endif
  150 continue

      if(lprint) write(6,201) taul
  201 format(' taul in solar5',/, 6(/,8e9.2))

c  if fractional cloud layers, do this first
c  otherwise this is total cloudiness
      ntop = min0(ntopt,ntopf)
      if(lprint) write(6,2050) ntop
 2050 format(' ntop = ',i4)

      call cloud5 (ilay,nlay,ntop,cosz,lprint)

c  cloud5 returns fractional absorption in layer in al array
c  and visual cloud albedo as rcloud

      acloud = fcloud*scosz
      do 260 n=ilay,nlay
      as(n) = as(n) + acloud*al(n)

  260 continue
      sg=sg+fcloud*sgcosz*al(nlay1)
      rmean=fcloud*rcloud

      if(.not.lprint) go to 257
      write(6,253)
  253 format(/,' after 1st call to cloud5')
      do 256 m=1,nlay1 
  256 write(6,255) m,pl(m),al(m),as(m)
  255 format(i4,f10.2,5e15.3)
      write(6,1010) fcloud,fclear,ntopf,ntopt,rsurf,rcloud
 1010 format(//,' fcloud = ',f3.1,/,' fclear = ',f3.1,/,
     * ' ntopf = ',i4,/,' ntopt = ',i4,/,' rsurf = ',f10.2,/,
     * ' rcloud = ',f10.2)
  257 continue

      if (fcloud.gt.0.99)  go to 300

c  if did fractional cloud layers, now do part of atmosphere
c  with total cloudiness

      ntop = ntopt
      if(lprint) write(6,2050) ntop      

      if(ntop.eq.nlay1) go to 300

      fcld1 = 1.0 - fcloud
      do 270 n=ilay,nlay
      if (cloud(n).lt.0.99)  taul(n)=0.0
  270 continue

      call cloud5 (ilay,nlay,ntop,cosz,lprint)

      acloud = fcld1*scosz

      do 280 n=ilay,nlay
      as(n) = as(n) + acloud*al(n)


  280 continue

      if(lprint) then
        write(6,271)
  271 format(/,' after 2nd call to cloud5')
        do 276 m=1,nlay1
  276   write(6,255) m,pl(m),al(m),as(m)
      endif

c  sg=sg+fcloud*(fscat*(1.0-rcloud)+sgcosz*al(nlay1))
      sg=sg+fcld1*sgcosz*al(nlay1)
      rmean=rmean + fcld1*rcloud
      if(lprint) write(6,1010) fcloud,fclear,ntopf,ntopt,rsurf,rcloud
c  return mean value of cloud albedo for o3 chappuis band absorption
      rcloud = rmean

  300 continue
c
c  convert as array to heating rate
      fac=3.6*24./10030.*.98
      do 345 n=ilay,nlay
c convert flux to ergs/cm**2/sec
      fxnet=as(n)*.48449074e3
      heat(n)= fxnet*fac/dpl(n)

  345 continue
      if(lprint) write(6,370)
  370 format(' heating rate: h2o',/)
      do 380 n=ilay,nlay
      as(n)=heat(n)

      if(lprint) write(6,375) n,pl(n),as(n)
  375 format(i4,f12.2,f12.3)
  380 continue
      return
      end
c
      subroutine solaer (nlay,jlat,fday,cosz,lprint)
c
c  compute solar heating due to water vapor and aerosols 
c  lacis and hansen parameterization
c  no sulfate aerosol with a cloud
c
      include 'PARAM.INC'
      include 'COMMONR1.INC'
      include 'aero.inc'
c      include 'psc.inc'
      real lat
      common/hght/fcloud,fclear,ntopf,
     * ntopt,rhostd,dlat(n$),long,rlat(n$)
c 
      common/cirrus/icir(45),hcir(45),sigcir(45),taucir(45),
     * epsa(45),eps11(45),tauvis(45),prmax

      common/pscld/ipsc(45),pscabs(45),pscext(45),psch(45),psc96(45),
     * psc15(45)

      logical lprint
      dimension heat(45)
      data gt0p0/120.1612/, delta/.0001/
c fraction of solar flux absorbed by an amount of h2o (cm)
ccelf      awater(x) = 2.9*x/((1. + 141.51*x)**.635 + 5.925*x)
c * * *
      if(lprint) write(6,8)
    8 format(' in solaer')
      nlay1 = nlay+1
      ilay=2
      do 41 n=1,nlay
      
      al(n)=0.0
      swil(n)=0.0
      swale(n)=0.0
   41 continue
      al(nlay1)=0.0
      swale(nlay1)=0.0
      ilay1=ilay+1
      rcloud=0.0
      rmean=0.0
c l&h magnification factor
      cosmag = 35.0/sqrt(1224.0*cosz*cosz + 1.0)
c * * *
c  partition of incident flux subject to scattering
      s0=1365./.48449074
      s0=s0*fday
      sgcosz= s0*cosz
      scosz=s0/cosmag
c
c  scaled water vapor content above each layer edge
      db = ple(ilay)*ple(ilay)
      swale(ilay) = db*shl(ilay-1)/(sqrt(tl(ilay-1))*gt0p0)
c      write(6,105) swale(ilay)
c  105 format(' swale(ilay):',e15.3)
               do 120 n=ilay,nlay
      m = n + 1
      da = ple(m)*ple(m)
      w = (da - db)*shl(n)/(sqrt(tl(n))*gt0p0)
      swil(n) = w
      swale(m) = swale(n) + w
      db = da
  120          continue
c      write(6,122)
c  122 format(' shl, scaled water vapor in layer and column')
c      do 125 n=ilay,nlay
c      write(6,124) n,pl(n),shl(n),swil(n),swale(n+1)
c  124 format(i4,f10.2,3e15.3)
c  125 continue

      do 130 n=ilay,nlay
  130 taul(n)=0.0
    
c  get scattering optical depth due to aerosol
      do 140 n=ilay,nlay
      if( iaer(n).eq.1 ) taul(n)=aerext(n)
  140 continue

      if (fclear.lt.0.01)  go to 250
c 
c 
c  absorption in atmosphere with no clouds
c  also, no sulfate
 
      ntop = ntaer
      if(lprint) write(6,205) ntop
  205 format(' ntop = ',i4)
      
      call cloud5 (ilay,nlay,ntop,cosz,lprint)

c  clouds returns fractional absorption in layer in al array
c  and visual cloud albedo as rcloud

      aclear = fclear*scosz
      do 260 n=ilay,nlay
      as(n) = as(n) + aclear*al(n)
  260 continue
c  sg=sg+fclear*(fscat*(1.0-rcloud)+sgcosz*al(nlay1))
      sg=sg+fclear*sgcosz*al(nlay1)
      rmean=fclear*rcloud

      if(.not.lprint) go to 257
      write(6,253)
  253 format(/,' after 1st call to cloud5 in solaer if no clouds')
c      do 256 m=1,nlay1 
c  256 write(6,255) m,pl(m),al(m),as(m)
c  255 format(i4,f10.2,5e15.3)
      write(6,1010) fcloud,fclear,ntopf,ntopt,rsurf,rcloud
 1010 format(//,' fcloud = ',f3.1,/,' fclear = ',f3.1,/,
     * ' ntopf = ',i4,/,' ntopt = ',i4,/,' rsurf = ',f10.2,/,
     * ' rcloud = ',f10.2)
  257 continue

      if(fclear .gt. 0.99) go to 300

  250 continue

c  absorption by water vapor in cloudy atmosphere

c  if fractional cloud layers, do this first
c  otherwise this is total cloudiness
c
c  set optical depth of layers
c
      do 150 n=ilay,nlay
      if( cloud(n).gt.0.01 )  then 
        taul(n)=tcond(n)
        if(icir(n).eq.1) taul(n)=tauvis(n)
        if(iaer(n).eq.1) taul(n)=aerext(n)
        if(ipsc(n).eq.1) taul(n)=pscext(n)
      endif
  150 continue

      ntop=ntaer
      if(lprint) write(6,151) ntop
  151 format(' ntop = ',i4)

      call cloud5 (ilay,nlay,ntop,cosz,lprint)

c  clouds returns fractional absorption in layer in al array
c  and visual cloud albedo as rcloud

      acloud = fcloud*scosz
      do 2600 n=ilay,nlay
      as(n) = as(n) + acloud*al(n)
 2600 continue
c  sg=sg+fcloud*(fscat*(1.0-rcloud)+sgcosz*al(nlay1))
      sg=sg+fcloud*sgcosz*al(nlay1)
      rmean=fcloud*rcloud

      if(.not.lprint) go to 258
      write(6,2530)
 2530 format(/,' after first call to cloud5 for clouds')
c      do 2570 m=1,nlay1 
c 2570 write(6,255) m,pl(m),al(m),as(m)
c  255 format(i4,f10.2,5e15.3)
      write(6,1010) fcloud,fclear,ntopf,ntopt,rsurf,rcloud
  258 continue

      if (fcloud.gt.0.99)  go to 300

c  if did fractional cloud layers, now do part of atmosphere
c  with total cloudiness

      if(ntopt .eq. nlay1) go to 300

      ntop=ntaer

      if(lprint) write(6,301) ntop
  301 format(' ntop = ',i4)

      fcld1 = 1.0 - fcloud
      do 270 n=ilay,nlay
      if (cloud(n).lt.0.99) then  
        taul(n)=0.0
        if(iaer(n) .eq. 1) taul(n)=aerext(n)
      endif
  270 continue

      call cloud5 (ilay,nlay,ntop,cosz,lprint)

      acloud = fcld1*scosz
      do 280 n=ilay,nlay
      as(n) = as(n) + acloud*al(n)
  280 continue

      if(lprint) then
        write(6,271) fcld1
  271 format(/,' after 2nd call to cloud5 for clouds',/,
     * ' fcld1 = ',f5.2)
        do 276 m=1,nlay1
  276   write(6,255) m,pl(m),al(m),as(m)
  255 format(i4,f10.2,5e15.3)
      endif

c  sg=sg+fcloud*(fscat*(1.0-rcloud)+sgcosz*al(nlay1))
      sg=sg+fcld1*sgcosz*al(nlay1)
      rmean=rmean + fcld1*rcloud
      if(lprint) write(6,1010) fcloud,fclear,ntopf,ntopt,rsurf,
     * rcloud
c  return mean value of cloud albedo for o3 chappuis band absorption
      rcloud = rmean

  300 continue
c 
c  convert as array to heating rate
      fac=3.6*24./10030.*.98
      do 345 n=ilay,nlay
c convert flux to ergs/cm**2/sec
      fxnet=as(n)*.48449074e3
  345 heat(n)= fxnet*fac/dpl(n)
      if(lprint) write(6,370)
  370 format(' heating rate: h2o',/)
      do 380 n=ilay,nlay
      as(n)=heat(n)
      if(lprint) write(6,375) n,pl(n),as(n)
  375 format(i4,f12.2,f12.3)
  380 continue

      return
      end
c
      subroutine cloud5 (ilay,nlay,ntop,zen,lprint)
      include 'aero.inc'
      include 'PARAM.INC'
      include 'COMMONR1.INC'
      common/hght/fcloud,fclear,ntopf,
     1 ntopt,rhostd,dlat(n$),long,rlat(n$)
 
      common/cirrus/icir(45),hcir(45),sigcir(45),taucir(45),epsa(45),
     * eps11(45),tauvis(45),prmax
 
      common/pscld/ipsc(45),pscabs(45),pscext(45),psch(45),psc96(45),
     * psc15(45)

      logical land, ocean, ice, snow, mixwi, frost
      logical lprint
c  radiation and source term fields
 
      awater(x) = 2.9*x/((1. + 141.51*x)**.635 + 5.925*x)
      rtop(tau,pi0,cosz) = 1. + .00001*tau - .0001*pi0 + .00001/cosz
      ttop(tau,pi0,cosz) = 1. - .00001*tau + .0001*pi0 - .00001/cosz
 
      if(lprint) write(6,10) ntop
   10 format(' in cloud5 with ntop = ',i4)

      nlay1 = nlay + 1
      nntop=ntop
      cosmag = 35.0/sqrt(1224.0*zen**2 + 1.)
      cosz=1./cosmag
 
c  absorption of incident flux above clouds
 
      nclear = ntop - 1
      w = swale(ilay)*cosmag
      db = awater(w)
      do 110 n=ilay,nclear
      w = swale(n+1)*cosmag
      da = awater(w)
      al(n) = da - db
      db = da
  110 continue
      do 120 n=nntop,nlay1
      al(n) = 0.0
  120 continue
c * * *
c  reflectivity of cloudy atmosphere for visual light
      tau = 0.0
      do 130 n=ilay,nlay
      tau = tau + taul(n)
  130 continue
      pi0 = .99999
      taup = .212132e-04*tau
      etp = exp(taup)
      d = .1500245e+09*etp - .1499755e+09/etp
      rnn = .1500e+09*(etp - 1./etp)/d
      srnn = (1. - rnn)*rsurf/(1. - rnn*rsurf)
      rnn = rnn*rtop(tau,pi0,cosz)
      rcloud = rnn + (1.0 - rnn)*srnn
               if (db.gt.0.999)  return
c * * *
c  absorption in clouds with multiple reflections for k-distribution

      do 250 k=1,nfk

      fkk = fk(k)
      xkk = xk(k)
      wk = w*xkk
               if (wk.gt.7.0)  return
      sk = exp(-wk)
      sb = exp(-1.66*wk)
 
c  sum over layers
      do 180 n=nntop,nlay
      if(icir(n).eq.1) go to 140
      if(ipsc(n).eq.1) go to 142
      if(iaer(n).eq.1) go to 143
c   liquid water clouds
      tausc = taul(n)
      tauab = amin1(swil(n)*xkk,20.)
      if (tausc.lt.0.01)  go to 160
      tau = tausc + tauab
c pi0 is single scattering albedo
      pi0 = amin1(tausc/tau,.9999)
      go to 145
  140 continue
c  ice water clouds
      pi0=0.95
      tau = amin1(swil(n)*xkk,20.) + taul(n)
      if(iaer(n).eq.1) tau = tau + aerext(n)
      tausc = pi0*tau
      tauab = tau - tausc
      go to 145
  142 continue
c polar stratospheric clouds
      tau = amin1(swil(n)*xkk,20.) + pscext(n)
      tausc = pscext(n) - pscabs(n)
      tauab = pscabs(n)
      pi0=tausc/tau
      go to 145      
  143 continue
c aerosols only
      tau = amin1(swil(n)*xkk,20.) + taul(n)
      tausc = taul(n) - aerabs(n)
      tauab = aerabs(n)
      pi0 = tausc/tau     
  145 continue

c      if(lprint) write(6,12) k,n,tausc,tauab,pi0
c   12 format(' k = ',i4,/,' n = ',i4,' tausc, tauab, pi0 = ',2e15.6,
c     * f10.4)
c      if(tauab .lt. 0.0) then
c       print *,' tauab < 0.0'
c       print *,' icir,ipsc,cloud,aerfl=',icir(n),ipsc(n),cloud(n),aerfl
c       print *,' iaer = ',iaer(n)
c       stop10 
c      endif


c  check for clear layer
cjer      if(pi0.lt.0.1) go to 160
      if(pi0.lt.1.e-5) go to 160
c * * *
c  individual cloud layer reflectivity and transmissivity
      u = sqrt((1. - 0.850*pi0)/(1. - pi0))
      taup = 1.732051*u*(1.0 - pi0)*tau
      etp = exp(taup)
      d = (u + 1.)**2*etp - (u - 1.)**2/etp
      rnn = (u**2 - 1.)*(etp - 1./etp)/d
      tnn = 4.*u/d
      rn(n) = rnn
      tn(n) = tnn
               if (n.eq.ntop)  go to 150
c * * *
c  sum reflectivity and transmissivity for top,bottom illumination
      denom = 1.0 - srsn*rnn
      srnn = srnn + stnn*rnn*stsn/denom
      stnn = stnn*tnn/denom
      srsn = rnn + tnn**2*srsn/denom
      stsn = tnn*stsn/denom
      tfk = fkk*stnn
               go to 170
c * * *
c  top cloud zenith angle dependent reflectivity and transmissivity
  150          continue
      rnk = rnn*rtop(tau,pi0,cosz)
      tnk = tnn*ttop(tau,pi0,cosz)
      srnn = sk*rnk*sb
      stnn = sk*tnk
      srsn = rnn
      stsn = tnn*sb
      tfk = fkk*stnn
               go to 170
c * * *
c  clear layer diffuse transmission
  160          continue
      sb = exp(-1.66*tauab)
      stnn = stnn*sb
      srsn = srsn*sb**2
      stsn = stsn*sb
      rn(n) = 0.
      tn(n) = sb
      tfk = fkk*stnn
  170          continue
               if (tfk.lt..001)  go to 190
      stn(n) = stnn
      srs(n) = srsn
  180          continue
c * * *
c  absorption at ground
      rnn = rsurf
      denom = 1. - srsn*rnn
      da = stnn*(1. - rnn)/denom*fkk
c
c
      al(n) = al(n) + da
      db = da
               go to 200
c * * *
c  distribution of absorption among individual cloud layers
  190          continue
      m = n - 1
               if (m.lt.ntop)  go to 220
      srsn = srs(m)
      stnn = stn(m)
      denom = 1. - srsn*rnn
      da = stnn*(1. - rnn)/denom*fkk
      al(n) = al(n) + da
      db = da
  200          continue
      n = n - 1
      m = n - 1
      srsn = rn(n)
      stnn = tn(n)
      denom = 1. - srsn*rnn
               if (m.lt.ntop)  go to 215
      n1 = m
               do 210 nn=nntop,n1
      rnn = srsn + stnn**2*rnn/denom
      srsn = srs(m)
      stnn = stn(m)
      denom = 1. - srsn*rnn
      da = stnn*(1. - rnn)/denom*fkk
c don't allow negative heating rate
      if(da.gt.db) al(n) = al(n) + da - db
c      if(lprint) write(6,198) nn,m,n,al(n),da,db
  198 format(' nn,m,n,al(n),da,db:',3i4,2x,e10.4,2e13.6)
      db = da
      n = n - 1
      m = n - 1
      srsn = rn(n)
      stnn = tn(n)
      denom = 1. - srsn*rnn
  210          continue
  215          continue
      rnn = rnk + tnk*stnn*rnn/denom
      da = sk*(1. - rnn)*fkk
c don't allow negative heating rate
      if(da .gt. db) al(n) = al(n) + da - db
c      if(lprint) write(6,198) nn,m,n,al(n),da,db
      db = da
               go to 230
c * * *
c  absorption of reflected flux above clouds
  220          continue
      da = sk*(1. - rnk)*fkk
      al(n) = al(n) + da
      rnn = rnk
  230          continue
      fkk = sk*rnn*fkk
               if (fkk.lt..001)  go to 250
      db = 0.
      ilay1=ilay+1
               do 240 nn=ilay1,nntop
      n = nntop + ilay - nn
      wk = 1.66*(w - swale(n))*xkk
      da = (1. - exp(-wk))*fkk
      al(n) = al(n) + da - db
      db = da
  240          continue

c end of k loop
  250          continue
               return
      end
c
      subroutine fit(pr,pr1,pr2,v,v1,v2)
c...linear in log p interpolation
      plog=alog(pr)
      plog1=alog(pr1)
      plog2=alog(pr2)
      v=v2 + (v1 - v2)*(plog - plog2)/(plog1 - plog2)
      return
      end
c
      subroutine irdriv(nvert,lprint)
ccelf      include 'test_value.inc'
ccelf      logical nantest
      logical flag,flag1
      logical lprint,co2var
      include 'aero.inc'
      include 'COMMONR1.INC'      
      real solflx(46),irflx(46)
      common /radflx/  ntrop, solflx, irflx
      common /temp/ ts
      common /init/ nup1,nsb,nc,npt,nw,nt,im
      common /amhn/ pu(48), ta(49), pa(48), wa(48), ps
      common /agios/ tui(49),rawi(48),
     +               ubar(49),wbar(49,2)
      common/height/ha(49)
      common /cool/ coolr(48),sflux
      common /fluxh/ fluxhu(48),fluxhd(48)
      common /fluxo/ fluxou(48),fluxod(48)
      common/cloud1/cfrac(3,48),pline(3,48,48)
      common/cirrus/icir(45),hcir(45),sigcir(45),taucir(45),epsa(45),
     *  eps11(45),tauvis(45),prmax
      common/pscld/ipsc(45),pscabs(45),pscext(45),psch(45),psc96(45),
     * psc15(45)
      dimension fxwav(48),fx15c(48)
      common/co2mr/co2vmr(45),co2lw(48),co2var
      
c npt=48!

      do 5 n=1,npt
      do 5 ib=1,3
    5 cfrac(ib,n)=0.0
c  initialize variables from radcom
      do 1 n=3,npt
      pu(n)=ple(n-2)
      tui(n)=tle(n-2)
      if(n.eq.3) go to 4
      pa(n)=pl(n-3)
      ta(n)=tl(n-3)
      wa(n)=shl(n-3)
      ha(n)=hi(n-3)
      co2lw(n)=co2vmr(n-3)
      do 3 ib=1,3
      cfrac(ib,n)=cfrac(ib,n) + cloud(n-3)
    3 continue

c cirrus can have fractional cloudiness
c can have either cirrus, psc, or sulfate in a grid box, not both
c subvisible cirrus are treated like psc

      if(icir(n-3).eq.1) cfrac(1,n)=cfrac(1,n)*eps11(n-3)
      if(icir(n-3).eq.1) cfrac(2,n)=cfrac(2,n)*eps11(n-3)
      if(icir(n-3).eq.1) cfrac(3,n)=cfrac(3,n)*eps11(n-3)

      if(ipsc(n-3).eq.1) cfrac(1,n)=psch(n-3)
      if(ipsc(n-3).eq.1) cfrac(2,n)=psc96(n-3)
      if(ipsc(n-3).eq.1) cfrac(3,n)=psc15(n-3)

      if(iaer(n-3).eq.1) cfrac(1,n)=aerh(n-3)
      if(iaer(n-3).eq.1) cfrac(2,n)=aer96(n-3)
      if(iaer(n-3).eq.1) cfrac(3,n)=aer15(n-3)
      
    4 continue
    1 continue
      ts=tgjr
      ps=ple(nvert+1)
      delog=alog(pu(3))-alog(pu(4))
      pu(2)= exp(alog(pu(3))+delog)
      pu(1)= exp(alog(pu(3))+2.*delog)
      pa(3)=(pu(2)*pu(3))**0.5
      pa(2)=(pu(1)*pu(2))**0.5
      ta(3)=extrap(tl(1),pl(1),tl(2),pl(2),pa(3))
      ta(2)=extrap(tl(1),pl(1),tl(2),pl(2),pa(2))
      tui(2)=extrap(tl(1),pl(1),tl(2),pl(2),pu(2))
      tui(1)=extrap(tl(1),pl(1),tl(2),pl(2),pu(1))
      wa(3)=wa(4)
      wa(2)=wa(3)
      pa(1)=pa(2)
      ta(1)=ta(2)
      wa(1)=wa(2)
      co2lw(3)=co2lw(4)
      co2lw(2)=co2lw(3)
      co2lw(1)=co2lw(2)
      ha(3)=ha(4)+(ta(3)+ta(4))*.5/34.17*
     +  alog(pa(4)/pa(3))
      ha(2)=ha(3)+(ta(2)+ta(3))*.5/34.17*
     +  alog(pa(3)/pa(2))
   11 continue
      do 12 ib=1,3
      cfrac(ib,1)=0.
      cfrac(ib,2)=0.
      cfrac(ib,3)=0.
   12 continue
c
c...get probabilities of clear line of sight
      call cloudy
c
c...do h2o calculation first (need scaled h2o amounts, used in iro3)
      call gamoto
c
c...get co2 and o3 transmissions
      if( .not. co2var) call radco2(nvert,pa,ta,pu,lprint)
c      if(flag .and. flag1) print *,'before trco2v'
      if( co2var ) call trco2v(nvert,pa,ta,pu,lprint)
c      if( co2var ) call trco2v(nvert,pa,ta,pu,lprint,flag,flag1)
c      if(flag .and. flag1) print *,'before rado3'
      call rado3(nvert,ta,pa,pu,lprint)
c
c...compute o3 and co2 fluxes
c      if(flag .and. flag1) print *,'before iro3'
      call iro3
c
      if(.not.lprint) go to 100
      write(6,10)
   10 format(13x,'npt',3x,'nc',2x,'ts',6x,'ps')
      write(6,120) npt,nc,ts,ps
  120 format(13x,i2,4x,i2,1x,f5.1,2x,f6.1)
      write(6,14)
   14 format(/,' n        pa(n)   ta(n)     wa(n)     cfrac(n)')
      do 88 ip=1,npt
   88 write(6,15) ip,pa(ip),ta(ip),wa(ip),(cfrac(ib,ip),ib=1,3)
   15 format(1x,i2,2x,f10.4,1x,f8.2,1x,e10.3,1x,3e12.2)
  100 continue
      if(lprint) write(6,30)
   30 format(//,3x,'pressure',4x,'fluxes in h2o bands',3x,'fluxes in 9.
     $6 & 15 micron bands')
      if(lprint) write(6,40)
   40 format(5x,'(mb)',7x,'(ergs/cm**2/sec)',11x,'(ergs/cm**2/sec)')
      if(lprint) write(6,50)
   50 format(15x,'down',8x,'up',14x,'down',8x,'up',7x,'net h2o',4x,
     +           'net 9.6/15',2x,'net total')
      do 200 ip=1,npt
      ipm1=ip-1
         diff1=fluxhu(ip)-fluxhd(ip)
         diff2=fluxou(ip)-fluxod(ip)
         tot=diff1+diff2
          if(lprint)  write(6,201) pu(ip),fluxhd(ip),fluxhu(ip),
     +   fluxod(ip),fluxou(ip),diff1,diff2,tot
  200    continue
  201 format(1x,f9.4,2f11.2,4x,2f12.2,3f12.2)
      fac=3.6*24./10030.*.98
      do 69 ip=1,npt
      fxwav(ip)=fluxhd(ip)-fluxhu(ip)
c      nantest=test_nan(fxwav(ip))
c      if(nantest) then 
c        print *,'nan in fxwav in irdriv:'
c        print *,'ip,fluxhd,fluxhu = ',ip,fluxhd(ip),fluxhu(ip)
c      endif
      fx15c(ip)=fluxod(ip)-fluxou(ip)
c      nantest=test_nan(fx15c(ip))
c      if(nantest) then
c        print *,'nan in fx15c in irdriv:'
c	print *,'ip,fluxod,fluxou = ',ip,fluxod(ip),fluxou(ip)
c      endif
   69 continue
      tir=( fluxhu(1)+fluxou(1) )*1.e-3
      do 90 ip=2,npt
      ipm1=ip-1
      dp=pu(ip)-pu(ipm1)
      x1=fx15c(ip)+fxwav(ip)
      if(ip.gt.2) irflx(ip-2) = -x1*1.e-3
      x2=fx15c(ipm1)+fxwav(ipm1)
      fxnet=x1-x2
      coolr(ip)=fxnet*fac/dp
c      nantest=test_nan(coolr(ip))
c      if(nantest) then
c        print *,'nan in coolr in irdriv: ip,fxnet,fac,dp='
c	print *,ip,fxnet,fac,dp
c	print *,'x1,x2 = ',x1,x2
c	print *,'fx15c(ip),fx15c(ipm1) = ',fx15c(ip),fx15c(ipm1)
c	print *,'fxwav(ip),fxwav(ipm1) = ',fxwav(ip),fxwav(ipm1)
c      endif
      if(ip.eq.npt)sflux=x1

   90 continue

      if(.not.lprint) go to 95
      write(6,85)
   85 format(//,' net upward ir flux (w/m**2)')
      do 87 ip=1,nvert+1
      write(6,89) ip,ple(ip),irflx(ip)
   87 continue
   89 format(i4,4x,f10.4,5x,f10.2)
   
      write(6,6)
    6 format(//16x,'pressure',10x,'cooling rate')
      write(6,7)
    7 format(16x,'(mbar)',7x,'(degree celcius/day)')
      write(6,8)
    8 format(12x,'from',6x,'to')
      do 92 ip=2,npt
      ipm1=ip-1
      if(lprint)write(6,91) ha(ip),pu(ipm1),pu(ip),coolr(ip)
   92 continue
   91 format(1x,3f10.4,4x,f10.2,5x,f10.2,2x,f10.2,2x,f10.2)
   95 continue

      return
      end
c
      subroutine gamoto
      dimension pre(2),rx(2),tre(2),dpi(48)
      common/cloud1/cfrac(3,48),pline(3,48,48)
      common /init/ nup1,nsb,nc,npt,nw,nt,im
      common /agios/ tui(49),rawi(48),
     +               ubar(49),wbar(49,2)
      common /amhn/ pu(48), ta(49), pa(48), wa(48), ps
      dimension pdop(2)
      data pre/275.,550./
      data tre/225.,256./,rx/.005,.016/
      data pdop/1.0, 15.0/
c
      ndm1=npt-1
      npt1=npt-1
      npt2=npt1+1
      do 1010 ix=2,npt
      dpi(ix)=(pu(ix)-pu(ix-1))*1.02
 1010 continue
c
c****** temperature and humidity interpolations *************
c
      do 1030 ix=2,npt
        rawi(ix)=wa(ix)
      if (rawi(ix) .lt. 0.0) rawi(ix)=0.0
      rawi(ix)=rawi(ix)*dpi(ix)
 1030 continue
c
c**********  compute scaled water vapor amounts **************
c
      ubar(1)=0.
      wbar(1,1)=0.0
      wbar(1,2)=0.0
 1031 continue
      do 14 ix=2,npt
      do 12 kb=1,2
      dt=ta(ix)-tre(kb)
      sz2=exp(rx(kb)*dt)
      peff=pa(ix)/pre(kb)
      if(pa(ix).lt.pdop(kb)) peff=pdop(kb)/pre(kb)
      wbar(ix,kb)=wbar(ix-1,kb)+rawi(ix)*sz2*peff
   12 continue
      xx=(pa(ix)/630.)*rawi(ix)*rawi(ix)/dpi(ix)
      xx=xx*exp(1800./ta(ix)-6.0811)
      ubar(ix)=ubar(ix-1)+xx
   14 continue
c
c...compute h2o transmissions and planck functions
c...irh2o calls fdiver
      call irh2o
      return
      end
c
      subroutine irh2o
c  h2o transmissions (see m.-d. chou, j. atm. sci., 41, 1775-1778,1984;
c     harshvardhan and corsetti, nasa tech memo 86072, 1984)        
      dimension sw(2),ss(49,2),tui1(49),tui2(49)
      dimension ubarm(49),wbarm(49,2)
      dimension a(49),b(49),c(49),ww(49),wv(49)
      dimension yv(49),yw(49),zv(49),zw(49)
      common /agios/ tui(49),rawi(48),
     +               ubar(49),wbar(49,2)
      common /h2o/ gl(30,20),b250(41,2),coeff(30,20,2)
      common /temp/ ts
      common /init/ nup1,nsb,nc,npt,nw,nt,im
      common /fdiv/ sh(49),sg(48,48),sht(48),shu(48)
      data temp1/190./,sw/-6., -5.4/
      data dw/0.3/, du/0.003/,dt/5./
      nptp1=npt+1
      nptm1=npt-1
      tui(nptp1)=ts
      do ip=1,nptp1
      end do
      do 43 ip=1,nptp1
      tui1(ip)=tui(ip) - 250.
      tui2(ip)=tui1(ip)*tui1(ip)
      sh(ip)=0.
   43 continue
      do 44 ip=1,npt
      sht(ip)=0.0
      shu(ip)=0.0
      do 44 ix=1,npt
      sg(ix,ip)=0.
   44 continue
c   planck function for level temperature
      do 35 ix=1,nptp1
      fh=(tui(ix)-temp1)/dt+1.5
      it=fh
      it=max0(min0(it,nt-1),1)
      f1=real(it-1)
      dh=tui(ix)-(temp1+f1*dt)
c kb=1 refers to band center, kb=2 to band wing
      do 30 kb=1,2
      ss(ix,kb)=b250(it,kb) + (b250(it+1,kb)-b250(it,kb))*dh/dt
   30 continue
      sh(ix)=ss(ix,1) + ss(ix,2)
   35 continue

      fluxus=ss(nptp1,1)+ss(nptp1,2)

      do 45 ix=2,npt
      ubarm(ix)=0.5*(ubar(ix) + ubar(ix-1))
      do 45 kb=1,2
      wbarm(ix,kb)=0.5*(wbar(ix,kb) + wbar(ix-1,kb))
   45 continue
c
c    downward flux exchange term sg(ix,ip)
c
      do 420 ip=2,npt
      do 400 ix=2,ip
c   scaled h2o for e-type absorption in path
      a(ix-1)=abs(ubar(ip)-ubarm(ix))
c
c     band center region
c
      b(ix-1)= abs(wbar(ip,1)-wbarm(ix,1))
c
c... band wing region
c
      c(ix-1)= abs(wbar(ip,2)-wbarm(ix,2))
c
c...compute wv,ww,yv,yw,zv,zw from functions
c
  400 continue
      jp=ip-1

      call gfunc(jp,a,b,c,wv,ww,yv,yw,zv,zw)
c
c  put everything together to compute sg and sht
c    x is tau(w,u;t)*b
c   sg is d ( tau(w,u;t)*b  )
c  v refers to band center region, w to band wing region
      do 410 ix=2,ip
      xv=wv(ix-1)*(1.+yv(ix-1)*tui1(ix-1)+
     1 zv(ix-1)*tui2(ix-1))*ss(ix-1,1)
      xw=ww(ix-1)*(1.+yw(ix-1)*tui1(ix-1)+
     1 zw(ix-1)*tui2(ix-1))*ss(ix-1,2)
      sg(ix-1,ip)=sg(ix-1,ip)+xv + xw
c
      xv=wv(ix-1)*(1.+yv(ix-1)*tui1(ix)+    
     1 zv(ix-1)*tui2(ix))*ss(ix,1)      
      xw=ww(ix-1)*(1.+yw(ix-1)*tui1(ix)+
     1 zw(ix-1)*tui2(ix))*ss(ix,2)
      sg(ix-1,ip)=sg(ix-1,ip)-xv - xw
c
  410 continue
  420 continue
c
c   upward flux exchange terms sg(ix,ip)
c
      do 42 ip=1,nptm1
      ip1=npt-ip
      
      do 40 ix=1,ip1
c   scaled h2o for e-type absorption in path
      a(ix)= abs(ubar(ip)-ubarm(ix+ip))
      b(ix)= abs(wbar(ip,1)-wbarm(ix+ip,1))
      c(ix)= abs(wbar(ip,2)-wbarm(ix+ip,2))
   40 continue
      jp=ip1

      call gfunc(jp,a,b,c,wv,ww,yv,yw,zv,zw)

      do 41 ix=1,ip1
      xv=wv(ix)*(1.+yv(ix)*tui1(ix+ip)+
     1 zv(ix)*tui2(ix+ip))*ss(ix+ip,1)
      xw=ww(ix)*(1.+yw(ix)*tui1(ix+ip)+
     1 zw(ix)*tui2(ix+ip))*ss(ix+ip,2)
      sg(ix+ip,ip)=sg(ix+ip,ip)+xv + xw
c
      xv=wv(ix)*(1.+yv(ix)*tui1(ix+ip-1)+
     1 zv(ix)*tui2(ix+ip-1))*ss(ix+ip-1,1)
      xw=ww(ix)*(1.+yw(ix)*tui1(ix+ip-1)+
     1 zw(ix)*tui2(ix+ip-1))*ss(ix+ip-1,2)
      sg(ix+ip,ip)=sg(ix+ip,ip)-xv - xw
c
   41 continue
   42 continue
c
c   downward flux from top sht(ip)
c
      do 140 ip=2,npt
      a(ip-1)= abs(ubar(ip))
      b(ip-1)= abs(wbar(ip,1))
      c(ip-1)= abs(wbar(ip,2))
  140 continue
      jp=nptm1

      call gfunc(jp,a,b,c,wv,ww,yv,yw,zv,zw)

      do 145 ip=2,npt
      xv=wv(ip-1)*(1.+yv(ip-1)*tui1(1)+
     1 zv(ip-1)*tui2(1))*ss(1,1)
      xw=ww(ip-1)*(1.+yw(ip-1)*tui1(1)+
     1 zw(ip-1)*tui2(1))*ss(1,2)
      sht(ip)=sht(ip)+xv + xw
  145 continue
c
c   upward flux from bottom shu(ip)
c
      do 150 ip=1,nptm1
      a(ip)= abs(ubar(npt) - ubar(ip))
      b(ip)= abs(wbar(npt,1) - wbar(ip,1))
      c(ip)= abs(wbar(npt,2) - wbar(ip,2))
  150 continue
      jp=nptm1

      call gfunc(jp,a,b,c,wv,ww,yv,yw,zv,zw)

      do 155 ip=1,nptm1
      xv=wv(ip)*(1.+yv(ip)*tui1(nptp1)+
     1 zv(ip)*tui2(nptp1))*ss(nptp1,1)
      xw=ww(ip)*(1.+yw(ip)*tui1(nptp1)+
     1 zw(ip)*tui2(nptp1))*ss(nptp1,2)

      shu(ip)=shu(ip)+xv + xw

      xv=wv(ip)*(1.+yv(ip)*tui1(npt)+
     1 zv(ip)*tui2(npt))*ss(npt,1)
      xw=ww(ip)*(1.+yw(ip)*tui1(npt)+
     1 zw(ip)*tui2(npt))*ss(npt,2)

      shu(ip)=shu(ip) - xv - xw

  155 continue
c
      call fdiver
c
      return
      end
 
      subroutine gfunc(n,a,b,c,wv,ww,yv,yw,zv,zw)
      dimension a(49),b(49),c(49),a1(49),a2(49),a3(49)
      dimension wv(49),ww(49),yv(49),yw(49),zv(49),zw(49)
      dimension b1(49),b2(49),b3(49),b4(49)
      dimension c1(49),c2(49),c3(49),c4(49)
      do 20 j=1,n
      a1(j)=(1.0+32.2095*a(j))/(1.0+52.85*a(j))
      a2(j)=(0.534874+199.0*a(j)-1990.63*a(j)*a(j))/(1.0+333.244*a(j))
      a3(j)=(1.0+74.144*a(j))/(0.43368+24.7442*a(j))
   20 continue
      do 30 j=1,n
      ww(j)=(a1(j)+a2(j)*sqrt(c(j)))/(1.0+a3(j)*sqrt(c(j)))
      ww(j)=amax1(ww(j),0.0)
      wv(j)=1.0/(1.0+9.22411*sqrt(b(j))+33.1236*b(j)+176.396*b(j)*b(j))
      wv(j)=amax1(wv(j),0.0)
   30 continue
      do 40 j=1,n
      a(j)=amin1(a(j),0.06)
      b(j)=amin1(b(j),2.0)
      c(j)=amin1(c(j),8.0)
   40 continue
      do 50 j=1,n
      yv(j)=0.1*(0.0851069*sqrt(b(j))+0.323105*b(j)
     1 -0.187096*b(j)*sqrt(b(j)))
      zv(j)=0.001*(0.239186*b(j)-0.0922289*b(j)*sqrt(b(j))
     1 -0.0167413*b(j)*b(j))
   50 continue
      do 60 j=1,n
      b1(j)=(1.34927*a(j)-49.2303*a(j)*a(j))/(1.0+257.949*a(j))
      b2(j)=(0.0779555+4.40720*a(j)+3.15851*a(j)*a(j))/
     1 (1.0+40.2298*a(j))
      b3(j)=(-0.0381305-3.63684*a(j)+7.98951*a(j)*a(j))/
     1 (1.0+62.5692*a(j))
      b4(j)=(0.00621039+0.710061*a(j)-2.85241*a(j)*a(j))/
     1 (1.0+70.2912*a(j))
      yw(j)=0.1*(b1(j)+b2(j)*sqrt(c(j))+b3(j)*c(j)+
     1 b4(j)*c(j)*sqrt(c(j)))
   60 continue
      do 70 j=1,n
      c1(j)=(0.201563*a(j)-1.47769*a(j)*a(j))/(1.0-4.79279*a(j))
      c2(j)=(-0.0291325-2.30007*a(j)+10.9460*a(j)*a(j))/
     1 (1.0+63.5190*a(j))
      c3(j)=(0.0143812+1.80265*a(j)-10.1311*a(j)*a(j))/
     1 (1.0+98.4758*a(j))
      c4(j)=(-0.00239016-0.371427*a(j)+2.35443*a(j)*a(j))/
     1 (1.0+120.228*a(j))
      zw(j)=0.001*(c1(j)+c2(j)*sqrt(c(j))+c3(j)*c(j)+
     1 c4(j)*c(j)*sqrt(c(j)))
   70 continue
      return
      end

      subroutine fdiver
c
c   ...computes fluxes in h2o bands
c
      common /fluxh/ fluxu(48),fluxd(48)
      common/cloud1/cfrac(3,48),pline(3,48,48)
      common /init/ nup1,nsb,nc,npt,nw,nt,im
      common/fdiv/sh(49),sg(48,48),sht(48),shu(48)
      npm1=npt-1
      npp1=npt+1
c***  compute downward fluxes
      fluxd(1)=0.
      do 200 ip=2,npt
      fluxd(ip)=0.0
      ipm1=ip-1
      do 130 ix=1,ipm1
      fluxd(ip)=fluxd(ip)+sg(ix,ip)*pline(1,ix,ip)
  130 continue
      fluxd(ip)=fluxd(ip)-sht(ip)*pline(1,1,ip)+sh(ip)
  200 continue
c
c***  compute upward fluxes
c
      do 400 ip=1,npm1
      fluxu(ip)=0.0
      ipp1=ip+1
      do 300 ix=ipp1,npt
      fluxu(ip)=fluxu(ip)+sg(ix,ip)*pline(1,ix,ip)
  300 continue
      fluxu(ip)=fluxu(ip)+sh(ip)+shu(ip)*pline(1,ip,npt)
  400 continue
      fluxu(npt)=sh(npp1)
      
      return
      end
c
      function extrap(x1,p1,x2,p2,pi)
      slope=(x2-x1)/(alog(p2/p1))
      extrap=slope*alog(pi/p2) + x2
      return
      end
c
      subroutine iro3
      dimension bai(49),buic(49),buiw1(49),buiw2(49)
      dimension bwi1(49),bwi2(49),blkwin(60,2),blkco2(60)
      dimension tauco2(48,48)
      dimension taustd(49,49),fsn(19600)
      common /co2/ tauco2,taustd,fsn
      common/tro3/tauf1(48,48),tauf2(48,48)
      dimension fluxw1(48), fluxw2(48)
      dimension fluxdc(48),fluxuc(48)
      dimension fluxdw(48),fluxuw(48)
      dimension fluxdw1(48),fluxuw1(48)
      dimension fluxdw2(48),fluxuw2(48)
      common /temp/ ts
      common /fluxo/ fluxu(48),fluxd(48)
      common/cloud1/cfrac(3,48),pline(3,48,48)
      common /agios/ tui(49),rawi(48),ubar(49),wbar(49,2)
      common /init/ nup1,nsb,nc,npt,nw,nt,im
      common /amhn/ pu(48), ta(49), pa(48), wa(48), ps
      common/height/ha(49)
c
c   ...blkwin is integrated planck function for window region at
c   ....60 temperatures
      data blkwin/38.,  57.,    85.,   122.,   171.,   236.,
     *       318.,     421.,     548.,     704.,     891.,    1114.,
     *      1376.,    1682.,    2036.,    2441.,    2902.,    3423.,
     *      4007.,    4658.,    5380.,    6176.,    7050.,    8005.,
     *      9043.,   10168.,   11382.,   12688.,   14087.,   15583.,
     *     17176.,   18868.,   20662.,   22557.,   24556.,   26660.,
     *     28868.,   31183.,   33604.,   36131.,   38766.,   41508.,
     *     44357.,   47314.,   50377.,   53548.,   56824.,   60207.,
     *     63695.,   67288.,   70984.,   74785.,   78688.,   82693.,
     *     86799.,   91005.,   95310.,   99713.,  104214.,  108811.,
     *        15.,      23.,      34.,      50.,      70.,      96.,
     *       130.,     173.,     225.,     290.,     367.,     459.,
     *       568.,     694.,     841.,    1009.,    1200.,    1415.,
     *      1657.,    1927.,    2227.,    2557.,    2919.,    3315.,
     *      3745.,    4212.,    4715.,    5256.,    5836.,    6456.,
     *      7116.,    7817.,    8559.,    9344.,   10172.,   11043.,
     *     11957.,   12915.,   13916.,   14962.,   16052.,   17185.,
     *     18364.,   19586.,   20852.,   22162.,   23516.,   24913.,
     *     26354.,   27838.,   29364.,   30933.,   32545.,   34198.,
     *     35892.,   37627.,   39403.,   41220.,   43076.,   44971./
      data temp1/130./,dt/5.0/,ntt/60/
c
c
c        new planck function (wide band)
      data blkco2/1946., 2526., 3222., 4044., 5002., 6106.,
     *  7367.,    8790.,   10386.,   12160.,   14118.,   16265.,
     * 18606.,   21145.,   23884.,   26826.,   29972.,   33324.,
     * 36881.,   40645.,   44614.,   48788.,   53166.,   57746.,
     * 62525.,   67503.,   72676.,   78042.,   83599.,   89344.,
     * 95274.,  101385.,  107676.,  114140.,  120779.,  127586.,
     * 134558.,  141692.,  148986.,  156435.,  164037.,  171788.,
     * 179685.,  187724.,  195904.,  204219.,  212668.,  221248.,
     * 229955.,  238786.,  247739.,  256812.,  266000.,  275302.,
     * 284714.,  294234.,  303861.,  313591.,  323422.,  333351./
      t(xw,a,b,xm,w)=xw*exp(-a*w/(1.+b*w**xm))
      nptp1=npt+1
      nptm1=npt-1
      ta(nptp1)=ts
      do 15 ip=1,nptp1
c...planck function for layer center
      dx=ta(ip)-temp1
      it=dx/dt+1.5
      if(it .lt. 2) it=2
      if(it .gt. ntt) it=ntt
      itm1=it-1
      x=temp1+dt*float(itm1)
      ds=ta(ip)-x
      bai(ip)=blkco2(it)+(blkco2(it)-blkco2(itm1))*ds/dt
      bwi1(ip)=blkwin(it,1)+(blkwin(it,1)-blkwin(itm1,1))*ds/dt
      bwi2(ip)=blkwin(it,2)+(blkwin(it,2)-blkwin(itm1,2))*ds/dt
   15 continue
c  surface
      bs=bai(nptp1)
      bss1=bwi1(nptp1)
      bss2=bwi2(nptp1)
      do 17 ip=1,npt
c  planck function for layer edge
      dx=tui(ip)-temp1
      it=dx/dt+1.5
      if(it .lt. 2) it=2
      if(it .gt. ntt) it=ntt
      itm1=it-1
      x=temp1+dt*float(itm1)
      ds=tui(ip)-x
      buic(ip)=blkco2(it)+(blkco2(it)-blkco2(itm1))*ds/dt
      buiw1(ip)=blkwin(it,1)+(blkwin(it,1)-blkwin(itm1,1))*ds/dt
      buiw2(ip)=blkwin(it,2)+(blkwin(it,2)-blkwin(itm1,2))*ds/dt
   17 continue
   18 continue
c  collect transmission functions
      do 28 ip=1,nptm1
      ipp1=ip+1
      do 32 ix=ipp1,npt
      if(ip.eq. ix) go to 32
c   .... 15 micron h2o line
      dxl=abs(wbar(ip,2)-wbar(ix,2))
      if (dxl.lt.0.00001) go to 7000
      x=t(1.,6.7,16.,0.6,dxl)
      go to 7010
 7000 continue
      x=1.0
 7010 continue
c... h2o e-type
      dxc=abs(ubar(ip)-ubar(ix))
      if (dxc.lt.0.00005) go to 38
c   ...15 micron h2o e-type
      y=27.*dxc**.83
      y=1./exp(y)
c   ...9.6 micron h2o e-type
      yy=1./exp(9.79*dxc)
      go to 40
   38 continue
      y=1.
      yy=1.
   40 continue
      tauco2(ip,ix)=tauco2(ip,ix)*x*y
      tauco2(ix,ip)=tauco2(ip,ix)
      tauf1(ip,ix)=tauf1(ip,ix)*yy
      tauf2(ip,ix)=tauf2(ip,ix)*yy
      tauf1(ix,ip)=tauf1(ip,ix)
      tauf2(ix,ip)=tauf2(ip,ix)
   32 continue
      tauco2(ip,ip)=1.0
      tauf1(ip,ip)=1.0
      tauf2(ip,ip)=1.0
   28 continue
      tauco2(npt,npt)=1.0
      tauf1(npt,npt)=1.0
      tauf2(npt,npt)=1.0
c  get 9.6 micron band fluxes
      ib=2
      call flux(bwi1,buiw1,bss1,tauf1,fluxuw1,fluxdw1,ib)
      call flux(bwi2,buiw2,bss2,tauf2,fluxuw2,fluxdw2,ib)
      do 95 ip=1,npt
      fluxuw(ip)=fluxuw1(ip) + fluxuw2(ip)
   95 fluxdw(ip)=fluxdw1(ip) + fluxdw2(ip)
c
c  get 15 micron band fluxes
      ib=3
      call flux(bai,buic,bs,tauco2,fluxuc,fluxdc,ib)
c
c
c********* total fluxes in the 9.6 and 15 micron bands *************
c
      do 97 ip=1,npt
      fluxu(ip)=fluxuc(ip)+fluxuw(ip)
      fluxd(ip)=fluxdc(ip)+fluxdw(ip)
   97 continue

      return
      end
c
      subroutine radco2(nvert,pa,ta,pu,lprint)
c
c    co2 15 micron transmissions from linear expansion model
c
      dimension pui(49),tai(49),pai(49)
      dimension puil(49),pail(49)
      dimension tstd(49),dt(49),fsn(19600)
      dimension tauco2(48,48),ta(48),pa(48),pu(48),palog(48),pulog(48)
      logical lprint
      dimension tauf(49,49),taustd(49,49),m1(49)
      common/co2/tauco2,taustd,fsn
c
c    fixed grid edge pressures
c
      data pui/      0.096000e-02, 0.132000e-02,
     + 0.180000e-02, 0.245944e-02, 0.336047e-02, 0.459160e-02,
     + 0.627376e-02, 0.857219e-02, 0.117126e-01, 0.160036e-01,
     + 0.218667e-01, 0.298776e-01, 0.408234e-01, 0.557793e-01,
     + 0.762144e-01, 0.104136e+00, 0.142287e+00, 0.194414e+00,
     + 0.265639e+00, 0.362957e+00, 0.495928e+00, 0.677613e+00,
     + 0.925861e+00, 0.126506e+01, 0.172851e+01, 0.236176e+01,
     + 0.322701e+01, 0.440924e+01, 0.602459e+01, 0.823173e+01,
     + 0.112475e+02, 0.153680e+02, 0.209982e+02, 0.286910e+02,
     + 0.392021e+02, 0.535639e+02, 0.731874e+02, 0.100000e+03,
     + 0.135000e+03, 0.180000e+03, 0.250000e+03, 0.340000e+03,
     + 0.435800e+03, 0.548700e+03, 0.661600e+03, 0.774500e+03,
     + 0.887400e+03, 0.100000e+04, 0.105000e+04/
c
c    fixed grid mid pressures
c
      data pai/
     +  0.480000e-03, 0.114000e-02, 0.156000e-02, 0.210404e-02,
     +  0.287487e-02, 0.392810e-02, 0.536718e-02, 0.733348e-02,
     +  0.100201e-01, 0.136910e-01, 0.187068e-01, 0.255602e-01,
     +  0.349243e-01, 0.477190e-01, 0.652011e-01, 0.890879e-01,
     +  0.121726e+00, 0.166321e+00, 0.227253e+00, 0.310508e+00,
     +  0.424265e+00, 0.579696e+00, 0.792070e+00, 0.108225e+01,
     +  0.147874e+01, 0.202048e+01, 0.276069e+01, 0.377209e+01,
     +  0.515401e+01, 0.704221e+01, 0.962218e+01, 0.131473e+02,
     +  0.179639e+02, 0.245450e+02, 0.335372e+02, 0.458238e+02,
     +  0.626115e+02, 0.855496e+02, 0.116189e+03, 0.155885e+03,
     +  0.212132e+03, 0.291547e+03, 0.384931e+03, 0.489002e+03,
     +  0.602511e+03, 0.715827e+03, 0.829030e+03, 0.942019e+03,
     +  0.102500e+04/
c
      data tstd/180.650, 180.650,
     + 180.650, 180.650, 180.650, 180.650, 180.650, 180.650,
     + 182.280, 186.648, 193.579, 200.769, 208.247, 216.002, 224.040,
     + 232.373, 241.023, 249.663, 255.948, 260.668, 265.473, 270.363,
     + 270.650, 269.535, 264.388, 257.706, 251.191, 244.852, 238.673,
     + 232.648, 227.975, 225.902, 223.848, 221.812, 219.795, 217.796,
     + 216.650, 216.650, 216.650, 216.650, 218.546, 227.426, 239.757,
     + 250.892, 261.043, 269.793, 277.400, 284.246, 288.753/
c
      data npt/49/
c
      data itime/1/
c
c  linear fitting function
      v(pr,pr1,pr2,v1,v2)=v2 + (v1-v2)*(pr-pr2)/(pr1-pr2)
c
      nptm1=npt-1
      ndata=nvert + 3
      ndata1=ndata-1
c if first time routine is called, get logs
      if(itime.eq.0) go to 20
      do 10 n=1,npt
      puil(n)=alog(pui(n))
      pail(n)=alog(pai(n))
   10 continue
      itime=0
   20 continue
      do 30 n=1,ndata
      palog(n)=alog(pa(n))
   30 pulog(n)=alog(pu(n))
c
c    interpolate temperatures from model grid to fixed grid
c  loop over longitudes
c
      do 50 n=2,npt
      pr=pai(n)
      do 45 k=3,ndata
      k1=k
      if(pr.le.pa(k)) go to 46
   45 continue
   46 continue
      pr=pail(n)
      pr1=palog(k1-1)
      pr2=palog(k1)
      v1=ta(k1-1)
      v2=ta(k1)
      tai(n)=v(pr,pr1,pr2,v1,v2)
   50 continue
c
c    compute transmissions on fixed pressure grid
c
      do 200 i=1,npt
      do 200 j=i,npt
      tauf(i,j)=taustd(i,j)
  200 continue
c
      is=0
      do 235 i=2,npt
  235 dt(i)=tai(i)-tstd(i)
c
      if(.not. lprint ) go to 239
      write(6,236)
      do 238 j=2,npt
      write(6,237) j,pai(j),tai(j),tstd(j),dt(j)
  238 continue
  239 continue
c
      do 240 i=1,nptm1
      i1=i+1
      do 245 j=i1,npt
      do 250 k=i1,j
      is=is+1
      tauf(i,j)=tauf(i,j)+fsn(is)*dt(k)
  250 continue
  245 continue
  240 continue
c
c      if(.not. lprint) go to 516
c      write(6,500)
c      do 515 i=1,npt
c      write(6,505) i
c      write(6,510)(tauf(i,j),j=i,npt)
c  515 continue
c  516 continue
c
c    interpolate taus from fixed to model grid
c
	do 691 jm=1,ndata
	do 661 l=1,npt
	m1(jm)=l
	if (pu(jm).le.pui(l)) go to 666
661	continue
666	continue
691	continue
 
 
      do 690 i=1,ndata1
      i1=i+1
      do 690 j=i1,ndata
      prj=pu(j)
      pri=pu(i)
	l1=m1(j)
	k1=m1(i)
c
c    before fitting, see if the model pressures currently being
c    considered lie above the fixed grid
c
      if(pri .ge. pui(1)) go to 680
      if(prj .lt. pui(1)) go to 675
c
c    only pri is out of range
c
      pr=pulog(j)
      pr1=puil(l1-1)
      pr2=puil(l1)
      v1=tauf(1,l1-1)
      v2=tauf(1,l1)
      tauco2(i,j)=v(pr,pr1,pr2,v1,v2)
      go to 689
c
c    for both out of range, set tauco2 to unity
c
  675 tauco2(i,j)=1.00
      go to 689
c
c    fully in range, do fitting
c
  680 continue
      pr=pulog(j)
      pr1=puil(l1-1)
      pr2=puil(l1)
      v1=tauf(k1-1,l1-1)
      v2=tauf(k1-1,l1)
      taukm1=v(pr,pr1,pr2,v1,v2)
      v1=tauf(k1,l1-1)
      v2=tauf(k1,l1)
      tauk=v(pr,pr1,pr2,v1,v2)
      pr=pulog(i)
      pr1=puil(k1-1)
      pr2=puil(k1)
      v1=taukm1
      v2=tauk
      tauco2(i,j)=v(pr,pr1,pr2,v1,v2)
c
  689 continue
c
  690 continue
c
c    next longitude
c
    1 continue
c
      do 692 i=1,ndata
      tauco2(i,i)=1.0
  692 continue
c
      do 695 i=1,ndata1
      i1=i+1
      do 695 j=i1,ndata
      tauco2(j,i)=tauco2(i,j)
  695 continue
c
c      if(.not. lprint) go to 717
c      write(6,700)
c      do 715 i=1,ndata
c      write(6,505) i
c      write(6,510) (tauco2(i,j),j=i,ndata)
c  715 continue
c  717 continue
c
      return
c
  236 format(1h0,'printout from radco2',/,
     +       14x,'pai ',5x,'tai',8x,'tstd',8x,'dt',/)
  237 format(5x,i4,f10.4,3f10.2)
  500 format(1h0,'model transmissions on fixed grid for first point',/)
  505 format(' i =',i4)
  510 format(10(1x,f8.5))
  700 format(1h0,'transmissions on variable grid for first point',/)
      end
c
      subroutine rado3(nvert,tai,pai,pui,lprint)
c...o3 9.6 micron transmissions
c...new model
      dimension dpt(48),pint(48),tint(48)
      logical lprint, wing
      dimension tauf1(48,48),tauf2(48,48)
      dimension taustd(49,49),fsn(19600)
      dimension pai(48),tai(48),rawi(48),dpi(48),pui(48),skw(48)
      dimension tauco2(48,48),raw(48)
      common/co2/tauco2,taustd,fsn
      common/tro3/tauf1,tauf2
      common/ozot/ vmro3(48)
 
      npt=nvert+3
      dpck=.0001
      nptm1=npt-1
      factor=48./28.964
      do 1010 l=2,npt
      dpi(l)=pui(l)-pui(l-1)
      dpi(l)=dpi(l)*1.02
 1010 continue
      pai(1)=.5*pui(1)
      do 4 i=2,npt
    4 raw(i)=vmro3(npt+2-i)
      if(.not.lprint) go to 9
      write(6,5)
    5 format(//,' printout from rado3')
      write(6,8) (i,pai(i),tai(i),raw(i),i=2,npt)
    8 format(i4,f10.4,f10.2,e15.3)
    9 continue
      do 1030 l=2,npt
      rawi(l)=raw(l)*dpi(l)*factor*1.e-6
 1030 continue
      skw(npt)=0.
      pint(npt)=0.
      tint(npt)=0.
      do 15 i=1,nptm1
      skw(i)=0.
      pint(i)=0.
      tint(i)=0.
      do 15 j=1,npt
      tauf1(i,j)=1.
      tauf2(i,j)=1.
   15 continue
c
c        compute column ozone and effective pressure
c
      do 50 l=2,npt
      skw(l)=skw(l-1)+rawi(l)
      pint(l)=pint(l-1)+pai(l)*rawi(l)
      tint(l)=tint(l-1)+tai(l)*rawi(l)
   50 continue
c
c       compute transmittances
c
      do 70 ip=2,npt
      ipm1=ip-1
      do 75 lp=1,ipm1
      ix=ip-lp
      dx=skw(ip)-skw(ix)
      if(dx.le.1.e-30) go to 75
      dpx=pint(ip)-pint(ix)
      peff=dpx/dx
      teff=(tint(ip)-tint(ix))/dx
c  get dx in cm-atm
      dx=dx*1.2553d22/2.69d19
c  band wing
      wing=.true.
      call trano3(dx,peff,teff,trans,wing)
      tauf1(ix,ip)=trans
c  band center
      wing=.false.
      call trano3(dx,peff,teff,trans,wing)
      tauf2(ix,ip)=trans
      if(tauf1(ix,ip).lt.dpck.and.tauf2(ix,ip).lt.dpck) go to 70
   75 continue
   70 continue
c
      do 31 i=1,npt
      tauf1(i,i)=1.0
      tauf2(i,i)=1.0
   31 continue
      do 32 i=1,nptm1
      j1=i+1
      do 32 j=j1,npt
      tauf1(j,i)=tauf1(i,j)
      tauf2(j,i)=tauf2(i,j)
   32 continue
c      if(.not.lprint) go to 711
c      write(6,710)
c  710 format(/,' *** o3 9.6 micron transmittance matrix',/)
c      write(6,712)
c  712 format('  band wing region ',/)
c      do 41 i=1,npt
c      write(6,43) i
c   43 format(/,' i = ',i4)
c      write(6,42) (tauf1(i,j),j=i,npt)
c   42 format(1x,10f7.5)
c   41 continue
c      write(6,713)
c  713 format(//,' band center region ',//)
c      do 741 i=1,npt
c      write(6,43) i
c      write(6,42) (tauf2(i,j),j=i,npt)
c  741 continue
c  711 continue
      return
      end
c
      subroutine trano3(dx,peff,teff,trans,wing)
      logical wing
      dimension absc(3,2),tmp(3)
      dimension alp(10,2),pra(10)
      dimension alp1(10)
      dimension delt(7),pr(7)
      data pi/3.1415927/
c  line width parameter
      data pra/0.0, 0.25, 1., 2.51, 10., 25.1, 100., 251., 1000., 1050./
      data alp/.00028, .00028, .00035, .0005, .0011, .0020, .0050,
     * 0.0070, 0.0090, 0.0090,
     *         0.0006, 0.0006, .00068, .00088, .0022, .0050, .016,
     *  0.035, 0.070, 0.070/
c  line strength parameter
      data tmp/ 200., 250., 300./
      data absc/ 0.1500 , 0.2050 , 0.2520 ,
     *          0.9200 , 0.8850 , 0.8050 /
c  line spacing parameter
      data pr/0.0, 0.25, 1., 10., 100., 1000., 1050./
      data delta/0.1/
c
      nb=2
      if(wing) nb=1
c  line spacing parameter is delta ( for center)
c  line width parameter
      ind=2
      do 20 j=2,10
      if(peff.lt.pra(j)) go to 20
      ind = j + 1
   20 continue
      alpha=alp(ind,nb)+(alp(ind-1,nb)-alp(ind,nb))*(peff-pra(ind))
     *  /(pra(ind-1)-pra(ind))
c  line strength parameter
      ind=3
      if(teff.lt.250.) ind=2
      ak=absc(ind,nb)+(absc(ind-1,nb)-absc(ind,nb))*(teff-tmp(ind))
     *  /(tmp(ind-1)-tmp(ind))
c
      u=dx*1.66
      xe=ak*u/(2.*pi*alpha)
      term1=2.*pi*alpha*xe/delta
      term2=1. + 8.*xe/pi
      tax=-term1/sqrt(term2)
      if(tax.lt.-50.) tax=-50.
      trans=exp(tax)
      return
      end
c
      subroutine cloudy
      logical aerfl
      common /init/nup1,nsb,nc,npt,nw,nt,im
      common/cloud1/cfrac(3,48),pline(3,48,48)
      common /amhn/pu(48),ta(49),pa(48),wa(48),ps
      common/cirrus/icir(45),hcir(45),sigcir(45),taucir(45),epsa(45),
     *  eps11(45),tauvis(45),prmax
      dimension fmax(48),fran(48)
      dimension pmax(48,48),pran(48,48)
c... returns probabilities of clear line of sight between levels
c...fmax, fran are maximum and random overlap cloud fractions
c...take clouds overlapping maximally below prmax mb,
c...randomly above prmax mb
      do 1000 ib=1,3
      do 5 ip=1,npt
      do 5 ix=1,npt
    5 pline(ib,ip,ix)=1.0
 1000 continue

c
      do 2000 ib=1,3

      do 10 n=1,npt
      fran(n)=0.
      fmax(n)=0.
      if(pu(n).lt.prmax) fran(n)=cfrac(ib,n)
      if(pu(n).ge.prmax) fmax(n)=cfrac(ib,n)
c     write(6,777) n,pu(n),cfrac(ib,n),fran(n),fmax(n)
c 777 format(' n,pu,cfrac,fran,fmax:',i4,f10.3,3f10.3)
   10 continue
      nptm1=npt-1
      do 15 ip=1,npt
      pmax(ip,ip)=1.0
   15 pran(ip,ip)=1.0
      do 20 ip=1,nptm1
      ipp1=ip+1
      do 20 ix=ipp1,npt
      pmax(ip,ix)=amin1(pmax(ip,ix-1),(1.-fmax(ix)))
   20 continue
      do 25 ip=1,nptm1
c     write(6,800) ip
c 800 format(' ip = ',i4)
      ipp1=ip+1
      do 25 ix=ipp1,npt
      pran(ip,ix)=pran(ip,ix-1)*(1.-fran(ix))
c     write(6,805) ix,pran(ip,ix)
c 805 format(' ix = ',i4,' pran = ',f10.3)
   25 continue
      do 30 ip=1,nptm1
      ipp1=ip+1
      do 30 ix=ipp1,npt
      pline(ib,ip,ix)=pmax(ip,ix)*pran(ip,ix)
      pline(ib,ix,ip)=pline(ib,ip,ix)
   30 continue
 2000 continue
      return
      end
c
      subroutine flux(bmid,bedge,bs,tau,fluxu,fluxd,ib)
c...routine to do flux computation
      dimension bmid(49),bedge(49),tau(48,48)
      dimension fluxu(48),fluxd(48)
      common/cloud1/cfrac(3,48),pline(3,48,48)
      common /init/ nup1,nsb,nc,npt,nw,nt,im
      nptm1=npt-1
c
c***  compute downward fluxes
c
c   from atmosphere
      fluxd(1)=0.0
      do 100 ip=2,npt
      fluxd(ip)=0.0
      do 50 ix=2,ip
      fluxd(ip)=fluxd(ip)+bmid(ix)*pline(ib,ip,ix-1)*
     *          (tau(ip,ix)-tau(ip,ix-1))
   50 continue
c   from cloud bottoms
      do 60 ix=2,ip
      fluxd(ip)=fluxd(ip)+bedge(ix)*tau(ix,ip)*
     *          amax1(0.,(pline(ib,ix,ip)-pline(ib,ix-1,ip)))
   60 continue
   65 continue
  100 continue
c
c
c***  compute upward fluxes
c
c   from atmosphere
      fluxu(npt)=bs
      do 300 ip=1,nptm1
      fluxu(ip)=0.0
      ipp1=ip+1
      do 200 ix=ipp1,npt
      fluxu(ip)=fluxu(ip)+bmid(ix)*pline(ib,ix,ip)*
     *          (tau(ix-1,ip)-tau(ix,ip))
  200 continue
c from ground
      fluxu(ip)=fluxu(ip)+bs*tau(ip,npt)*pline(ib,ip,npt)
c  from cloud tops
      do 250 ix=ip,nptm1
      fluxu(ip)=fluxu(ip)+bedge(ix)*tau(ix,ip)*
     *          amax1(0.,(pline(ib,ix,ip)-pline(ib,ix+1,ip)))
  250 continue
  255 continue
  300 continue

      write (776, 1333) 
 1333 format(1x)

      return
      end

      subroutine aeros1(nz,ideltyr,month,nday,j,lprint)
c
c gets optical depths for sulfate aerosols
c used in radiative heating code
c
      include 'PARAM.INC'
      include 'aero.inc'
      include 'COMMONR1.INC'
      logical aeread, back,volgdh
      integer ideltyr,ixsrd
      real lataer,irflx

      common/aerheat/nzaer,nltaer,paer(45),lataer(36),backext(12,36,45),
     * volcext(96,36,45),aeread

      dimension aerden(45),volgdh(45),ext(45)

      common/hght/fcloud,fclear,ntopf,ntopt,rhostd,dlat(n$),
     * long,rlat(n$)
      logical lprint 
      real los

      common/cyear/iyr,dyt,dayst,dayin
      COMMON /YEAR360/ DAY360,DYT360, IDAY360

      common/aerxs/vxs_15m(6,45),vxs_96m(6,45),vxs_h2o(6,45),
     *   vxs_1m(6,45),vxs_nire(6,45),vxs_nira(6,45)
     
      data los/2.69e19/, pi/3.1415927/, ixsrd/0/
c
      
      if(lprint) write(6,2)
    2 format('  printout from aeros') 


  
c read aerosol cross sections on first call
      if(ixsrd .eq. 0) then
         open(unit=200,file='ext.1micron.2dn',status='old',
     *    convert="big_endian", form='unformatted')
	 read(200) vxs_1m
	 close (unit=200)
	 
         open(unit=200,file='abs.9.6micron.2dn',status='old',
     *    convert="big_endian", form='unformatted')
	 read(200) vxs_96m
	 close (unit=200)
	 
         open(unit=200,file='abs.15micron.2dn',status='old',
     *    convert="big_endian", form='unformatted')
	 read(200) vxs_15m
	 close (unit=200)
	 
         open(unit=200,file='abs.h2o.2dn',status='old',
     *    convert="big_endian", form='unformatted')
	 read(200) vxs_h2o
	 close (unit=200)
	 
         open(unit=200,file='ext.nearir.2dn',status='old',
     *    convert="big_endian", form='unformatted')
	 read(200) vxs_nire
	 close (unit=200)
	 
         open(unit=200,file='abs.nearir.2dn',status='old',
     *    convert="big_endian", form='unformatted')
	 read(200) vxs_nira
	 close (unit=200)
	 
	 ixsrd=1
	 
      endif	 
	 	  
c  read aerosol extinctions if first call

      if (.not. aeread) call rdaer1
c
c choose whether background or volcanic
     
      if(lprint) write(6,3) volcano
    3 format(' volcano = ',l4)

      back=.true.
      if (volcano .and. iyr.eq.56 .and. month.ge.6) back=.false.
      if (volcano .and. iyr.ge.57 .and. iyr .le. 62) back=.false.
   
      if( iyr .eq. 56 ) mon=month+12
      if( iyr .eq. 57 ) mon=month+24
      if( iyr .eq. 58 ) mon=month+36
      if( iyr .eq. 59 ) mon=month+48
      if( iyr .eq. 60 ) mon=month+60
      if( iyr .eq. 61 ) mon=month+72
      if( iyr .eq. 62 ) mon=month+84

        if(back) mon=month

      if(lprint) write(6,17) back,ideltyr,month,mon
   17 format(' AEROS1: back, ideltyr, month, mon = ',l4,3i4)

      do 400 n=1,nz

      iaer(n)=0
      
      volgdh(n)=.false.
      if( .not. back .and. (volcext(mon,j,n) .gt. backext(month,j,n)))
     * volgdh(n)=.true.
      ext(n)=backext(month,j,n)
      if(volgdh(n)) ext(n)=volcext(mon,j,n)
  400 continue
  
c  find topmost layer with aerosol (for solaer)
      do 700 n=1,nz
      ntaer = n
      if(ext(n).gt.0.0 .and. ext(n).lt.1.e4) go to 705
  700 continue
  705 continue
      

      do 500 n=ntaer,nz

      iaer(n)=1
      aerabs(n)=0.
      aerext(n)=0.
      aerh(n)=0.
      aer96(n)=0.
      aer15(n)=0.
      
      if(pl(n) .gt. 200.) go to 500

c convert to number density
c      if(back) sig1m=1.9758e-10
c      if(.not. back) sig1m=2.0329e-9

       if(back) sig1m=vxs_1m(1,n)
       if(.not. back .and. iyr .eq. 56) sig1m=vxs_1m(2,n)
       if(.not. back .and. iyr .eq. 57) sig1m=vxs_1m(3,n)
       if(.not. back .and. iyr .eq. 58) sig1m=vxs_1m(4,n)
       if(.not. back .and. iyr .eq. 59) sig1m=vxs_1m(5,n)
       if(.not. back .and. iyr .eq. 60) sig1m=vxs_1m(6,n)
       if(.not. back .and. iyr .eq. 61) sig1m=vxs_1m(6,n)
       if(.not. back .and. iyr .ge. 62) sig1m=vxs_1m(1,n)
        
       aerden(n)=ext(n)*1.e-5/sig1m
c
c  sulfate aerosol (75% h2so4) cross sections
c
      if ( volgdh(n) ) then

c Deshler, 8/2/91; bimodal
c
c  solar
c     values with no solar weighting
c         solext=1.1261e-09
c         solabs=6.0101e-11
c     values weighted by solar radiance   (PEM, 8/5/95, as per JER)
c         solext = 1.9928e-09
c         solabs = 3.1224e-12
cjer   values when weighted by solar radiance and with additional near-ir
c   absorption according to Pollack, GRL 8, 26-28, 1981.
c         solext = 1.9902e-09
c         solabs = 2.1716e-11
c  derived from Laramie profiles
c  values when weighted by solar radiance and with additional near-ir
         if(iyr .eq.56) then
	   solext=vxs_nire(2,n)
	   solabs = vxs_nira(2,n)
         endif
         if(iyr .eq.57) then
	   solext=vxs_nire(3,n)
	   solabs = vxs_nira(3,n)
         endif
         if(iyr .eq.58) then
	   solext=vxs_nire(4,n)
	   solabs = vxs_nira(4,n)
         endif
         if(iyr .eq.59) then
	   solext=vxs_nire(5,n)
	   solabs = vxs_nira(5,n)
         endif
         if(iyr .eq.60) then
	   solext=vxs_nire(6,n)
	   solabs = vxs_nira(6,n)
         endif
         if(iyr .eq.61) then
	   solext=vxs_nire(6,n)
	   solabs = vxs_nira(6,n)
         endif
         if(iyr .ge.62) then
	   solext=vxs_nire(1,n)
	   solabs = vxs_nira(1,n)
         endif
     
c
c  ir
c         sigh = 2.0895e-10
c         sig96= 2.8435e-10
c         sig15= 9.0882e-11

         if(iyr .eq. 56) then
           sigh = vxs_h2o(2,n)
	   sig96 = vxs_96m(2,n)
	   sig15 = vxs_15m(2,n)
         endif
         if(iyr .eq. 57) then
           sigh = vxs_h2o(3,n)
	   sig96 = vxs_96m(3,n)
	   sig15 = vxs_15m(3,n)
         endif
         if(iyr .eq. 58) then
           sigh = vxs_h2o(4,n)
	   sig96 = vxs_96m(4,n)
	   sig15 = vxs_15m(4,n)
         endif
         if(iyr .eq. 59) then
           sigh = vxs_h2o(5,n)
	   sig96 = vxs_96m(5,n)
	   sig15 = vxs_15m(5,n)
         endif
         if(iyr .eq. 60) then
           sigh = vxs_h2o(6,n)
	   sig96 = vxs_96m(6,n)
	   sig15 = vxs_15m(6,n)
         endif
         if(iyr .eq. 61) then
           sigh = vxs_h2o(6,n)
	   sig96 = vxs_96m(6,n)
	   sig15 = vxs_15m(6,n)
         endif
         if(iyr .ge. 62) then
           sigh = vxs_h2o(1,n)
	   sig96 = vxs_96m(1,n)
	   sig15 = vxs_15m(1,n)
         endif

      else

c
c  Hofmann background, 1990; rmod=0.05, sigma=2.2
c
c  solar
c     unweighted:
c         solext=1.1664e-10
c         solabs=8.3760e-12
c         sigbak=7.5371e-11
c
c     solar radiance weighted:
c         solext = 2.1610e-10
c         solabs = 4.3163e-13
c   values when weighted by solar radiance and additional near-ir absorption
c         solext = 2.1663e-10
c         solabs = 2.8561e-12

c  derived from Laramie profiles
          solext = vxs_nire(1,n)
	  solabs = vxs_nira(1,n)


c
c  ir
c         sigh = 2.9811e-11
c         sig96= 4.0506e-11
c         sig15= 1.3247e-11

         sigh = vxs_h2o(1,n)
	 sig96 = vxs_96m(1,n)
	 sig15 = vxs_15m(1,n)


      endif
c
c  column density in grid box
      htkm=hie(n) - hie(n+1)
      colden=aerden(n) * htkm *1.e5
c
c  optical depths
c
c  for near ir solar
c  don't allow tropospheric aerosols
      aerabs(n)=solabs*colden
      aerext(n)=solext*colden
      
c  for ir
      oaerh=sigh*colden
      oaer96=sig96*colden
      oaer15=sig15*colden
c backscatter coefficient
c      if(back) backs=sigbak*den(n)
c      if(back) backpi=sigbak*den(n)/(4.*pi)
c
      if(.not.lprint .or. pl(n).lt.20.) go to 121
      write(6,*) ' n, press = ',n,pl(n)
      write(6,110)
  110 format(/,'  ir optical depths',/)
      write(6,116)
  116 format('   oaerh    oaer96    oaer15')
      write(6,117) oaerh,oaer96,oaer15
  117 format(1pe12.4,1pe12.4,1pe12.4)
  121 continue

      aerh(n)=1.0 - exp(-1.66*oaerh)
      aer96(n)=1.0 - exp(-1.66*oaer96)
      aer15(n)=1.0 - exp(-1.66*oaer15)

      if(.not.lprint .or. pl(n).lt.20.) go to 131
      write(6,*) 'solext,solabs = ',solext,solabs
      write(6,120) 
  120 format(/,'  ir emissivities',/)
      write(6,118)
  118 format('    aerh   aer96   aer15')
      write(6,117) aerh(n),aer96(n),aer15(n)

      write(6,125)
  125 format(/,'  solar near infrared optical depths',/)
      write(6,130)
  130 format('  aerext   aerabs')
      write(6,117) aerext(n),aerabs(n)
  
c        if(back) write(6,140) backs,backpi
c  140 format(/,' backscatter coeff. = ',1pe12.4,/,
c     * 'backscatter coeff./4pi = ',1pe12.4)       

  131 continue
  
c end of loop over vertical levels
  500 continue
  
      return
      end

      integer function nearidx(x,a,dim)

c Finds the index of the nearest value in an array
c T.Atwater 7/25/91

c Input:
c    x    -- real value to figure index for
c    a    -- real array to find index in; must be montonically increasing
c    dim  -- dimension of a
c Output:
c    returns index of array value closest to x
c        returns 1 if x is less than a(1)
c        returns dim if x is greater than a(dim)

      real x
      integer dim
      real a(dim)

      nearidx = 1
      if (x.le.a(1)) return
      do 10 i=1,dim-1
         if (x.le.a(i+1)) then
            if (abs(x-a(i)).lt.abs(x-a(i+1))) then
               nearidx=i
            else
               nearidx=i+1
            end if
            return
         end if
10    continue
      nearidx = dim

      return
      end
