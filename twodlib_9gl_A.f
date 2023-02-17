      SUBROUTINE SETCON(IHEATS)

C   SET UP CONSTANTS AND ARRAYS FOR TWOD MODEL

C   TRANSPORTED CONSTITUENTS ARE SET UP IN SOUBROUTINE 


      PARAMETER (NCH=5)
      PARAMETER (NTZ=18)
      PARAMETER (NTY=73)
      PARAMETER (NMON=13)

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
	include  'COMMONW.INC'

      INCLUDE 'GWRAY.INC'

      COMMON /SLAKWK/ PWAT(23),WAT(23,18,12),
     1 PS(25),VNC0(25),goz(NCH),TOZ(NCH),ZOZ(NCH),
     2 TTROP(NTY,NTZ,NMON),ZTROP(NTZ),YTROP(NTY),TTLAT(N$,NTZ),
     3 ZTROP2(NTZ),CSBUV(52,18,30),HEATW(12,18,30),XLATDJ(18),
     4 ZDJ(8),TEHE2(18,30)
     
       REAL MWBRK0(90),YMB(90)

c     This REAL statement added 10/07/94 by PEM to dimension arrays and
c        avoid compiler errors under IRIX 5.2.  Variables may not actually
c        get used.
      REAL AY(N$,M$), CY(N$,M$), B(N$,M$), AZ(N$,M$), CZ(N$,M$)

       DATA ZTROP/1000.,850.,700.,500.,400.,300.,250.,200.,150.,
     1100.,70.,50.,30.,10.,5.,2.,1.,0.4/

      DATA PS/-6.90,-6.80,-6.55,-6.214,-5.85,-5.48,-4.67,-4.05,
     1-3.33,-2.1,-0.8,0.0,0.5,1.0,1.8,2.2,3.125,4.5,5.1,6.43,
     2 8.11,9.51,10.58,11.5,12.5/

      DATA VNC0/0.25,      0.15,      0.12,       0.09,      0.07,
     1 0.06,      0.06,      0.06,    0.061,     0.08,
     2 0.135,    0.212,     0.220,     0.20,      0.172,
     3 0.15,      0.1 ,      0.1,      0.2,       0.3,
     4 0.8,      1.2 ,      1.4 ,      1.4,       1.4/


      DATA GOZ/6.E-7,1.2E-5,1.0E-4,2.0E-3,1.0E-3/
      DATA TOZ/0.,1.2E-12,2.0E-12,1.0E-12,-1.5E-12/
      DATA ZOZ/30.,40.,50.,60.,70./


      PI=3.14159
C       DRY AIR GAS CONSTANT
      R=2.87E6
C       LENGTH OF DAY IN SECONDS
      DAYL=86164.
C       TWO OMEGA
      TOMEG=2.*2.*PI/DAYL
cc      TOMEG=10.*2.*PI/DAYL
C       SCALE HEIGHT (AVERAGE) IN CM
      H=7.0E5
C       GRAVITY
      GZ=980.
C       VERTICAL GRID SPACING (CM)
C     DZ=0.38*H
      DZ=ZRES*H
      DZP=DZ*1.0E-5
C       PRANDTL NUMBER
      PRNDTL=1.4
c       time step
      DT=DAYL/ISW(6)
      FRACD=DT/DAYL
C       DELTA THETA (LATITUDE)
      DTH=180./(N$+1)
        print 800,dTH,dzP,FRACD
800     FORMAT('0 DY = ',F10.4,' DEG  DZ = ',F10.4,' KM',
     1' DT = ',F10.4,' DAYS')
      CVT=PI/180.
C        PLANETARY RADIUS
      A=6.37E8
C        R/CP
      XKAPPA=0.288
C        DELTA Y
      DY= A*DTH*CVT
      DZ1=0.5/DZ
      DY1=0.5/DY
      DZ2=2.*DZ1
      DY2=2.*DY1
C        SURFACE PRESSURE
      PS0=ALOG(1000.)
      XNORM=0.
      DO 1 J=1,N$
c     y in latitude degs
      YP(J)=-90.+J*DTH
c     y in radians
      YY(J)=YP(J)*CVT
      YPP(J)= -90. + (J - 0.5)*DTH
      S(J)=SIN(YY(J))
      C(J)=COS(YY(J))
      CST(J)=COS(YPP(J)*CVT)      
      SST(J)=SIN(YPP(J)*CVT)      
      XNORM=XNORM+C(J)
      CF(J)=S(J)*TOMEG
      CFST(J)=SST(J)*TOMEG
      BETA(J)=C(J)*TOMEG/A
1     CONTINUE

      YPP(NP$)=YP(N$)+0.5*DTH
      CST(NP$)=COS(YPP(NP$)*CVT)
      SST(NP$)=SIN(YPP(NP$)*CVT)
      CFST(NP$)=SST(NP$)*TOMEG
	zspong=50.
	ral0=50.

      DO 2 K=1,M$
      ZZ(K)=K*DZ
      ZP(K)=ZZ(K)*1.0E-5
      ZPP(K)=(ZZ(K)-DZ*0.5)*1.0E-5
      RHO(K)=EXP(-ZZ(K)/H) !NOT REALLY DENSITY
      EZ2H(K)=EXP(ZZ(K)*0.5/H)
      TTOTH(K)=EXP(ZZ(K)*XKAPPA/H)
      RSTAR(K)=R/(TTOTH(K)*H)
      
C      RAYLEIGH FRICTION TERM 
c      RAL(K)=ral0
c      if (zp(k).ge.zspong) ral(k) = ral0*exp(-(zp(k)-zspong)/12.)

C      RAL(K)=40./SQRT(EZ2H(K))
c      IF ((ZP(K).GE.10.).AND.(ZP(K).LE.30.)) RAL(K)=15.
c      IF (RAL(K).LT.1.5) RAL(K)=1.5
2     CONTINUE


c NEW RAYLEIGH FRICTION (JTB 4/28/92)

      read(101,*) rayshrt, raylong, rayz1, rayz2, thermyn, thermcr
      write (*,*)  rayshrt,raylong,rayz1,rayz2,thermyn,thermcr
      close(unit = 101)

c      rayz1=30.    !km
c      rayz2=80.    !km
c      rayshrt=1.5  !days
c      raylong=30.  !days
cjer these are currently
c      rayz1=10.    !km
c      rayz2=50.    !km
c      rayshrt=5.0  !days
c      raylong=80.  !days


c      write(6,*) 'new rayleigh damping '

      DO K=1,M$
         IF( ZP(K) .GE. rayz2 ) RAL(K)=rayshrt
         IF((ZP(K) .GE. rayz1 ) .AND. (ZP(K) .LT. rayz2 )) 
     >         RAL(K)=raylong-( ZP(K)-rayz1 )*(raylong-rayshrt)
     >              /(rayz2-rayz1)
         IF( ZP(K) .LT. rayz1 ) RAL(K)=raylong
c         write(6,708) zp(k), ral(k)
      END DO

c         write(6,*) '******'

C     SMOOTH OUT RAL PROFILE
      DO 331 KS=1,4
      DO 333 K=2,M1$
333      THG(K)=0.25*(RAL(K-1)+RAL(K+1))+0.5*RAL(K)
      DO 334 K=2,M1$
334      RAL(K)=THG(K)     
331   CONTINUE

c      print 708,(zp(k),ral(k),K=1,M$)
708   format(' alt = ',f10.2,' km   ral = ',f10.2,' days')

C     CONVERT RAL TO (SEC)^-1
         DO 332 K=1,M$
332      RAL(K)=1./(DAYL*RAL(K))

       ZPP(MP$)=(ZZ(M$)+DZ*0.5)*1.0E-5

       DO K=1,M$
          RALP(K)=RAL(K)
       END DO
       RALP(MP$)=RAL(M$)
       
       DO K=1,MP$
          Z_IN_CM=ZPP(K)*1.0E5
          RHOST(K)=EXP( -Z_IN_CM / H )
       END DO
       

cjer get global mean temperature, tg
      CALL GETT(TG,ZP,M$)

      DO 3 K=1,M$
      THG(K)=TG(K)*TTOTH(K)
3     CONTINUE

      DO 31 K=2,M1$
      DTHDZ(K)=(THG(K+1)-THG(K-1))*DZ1
31    CONTINUE
      DTHDZ(1)=DTHDZ(2)
      DTHDZ(M$)=DTHDZ(M1$)


c  get global mean temperature for MP$ grid, ZPP(MP$), TGX(MP$), TTOTHX(MP$), THGX(MP$),- EF, 4/09
      CALL USTA62(ZPP,TGX,MP$)

      DO 776 K=1,MP$
        TTOTHX(K) = EXP(ZPP(K)*XKAPPA/7.)
        THGX(K) = TGX(K)*TTOTHX(K)
 776  CONTINUE


C     SET UP TRANSPORT ARRAYS AND CONSTANTS FOR SCHEME

      DO 301 K = 1,MP$
      XEP(K)=EXP(-(DZ*(K-.5))/H)
      DO 301 J = 1,NP$   
      RHO0(J,K)=CST(J)*XEP(K)
      if (k.le.m$) RHOY(J,K)=CST(J)*RHO(K)
      IF (J.LE.N$) RHOX(J,K)=C(J)*XEP(K)
301   CONTINUE

      DO K=1,MP$
      DO J=1,NP$
         RHO00(J,K)=RHO0(J,K)
         if (k.le.m$) RHOY0(J,K)=RHOY(j,K)
         IF (J.LE.N$) RHOX0(J,K)=RHOX(J,K)
C                                              initialize TEMPBX(NP$,MP$) in COMMOND.INC
         TEMPBX(J,K) = TGX(K)
      END DO
      END DO




      DO 605 J=1,N$
605   YPV(J+1)=YP(J)
      YPV(1)=-90.
      YPV(NPP$)=90.
      DO 606 J=1,M$
      XN2(J)=GZ*DTHDZ(J)/THG(J)
606   ZPW(J+1)=ZP(J)
      ZPW(1)=0.
      ZPW(MPP$)=ZPW(MP$)+DZP

      DO J=2,M$
         XN2ST(J) = XN2(J)*0.5 + XN2(J-1)*0.5
      END DO
 
      XN2ST(1)  =XN2(1)
      XN2ST(MP$)=XN2(M$)



C     SET UP SLAK ARRAYS

      DYX=1./(DY*DY)
      DZX=1./(DZ*DZ)
      H25=0.25/(H*H)

      DO 4 K=1,M$
      TEMP= DTHDZ(K)
      DO 4 J=1,N$
	PWMOM(J,K)=0.
	GWMOM(J,K)=0.
      AY(J,K) = TEMP*(DYX+S(J)*DY1/(A*C(J)))/(CF(J)*CF(J))
      CY(J,K) = TEMP*(DYX-S(J)*DY1/(A*C(J)))/(CF(J)*CF(J))
      B(J,K)  = (TEMP*(-2.*DYX)/(CF(J)*CF(J)) - 2.*DZX -H25)
      AZ(J,K) = DZX
      CZ(J,K) = DZX
4     CONTINUE

c	Zero the planetary wave fields

	do 420 m=1,mmax$
	do 420 k=1,m$
	do 420 j=1,n$
420	psic(j,k,m)=0.


C      GET HEATING RATES

      CALL PRODUC

      DO 5 J=1,25
      PS(J)=H*(PS(J)+PS0)*1.0E-5
      VNC0(J)=VNC0(J)/DAYL
5      CONTINUE

        CALL INTER (VNC,ZP,M$,VNC0,PS,25)
c        PRINT 52
52      FORMAT('    J       Z         NEWT COOLING    TEMP      
     1 RALEIGH F.  PRESSURE')
      ZP00=1.2*H
      DO 51 K=1,M$
      WGT=1.0
      IF (ZZ(K).GT.ZP00) WGT=EXP(-(ZZ(K)-ZP00)/H)
      vnc(k)=0.5e-5*wgt+vnc(k)*(1.-wgt)
      PRESS=1013.25*RHO(K)
      DUM1(1,K)=1./(VNC(K)*DAYL)
      DUM1(2,K)=1./(RAL(K)*DAYL)
c      PRINT 40,K,ZP(K),DUM1(1,K),TG(K),DUM1(2,K),PRESS
40      FORMAT(' ',I4,5F14.4)
51      CONTINUE

C   RADSTART no longer needed with GEOS5 heating rates
C
      IF (IHEATS .NE. 1  .and.  IHEATS .LE. 4) CALL RADSTART

C      SET UP TROPOSPHERIC HEATING FUNCTIONS

       DO 400 K=1,NTY
400       YTROP(K)=-90.+2.5*(K-1)
       DO 399 K=1,NTZ
399       ZTROP(K)=-H*ALOG(ZTROP(K)/1000.)

c       READ IN THE TROPOSPHERIC TEMPERATURES, UNIT 90
	READ(90,81) (((ttrop(i,j,k),i=1,nty),j=1,ntz),k=1,13) 
81      format(18f7.1)
        close(unit=90)
       DO 4015 I=1,NMON
       CALL Regrid(DUM1,YP,ZZ,N$,M$,TTROP(1,1,I),YTROP,ZTROP,NTY,NTZ,0)
       DO 401 J=1,N$
       DO 401 K=1,M$
  401  TFROP(J,K,I)=DUM1(J,K)


4015    CONTINUE


C    INTERPOLATE CHEMISTRY COEFFICENTS FROM STOLARSKI AND DOUGLASS(1985)

      DO 111 K=1,NCH
      GOZ(K)=ALOG10(GOZ(K))
111   CONTINUE

      CALL INTER(CHGAMA,ZPP,MP$,GOZ,ZOZ,NCH)
      CALL INTER(CHTHET,ZPP,MP$,TOZ,ZOZ,NCH)

      IF(IPRINT.EQ.1) PRINT 113
113   FORMAT('0  HEIGHT (KM)    GAMMA (SEC -1)  THETA (SEC-1,DEG-1)')
      DO 112 K=1,MP$
      CHGAMA(K)=10.**CHGAMA(K)
      IF(IPRINT.EQ.1) PRINT 114,ZPP(K),CHGAMA(K),CHTHET(K)
112      CONTINUE
114      FORMAT(' ',F10.2,2G14.4)

C     SET UP WATER VAPOR DISTRIBUTION, UNIT 95

C     READ IN TABLES

      PWat(1)=0.01
      PWAT(2)=0.05
      PWAT(3)=0.2
      READ(595,500) (PWAT(N),N=4,23)
  500 FORMAT(6X,6F10.1)

      DO 550 M=1,12
      READ(595,505) MON
  505 FORMAT(I4)
      DO 540 L=1,18
      READ(595,510) WLAT(L),(WAT(N,L,M),N=4,10)
  510 FORMAT(F5.0,7F8.1)
      READ(595,515) (WAT(N,L,M),N=11,17)
      READ(595,515) (WAT(N,L,M),N=18,23)

C     ADD MESOSPHERIC H2O

      WAT(1,L,M)=1.6
      WAT(2,L,M)=3.8
      WAT(3,L,M)=6.7

  515 FORMAT(5X,7F8.1)
  540 CONTINUE
  550 CONTINUE
        close(unit=95)

C     AT THIS POINT THE TABLES ARE UPSIDE DOWN IN PRESSURE

C     REORGANIZE TABLES
  
      DO 620 K=1,23
620   ZWAT(K)=H*ALOG(1000./PWAT(24-K))

      DO 621 K=1,23
      DO 621 M=1,12
      DO 621 L=1,18
621   WATZ(L,K,M)=WAT(24-K,L,M)

c	READ IN THE SBUV OZONE TABLES

	nzsbuv=17
	nolat=18
	do 38 il=1,nolat
38	xlsbuv(il)=-85.+(il-1)*10.
	read(91,405)(zsbuv(nzsbuv-n+1),n=1,nzsbuv)
405	format(1x,f9.3)
	do 39 i=1,nzsbuv
39      zsbuv(i)=H*alog(1023./zsbuv(i))
	do 85 mon=1,12
	do 83 il=1,nolat
	read(91,77)(sbuv(il,nzsbuv-ip+1,mon),ip=1,10)
	read(91,77)(sbuv(il,nzsbuv-ip+1,mon),ip=11,nzsbuv)
77	format(6x,10f6.2)
83	continue
c             write(*,*)'reading ozone mon =  ', mon      
85	continue

C	READ IN THE WEAK FIXED HEATING RATES FOR THE TROPOSPHERE

       DO 9499 II=1,12
       DO 9388 IJ=1,18
	XLATDJ(IJ)=-95.+IJ*10.

C     THE FIRST PART OF THE DATA IS THROWN AWAY

      READ(516,91)(HEATW(II,IJ,IK),IK=1,30)
91     FORMAT(1X,15F5.1)

9388       CONTINUE
        close(unit=91)

       DO 950 IJ=1,18
      READ(516,911)(HEATW(II,IJ,IK),IK=1,30)
911   FORMAT(1X,15F5.2)
950     CONTINUE
9499    CONTINUE
        close(unit=16)
      

C     REORGANIZE THE TABLES

	DO 832 II=1,12
	DO 830 K=1,8
	ZDJ(K)=H*1.0E-5*0.2844*(K-0.5)
	DO 830 J=1,18
	TEHE2(J,K)=HEATW(II,J,K)
830	CONTINUE

       CALL Regrid(DUM1,YP,ZP,N$,M$,TEHE2,XLATDJ,ZDJ,18,8,0)

	DO 842 K=1,10
	DO 842 J=1,N$
	TEHE(J,K,II)=DUM1(J,K)
842	CONTINUE

832	CONTINUE

c	get planetary wave heights

	call heights


c       read in gravity wave forcing

      do 555 j=1,181
      jj=182-j
      read(47,3011) glat(jj),gavg(jj),gvar(jj)
c      print 445,glat(jj),gavg(jj),gvar(jj)
c445   format(' '3f14.1)
C3011   format(' ',3f10.1)
3011   format(1X,3f10.1)
      gvar(jj)=100.*gvar(jj)
555   continue
        close(unit=47)
      call inter(gvint,ypp,np$,gvar,glat,181)
c      print 444,(ypp(j),gvint(j),j=1,np$)
444   format('  lat=',f10.2,' var=',g14.4)


C     Read in NEW Gravity Wave parameters.  Added 1/30/96 by PEM,
C     as per JTB.

      CALL GWDRAG(0)

C     End added 1/30/96 by PEM.


c ---- READ IN MOUNTAIN WAVE BREAKING LEVELS ------

      DO J=1,90
         YMB(J)=2.0*(J-1)-89.
      END DO

      READ(98,445 )  mwrmp 
      READ(98,445 )  mwbrk0
445   format( f8.0  )

      CLOSE(UNIT = 98)

      CALL INTER(MWBRK,YPP,NP$,MWBRK0,YMB,90)

c      write(*,*) '>>>>>>>>>>>> MWBRK <<<<<<<<<<<<<< '
c      write(*,*) '  RAMP UP FOR MTN WAVES --', mwrmp
c      write(*,*) mwbrk
c      write(*,*) '>>>>>>>>>>>> MWBRK <<<<<<<<<<<<<< '

      RETURN
      END


	SUBROUTINE HEIGHTS

C       READ IN THE PLANETARY WAVE HEIGHT FIELD DATA

	INCLUDE 'PARAM.INC'
	INCLUDE 'COMMONC.INC'
	INCLUDE 'COMMOND.INC'
	INCLUDE 'COMMONW.INC'
	COMPLEX ONE/(1.0,0.0)/,EYE/(0.0,1.0)/
	DIMENSION DATA(2,11,180),XLAT(180),HIN(180),HOUT(N$)

c      	OPEN(89,FILE='HEIGHTS.DAT',FORM='UNFORMATTED',STATUS='OLD')

	DO 10 J=1,180
	READ(89,81) XLAT(J),((DATA(K,I,J),K=1,2),I=1,11)
ccc	WRITE(6,81) XLAT(J),((DATA(K,I,J),K=1,2),I=1,11)
  10	CONTINUE
        close(unit=89)

c81      format(18f7.1)
81      format(23f7.1)

	DO 20 M=1,MMAX$
	DO 20 I=1,N$
	  HEIGHT(I,M) = (0.0,0.0)
  20	CONTINUE

C	NOW INTEROPLATE ONTO TWOD'S GRID FOR EACH WAVENUMBER
	DO 100 M=1,MMAX$

C	REAL PART (COSINE):
	DO 25 I=1,180
	  HIN(I) = DATA(1,M,I)
  25	CONTINUE

	CALL INTER(HOUT,Yp,N$,HIN,XLAT,180)

	DO 30 I=1,N$
	  HEIGHT(I,M) = ONE*HOUT(I)
  30	CONTINUE

C	IMAGINARY PART (SINE):
	DO 40 I=1,180
	  HIN(I) = DATA(2,M,I)
  40	CONTINUE

	CALL INTER(HOUT,Yp,N$,HIN,XLAT,180)

C	CONVERT HEIGHTS FROM METERS TO CENTIMETERS:
	DO 50 I=1,N$
	  HEIGHT(I,M) = 100.*(HEIGHT(I,M)+EYE*HOUT(I))
  50	CONTINUE


C       Modified 1/29/96 by PEM.  Forcing topography is now zeroed out
C       in JPWAVE.

C        WRITE(6,*) ' ZEROING PW FORCING POLEWARD OF ',ISW(63)
C        ZERLAT=1.00*ISW(63)

C	SET THE HEIGHTS OVER ANTARCTICA = 0.

C        DO I=1,N$
C           IF(YP(I).LT.-ZERLAT) HEIGHT(I,M) = 0.*HEIGHT(I,M)
C        END DO

C	SET THE HEIGHTS OVER ARCTICA (NORTH OF ZERLAT) = 0.

C        DO I=1,N$
C           IF(YP(I).GT. ZERLAT) HEIGHT(I,M) = 0.*HEIGHT(I,M)
C        END DO

C       End Modified 1/29/96 by PEM.


c	HEIGHT(1,M)=0.0
c	HEIGHT(2,M)=0.0
c	print 99,m
99	format(' surface height for m = ',i2)
c	print 88,(yp(i),height(i,m),i=1,n$)
88	format(' ',3g14.4)
 100	CONTINUE

	RETURN
	END


      SUBROUTINE GETO3

C     LOADS UP THE OZONE AND H2O PROFILES FOR USE BY THE
C     RADIATION ROUTINES -- LATER, IF OZONE OR H2O IS A TRANSPORTED CONSTITUENT
C     THIS ROUTINE IS THE TRANSFER POINT BETWEEN THE TRANSPORTED OZONE AND 
C     THE OZONE USED BY THE HEATING ROUTINES - MRS

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'

      include 'comcfg.h'

C     COMMON block renamed from /SLAKWK/ 7/14/95 by PEM to avoid name
C     conflict.

      COMMON/SLAKNEW/ TEP(M$),O3T(M$),WATR(18,23),OZZ(18,17)

cjer added 4/14/97 to use model h2o in radiation
      dimension h2otmp(n$,m$)

      MON=month

cjer iswo3, iswh2o switches (0 noninteractive) set in
c iswo3=1, interactive in ozone
c iswh2o=1, interactive in h2o
c f60_switches.dat, carried in comcfg.h


      if (iswo3.eq.0) then

c     if don't use model ozone
c     o3m is the ozone mixing ratio in ppm by volume

C     SBUV TABLES

	DO 88 K=1,17
	DO 88 J=1,18
88      OZZ(J,K)=alog(SBUV(J,K,MON))

      CALL Regrid (O3M,YP,ZZ,N$,M$,OZZ,XLSBUV,ZSBUV,18,17,0)

C ABOVE THE SBUV FIELDS ATTACH THE KRUEGER MINZER PROFILE

      DO 30 J=1,N$
      DO 20 K=1,M$
20    TEP(K)=T(J,K)
      CALL OZONEX(TEP,M$,.FALSE.,RHO,O3T)
      do 21 k=1,m$
      IF (ZSBUV(17).LT.ZZ(K)) O3M(J,K)=alog(O3T(K))
21    continue

c     smooth profile

      DO 96 II=1,2
      do 85 k=1,M$
85    TEP(K)=O3M(J,K)
      DO 86 K=2,M1$
86    O3M(J,K)=0.25*TEP(K-1)+0.5*TEP(K)+0.25*TEP(K+1)
96    CONTINUE

      DO 97 K=1,M$
97    O3M(J,K)=EXP(O3M(J,K))
      
30    CONTINUE

C      print 899,mon
C899	format(' month= ',i5)

      endif

cjer      CALL OUTA(25,NSTEP,NDAY,NYEAR,YP,ZP,O3M,N$,M$,0)
cjer      CALL OUTA(26,NSTEP,NDAY,NYEAR,WLAT,ZWAT,WATZ(1,1,mon),18,20,0)


C     GENERATE THE WATER VAPOR MIXING RATIO  (observed)

cjer this is an error      DO 180 K=1,20
      do 180 k=1,23
      DO 180 J=1,18
180   WATR(J,K)=Alog( WATZ(J,K,MON) )

      CALL Regrid (h2otmp,YP,ZZ,N$,M$,WATR,WLAT,ZWAT,18,23,0)

cjer added 4/14/97
      do j=1,n$
        do k=1,m$
          h2otmp(j,k)=exp(h2otmp(j,k))
        end do
      end do

c if not interactive, use observed h2o everywhere
      if(iswh2o.eq.0) then
          do j=1,n$
            do k=1,m$
              h2om(j,k) = h2otmp(j,k)
            end do
          end do
      endif
         
cjer don't use model h2o above 50 km

      DO K=1,M$
      DO J=1,N$
         IF ( ZP(K) .GT. 50.) THEN
            H2OM(J,K)=h2otmp(J,K)
         END IF
      END DO
      END DO


c      CALL OUTA(26,NSTEP,NDAY,NYEAR,YP,ZP,H2OM,N$,M$,0)


      RETURN
      END

      SUBROUTINE SWIT(REST)

C     SETS THE SWITCHES FOR A RUN
C     THE SWITCHES ARE NOT REMEMBERED FOR RESTARTING

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      LOGICAL REST
      COMMON /DUMP_1/ IDUMP,IFILE,ISW2(ISW$)
      CHARACTER*20 JUNKY

C      NSTEP IS STEP COUNTER
C      IPRINT IS PRINT DUMP COUNTER


C     READ IN THE SWITCHES AND PRINT OUT THEIR VALUES

      DO 1 J=1,NSWTCHS
        READ(9,100) ISW(J)
1     CONTINUE
      close(unit=9)

      DO J=1,12
         READ(40,100) LH_SW(J)
      END DO
      close(unit=40)

      DO J=1,15
         READ(44,100) KZZ_SW(J)
      END DO
      close(unit=44)

      READ(45,102) JUNKY
      READ(45,102) JUNKY

      DO J=1,15
         READ(45, 101 ) II , KYYPR(J), IYYPR(J), KZZPR(J), IZZPR(J)
     >                      , DMOLPR(J)
      END DO

100   FORMAT(I5)
101   FORMAT(I3,1X, F10.0,5X, F10.0,5X, F10.0,5X, F10.0,5X,F10.0)
102   FORMAT(A20)


      REST=.FALSE. !RESTART SWITCH
      IF (ISW(5).GT.0) REST=.TRUE.
      IF (ISW(7).GT.NCON) ISW(7)=NCON
      PRINT 1020,ISW(7)
1020  FORMAT(' NUMBER OF DYNAMICAL CONSTITUENTS TRANSPORTED = ',I5)
      
	PRINT 1000,ISW(1)
1000  FORMAT(' NUMBER OF STEPS = ',I4)

        PRINT 1001,ISW(2)
1001  FORMAT(' DUMP EVERY ',I5,' STEPS')

        IF (REST) PRINT 1002
1002  FORMAT(' CODE RESTARTED ')

      IF(ISW(3).EQ.1) PRINT 7000
 7000 FORMAT(' COOLING -> NEWTONIAN COOLING, HEATING -> STROBEL') 
      IF(ISW(3).EQ.0)  PRINT 7001
 7001 FORMAT(' COOLING AND HEATING -> ROSENFIELD RADIATION CODE') 

c      IF(ISW(4).EQ.1) PRINT 7100
 7100 FORMAT(' KZZ NOT ZERO') 
c      IF(ISW(4).EQ.0) PRINT 7101
 7101 FORMAT(' KZZ set to zero') 

c      if(isw(8).EQ.1) PRINT 7200
7200  FORMAT(' KYY NOT ZERO')
c      IF(ISW(8).EQ.0) PRINT 7210
7210  FORMAT(' KYY SET TO ZERO')

c      if(isw(9).EQ.1) PRINT 7201
7201  FORMAT(' KYZ NON ZERO')
c      IF(ISW(9).EQ.0) PRINT 7211
7211  FORMAT(' KYZ SET TO ZERO')

      if(isw(10).EQ.1) PRINT 7221
7221  FORMAT(' PERPETUAL SOLSTICE')
      IF(ISW(10).EQ.0) PRINT 7222
7222  FORMAT(' ANNUAL CYCLE')
      PRINT 7223,ISW(11)
7223  FORMAT(' W1*10 = ',I5)

c       read in the output switch array
        i26=26+ncon
        i27=26+ncon+mmax$
        read(81,231) (isw2(i),i=1,isw$)

ccc        if(ncon.gt.0) read(81,231) (isw2(i),i=27,i26)
ccc        if(mmax$.gt.0) read(81,231) (isw2(i),i=i26+1,i27)
231     format(i3)
c        print 442,(i,isw2(i),i=1,isw$)
442     format(' isw2 ',i5,i5)
        close(unit=81)

      RETURN
      END


      SUBROUTINE DERV4B(TYP,XIN,OUT,DEL,N,M,SW)

C      DO A FOURTH OR SECOND ORDER FIRST DERIVATIVE

C
C       IF TYP = 1 THEN SECOND ORDER
C       IF TYP = 2 THEN FOURTH ORDER
C

C      IF SW=1 THEN INNER INDEX
C      IF SW=2 THEN OUTER INDEX

C      BOUNDARIES ARE SET TO ZERO IF OUTER INDEX (Z)
C      BOUNDARIES ARE EXTRAPOLATED IF INNER INDEX (Y) (JTB, 9/2/93)

      DIMENSION XIN(N,M),OUT(N,M)
      INTEGER TYP,SW

      GO TO (1,2),TYP
      STOP ' WRONG TYP IN DERV4B'

2      AY1=1./(12.*DEL)
       AY2=2./(3.*DEL)
       AY3=1./(2.*DEL)
      GOTO 8
1       AY1=0.
      AY2=1./(2.*DEL)
      AY3=AY2
8      GO TO (10,20),SW
      STOP 'ILLEGAL CALL TO DERV4B'      

10      N1=N-1
      N2=N-2
      DO 51 K=1,M
      DO 50 J=3,N2
      OUT(J,K)=(-AY1*XIN(J+2,K)+AY2*XIN(J+1,K)-AY2*XIN(J-1,K)
     1+AY1*XIN(J-2,K))
50      CONTINUE
      OUT(2,K)=(XIN(3,K)-XIN(1,K))*AY3
      OUT(n1,K)=(XIN(n,K)-XIN(n2,K))*AY3
c      OUT(1,K)=0.
c      OUT(N,K)=0.
      OUT(1,K)=(XIN(2,K)-XIN(1,K))/DEL
      OUT(N,K)=(XIN(N,K)-XIN(N-1,K))/DEL 
51      CONTINUE
      RETURN

20    M1=M-1
      M2=M-2
      DO 61 J=1,N
      DO 60 K=3,M2
      OUT(J,K)=(-AY1*XIN(J,K+2)+AY2*XIN(J,K+1)-AY2*XIN(J,K-1)
     1+AY1*XIN(J,K-2))
60    CONTINUE
      OUT(J,2)=(XIN(J,3)-XIN(J,1))*AY3
      OUT(J,m1)=(XIN(J,m)-XIN(J,m2))*AY3
      OUT(J,1)=0.
      OUT(J,M)=0.
61    CONTINUE
      RETURN
      END

      SUBROUTINE CDERV4(TYP,XIN,OUT,DEL,N,M,SW)

C      DO A FOURTH OR SECOND ORDER FIRST COMPLEX DERIVATIVE

C
C       IF TYP = 1 THEN SECOND ORDER
C       IF TYP = 2 THEN FOURTH ORDER
C

C      IF SW=1 THEN INNER INDEX
C      IF SW=2 THEN OUTER INDEX

C      BOUNDARIES ARE SET TO ZERO 

      COMPLEX XIN(N,M),OUT(N,M)
      INTEGER TYP,SW

      GO TO (1,2),TYP
      STOP ' WRONG TYP IN CDERV4'

2      AY1=1./(12.*DEL)
       AY2=2./(3.*DEL)
      GOTO 8
1       AY1=0.
      AY2=1./(2.*DEL)

8      GO TO (10,20),SW
      STOP 'ILLEGAL CALL TO CDERV4'      

10      N1=N-1
      N2=N-2
      DO 51 K=1,M
      DO 50 J=3,N2
      OUT(J,K)=(-AY1*XIN(J+2,K)+AY2*XIN(J+1,K)-AY2*XIN(J-1,K)
     1+AY1*XIN(J-2,K))
50      CONTINUE
      J=2
      OUT(J,K)=(XIN(J+1,K)-XIN(J-1,K))/(2.*DEL)
      J=N1
      OUT(J,K)=(XIN(J+1,K)-XIN(J-1,K))/(2.*DEL)
      OUT(1,K)=0.
      OUT(N,K)=0.
51      CONTINUE
      RETURN
20      M1=M-1
      M2=M-2
      DO 61 J=1,N
      DO 60 K=3,M2
      OUT(J,K)=(-AY1*XIN(J,K+2)+AY2*XIN(J,K+1)-AY2*XIN(J,K-1)
     1+AY1*XIN(J,K-2))
60      CONTINUE
      K=2
      OUT(J,K)=(XIN(J,K+1)-XIN(J,K-1))/(2.*DEL)
      K=M1
      OUT(J,K)=(XIN(J,K+1)-XIN(J,K-1))/(2.*DEL)
      OUT(J,1)=0.
      OUT(J,M)=0.
61      CONTINUE
      RETURN
      END


cdbc 6/28/93 changed name of routine from zmean to zmean_sub to
cdbc distinguish it from function zmean.

	SUBROUTINE ZMEAN_sub(CIN1,CIN2,COUT,N,M)

	COMPLEX CIN1(N,M),CIN2(N,M)
	DIMENSION COUT(N,M)
	DO 200 J=1,M
	DO 200 I=1,N
        COUT(I,J) = 0.5*real(cin1(i,j)*conjg(cin2(i,j)))
 200	CONTINUE

	RETURN
	END





      SUBROUTINE GETALT2

C      GET ALTITUDE AND TEMPERATURES

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'

      TMAXx=0.
      TMINx=350.
       thmax=0

C     write(*,*) ' IN GETALT : IT2 , IT3  ',IT2,IT3
      DO 10 K=1,M$
      R7=1./TTOTH(K)
      DO 1 J=1,N$
      T(J,K)=(TH(J,K,IT3)+THG(K))*R7
      TMAXx=AMAX1(TMAXx,T(J,K))
      TMINx=AMIN1(TMINx,T(J,K))
      THMAX=AMIN1(THMAX,TH(J,K,1))
1     CONTINUE
c        write(*,*) '  in GETALT : ',  TH(3,k,it2), TH(3,k,it3), T(3,K)
10    CONTINUE
c      write(*,*) 'WITHIN GETALT: thmax, tmax ',thmax , tmaxx 

      DO 3 J=1,N$
      HS=R*T(J,1)/GZ
      ALT(J,1)=DZ/H*HS
3      CONTINUE

      DO 2 K=2,M$
      DO 2 J=1,N$
      ALT(J,K)=ALT(J,K-1)+R*DZ*(T(J,K-1)+T(J,K))*0.5/H/GZ
2     CONTINUE
      RETURN
      END


      SUBROUTINE GETHET

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INTEGER IDAYS(13)
	DATA IDAYS/1,32,60,91,121,152,182,213,244,274,305,335,365/

      CALL SOLHET 
C
C     ADD TROPOSPHERIC HEAT SOURCE FORM PRINN ET AL.

c     a duplicate section of this code exists in Radlib2 in getrad
c     The getrad section is called if the full radiation routines are 
c     used.
C
	DX=IDAYX/(IDAYS(MONTH+1)-IDAYS(MONTH))
      DO 111 J=1,N$
      DO 111 K=1,M$
  111 DUM1(J,K)=TFROP(J,K,MONTH)+(TFROP(J,K,MONTH+1)-TFROP(J,K,MONTH))
     1*DX

C      CALL GETGLM(DUM1,GHET,C,N$,M$)

      ZP00=1.5*H

	F1=1.0
	F2=1.5

      DO 210 K=1,M$
      IF (ZZ(K).GT.ZP00+1.5*H) GO TO 210
      WGT=1.0
      IF (ZZ(K).GT.ZP00) WGT=EXP(-(ZZ(K)-ZP00)/H)
      VNC2=VNC(K)
      DO 200 J=1,N$
      HEAT(J,K)=WGT*VNC2*DUM1(J,K)*F2 + (1.-WGT)*HEAT(J,K)
200   CONTINUE
 
210   CONTINUE


       CALL GETGLM(HEAT,GHET,C,N$,M$)

c     CALL OUTA(1,NSTEP,NDAY,NYEAR,YP,ZP,HEAT,N$,M$,0)


      DO 1 K=1,M$
      DO 1 J=1,N$
      HEAT(J,K)=HEAT(J,K)*TTOTH(K)
1     CONTINUE

      RETURN
      END

      SUBROUTINE GETCOL

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
C
C      NEWTONIAN COOLING ROUTINE, NOTE THAT COOLING IS POSITIVE
C


          write(*,*) '  using strobel rad. '    
  


      DO 1 K=1,M$
      DO 1 J=1,N$
      COOL(J,K)=VNC(K)*(T(J,K)-THG(K)/TTOTH(K))
1     CONTINUE

C     CALL OUTA(2,NSTEP,NDAY,NYEAR,YP,ZP,COOL,N$,M$,0)

      DO 2 K=1,M$
      DO 2 J=1,N$
      COOL(J,K)=COOL(J,K)*TTOTH(K)
2     CONTINUE


cc      CALL PRMINMAX(AX ,'AX-JSOR --',N$, M$,NDBUG)

      CALL PRMINMAX(cool ,'NEWT.-COOL',N$, M$,NDBUG)


      RETURN
      END



      SUBROUTINE START(REST,IHEATS)

C     RESTART ROUTINE

C     THIS ROUTINE ALSO SETS UP THE CONSTITUENT DISTRIBUTIONS

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
      DIMENSION OZON(M$),OZONP(MP$),XN2OIN(18,30),XLATS(18),
     1PXIN(30),ZXIN(30),XN2OPH(18,30)

      LOGICAL REST
      IF (REST) GOTO 50

C      NO RESTART

      IT1=3
      IT2=2
      IT3=1

      DO 1110 K=1,M$
         DO 1110 J=1,N$
            UB(J,K,IT1)=0.0
            UB(J,K,IT2)=0.0
            UB(J,K,IT3)=0.0
            TH(J,K,IT1)=0.0
            TH(J,K,IT2)=0.0
            TH(J,K,IT3)=0.0
            PSI(J,K)=0.0
 1110 CONTINUE


      CALL OZONEX (TG,M$,.FALSE.,RHO,OZON)

      CALL INTER (OZONP,ZPP,MP$,OZON,ZP,M$)
      IF(IPRINT.EQ.1) PRINT 102
102   FORMAT(' HEIGHT (KM)  OZONE MIXING RATIO (PPM)')
      IF(IPRINT.EQ.1) PRINT 101,(ZPP(K),OZONP(K),K=1,MP$)
101   FORMAT(' ',F10.2,G14.4)

      RETURN

50      CONTINUE
      
C     READ RESTART FILE

      CALL DUMPA(.TRUE.)
      IF (IHEATS .NE. 1  .and.  IHEATS .LE. 4) CALL RADSTART

C      **** A BIG BLOB IN THE CENTER *****
C
C      DO 20 K=1,MP$!2,M$
C      DO 20 J=1,NP$!8,12
C      X0(J,K,N)=10.*EXP(-((FLOAT(J)-10.)/3.)**2.)*
C     +          EXP(-((FLOAT(K)-20.)/8.)**2.)  !1.0
C      XC(J,K,N)=X0(J,K,N)
C      X(J,K,N)=X0(J,K,N)*RHO0(J,K)
C      X0(J,K,N)=X(J,K,N)
C20    CONTINUE

      RETURN
      END

      SUBROUTINE DUMPA(READIT)

      PARAMETER(NSWOLD=70)

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONR.INC'
      INCLUDE 'COMMONT.INC'
      INCLUDE 'COMMONW.INC'

      COMMON /SLAKO/ AY(N$,M$),CY(N$,M$),B(N$,M$),AZ(N$,M$),CZ(N$,M$)
     1,F(N$,M$)

      INTEGER ISWOLD(NSWOLD)

      LOGICAL READIT

      IF (READIT) GO TO 100

C     OVERWRITE FILE UNLESS FILE IS NOT THERE

      open(unit=13, status='OLD',ERR=99, 
     >     convert="big_endian", FORM='UNFORMATTED')

      GOTO 47
99    CONTINUE
      CLOSE(UNIT=13)

       open(unit=13, file='fort.13',status='NEW',  
     >      convert="big_endian", FORM='UNFORMATTED')

47    CONTINUE          

C     COMMONC

      WRITE (13) DZ,DY,DT,DTH,DY1,DY2,DZ1,DZ2,DELT,C00,GZ
     1 ,DYX, DZX, XNORM, DAYL, YP, YY, ZP, ZZ, RHO
     2 ,RSTAR ,h25, YPP, ZPP, XEP,XN2ST
     3 ,EZ2H ,S ,C, DTHDZ, THG, TG, TTOTH
     4 ,CF, BETA, XMM, PRNDTL, CST, XN2, zbreak
     5 ,W1,W2,ZRES,SST,CFST,THERMYN,THERMCR,RHOST



      WRITE (13) PI,R,H,CP,TOMEG,A,XKAPPA,IPRINT,NDBUG

      WRITE (13) ISW,LH_SW,KZZ_SW

      WRITE (13) KYYPR, IYYPR, KZZPR, IZZPR, DMOLPR

C     COMMOND

      WRITE (13) PSI 

      WRITE (13) AY ,CY ,B ,AZ ,CZ ,F 

      WRITE (13) UB,TH

      WRITE (13) UBY ,UBZ ,DTHY ,DTHZ 

      WRITE (13) FLX3

      WRITE (13) NSTEP,NDAY,NYEAR,DECL,IT1,IT2,IT3,MONTH,IDAYX

      WRITE (13) GCOL ,T ,COOL ,VNC 

      WRITE (13) GHET ,HEAT ,ALT, GHEATX 

      WRITE (13) PWMOM ,GWMOM ,RAL, 
     1GWMOMX ,UTHRM ,THRMW ,
     2UTHRMX,THRMWX,QUBX

      WRITE (13) KYY ,KYZ ,KZZ, PWKYY, KYY2, ZIKYY,
     1           ZIKZZ, DMOL

      WRITE (13) VERF ,HORF, EPFLX 

      WRITE (13) DUM1 ,DUM2 ,DUM3 ,D3,DUM4 ,DUMV,YPV,DUM1X,
     2DUM2X,D3X,ZPW,DUM3X,
     3DUM4X

      WRITE (13) TFROP,TEHE

      WRITE (13) IDUMP,IFILE,isw2

      WRITE (13) CHGAMA,CHTHET

      WRITE (13) UMAX,UMIN,TMAXx,TMINx,VMAX,VMIN,WMAX,WMIN

      WRITE (13) DIFZU ,DIFZT ,TMIX,
     1               UMIX         

      WRITE (13)O3M ,H2OM ,WLAT,ZWAT,WATZ,XLSBUV,ZSBUV,SBUV

      WRITE (13)QBAR, QBARY ,ubw 

      WRITE (13) glat,gavg,gvar,gvint 

C     COMMONR DATA

      WRITE (13) MZ,MZ1,MZTOT,ILBEDO

      WRITE (13) PTOT ,RHOTOT ,TTOT ,LAT ,ZPREV ,
     +              ZJRREV ,ZJRMID ,ZJR 
      WRITE (13) HEATAV ,COOLJR ,VNCJR 

C     COMMONR1 DATA

      WRITE (13)   AS ,LWCR ,PL ,PLE ,DPL ,
     1 HI ,HIS,HIE ,TL ,TLE ,TGJR,SHL ,SHLE ,SHG,
     2 CLOUD, SWALE ,SWIL ,AL ,TAUL ,RN ,TN , SRS ,STN ,
     3 TCOND ,TPENE ,TLOWL,TMIDL,THIGH,SG,SP,RSURF,RCLOUD,JALB,
     4 RHOJR ,FK,XK,NFK,MZRAD

C     COMMONT DATA

      WRITE (13)  WS ,VS 
      WRITE (13)  WRBR ,VRBR 
      WRITE (13)  WRB0 ,VRB0 
      WRITE (13)  GZPW ,GYPW 
   
      WRITE (13)  RALP ,UBX ,VFRC ,
     1THX , HEATX, DRAGX,
     2PSIX


      WRITE (13) X_D,SOR_D,X0_D,RHO0_D,RHOX_D,RHOY_D,XX_D,XY_D,XXX_D,
     .XXY_D,XYY_D,XC_D,SM_D

      WRITE (13) XC00_D,XCT_D,XMIX_D,RHO00_D ,RHOX0_D ,RHOY0_D
     .         , XN2OPHX

C    COMMONW DATA

      WRITE (13) psic ,psa , psb 

      WRITE (13) QBRYWX, MWBRK, MWRMP, QPRMY, QPRIME, QBRYW,
     1           MWVNS, GDISS, EQDMP


      WRITE (13) HEIGHT

      CLOSE(UNIT=13)
      PRINT 44,NSTEP
44    FORMAT(' RESTART FILE WRITTEN, NSTEP=',I6)
      RETURN
     
100      CONTINUE

C     COMMONC

      READ(11) DZ,DY,DT,DTH,DY1,DY2,DZ1,DZ2,DELT,C00,GZ
     1 ,DYX, DZX, XNORM, DAYL, YP, YY, ZP, ZZ, RHO
     2 ,RSTAR ,h25, YPP, ZPP, XEP,XN2ST
     3 ,EZ2H ,S ,C, DTHDZ, THG, TG, TTOTH
     4 ,CF, BETA, XMM, PRNDTL, CST, XN2, zbreak
     5 ,W1,W2,ZRES,SST,CFST,THERMYN,THERMCR,RHOST

      READ(11) PI,R,H,CP,TOMEG,A,XKAPPA,IPRINT,NDBUG


      WRITE(*,*)' !!!!! !!!!!! ???? !!!! ??  !!!! '
      WRITE(*,*)'<><><><><><><><><><><><><><><><><>'
      WRITE(*,*)'<><><><><><><><><><><><><><><><><>'
      WRITE(*,*)'<>^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^<> '
      WRITE(*,*)'<>| DID YOU                    |<>' 
      WRITE(*,*)'<>|    CHANGE # OF SWITCHES ?? |<>'
      WRITE(*,*)'<>vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv<> '
      WRITE(*,*)'<><><><><><><><><><><><><><><><><>'
      WRITE(*,*)'<><><><><><><><><><><><><><><><><>'

      READ(11) ISWOLD, LH_SW, KZZ_SW

      READ(11) KYYPR, IYYPR, KZZPR, IZZPR, DMOLPR

      WRITE(*,*) ' *---------------- >>> READ COMMONC DATA '



C    COMMOND

      READ(11) PSI 

      READ(11) AY ,CY ,B ,AZ ,CZ ,F 

      READ(11) UB,TH

      READ(11) UBY ,UBZ ,DTHY ,DTHZ 


      READ(11) FLX3

      READ(11) NSTEP,NDAY,NYEAR,DECL,IT1,IT2,IT3,MONTH,IDAYX

      READ(11) GCOL ,T ,COOL ,VNC 

      READ(11) GHET ,HEAT ,ALT, GHEATX 


      READ(11)  PWMOM ,GWMOM ,RAL, 
     1GWMOMX ,UTHRM ,THRMW ,
     2UTHRMX,THRMWX,QUBX

      READ(11) KYY ,KYZ ,KZZ, PWKYY, KYY2, ZIKYY,
     1          ZIKZZ, DMOL

      READ(11) VERF ,HORF, EPFLX

      READ(11) DUM1 ,DUM2 ,DUM3 ,D3,DUM4 ,DUMV,YPV,DUM1X,
     2DUM2X,D3X,ZPW,DUM3X,
     3DUM4X


      READ(11) TFROP,TEHE

      READ(11) IDUMP,IFILE,isw2

      READ(11) CHGAMA,CHTHET

      READ(11) UMAX,UMIN,TMAXx,TMINx,VMAX,VMIN,WMAX,WMIN

      READ(11) DIFZU ,DIFZT ,TMIX,
     1               UMIX         


      READ(11)O3M ,H2OM ,WLAT,ZWAT,WATZ,XLSBUV,ZSBUV,SBUV

      READ(11)QBAR, QBARY ,ubw 

      READ(11) glat,gavg,gvar,gvint 


      WRITE(*,*) ' *---------------- >>> READ COMMOND DATA '


C     COMMONR DATA
      READ(11) MZ,MZ1,MZTOT,ILBEDO

      READ(11) PTOT ,RHOTOT ,TTOT ,LAT ,ZPREV ,
     +              ZJRREV ,ZJRMID ,ZJR 
      READ(11) HEATAV ,COOLJR ,VNCJR 

      WRITE(*,*) ' *---------------- >>> READ COMMONR DATA '


C     COMMONR1 DATA

      READ(11)   AS ,LWCR ,PL ,PLE ,DPL ,
     1 HI ,HIS,HIE ,TL ,TLE ,TGJR,SHL ,SHLE ,SHG,
     2 CLOUD, SWALE ,SWIL ,AL ,TAUL ,RN ,TN ,
     3 SRS ,STN ,
     3 TCOND ,TPENE ,TLOWL,TMIDL,THIGH,SG,SP,RSURF,RCLOUD,JALB,
     4 RHOJR ,FK,XK,NFK,MZRAD

      WRITE(*,*) ' *---------------- >>> READ COMMONR1 DATA '


C     COMMONT DATA

      READ(11)  WS ,VS 
      READ(11)  WRBR ,VRBR 
      READ(11)  WRB0 ,VRB0 
      READ(11)  GZPW ,GYPW 

      READ(11)  RALP ,UBX ,VFRC ,
     1THX , HEATX, DRAGX,
     2PSIX

      READ(11) X_D,SOR_D,X0_D,RHO0_D,RHOX_D,RHOY_D,XX_D,XY_D,XXX_D,
     .XXY_D,XYY_D,XC_D,SM_D

      READ(11) XC00_D,XCT_D,XMIX_D,RHO00_D ,RHOX0_D ,RHOY0_D
     .         , XN2OPHX


      WRITE(*,*) ' *---------------- >>> READ COMMONT DATA '



C    COMMONW DATA

      READ(11) psic ,psa , psb 

      READ(11) QBRYWX, MWBRK, MWRMP, QPRMY, QPRIME, QBRYW,
     1          MWVNS, GDISS, EQDMP

      READ(11) HEIGHT

      WRITE(*,*) ' *---------------- >>> READ COMMONW DATA '


C     Added 8/18/95 by PEM to check restart of MWVNS

      PRINT*, " "
      PRINT*, "Planetary Wave Forcing Heights (from restart):"
      PRINT*, " "
      PRINT*, "MWVNS:"
      DO 321 IDEX=1, MMAX$
         PRINT*, MWVNS(1,IDEX), MWVNS(2,IDEX),
     1           MWVNS(3,IDEX), MWVNS(4,IDEX)
 321     CONTINUE
      PRINT*, " "

C     End Added 8/18/95




        close(unit=11)
      
      RETURN
      END            


      SUBROUTINE COMFORT
      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'

C     TEST ON MAX AND MIN VALUES FOR PRINTOUT

	gwmax=-9999.e10
	pwmax=-9999.e10
        gwmin=9999.e10
        pwmin=9999.e10
      UMAX=-10.E5
      VMAX=-10.E5
      WMAX=-10.E5
      vmin=10.e5
      WMIN=10.E5
      UMIN=10.E5

      DO 5 K=1,M$
      DO 5 J=1,N$
	gwmax=amax1(gwmax,gwmom(j,k))
	pwmax=amax1(pwmax,pwmom(j,k))
	gwmin=amin1(gwmin,gwmom(j,k))
	pwmin=amin1(pwmin,pwmom(j,k))
        UMAX=AMAX1(UMAX,UB(J,K,2))
        VMAX=AMAX1(VMAX,VS(J,K))
        WMAX=AMAX1(WMAX,WS(J,K))
        UMIN=AMIN1(UMIN,UB(J,K,2))
        VMIN=AMIN1(VMIN,VS(J,K))
        WMIN=AMIN1(WMIN,WS(J,K))
5     CONTINUE


C      PRINT 505,NSTEP,MONTH,IDAYX,NYEAR,DECL
505   FORMAT(' COMPLETE TIME STEP ',I5,' MONTH,DAY,YEAR = ',I2,'/',
     1I2,'/',I2,' DEC. ANGLE = ',F8.1)
      umax=umax*0.01
      umin=umin*0.01
C      PRINT 506,UMAX,UMIN,TMAXx,TMINx,VMAX,VMIN,WMAX,WMIN,gwmax,gwmin,
C     1pwmin,pwmax
506   FORMAT(' UMAX,UMIN = ',2f8.1,' TMAX,TMIN = ',2F8.1,
     1/,' VMAX,VMIN = ',2F8.1,' WMAX,WMIN = ',2F8.1,/,
     2' gwmin,max = ',2g12.3,' pwmin,max =',2g12.3)
C      PRINT 505,NSTEP,MONTH,IDAYX,NYEAR,DECL
      
C ------------------ FINISH TIME STEP -----------------------------

C     CALL OUTA(99,NSTEP,NDAY,NYEAR,YP,ZP,DUM2,1,1,0)

C      WRITE(6,*)'IT3, IT2, IT1  ',IT3,IT2,IT1


      RETURN
      END


      SUBROUTINE GETGLM(A,B,C,N$,M$)

      DIMENSION A(N$,M$),B(M$),C(N$)
C
C     THIS ROUTINE REMOVES THE GLOBAL MEAN FROM A AND PUTS IT IN B

C     C IS COS OF LATITUDE

      C0=0.
      DO 5 J=1,N$
5     C0=C0+C(J)
      C0=1./C0

      DO 1 K=1,M$
      B(K)=0.      
      DO 12 J=1,N$
12    B(K)=B(K)+A(J,K)*C(J)
      B(K)=B(K)*C0 
      DO 2 J=1,N$
2     A(J,K)=A(J,K)-B(K)
1     CONTINUE
      RETURN
      END    

      SUBROUTINE GETT(TG,ZP,N)
C
C     GET GLOBAL MEAN TEMPERATURE
C
      DIMENSION TG(N),ZP(N)
      REAL ttinit(117), zzinit(117)

      CALL USTA62(ZP,TG,N)

C
C  read in 1979-1980 Jan. 1 avg from NCEP data, interpolate to coupled model vertical grid
C                                    - this was written out in IDL on Linux, so use "little_endian"
C
        OPEN (78, FILE = 'tgl_jan7980.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
           READ (78) ttinit
           READ (78) zzinit
        CLOSE (78)

ccc        CALL LINTERP(ZZINIT, 117, TTINIT, ZP, N, 0, TG)

c        print *, '   '
c        print *, 'IN GETT: '
c        print *, zzinit
c        print *, ttinit
c        print *, '    '
c        print *, zp
c        print *, tg
c        print *, '    '

      RETURN
      END      


      SUBROUTINE PRODUC
C---------------------------------------------------------------------
C
C     THIS ROUTINE SETS UP CONSTANTS FOR STROBELS BAND MODEL TO COMPUTE THE
C      SOLAR HEATING. ROUTINE CALLED FROM SETCON.
C
C---------------------------------------------------------------------
      PARAMETER (NC=2)
      parameter (NW=4)

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'

          COMMON /PARM/ ABSORP(NC,NW),FSOL(NW),WAVE(NW),ENERGY(NC)
          COMMON /PARM2/ HUI1,HUI2,HUX,HUM,YLONG,YSHOR,SR1,SR2,XL,XM,XS,
     1     SRE1,SRE2,XLE,XME,XSE,SRBA,SRBB,SRBJ,SRBAE,SRBBE,SRBJE,ALBEDO

      DATA WAVE/6000.,2500.,2290.,1430./
      DATA FSOL/3.7E5,4.58E3,1.2E3,1.1/
      DATA ABSORP/0.,2.85E-21,0.,8.8E-18,6.6E-24,4.9E-18,1.0E-17,0./
      DATA ENERGY/5.12,1.050/
      DATA HUI1/.592E2/,HUI2/.401E2/,HUX/.013/,HUM/.013/YLONG/3055./
     1,YSHOR/.28E4/
      DATA SR1/3.43/,SR2/1.35/,XL/2.9E-19/,XM/1.54E-18/,XS/1.1E-17/
      DATA SRE1/.984/,SRE2/.43/,XLE/2.9E-19/,XME/1.7E-18/,
     1XSE/1.15E-17/
      DATA SRBA/.14/,SRBB/.96E9/,SRBJ/.90E-18/
      DATA SRBAE/.67/,SRBBE/3.4/,SRBJE/.24E-18/
          IF(IPRINT.EQ.1) PRINT 200
200    FORMAT (' WAVELENGTHS',5X,'FLUXES',5X,'CROSS SECTIONS FOR O,O2')
          IF(IPRINT.EQ.1)
     +    PRINT 190,(WAVE(I),FSOL(I),(ABSORP(K,I),K=1,NC),I=1,NW)
190       FORMAT (' ',F7.0,10X,3E10.3)
          IF(IPRINT.EQ.1) PRINT 170
  170      FORMAT (' DISSOCIATION ENERGIES')
          IF(IPRINT.EQ.1) PRINT 150,(ENERGY(I),I=1,NC)
  150      FORMAT (' ',1P10E12.3)
          IF(IPRINT.EQ.1) PRINT 230,HUI1,HUI2,HUX,HUM,YLONG,YSHOR
          IF(IPRINT.EQ.1) PRINT 250,SR1,SR2,XL,XM,XS

           IF(IPRINT.EQ.1) PRINT 250,SRE1,SRE2,XLE,XME,XSE
           IF(IPRINT.EQ.1) PRINT 260,SRBA,SRBB,SRBJ
           IF(IPRINT.EQ.1) PRINT 260,SRBAE,SRBBE,SRBJE
  210      FORMAT (F10.0,5E10.3)
  220      FORMAT (6E10.2)
  230      FORMAT ('0HUGGINS BAND INTENSITIES = ',2E10.3,2X,'CROSS SECTI
     1 ON =',E9.2,' WITH M= ',E9.2,',YLONG = ',F7.1,/,'YSHORT = ',G12.2
     2)
  250      FORMAT ('0SCHUMANN RUNGE CONT. INTEN = ',1PE9.2,2X,E9.2,2X,
     1   'CROSS   SECTIONS = ',E9.2,2X,E9.2,2X,E9.2)
260   FORMAT('0 SCHUMANN RUNGE BAND COEF =',2E10.2,'OPTICALLY THIN RATE
     1 = ',E9.2)
           RETURN
      END
C
           SUBROUTINE SOLHET
C---------------------------------------------------------------------
C
C     THIS ROUTINE USES STROBELS BAND MODEL TO COMPUTE THE
C      SOLAR HEATING. ROUTINE CALLED FROM GETHET.
C
C---------------------------------------------------------------------
      PARAMETER (NC=2)
      parameter (NW=4)

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'

      COMMON /SOLHWK/Y1(M$),DTAB(M$,NC),
     1 PS(NC,NW),PT(NC),CSOL(NC),THETA,TAU,SIG,HSOL(NC),HE(M$),HB(M$)
     2 ,DNZ(NC)
           COMMON /PARM/ ABSORP(NC,NW),FSOL(NW),WAVE(NW),ENERGY(NC)
           COMMON /PARM2/HUI1,HUI2,HUX,HUM,YLONG,YSHOR,SR1,SR2,XL,XM,XS,
     1     SRE1,SRE2,XLE,XME,XSE,SRBA,SRBB,SRBJ,SRBAE,SRBBE,SRBJE,ALBEDO

      DATA RMA/4.79E-23/,BC/1.3803E-16/,CP0/1.1643E2/
     1,OOP/1.74533E-2/,OME/7.2722E-5/,PCGS/1.01325E6/

      INET=1
      ITOT=0

C     LOOP OVER LATITUDE POINTS

      DO 60 LL =1,N$
      HEAT(LL,1) = 0.0
      DO 65 LU=1,M$
C      Y1(LU) = T(LL,LU)
      DTAB(LU,2) = O3M(LL,LU)*0.2843e-6*2.68E19*1013.25*RHO(LU)/T(LL,LU)
65    DTAB(LU,1) = 0.21*PCGS *RHO(LU)/(BC*T(LL,LU))  

C      CALL OZONEX(Y1,M$,.TRUE.,RHO,DTAB(1,2))

      XLAT = YP(LL)
           XA=SIN(XLAT*OOP)*SIN(DECL*OOP)
           XB=COS(XLAT*OOP)*COS(DECL*OOP)
           TEST=-XA/XB
           IF (ABS(TEST).GE.0.999) GO TO 11
           TX=ACOS(TEST)/OME
           GO TO 21
   11      IF (-TEST.GE.0.999) TX=4.32E4
           IF (-TEST.LE.-0.999) GO TO 63
   21      XMU=2.*(XA*TX+XB*SIN(TX*OME)/OME)
           XMU2=2.*(XA*XA*TX+2.*XA*XB*SIN(OME*TX)/OME+XB*XB
     1*(TX/2.+SIN(OME*TX*2.)/(4.*OME)))
      AVE = XMU*SQRT(2.*TX/XMU2)/8.64E4
           THETA=ACOS(SQRT(XMU2/(2.*TX)))/OOP
           ZANG=THETA/57.2957
      IF(ABS(THETA).GT.89.2) GO TO 63
      SECT = 1.0/COS(ZANG)
      CSOL(1) = SECT*DTAB(M$,1)*BC*T(LL,M$)/(GZ*5.31E-23)
      CSOL(2) = SECT*DTAB(M$,2)*BC*T(LL,M$)/(GZ*7.97E-23)
      DO 61 LZ =1, M1$
      LU = M$-LZ+1
      J=1
           DO 30 I=1,NC
      DNZ(I) = DTAB(LU,I)
           PT(I)=0.0
   30      PS(I,J)=0.0
           DO 90 I=1,NC
           DO 80 J=1,NW
           D=0.0
           WE=1.240E4/WAVE(J)
           DELE=(WE-ENERGY(I))/WE
           DO 70 L=1,NC
   70      D=D+ABSORP(L,J)*CSOL(L)
           FZ=FSOL(J)*EXP(-AMIN1(D,100.))
           PW=FZ*ABSORP(I,J)*DNZ(I)
           PS(I,1)=PS(I,1)+DELE*PW
   80      PT(I)=PT(I)+PW
   90      CONTINUE
           TLONG=HUX*CSOL(2)*EXP(-HUM*YLONG)
           TSHOR=HUX*CSOL(2)*EXP(-HUM*YSHOR)
      EDEL=0.98/(1.+1.2E14/((t(ll,lu)/300.)**0.8*DNZ(1)))
      EHA=(2.39 + 1.63*(7.5E-11*EXP(107./t(ll,lu)) + 2.9E-11*EXP(67./
     1t(ll,lu))/(1.+1.06E13/DNZ(1)))/(7.5E-11*EXP(107./t(ll,lu)) +
     1 2.9E-11*EXP(67./t(ll,lu))) +EDEL)/5.
       QHU=DNZ(2)*(HUI1+(HUI2-HUI1)*EXP(-AMIN1(TLONG,100.))-HUI2*EXP
     1     (-AMIN1(TSHOR,100.)))/(HUM*CSOL(2))
       QSR=(SR1*EXP(-AMIN1(XL*CSOL(1),100.))+(SR2-SR1)*
     1     EXP(-AMIN1(XM*CSOL(1),100.))
     2     -SR2*EXP(-AMIN1(XS*CSOL(1),100.)))*DNZ(1)/CSOL(1)
          QSRE=(SRE1*EXP(-AMIN1(XLE*CSOL(1),100.))+(SRE2-SRE1)*
     1     EXP(-AMIN1(XME*CSOL(1),100.))-
     2     SRE2*EXP(-AMIN1(XSE*CSOL(1),100.)))*DNZ(1)/CSOL(1)
          QSRB=DNZ(1)/(SRBA*CSOL(1)+SRBB*SQRT(CSOL(1)))
           IF (CSOL(1).LT.1.0E18) QSRB=DNZ(1)*SRBJ
          QSRBE=DNZ(1)/(SRBAE*CSOL(1)+SRBBE*SQRT(CSOL(1)))
           IF (CSOL(1).LT.1.0E18) QSRBE=DNZ(1)*SRBJE
      QSCAT=0.3*FSOL(1)*ABSORP(2,1)*DNZ(2)*EXP(-AMIN1(ABSORP(2,1)*
     1CSOL(2),100.))
      QHA=FSOL(2)*ABSORP(2,2)*DNZ(2)*EXP(-AMIN1(ABSORP(2,2)
     1     *CSOL(2),100.))
           PT(1)=PT(1)+QSR+QSRB
           PS(1,1)=PS(1,1)+QSRE+QSRB*ITOT + QSRBE*INET
      PT(2)=PT(2)+QHU*(3.02 +EDEL)/4. + QSCAT + QHA*(EHA-1.)
           PS(2,1)=PS(2,1)+QHU*0.75+QSCAT*0.50
      CPD = CP0*RHO(LU)*PCGS/(BC*T(LL,LU))*RMA
      HEAT(LL,LU) =(HB(LU) + (PS(1,1)+PT(2))*AVE/CPD)/DAYL
      DO 68 I=1,NC
      HSOL(I) = (ALT(LL,LU)-ALT(LL,LU-1))/(ALOG(DTAB(LU-1,I)/
     1DTAB(LU,I))+1.0E-5)
68    CSOL(I) = CSOL(I)+SECT*DTAB(LU-1,I)*HSOL(I)*
     1      (1.0-EXP(-(ALT(LL,LU)-ALT(LL,LU-1))/HSOL(I)))
61    CONTINUE
      GO TO 60
63    CONTINUE
      SECT=72.0
      DO 62 LU=1,M$
62    HEAT(LL,LU) = 0.0
60    CONTINUE
           RETURN
           END


      SUBROUTINE OZONEX(T,I,ILP,PTOT,O3MIX)

C THIS ROUTINE COMPUTES OZONE MIXING RATIOS
C THE STANDARD KRUEGER-MINZNER MODEL IS USED THROUGHOUT
C EXCEPT SLIGHTLY LOWER VALUES USED NEAR THE STRATOPAUSE
C
C T IS TEMPERATURE ARRAY
C I IS SIZE OF ARRAY
C PTOT IS PRESSURE ARRAY IN ATMOSPHERES
C IF ILP IS TRUE, O3MIX IS NO.DENS. IF FALSE IT IS mixing ratio in ppm.

      PARAMETER (l$= 50)
      LOGICAL ILP
ccelf       DIMENSION T(I),PTOT(I),O3MIX(I),PMB(l$),PMBPT(38),
       DIMENSION T(I),PTOT(I),O3MIX(I),PMB(I),PMBPT(38),
     *O3RAT(38)
      DATA PLAW/0.565/
      DATA PMBPT/1013.,795.,616.,472.,356.,264.,193.,141.,103.,75.,
     *54.7,40.,29.3,21.5,15.9,11.7,8.67,6.46,4.84,3.65,2.77,2.12,1.63,
     *1.26,.976,.758,.589,.457,.353,.272,.208,.159,.12,.0899,.0667,
     *.0489,.0354,.0254/
      DATA O3RAT/.032,.033,.034,.041,.06,.132,.31,.50,.85,1.44,2.70,
     14.85,7.50,9.35,10.2,10.5,10.22,9.33,8.09,7.84,7.3,5.34,5.08,
     14.12,3.38,
     *3.11,2.29,1.92,1.56,1.36,1.13,.96,.82,.70,.59,.50,.41,.34/
      NCUR=1
      DO  1  N=1,I
  1   PMB(N)=1013.25*PTOT(N)
      DO  2  N=1,I
      IF (PMB(N) .GE. PMBPT(1)) GO TO 3
      IF (PMB(N) .LE. PMBPT(38)) GO TO 4
       DO 6 J=1,37
       IF(PMBPT(J).GE.PMB(N).AND.PMBPT(J+1).LT.PMB(N)) NCUR=J
6     CONTINUE

C LINEAR INTERPOLATION

      O3MIX(N)=O3RAT(NCUR)+((O3RAT(NCUR+1)-O3RAT(NCUR))/(PMBPT(NCUR+1)-
     *PMBPT(NCUR)))*(PMB(N)-PMBPT(NCUR))
      GO TO 2
  3   O3MIX(N)=O3RAT(1)
      GO TO 2
  4   O3MIX(N)=O3RAT(38)*(PMB(N)/PMBPT(38))**PLAW
  2   CONTINUE
      IF (.not.ILP) return
      DO 66 N=1,I
66    O3MIX(N)=2.843E-7*O3MIX(N)*2.68E19*PMB(N)/T(N)
      RETURN
      END

           SUBROUTINE INTER (A,B,M,Y,X,N)
C
C     THIS ROUTINE INTERPOLATES BETWEEN  TWO LINEAR ARRAYS.
C     THE INPUT ARRAY IS X,Y AND B. THE OUTPUT ARRAY IS A. A THIRD ORDER
C     POLYNOMIAL FIT IS USED BETWEEN ARRAY POINTS. OUTSIDE REGION
C     v:A LINEAR INTERPOLATION IS USED.
C
           DIMENSION A(M),B(M),X(N),Y(N),Q(4),R(4),S(3),T(2),E(3)
           DO 90 J=1,M
           IF (B(J).GE.X(N)) GO TO 70
           IF (B(J).LE.X(1)) GO TO 80
           DO 10 K=1,N
           IF (B(J).LE.X(K)) GO TO 20
   10      CONTINUE
   20      LQ=K-2
           IF (LQ.LE.1) LQ=1
           IR=N-3
           IF (LQ.GE.IR) LQ=N-3
           DO 30 L=1,4
           LN=L+LQ-1
           Q(L)=X(LN)
   30      R(L)=Y(LN)
           DO 40 L=1,3
   40      S(L)=(R(L+1)-R(L))/(Q(L+1)-Q(L))
           DO 50 L=1,2
   50      T(L)=(S(L+1)-S(L))/(Q(L+2)-Q(L))
           D=(T(2)-T(1))/(Q(4)-Q(1))
           DO 60 L=1,3
   60      E(L)=B(J)-Q(L)
           A(J)=R(1)+S(1)*E(1)+T(1)*E(1)*E(2)+D*E(1)*E(2)*E(3)
           GO TO 90
   70      A(J)=(Y(N)-Y(N-1))/(X(N)-X(N-1))*(B(J)-X(N-1))+Y(N-1)
           GO TO 90
   80      A(J)=(Y(2)-Y(1))/(X(2)-X(1))*(B(J)-X(1))+Y(1)
   90      CONTINUE
           RETURN
           END

           Subroutine Regrid (AOUT,XNEW,YNEW,NXNEW,NYNEW,AIN,
     1XOLD,YOLD,NXOLD,NYOLD,SW)

        PARAMETER(MAXDIM=100)

C       2-D/ 1-D INTERPOLATOR

C	THE INPUT ARRAY IS AIN(NXOLD,NYOLD) WITH CORRESPONDING 
C       VECTORS X(NXOLD) AND Y(NYOLD). THE INTERPOLATED ARRAY IS
C       AOUT(NXNEW,NYNEW) WITH VECTORS XNEW(NXNEW) YNEW(NYNEW)
C       SW = 0 DO FULL INTERPOLATION
C       SW = 1 DO X ONLY INTERPOLATION (YNEW=YOLD)
C       SW = 2 DO Y ONLY INTERPOLATION (XNEW=XOLD)


C	ARRAY SIZES MUST BE LESS THAN MAXDIM

C       M. R. SCHOEBERL



           DIMENSION AOUT(NXNEW,NYNEW),XNEW(NXNEW),YNEW(NYNEW),
     1AIN(NXOLD,NYOLD),XOLD(NXOLD),YOLD(NYOLD)
           DIMENSION TEMP(MAXDIM,MAXDIM),
     2         DX3(MAXDIM),DX4(MAXDIM)

C     TEMP MUST BE AS BIG AS NXNEW,NYOLD

C      DO THE X INTERPOLATION FIRST

       IF (NXNEW.GT.MAXDIM) STOP ' INSUFFICENT X MEMORY IN REGRID'
       IF (NYold.GT.MAXDIM) STOP ' INSUFFICENT Y MEMORY IN REGRID'
       IF (SW.EQ.1.AND.NYNEW.NE.NYOLD) STOP 'Y LENGTH ERROR IN REGRID'
       IF (SW.EQ.2.AND.NXNEW.NE.NXOLD) STOP 'X LENGTH ERROR IN REGRID'

       IF (SW.EQ.2) GO TO 65

C      DO THE X INTERPOLATION FIRST

       DO 401 K=1,NYold
      IF (SW.EQ.0) CALL INTER(TEMP(1,K),XNEW,NXNEW,AIN(1,K),XOLD,NXOLD)
      IF (SW.EQ.1) CALL INTER(AOUT(1,K),XNEW,NXNEW,AIN(1,K),XOLD,NXOLD)
401    CONTINUE
       IF (SW.EQ.1) RETURN

65     CONTINUE

C      DO THE Y INTERPOLATION

       DO 402 J=1,NXNEW
       DO 4020 K=1,NYOLD
       IF (SW.EQ.0) DX3(K)=TEMP(J,K)
       IF (SW.EQ.1) DX3(K)=AIN(J,K)
4020   CONTINUE

       CALL INTER(DX4,YNEW,NYNEW,DX3,YOLD,NYOLD)

       DO 403 K=1,NYNEW
403    AOUT(J,K)=DX4(K)
402    CONTINUE

       RETURN
       END

	SUBROUTINE USTA62(ZZ,T,N)
	DIMENSION ZZ(N),T(N)
	COMMON /USTAWK/ SM(100)
C
C     COMPUTES THE US STANDARD ATMOS TEMP PROFILE
C
C     Z IS A HEIGHT ARRAY IN KM
C     T IS OUTPUT TEMPERATURE ARRAY IN K
C
C
	IF(N.GT.401) STOP 'INCREASE SM IN USTA62'
	DO 1 K=1,N
	Z=ZZ(K)
      IF(Z.LE.11.0) T(K)=288.15-6.5*Z
      IF(Z.LE.20.0.AND.Z.GE.11.0) T(K)=216.65
      IF(Z.LE.32.0.AND.Z.GE.20.0) T(K)=216.65+1.0000*(Z-20.)
      IF(Z.LE.47.0.AND.Z.GE.32.0) T(K)=228.65+2.800*(Z-32.0)
      IF(Z.LE.52.0.AND.Z.GE.47.0) T(K)=270.65
      IF(Z.LE.61.0.AND.Z.GE.52.0) T(K)=270.65-2.0*(Z-52.0)
      IF(Z.LE.79.0.AND.Z.GE.61.0) T(K)=252.65-4.0*(Z-61.0)
      IF (Z .LE. 90. .AND. Z .GE. 79.0)  T(K)=180.65
      IF (Z .LE. 100. .AND. Z .GE. 90.)  T(K)=180.65+2.937*(Z-90.)
      IF (Z .LE. 110. .AND. Z .GE. 100.) T(K)=210.02+4.698*(Z-100.)
      IF (Z .LE. 140. .AND. Z .GE. 110.) T(K)=257.00+15.00*(Z-110.)
1	CONTINUE
	N1=N-1
	DO 4 I=1,4 !NUMBER OF SMOOTH PASSES
	DO 2 K=2,N1
	SM(K)=0.25*(T(K-1)+T(K+1))+0.5*T(K)
2	CONTINUE
	DO 3 K=2,N1
3	T(K)=SM(K)
4	CONTINUE
      RETURN
      END

      SUBROUTINE SMTHR3( XIN ,N,M )

      INCLUDE 'PARAM.INC'
      
      REAL WK(NPP$,MPP$),XIN(N,M)


      DO 1 K=2,M-1
      DO 1 J=2,N-1          
          WK(J,K)=
     >    (XIN(J+1,K)+XIN(J-1,K)+XIN(J,K+1)+XIN(J,K-1))*0.125
     >   + XIN(J,K)*0.5
 1    CONTINUE

      DO K=2,M-1
           WK(1,K)=
     >     (XIN(1,K+1)+XIN(1,K-1)+XIN(2,K))/6.
     >    + XIN(1,K)*0.5      

           WK(N,K)=
     >     (XIN(N,K+1)+XIN(N,K-1)+XIN(N-1,K))/6.
     >    + XIN(N,K)*0.5      
      END DO

      DO J=2,N-1
           WK(J,1)=
     >     (XIN(J+1,1)+XIN(J-1,1)+XIN(J,2))/6.
     >    + XIN(J,1)*0.5     

           WK(J,M)=
     >     (XIN(J+1,M)+XIN(J-1,M)+XIN(J,M-1))/6.
     >    + XIN(J,M)*0.5     
      END DO

      WK(1,1)=(XIN(1,2)+XIN(2,1))/4.+XIN(1,1)*0.5
      WK(1,M)=(XIN(1,M-1)+XIN(2,M))/4.+XIN(1,M)*0.5
      WK(N,M)=(XIN(N,M-1)+XIN(N-1,M))/4.+XIN(N,M)*0.5
      WK(N,1)=(XIN(N,2)+XIN(N-1,1))/4.+XIN(N,1)*0.5

      DO 2 K=1,M
      DO 2 J=1,N
         XIN(J,K)=WK(J,K)
 2    CONTINUE
     

       
      RETURN
      END

      
      SUBROUTINE UTFRMX     
C                         !(XC,XC00,UB,TH)

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
      
ccelf      include 'test_value.inc'
ccelf      logical nantest

C     DUM1X already declared in COMMOND.INC
C      DIMENSION DUM1X(NP$,MP$)



C--------------------------- TH FROM TOTAL POTENTIAL ADV.
      DO K=1,MP$
      DO J=1,NP$
          THX(J,K)=XC(J,K,1)-XC00(J,K,1)
      END DO
      END DO

c      CALL OUTA(191,NSTEP,NDAY,NYEAR,YPP,ZPP,THX,NP$,MP$,0)

      DO K=1,M$
      DO J=1,N$
          TH(J,K,IT1)= (
     >      THX(J,K)   + THX(J+1,K) 
     >    + THX(J,K+1) + THX(J+1,K+1)
     >                  )/4.
      END DO
      END DO
      

C---------------------------- UB FROM ZONAL MOMENTUM 
      DO K=1,MP$
      DO J=1,NP$
c          DUM1X(J,K)=XC(J,K,2)                ! UB
          DUM1X(J,K)=XC(J,K,2)-XC00(J,K,2)    ! M

c      IF (DUM1X(J, K) .NE. DUM1X(J, K)) THEN
c       nantest=test_nan( dum1x(j,k) )
c       if(nantest) then
c         PRINT*, "NAN in UTFRMX:"
c         PRINT*, "     XC = ",   XC(J, K, 2)
c         PRINT*, "   XC00 = ", XC00(J, K, 2)
c         STOP
c         ENDIF

      IF (DUM1X(J,K) .GT. 1.0E+14) THEN
         PRINT*, "DUM1X blowup in UTFRMX:"
         PRINT*, "     XC(1) = ",   XC(J, K, 1)
         PRINT*, "   XC00(1) = ", XC00(J, K, 1)
         PRINT*, "     XC(2) = ",   XC(J, K, 2)
         PRINT*, "   XC00(2) = ", XC00(J, K, 2)
         STOP
         ENDIF

      END DO
      END DO

      DO K=1,MP$
      DO J=1,NP$
c          UBX(J,K)=DUM1X(J,K)               ! UB
          UBX(J,K)=DUM1X(J,K)/CST(J)        ! M

c      IF (UBX(J, K) .NE. UBX(J, K)) THEN
c       nantest=test_nan( ubx(j,k) )
c       if(nantest) then
c         PRINT*, "NAN in UTFRMX:"
c         PRINT*, "   DUM1X = ", DUM1X(J, K)
c         PRINT*, "     CST = ", CST(J)
c         STOP
c         ENDIF

      IF (ABS(UBX(J, K)) .GE. 1.0E+15) THEN
         PRINT*, "UBX blowup in UTFRMX:"
         PRINT*, "   DUM1X = ", DUM1X(J, K)
         PRINT*, "     CST = ", CST(J)
         STOP
         ENDIF

      END DO
      END DO
C
C  divide by cos(+-85) instead of cos(+-87.5) above -  DON'T DO w/ adjustment below
C
ccpw17      DO K=1,MP$
ccpw17          J=1
ccpw17          UBX(J,K)=DUM1X(J,K)/C(J) 
ccpw17          J=NP$
ccpw17          UBX(J,K)=DUM1X(J,K)/C(J-1) 
ccpw17      END DO


C
C  Now adjust high latitude UBAR is NOT LARGER than the value required to ensure the
C  relative angular velocity is constant, so that UBAR -> 0 at poles  (EF, April 2008)
C  if the value is LESS than that of the relative angular velocity, DO NOT ADJUST
C  do the same in both hemispheres, THIS IS NEEDED - EF, April, 2009
C     UBX(NP$,MP$);  CST(NP$) -   do poleward of 80N, 80S (Oct 2010) - wconv89
C
      DO 770 K=1,MP$
         DO 771 J=1,NP$
            if (ypp(J) .gt. 80.) then
               ubang = UBX(J-1,K)*CST(J)/CST(J-1)
ccpw9               UBX(J,K) = ubang
               if (ABS(UBX(J,K)) .gt. ABS(ubang)) UBX(J,K) = ubang
            endif
 771     CONTINUE

         DO 772 J=NP$,1,-1
            if (ypp(J) .lt. -80.) then
               ubang = UBX(J+1,K)*CST(J)/CST(J+1)
ccpw9               UBX(J,K) = ubang
               if (ABS(UBX(J,K)) .gt. ABS(ubang)) UBX(J,K) = ubang
            endif
 772     CONTINUE
 770   CONTINUE



      DO K=1,MP$
      DO J=1,NP$
          QUBX(J,K)=UBX(J,K)         
      END DO
      END DO

c      CALL OUTA(192,NSTEP,NDAY,NYEAR,YPP,ZPP,UBX,NP$,MP$,0)

      DO K=1,M$
      DO J=1,N$
          UB(J,K,IT1)= (                 ! PUT ADVECTED ZONAL
     >      UBX(J,K)   + UBX(J+1,K)      ! MOMENTUM INTO ARRAY
     >    + UBX(J,K+1) + UBX(J+1,K+1)    ! UB ON MRS GRID
     >                  )/4.             !
      END DO
      END DO


C---------------------------- VSTAR FROM XC(J,K,3) 
      DO K=1,M$
      DO J=1,N$
c          VS(J,K)= (
c     >      XC(J,K,3)   + XC(J+1,K,3) 
c     >    + XC(J,K+1,3) + XC(J+1,K+1,3)
c     >                  )/4.
      END DO
      END DO


cccc          write (1561) N$, M$
cccc          write (1561) yp, ypp, zp, zpp
cccc          write (1561) dum1x, ubx, ub, cst, c


      RETURN
      END


      SUBROUTINE MRS2RBR(AMRS,ARBR)
      INCLUDE 'PARAM.INC'
      DIMENSION AMRS(N$,M$) , ARBR(NP$,MP$)


      DO K=2,M$
      DO J=2,N$    
         ARBR(J,K)=
     >            ( AMRS(J-1,K-1)+AMRS(J,K-1)+AMRS(J-1,K)
     >             +AMRS(J,K) )/4.
      END DO
      END DO
      DO J=2,N$ 
         ARBR(J,MP$)=(AMRS(J,M$)+AMRS(J-1,M$))/2.
         ARBR(J,1)  =(AMRS(J,1)+AMRS(J-1,1))/2.
      END DO
      DO K=2,M$ 
         ARBR(NP$,K)=(AMRS(N$,K-1)+AMRS(N$,K))/2.
         ARBR(1,K)  =(AMRS(1,K-1)+AMRS(1,K))/2.
      END DO
      ARBR(1,1)=AMRS(1,1)
      ARBR(1,MP$)=AMRS(1,M$)
      ARBR(NP$,1)=AMRS(N$,1)
      ARBR(NP$,MP$)=AMRS(N$,M$)

      RETURN
      END

      SUBROUTINE GETQB

C     COMPUTE Q, THE QG POTENTIAL VORTICITY
C     FROM WINDS

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
      INCLUDE 'COMMONW.INC'

      DO K=1,MP$     
      DO J=1,NP$     
         IF(ISW(30).EQ.0) DUM1X(J,K)= UTHRMX(J,K)
         IF(ISW(30).EQ.1) DUM1X(J,K)= QUBX(J,K)
      END DO
      END DO


CC      CALL DERV4B(1,UTHRMX,DUM2X,DY,NP$,MP$,1)
 
      DO K=1,MP$     
      DO J=2,N$     
         DUM2X(J,K)= ( DUM1X(J+1,K)-DUM1X(J-1,K) )*DY1
      END DO
      END DO
      DO K=1,MP$     
         J=1
         DUM2X(J,K)= ( DUM1X(J+1,K)-DUM1X(J,K) )*DY1
         J=NP$
         DUM2X(J,K)= ( DUM1X(J,K)-DUM1X(J-1,K) )*DY1
      END DO



      CALL PRMINMAX(UBX,  'UBX-QBY---',NP$, MP$,1)
      CALL PRMINMAX(DUM2X,'DUM-QBY---',NP$, MP$,1)
      CALL PRMINMAX(CFST ,'CFS-QBY---',NP$,  1 ,1)

 
      DO K=1,MP$     
      DO J=1,NP$     
    
         QBAR(J,K)=  CFST(J) - DUM2X(J,K)

      END DO
      END DO



C***************************************************************
C***   MATSUNO'S MODIFIED QG PLANETARY POTENTIAL VORT.     *****
C***      (see Andrews, Holton, Leovy [Eq. 5.3.4] )        *****
C***************************************************************

      CALL DERV4B(1,DUM1X,DUM2X,DZ,NP$,MP$,2)   ! D( U )/Dz

      DO K=1,MP$
      DO J=1,NP$
         
         DUM2X(J,K)=( RHOST(K)*CFST(J)*CFST(J)/XN2ST(K)  )*
     2               DUM2X(J,K)

      END DO
      END DO

      CALL DERV4B(1,DUM2X,DUM3X,DZ,NP$,MP$,2)   !  D( f^2*Rho* (DU/Dz) / N^2 )/Dz


      DO K=1,MP$
      DO J=1,NP$
  
         DUM3X(J,K)=DUM3X(J,K)/RHOST(K)

      END DO
      END DO

      DO K=1,MP$
      DO J=1,NP$
  
         DUM2X(J,K)=DUM1X(J,K)*CST(J)          !  U*COS(y)

      END DO
      END DO

      
      CALL DERV4B(1,DUM2X,DUM4X,DY,NP$,MP$,1)   !  D( U*COS(y) )/Dy
      DO K=1,MP$
         J=1                                                !
         DUM4X(J,K)= ( DUM2X(J+1,K) - 0.0000 ) / ( 2.*DY )  !  NO ZONAL WIND
         J=NP$                                              !  AT EDGES 
         DUM4X(J,K)= ( 0.0000 - DUM2X(J-1,K) ) / ( 2.*DY )  !
      END DO






      DO K=1,MP$
      DO J=1,NP$
  
         DUM2X(J,K)=DUM4X(J,K)/CST(J)

      END DO
      END DO

      CALL DERV4B(1,DUM2X,DUM4X,DY,NP$,MP$,1)   !  D( D( U*COS(y) )/Dy/COS(y) )/Dy
      DO K=1,MP$
         J=1                                   !
         DUM4X(J,K)= 0.0000                    !  NO CURVATURE AT -/+DY FROM EDGE
         J=NP$                                 !
         DUM4X(J,K)= 0.0000                    !
      END DO
      

      DO K=1,MP$
      DO J=1,NP$

         QBRYWX(J,K)=
     2                TOMEG*CST(J)/A   
     3           -    DUM4X(J,K)   
     4           -    DUM3X(J,K)
      END DO
      END DO     

c      DO K=1,MP$
c         J=1
c         QBRYWX(J,K)=
c     2                TOMEG*CST(J)/A   
c     3           -    DUM4X(J+1,K)   
c     4           -    DUM3X(J+1,K)
c         J=NP$
c         QBRYWX(J,K)=
c     2                TOMEG*CST(J)/A   
c     3           -    DUM4X(J-1,K)   
c     4           -    DUM3X(J-1,K)
c      END DO

C*************************************************
C**  ------- QBRYWX IS MATSUNO QG-PV  -------   **
C*************************************************



      DO K=1,MP$     
      DO J=2,N$     
         QBARY(J,K)= ( QBAR(J+1,K)-QBAR(J-1,K) )*DY1
      END DO
      END DO
      DO K=1,MP$     
         J=1
         QBARY(J,K)= ( QBAR(J+1,K)-QBAR(J,K) )/DY
         J=NP$
         QBARY(J,K)= ( QBAR(J,K)-QBAR(J-1,K) )/DY
      END DO



      RETURN
      END      





