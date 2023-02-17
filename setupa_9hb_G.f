	SUBROUTINE SETUPA(DTIMEB)

        include "com2d.h"
        include "com_aerd.h"
        include "com_aerg.h"


      DOUBLE PRECISION DTIMEB(L$,NTIME)
      REAL PMTIME(NTIME),DLITE(NHT)
	COMMON/CHIDSFAC/CHIDA(18,L$),SFACA(18,L$)
      DIMENSION gcrout(L$,Z$)

      REAL*8 sjfluxday, hjfluxday, flxa


c   solar flux data from Judith Lean:
C   index 1 = time;  2-40 are the 39 model wavelength bins (in ph/cm^2/sec); 41=TOA flux (W/m2)
C      SJFLUXDAY(41) are TD values interpolated to the current day 
C
C    HJFLUXDAY(12) are for the heating rate intervals (in mW/m2), interpolated to current day 
C
      COMMON/CSJFLUXD/ SJFLUXDAY(41), HJFLUXDAY(12), FLXA(41), ISOLCYC


      SAVE

C     DECLINATION, ZENITH ANGLE, 
C     CHAPMAN FUNCTION AT LATITUDE LATP AND ELAPSED TIME T.
C         RESET DECLINATION WITH TIME IF REQUESTED. FORMULA AND
C         CONSTANTS ARE FROM THE ASTRONOMICAL ALMANAC, WITH
C         MEANL=MEAN SOLAR LONGITUDE (CORRECTED FOR ABERRATION),
C         MEANA=MEAN ANOMALY, ECLIPL=ECLIPTIC LONGITUDE.
          MEANL=DEC1+DEC2*DAY
          MEANA=DEC3+DEC4*DAY
          ECLIPL=MEANL+DEC5*SIN(MEANA)+DEC6*SIN(2.*MEANA)
          SIND=DEC7*SIN(ECLIPL)

          R=ASIN(SIND)
c R is the solar Declination 
	  DECLIN=R

          COSD=COS(R)
          DECD=R*RTD
c          print *,' dec1,dec2,dec3,dec4,dec5,dec6,meanl,meana,
c     *  eclipl,dec7,dec8,rtd,r',dec1,dec2,dec3,dec4,dec5,
c     *  dec6,meanl,meana,eclipl,dec7,dec8,rtd,r
c          print *,' day,decd,sind,cosd=',day,decd,sind,cosd

C           CORRECT SOLAR FLUXES AT THE TOP OF THE ATMOSPHERE FOR
C           THE CURRENT ORBITAL POSITION OF THE EARTH.
            CORB=1.0/(1.0-DEC8*COS(MEANA))**2


C   use solar flux from J. Lean here;  SJFLUXDAY(41) are for the current day
C                                
            DO 200 JL=1,IL$
             SOLCYC=1.0
            FLUX(JL)=CORB*SOLCYC*SJFLUXDAY(JL+1)
ccccccc            FLUX(JL)=CORB*SOLCYC*FLUX0(JL)
200	continue

C
C
c  select solar max gcr's (or solar min), interpolate to new model grid, GCR(L$,Z$), gcrout(L$,Z$)
C        LATIN(18), ZZ46(46),  LAT4(L$), ZALT90(Z$) are REAL*4 in COMMON, as are GCRMIN(18,46),GCRMAX(18,46)
C 
            CALL BINTERP(LATIN, 18, zz46, 46, GCRMAX, 
     >                   LAT4, L$, ZALT90, Z$, 0, 0, GCROUT)

	DO 2003 IK=1,Z$
	DO 2003 IJ=1,L$
		GCR(IJ,IK)=GCROUT(IJ,IK)
2003	CONTINUE

c  zero values above 60 km

      DO 110 IK=1,Z$
         IF (zalt90(ik) .GT. 60.) THEN
                DO 111 IJ=1,L$
 111               GCR(ij,ik) = 0.0
         ENDIF
 110  CONTINUE


	DO 3000 IJ=1,L$
	CALL FILL(IJ,Z$)

      B=TAN(DECLIN)*TAN(LATP*DTR)
      TEST1= -B

c	print *,' b declin ij latp dtr tanr tanl test1',
c     *  b,declin,ij,latp,dtr,tan(r),tan(latp*dtr),test1

      IF(TEST1 .GE. 0.9998) THEN
         DL90=0.
      ELSE IF(TEST1 .LE. -0.9998) THEN
         DL90=1.0
      ELSE
         DL90=ACOS(-B)/PI
      END IF

c	print *,' dl90=',dl90

      SINL=SIN(LATP*DTR)
      COSL=COS(LATP*DTR)+1.e-20

C     ADJUST THE LENGTH OF DAY AND LENGTH OF NIGHT FOR THE
C     CURRENT LATITUDE. TOTAL DAYLIGHT=24*ACOS(-TAN(L)*TAN(D))/PI
C     FROM RUNDEL (1977), WITH L=LATITUDE, D=DECLINATION.
        R=SINL*SIND/(COSL*COSD)
        IF(R.GE.1.e0)THEN
          TDF=1.
        ELSE IF(R.LE.-1.)THEN
          TDF=0.
        ELSE
          TDF=ACOS(-R)*RPI
        END IF
        TNF=1.-TDF
       TFD(IJ)=TDF

c	print *,' nht==',nht

	DO 3500 INDT=1,NHT
	   DLITE(INDT)=TFD(IJ)
 3500	   CONTINUE

c	   print *,' dlite=',dlite

         CALL CLOCKAER(IJ,DTIMEB,PMTIME,DLITE,DL90,L$)


c      print *,' ij latp cosl sinl',ij,latp,cosl,sinl

C  SOLAR ZENITH ANGLE CHID. - problem with directly overhead sun and COSBC=1, it may actually be 
C   numerically slightly greater than 1 (eg, 1.0000001), so that ACOS(COSBC) = NaN, so reset COSBC = 1.
C
C   COSBC=SIND*SINL+COSD*COSL*COS(OMEGAD*DTR*(CLOCK+12.0))   ! omegad no longer used w/ diurnal model
C    also, the 15 factor converts from hours (CLOCK) to degrees, and DTR=pi/180., converts from degs to rads
C

       DO 500 ICLOCK=1,18    ! SUNSET
	  CLOCK=((PMTIME(ICLOCK)/PI)*12.) + 12.

c	  print *,' iclock clock=',iclock,clock

      COSBC=SIND*SINL+COSD*COSL*COS(15.*(CLOCK+12.0)*DTR)
        if (COSBC .gt. 1.) COSBC=1.
        if (COSBC .lt. -1.) COSBC=-1.

      R=ACOS(COSBC)
      SINC=SIN(R)
      CHIDA(ICLOCK,IJ)=R*RTD

	COSAER=SIND*SINL+COSD*COSL*COS(PMTIME(ICLOCK))
	SZAAER=ACOS(COSAER) * RTD

c	print *,' chida iclock ij cosaer szaaer=',chida(iclock,ij),
c     *  iclock,ij,cosaer,szaaer

C     GRAZING ALTITUDE ZGRZ. AT ALTITUDES FOR WHICH IT IS ZERO OR
C     LESS, THE SUN IS NOT UP; THE J COEFFICIENTS SHOULD BE SET TO
C     ZERO, AND THE CHAPMAN FUNCTION IS NOT REQUIRED.
      IF(CHIDA(ICLOCK,IJ).LE.90.0)THEN
        ZGRZ(IJ)=1000.0
      ELSE
        ZGRZ(IJ)=SINC*(ZP+R0)-R0
      END IF
c	type *,'ij,zgrz(ij),zp,r0,chid(ij)',ij,zgrz(ij),zp,r0,chid(ij)
C     CHAPMAN FUNCTION. APPROXIMATE IT WITH THE SECANT(CHI) FOR
C     CHI LESS THAN 70 DEGREES.
      IF(CHIDA(ICLOCK,IJ).LT.70.0)THEN       
        SFACA(ICLOCK,IJ)=1.0/COSBC
c	print *,' ij=',ij,' cosbc=',cosbc,' sfac=',sfac(ij)
      ELSE IF(ZGRZ(IJ).GT.0.0)THEN
C       CHAPMAN FUNCTION SFAC FOR CHID GREATER THAN OR EQUAL TO 70
C       DEGREES WHEN ZGRZ IS POSITIVE.
        R=SQRT(0.50*R0/HBAR)
        S=R*ABS(COSBC)

        IF(S.LE.8.0)THEN
          S=(D1+D2*S)/(D3+D4*S+S**2)
        ELSE
          S=D5/(D6+S)
        END IF

        R=R*SQRT(PI)
        SFACA(ICLOCK,IJ)=R*S
        IF(CHIDA(ICLOCK,IJ).GT.90.0)SFACA(ICLOCK,IJ)=
     .  2.0*R*EXP((R0+ZBAR)*(1.0-SINC)/HBAR) - SFACA(ICLOCK,IJ)
      ELSE
        SFACA(ICLOCK,IJ)=-9999.0
      END IF

  500 CONTINUE

c	 if(ij.eq.4)stop
c	print *,' ij,sfac, zgrz tfd=',ij,sfac(ij),zgrz(ij),
c     *  tfd(ij)
3000	CONTINUE

c	 if(1.gt.0)stop

C	print *,' sfac=',sfac
C	print *,' zgrz=',zgrz

      RETURN
      END
