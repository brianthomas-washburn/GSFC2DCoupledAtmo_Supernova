C **********************************************************************
C
      SUBROUTINE CLOCKAER(IJ,DTIMEB,PMTIME,DLITE,DL90,L$)

c      SUBROUTINE CLOCK(DTIME,PMTIME,DLITE,DL90,JHT,LOLITA)
C
C **********************************************************************
C    Taken from AER routine "clock.f"
C
C         THIS SUBROUTINE SETS UP THE TIME GRID POINTS FOR 24 HOURS
C
C    DLITE IS HOURS OF DAYLIGHT (FUNCTION OF ALTITUDE).
C    DL90 IS HOURS OF DAYLIGHT BETWEEN 90 DEGREE ZENITH ANGLES.
C    HOURS FROM NOON TO 90 DEGREE ZENITH ANGLE DIVIDED INTO NDIVDAY INTERVALS
C    HOURS FROM 90 DEGREE ZENITH ANGLE TO SUNSET DIVIDED INTO NDIVTWI ITERVALS
C    HOURS FROM SUNSET TO SUNRISE DIVIDED INTO NDIVNIT INTERVALS
C    TOTAL NUMBER OF TIME INTERVALS IS 2*(NDIVDAY+NDIVTWI)+NDIVNIT
C    TOTAL NUMBER OF TIME GRID POINTS IS 2*(NDIVDAY+NDIVTWI)+NDIVNIT=NTIME
C    DTIMEA(I)=LENGTH OF I-TH INTERVAL IN SECONDS
C    PMTIME(I)=TIME PAST NOON OF THE LEFT END OF I-TH INTERVAL (RADIANS)
C    LOLITA=PRINTOUT CONTROL PARAMETER
C

       include "com_aerd.h"
       include "com_aerg.h"

c      include 'diurnal.i'
c      include 'grid.i'

      INTEGER L$

      DOUBLE PRECISION DTIMEB(L$,NTIME)
      REAL PMTIME(NTIME),DLITE(NHT)
C        CONVERSION FACTOR FROM FRACTIONS OF DAY TO RADIANS = 2*PI
      PARAMETER (FACTOR=2.*3.14159265)
C
c      COMMON/FLAGS/IFL(90)
C
C         WRITE THE SUBROUTINE ID THE FIRST TIME AROUND
c      IF(IFL(45).EQ.0) THEN
c         WRITE(16,*) '@(#)clock.f	1.6  06/28/99'
c         IFL(45)=1
c      END IF
c      IF(DL90.GT.DLITE(JHT)) THEN
c         WRITE(16,*) 'ERROR IN CLOCK:  DL90 > DLITE'
c         STOP 'ERROR IN CLOCK:  DL90 > DLITE'
c      END IF
C
C          TIME OF 90 DEGREE ZENITH ANGLE IN RADIANS PAST NOON
      TTWI=DL90/2.*86400.
C          LOCAL SUNSET AND SUNRISE TIMES IN RADIANS PAST NOON

      jht=nht
c      print *,' jht nht',jht,nht

      TSET=DLITE(JHT)/2.*86400.
      TRIS=(1.-DLITE(JHT)/2.)*86400.
      TTOP=DLITE(NHT)/2.*86400.

c       print *,' tset dlite(jht) tris dlite(nht) ttop jht nht',
c     * tset,dlite(jht),tris,dlite(nht),ttop,jht,nht

C
C      DAYTIME TIME STEPS MUST BE THE SAME FOR ALL HEIGHTS AT ONE LATITUDE
C      SO THAT UCI PHOTOLYSIS CODE WILL RUN EFFICIENTLY
C
C          TIME STEP IN THE DAYTIME (NOON TO TWILIGHT) IN RADIANS
      nday=ndivday

      HDAY=TTWI/NDAY
c      print *,' nday hday',nday,hday

C          TIME STEP IN TWILIGHT (90 DEGREES TO LOCAL SUNSET) IN RADIANS
      HTWI1=(TTOP-TTWI)/NDIVTWI
      ntwi1=0
      do i=1,ndivtwi
         if(i*htwi1.le.(tset-ttwi)) ntwi1=i
      end do
      htwi2=(tset-ttwi-ntwi1*htwi1)
      ntwi2=0
      if(htwi2.gt.0.0001) ntwi2=1
      if(jht.eq.nht) then
         ntwi1=ndivtwi
         ntwi2=0
         htwi2=0
      end if
C          TIME STEP AT NIGHT (SUNSET TO SUNRISE) IN RADIANS
      nnit=ndivnit+(ndivtwi-(ntwi1+ntwi2))*2
      HNIT=(TRIS-TSET)/NNIT
C
C         DTIMEB IN SECONDS
CVD$ NOVECTOR
      JCLOCK=1
      DO I=1,NDAY
         DTIMEB(IJ,JCLOCK)=HDAY
         IF(DTIMEB(IJ,JCLOCK).lt.0.0E0)DTIMEB(IJ,JCLOCK)=0.0E0

c         print *,' jclock dtimea',jclock,dtimea(jclock)

         JCLOCK=JCLOCK+1
      END DO
      SET90=JCLOCK
      DO I=1,ntwi1
         DTIMEB(IJ,JCLOCK)=Htwi1
         IF(DTIMEB(IJ,JCLOCK).lt.0.0E0)DTIMEB(IJ,JCLOCK)=0.0E0

c         print *,'  jclock dtimea',jclock,dtimea(jclock)

         JCLOCK=JCLOCK+1
      END DO
      DO I=1,Ntwi2
         DTIMEB(IJ,JCLOCK)=HTWI2
         IF(DTIMEB(IJ,JCLOCK).lt.0.0E0)DTIMEB(IJ,JCLOCK)=0.0E0

c         print *,'   jclock dtimea',jclock,dtimea(jclock)

         JCLOCK=JCLOCK+1
      END DO
      SUNSET=JCLOCK

c      print *,' SUNSET=',sunset

      DO I=1,NNIT
         DTIMEB(IJ,JCLOCK)=HNIT
         IF(DTIMEB(IJ,JCLOCK).lt.0.0E0)DTIMEB(IJ,JCLOCK)=0.0E0

c         print *,'    jclock dtimea',jclock,dtimea(jclock)

         JCLOCK=JCLOCK+1
      END DO
      SUNRIS=JCLOCK
      DO I=1,Ntwi2
         DTIMEB(IJ,JCLOCK)=HTWI2
         IF(DTIMEB(IJ,JCLOCK).lt.0.0E0)DTIMEB(IJ,JCLOCK)=0.0E0

c         print *,'     jclock dtimea',jclock,dtimea(jclock)

         JCLOCK=JCLOCK+1
      END DO
      DO I=1,ntwi1
         DTIMEB(IJ,JCLOCK)=Htwi1
         IF(DTIMEB(IJ,JCLOCK).lt.0.0E0)DTIMEB(IJ,JCLOCK)=0.0E0

c         print *,'      jclock dtimea',jclock,dtimea(jclock)

         JCLOCK=JCLOCK+1
      END DO
      RIS90=JCLOCK
      DO I=1,NDAY
         DTIMEB(IJ,JCLOCK)=HDAY
         IF(DTIMEB(IJ,JCLOCK).lt.0.0E0)DTIMEB(IJ,JCLOCK)=0.0E0

c         print *,'       jclock dtimea',jclock,dtimea(jclock)

         JCLOCK=JCLOCK+1
      END DO
C
C         PMTIME IN RADIANS FROM NOON
      PMTIME(1)=0.
      DO 300 I=2,NTIME
         PMTIME(I)=PMTIME(I-1)+DTIMEB(IJ,I-1)/86400*FACTOR

c         print *,' i pmtime factor',i,pmtime(i),factor

  300 CONTINUE

c      print *,' pmtime=',pmtime

C
C         PRINT OUT IF REQUIRED
c      IF(LOLITA.EQ.0)RETURN
C
c      WRITE(16,6101) (PMTIME(I)/FACTOR*24.,I=1,NTIME)
c 6101 FORMAT(/10X,'TIME GRID POINTS (HOURS FROM NOON) ARE'/20(1X,F5.2))
c      WRITE(16,6102) (DTIMEB(IJ,I)/3600.,I=1,NTIME-1)
c 6102 FORMAT(10X,'TIME INTERVALS (HOURS) ARE'/20(1X,F5.2))
      RETURN
      END



