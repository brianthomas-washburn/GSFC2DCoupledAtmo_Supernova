C **********************************************************************
C
      SUBROUTINE TIMEAVJ(JRA,JEXR,DTIMEA,I,J,TFD,L$,Z$,JAVAER,JAER, 
     >                   PH$, IL$, JZ39, J39)
C 
C
C  this computes the diurnally averaged J-coeffs, 
C     adapted from the AER TIMEAV routine 
C     also load in diurnal Js into JAER
C
C **********************************************************************
C

      include  'com_aerd.h'
      include  'com_aers.h'


      DOUBLE PRECISION JRA(NJR,NTIME), DTIMEA(NTIME), JEXR(3,NTIME)

      INTEGER L$,Z$, PH$, IL$
      REAL TFD(L$)
      REAL*8 JAVAER(NJR+3,L$,Z$), SPDAY, JAER(NJR+3,NTIME,L$,Z$)


C  wavelength dependent JZ39 array (for diagnostics ONLY), J39 is the output

      REAL*8 JZ39(18,PH$,5,L$,Z$), J39(PH$,5,L$,Z$)


      COMMON/FLAGS/IFL(90)

C
C         WRITE THE SUBROUTINE ID THE FIRST TIME AROUND
CCC      IF(IFL(44).EQ.0) THEN
CCC         WRITE(16,*) '@(#)timeav.f	1.5  06/29/99'
CCC         IFL(44)=1
CCC      END IF
C

      HDAY = TFD(I)*86400.
      HNITE = (1.-TFD(I))*86400.

C      HDAY=DLITE*86400.
C      HNITE=(1.-DLITE)*86400.


C    Compute diurnal avg Js - DO THIS AS IN AER's TIMEAV.F, 
C           ie, - only compute daytime values, assume nitetime values are 0
CCCCC -             DTIMEA(NTIME),  JRA(NJR,NTIME), JAVAER(NJR+3,L$,Z$)

      DO 200 IJR=1,NJR
        SPDAY=0.D0

         do 410 IT=1,SUNSET-1
 410      SPDAY = SPDAY + (JRA(IJR,IT) + JRA(IJR,IT+1))*0.5D0*DTIMEA(IT)

         do 420 IT=SUNRIS,NTIME-1
 420      SPDAY = SPDAY + (JRA(IJR,IT) + JRA(IJR,IT+1))*0.5D0*DTIMEA(IT)
 
        JAVAER(ijr,I,J) = SPDAY/86400.0D0
 200  CONTINUE


C    Compute diurnal avg Js for extras JEXR(3,NTIME), JAVAER(NJR+3,L$,Z$)

      DO 202 IJR=1,3
        SPDAY=0.D0

        do 402 IT=1,SUNSET-1
 402    SPDAY = SPDAY + (JEXR(IJR,IT) + JEXR(IJR,IT+1))*0.5D0*DTIMEA(IT)

        do 403 IT=SUNRIS,NTIME-1
 403    SPDAY = SPDAY + (JEXR(IJR,IT) + JEXR(IJR,IT+1))*0.5D0*DTIMEA(IT)

        JAVAER(NJR+ijr,I,J) = SPDAY/86400.0D0
 202  CONTINUE



C  also load in Diurnally varying J's into JAER(NJR+3,NTIME,L$,Z$)

      DO 700 IT=1,NTIME
      DO 700 IJR=1,NJR
 700     jaer(ijr,it,I,J) = JRA(IJR,IT)

      DO 800 IT=1,NTIME
      DO 800 IJR=1,3
 800     jaer(NJR+ijr,it,I,J) = JEXR(IJR,IT)



C    Compute diurnal avg Js for wavelength dependent JZ39(18,PH$,5,L$,Z$) (for diagnostics ONLY)
C                                 load into J39(PH$,5,L$,Z$)
      DO 277 ijj0=1,5
      DO 277 IJR=1,PH$
         SPDAY=0.D0

         do 510 IT=1,SUNSET-1
         SPDAY = SPDAY + (JZ39(IT,ijr,ijj0,i,j)+JZ39(IT+1,ijr,ijj0,i,j))
     >                    *0.5D0*DTIMEA(IT)
 510     CONTINUE


         do 520 IT=SUNRIS,NTIME-1
         SPDAY = SPDAY + (JZ39(IT,ijr,ijj0,i,j)+JZ39(IT+1,ijr,ijj0,i,j))
     >                    *0.5D0*DTIMEA(IT)
 520     CONTINUE

         J39(ijr,ijj0,I,J) = SPDAY/86400.0D0
 277  CONTINUE



Cj
cj      COMMON/JAER1/JAVAER(NJR,18,46), AJOUT(2,ntime,18),
cj     >             AJAVDAY(3,NJR,18,46)
cj
Cj   DTIMEA(NTIME),  JRA(NJR,NTIME), AJAVDAY(3,NJR,L$,Z$);   Note: DTIMEA(NTIME=18) is always 0.0
cj
cj   THIS IS CODE FOR COMPUTING DAYTIME and NITETIME AVERAGED Js (of course the NITETIME avgs are all ZERO!)
Cj       and also the diurnal avg computed from the day and nite avgs, which is identical to 
Cj       the computation above -    use this code only if you need the DAYTIME AVG J-COEFFS.
Cj
Cj      DO 150 IJR=1,NJR
Cj        SPDAY=0.
Cj        SPNITE=0.
Cj
Cj      DO IT=1,SUNSET-1
Cj         SPDAY = SPDAY + (JRA(IJR,IT) + JRA(IJR,IT+1))*0.5*DTIMEA(IT)
Cj      END DO
Cj
Cj      DO IT=SUNSET,SUNRIS-1
Cj         SPNITE = SPNITE + (JRA(IJR,IT) + JRA(IJR,IT+1))*0.5*DTIMEA(IT)
Cj      END DO
Cj
Cj      DO IT=SUNRIS,NTIME-1
Cj         SPDAY = SPDAY + (JRA(IJR,IT) + JRA(IJR,IT+1))*0.5*DTIMEA(IT)
Cj      END DO
Cj
Cj
Cj      IF(HDAY.GT.0.) THEN
Cj         ajavday(1,ijr,I,J) = SPDAY/HDAY
Cj      ELSE
Cj         ajavday(1,ijr,I,J) = 0.0
Cj      END IF
Cj
Cj      IF(HNITE.GT.0.) THEN
Cj         ajavday(2,ijr,I,J) = SPNITE/HNITE
Cj      ELSE
Cj         ajavday(2,ijr,I,J) = 0.0
Cj      END IF
cc                                                         diurnal (24-hour) avg
Cj       ajavday(3,ijr,I,J) = (SPDAY + SPNITE)/86400.
Cj  150   CONTINUE
Cj 

      RETURN
      END

