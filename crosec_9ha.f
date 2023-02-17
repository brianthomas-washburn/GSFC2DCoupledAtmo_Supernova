	SUBROUTINE CROSEC(IJC,IKC)

        include "com2d.h"


C  common for NO2 cross section for photolysis heating rates:

        COMMON/CXNO2/xsectno2t(IL$,201), xsectno2(IL$)

        REAL ss1(157), frpr(1), frout(1)


	SAVE

C
C  FIRST, FIND TEMP INDEX FOR TEMP DEPENDENT XSECTIONS FROM JPL-10
C  GET ALL Xsections now from the temp dependent array: XSECTTD(PH$,IL$,201) in COMMON
C
        IJTS = IFIX(TS + 0.5) - 119
          IF (TS .LT. 120.) IJTS = 1
          IF (TS .GT. 320.) IJTS = 201
C
C  Load in all X-sections into XSECT(PH$,IL$) (in COMMON)
C
        DO 701 I = 1,IL$
        DO 701 iph = 1,PH$
           XSECT(iph,I) = XSECTTD(iph,I,IJTS)
 701    CONTINUE

C                              load in NO2 cross section for photolysis heating rates
        DO 705 I = 1,IL$
 705	   xsectno2(i) = xsectno2t(i,IJTS)


C
C O2 DISSOCIATION CROSS SECTIONS
	DO 400 I=5,18
	XSECT(1,I)=XSCHRUN(I,IJC,IKC)
C NEED TO ASSIGN CROSS SECTIONS ONLY IN SCHUMANN-RUNGE BANDS.
C THE CROSS SECTIONS OUTSIDE THE SCHUMANN-RUNGE BANDS HAVE ALREADY
C BEEN ASSIGNED.
400	CONTINUE


C NO DISSOCIATION CROSS SECTIONS 
C  -  NOT Sure if this is used, J(16) is overwritten from NOPHOT in SPDR
C
        DO 520 I=8,9
	XSECT(16,I)=NOSIG(1,IJC,IKC)
520	CONTINUE
        DO 500 I=13,14
        XSECT(16,I)=NOSIG(2,IJC,IKC)
500	CONTINUE

C
C  J COEFFICIENTS FOR F113, F114, F115 FROM SIMON ET AL. (1988), DONE 12/30/88 
C     now use JPL-2010, done above
C
C
C   PRESSURE-DEPENDENT QUANTUM YIELD FOR CH2O => H2+CO FOR WAVELENGTHS > 330nm
C
C     this has been updated for JPL-2010, using look-up table for model
C     wavelength intervals 34, 35, 36 as a function of pressure (fraction of ATM)
C     overwrite XSECT array; first find pressure level; Surface pressure = 1013 mb
C  
C    XSECTTD(PH$=81,IL$=39,201); xspd(3,201,157), aprpd(157); PRES90(Z$), ss1(157)
C      XSECT(PH$,IL$) are all in COMMON;  frpr(1), frout(1)
C      xstd(11,*,*) is J[CH2O] => H2+CO
C
C   NOTE: in LINTERP, APRPD array must be monotonically INCREASING, 
C         ie, APRPD(157) = -ALOG(prpd(157))
C

        frpr(1) = -ALOG(pres90(ikc)/1013.)

        do 770 ijw=34,36

            do 771 i=1,157
 771	       ss1(i) = xspd(ijw-33,ijts,i)

            CALL LINTERP(aprpd, 157, ss1, frpr, 1, 1, frout)

            XSECT(11,ijw) = frout(1)

 770	CONTINUE



C
C  NOW  WRITE  OUT  CROSS SECTIONS TO FORT.81
C
C          IF (IYR .GE. 19 .AND. DAY360 .GE. 350. 
C     *         .AND. IJC .EQ. 8 .AND. IKC .EQ. 17) THEN 
C
C	DO 140 I=1,PH$
C	WRITE (81,210)
C	WRITE (81,20) TS, IYR, DAY360, IPHOT(I), PHOTCHAR(I)
C
C	WRITE (81,13) (XSECT(I,JL),JL=1,IL$)
C140	CONTINUE
C
C        END IF

13      FORMAT(1P8E10.3)
20	FORMAT(F7.1, 3X, I3, 2X, F5.1, 5X, I2,3X,A37)
 210    FORMAT(A10)

	RETURN
         END
