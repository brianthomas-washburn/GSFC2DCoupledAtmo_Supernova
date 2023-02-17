	SUBROUTINE CROSSIN

       include "com2d.h"


C  common for NO2 cross section for photolysis heating rates:

      COMMON/CXNO2/xsectno2t(IL$,201), xsectno2(IL$)


      DATA SRL/1762.0,1990.0/,
     .  PRA/55.2930,11.9700,2.87140,1.8359e-3/,
     .  SIGNO/1810.0,1896.0,1826.50,1914.0/,
     .  SIGNOUSE/1810.0,1896.0,1826.50,1914.0/,
     .  CORSR/0.0,0.0,0.0,-3.84E-24,-3.89E-24,-3.87E-24,
     .  -3.85E-24,-3.77E-24,-3.70E-24,-3.62E-24,-3.54E-24,
     .  -3.38E-24,-3.27E-24,-3.24E-24,0.0,0.0,0.0/,
     .  TLIM/300.0/,SIGLM/1.e-15/,ESIGLM/-15.0/
      DATA JBL/4,14,17/,JBU/13,16,17/,
     .  ILLB/17/

       DATA WO1D/2985.0,3250.0/, WN2O5/2850.0,3800.0/
C
      DATA TS320/320./,TS300/300.0/,TS230/230.0/, TS220/220.0/,
     .  TS225/225.0/,TS180/180.0/, WCH2O/3300.0/
C
C  COEFFS FOR NO2 TEMP CORRECTION FOR 2667-5475A  JPL-94 (no change for JPL-00)  
C
       DATA NO2T/25*0.0, 0.067, -0.043, -0.082, -0.280, -0.364, -0.534,
     . -0.681, -0.800, -1.246, -1.640, -1.100, -0.866, -0.696, 0.0/
C
C  PUT IN PARAMETER VALUES FOR CROSS SECTIONS FOR F113, F114, F115
C  12/30/88  FROM SIMON ET AL. (1988)
       DATA A0/-1087.9,-160.5,5.8281/
       DATA A1/20.004,2.4807,-0.299/
       DATA A2/-0.1392,-1.5202E-2,1.3525E-3/
       DATA A3/4.2828E-4,3.8412E-5,-2.6851E-6/
       DATA A4/-4.9384E-7,-3.4373E-8,0.0E0/
       DATA B0/12.493,-1.5296,0.0E0/
       DATA B1/-0.23937,3.5248E-2,0.0E0/
       DATA B2/1.7142E-3,-2.9951E-4,0.0E0/
       DATA B3/-5.4393E-6,1.1129E-6,0.0E0/
       DATA B4/6.4548E-9,-1.5259E-9,0.0E0/
c
c put in parameter values for HOBr cross section, from Burkholder
c measurements. David B. Considine, 950216  -- HOBr included as usual in JPL-00 ASCII file
cc       data hobrlam/268.,286.,350./
cc       data hobra/0.113,0.262,0.062/
cc       data hobrb/0.00349,0.00279,0.000822/
C
	SAVE


C   read in  NO2 cross section for photolysis heating rates:  xsectno2t(IL$=39,201)
C   - this was written out in IDL on Linux, so use "little_endian"  (NO change for JPL10)
C
      OPEN (1789, FILE = 'xsect_NO2_heat.dat',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='stream', STATUS = 'OLD')
           READ (1789) xsectno2t
        CLOSE (1789)


C INPUT PHOTODISSOCIATION CROSS-SECTIONS, first the non-temp dependency from JPL-10 (ASCII file)
C INPUT FIT PARAMETERS FOR THE SCHUMANN-RUNGE BAND X-SECTIONS FROM ALLEN AND FREDERICK (1982).

	WRITE(6,3103)
3103   FORMAT(' READ XSECT10_JPL10.DAT')

	DO 140 I=1,PH$
	READ(93,21)
	READ(93,20) IPHOT(I), PHOTCHAR(I)
20	FORMAT(I2,7X,A37)
C	WRITE(6,20) PHOTCHAR(I)
	READ(93,13)(XSECT(I,JL),JL=1,IL$)
140	CONTINUE
C INPUT CROSS SECTIONS
C
C  Now READ in SCHUMANN-RUNGE BANDS
C
	READ(93,21)
	READ(93,21)

	DO 50 JO=1,17
	READ(93,14)(SR1(I,JO),I=1,9)
50	CONTINUE
	DO 52 JO=1,17
	READ(93,15)(SR2(I,JO),I=1,6)
52	CONTINUE
	DO 54 JO=1,2
	READ(93,14)(SRNO1(I,JO),I=1,9)
54	CONTINUE
	DO 56 JO=1,2
	READ(93,16)(SRNO2(I,JO),I=1,5)
56	CONTINUE
C
C     FIT PARAMETERS FOR SCHUMANN-RUNGE BAND CALCULATION FROM
C     ALLEN AND FREDERICK. 
C
C
C  NOW READ IN TEMP DEPENDENT X-SECTIONS FROM JPL-10
C  NOTE: ALL X-sections are now in the temp dependent array, 
C  so the ASCII file  XSECT10_JPL10.DAT(PH$,IL$) is only used for the SR bands etc at the end
C
C  XSECTTD(PH$=81,IL$=39,201); xspd(3,201,157), prpd(157) are in COMMON
C     - this was written out in IDL on Linux, so use "little_endian"
C
        OPEN (23, FILE = 'xsecttd10_jpl10.dat', convert="little_endian",
     >   FORM = 'unformatted', ACCESS='stream', STATUS = 'old')
            READ (23) xsecttd
            READ (23) xspd, prpd
        CLOSE (23)


C  convert prpd to -ALOG for proper interpolation in CROSEC;  aprpd(157) is in COMMON

        do 177 ihj = 1,157
 177	   aprpd(ihj) = -ALOG(prpd(ihj))

C
C SET UP ARRAYS AND add the two O2 and O3 x-sections - JPL-06/JPL-10
C   use the temp dependent array for 240K - index 121  (index 1 = 120K)
C
	DO 142 JL=1,IL$
cccccc           XNO2(JL) = xsecttd(7,JL,121) - don't need this any longer
	   XO2(JL)  = xsecttd(1,JL,121) + xsecttd(46,JL,121)
           XO3(JL)  = xsecttd(2,JL,121) + xsecttd(3,JL,121)
142	CONTINUE
C
C       INDICES OF WAVELENGTHS DEFINING THE S.R. BANDS.
C
        I=1
        DO 110 JJ=1,2
        LAM=SRL(JJ)
        DO 120 JI=I,IL$
        IF(WVL(JI).GE.LAM)GO TO 130
  120   CONTINUE
  130   ISR(JJ)=JI
        IF(ABS(LAM-WVL(JI-1)).LT.1.e-5)ISR(JJ)=JI-1
  110   I=JI
        ILL=ISR(1)
        ILU=ISR(2)
	WRITE(6,1122)ISR(1),ISR(2)
1122	FORMAT(' ISR1=',I5,' ISR2=',I5)
C       INDICES OF WAVELENGTHS DEFINING NO DEL BANDS.
        I=1
        DO 1140 JJ=1,4
        LAM=SIGNOUSE(JJ)
        DO 150 JI=I,IL$
        IF(WVL(JI).GE.LAM)GO TO 160
  150   CONTINUE
  160   INOUSE(JJ)=JI
        IF(ABS(WVL(JI-1)-LAM).LT.1.e-5)INOUSE(JJ)=JI-1
        IF(JJ.EQ.1)I=JI
 1140   CONTINUE
        INOL(1)=INOUSE(1)
        INOL(2)=INOUSE(2)
        INOU(1)=INOUSE(3)
        INOU(2)=INOUSE(4)
	WRITE(6,1123)INOL(1),INOL(2),INOU(1),INOU(2)
1123	FORMAT(' INOL1=',I5,' INOL2=',I5,' INOU1=',I5,' INOU2=',I5)
C       COMMON LOGS OF U.S. STANDARD ATMOSPHERE PRESSURES
C       (MB) AT 20,30,40,90 KM, USED IN EXPANSIONS.
        DO 170 JJ=1,4
  170   PRAL(JJ)=aLOG10(PRA(JJ))
C
C       IDENTIFY LOWER WAVELENGTH FOR PRESSURE-DEPENDENT CH2O=H2+CO
C       CROSS-SECTIONS.
        DO 5160 JI=1,IL$
        IF(WCH2O.LT.WVL(JI).OR.ABS(WCH2O-WVL(JI)).LT.1.e-10)GO TO 5170
 5160   CONTINUE
 5170   ICH2O=JI
	WRITE(6,1126)ICH2O
1126	FORMAT(' ICH2O=',I5)

	RETURN

   13 FORMAT(1P8E10.3)
   14 FORMAT(1P4E14.6/1P4E14.6/1PE14.6)
   15 FORMAT(1P4E14.6/1P2E14.6)
   16 FORMAT(1P4E14.6/1PE14.6)
   21 FORMAT(A10)
      END


