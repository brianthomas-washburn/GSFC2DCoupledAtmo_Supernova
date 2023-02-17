	SUBROUTINE SRBAND(IJC,IKC)

        include "com2d.h"

	DIMENSION SUM(17),SRSIG(17),SRSIG0(17),NOSIG0(2)

        SAVE


C     CALCULATION OF THE
C     O2 DISSOCIATION CROSS-SECTIONS IN THE SCHUMANN-RUNGE
C     BAND INTERVALS 1 TO 17 BASED ON THE WORK OF ALLEN
C     AND FREDERICK (1981). CROSS-SECTIONS ARE STORED ON
C     ACKERMAN INTERVAL NUMBER (INCREASING WAVELENGTH);
C     CALCULATIONS ARE MADE ON INTERVAL NUMBER (DECREASING
C     WAVELENGTH). NO DISSOCIATION CROSS-SECTIONS IN THE
C     DEL(0,0) AND DEL(1,0) GAMMA BANDS ARE ALSO COMPUTED.
C
C     NOTE: CODE MODIFIED TO USE BANDS 4-17 (1762-1990A)
C           WITH THE HERMAN AND MENTALL O2 CROSS-SECTIONS
C           IN THE HERZBERG CONTINUUM. BAND 1=2051.5A.


	CALL FILL(IJC,IKC)

C     PRESSURE (MB) AT CURRENT POSITION AND OTHER PARAMETERS.
      PRIN=PRP
      IF(PRIN.LT.PRA(4))PRIN=PRA(4)
      PRINL=aLOG10(PRIN)
      NCL=aLOG10(NCOL(ISCO2))

C     O2 CROSS-SECTIONS AT SOLAR ZENITH ANGLE=0 DEGREES.
C     INITIALIZE EXPANSION SUMMATIONS.
      DO 200 JI=1,17
  200 SUM(JI)=SR1(1,JI)

C     BANDS 1-3 USE A TEMPERATURE EXPANSION VALID FOR T<300 DEG K.
      T0=TS
      IF(T0.GT.TLIM)T0=TLIM
      TI=1.0
      DO 220 I=2,9
      TI=TI*T0
      DO 220 JI=1,3
  220 SUM(JI)=SUM(JI)+SR1(I,JI)*TI

C     BANDS 4-17 USE A PRESSURE EXPANSION. BANDS 4-13 (JJ=1) USE
C     PRIN IF PRIN<PR20=PRESSURE(MB) AT 20 KM AND PR20 OTHERWISE.
C     SIMILARLY, BANDS 14-16 (JJ=2) CUT OFF AT 30 KM; BAND 17
C     (JJ=3) AT 40 KM.
      DO 240 JJ=1,3
      JL=JBL(JJ)
      JU=JBU(JJ)
      IF(PRIN.LT.PRA(JJ))THEN
        P0=PRINL
      ELSE
        P0=PRAL(JJ)
      END IF
      PRI=1.0
      DO 240 I=2,9
      PRI=PRI*P0
      DO 240 JI=JL,JU
  240 SUM(JI)=SUM(JI)+SR1(I,JI)*PRI

C     STORE CROSS-SECTIONS AT CHI=0 ON ACKERMAN INTERVAL
C     NUMBER, AND INITIALIZE EXPANSION AT CURRENT SOLAR
C     ZENITH ANGLE.
      DO 260 JI=1,17
      SRSIG0(18-JI)=10.**SUM(JI)

       IF(YSCORSR)SRSIG0(18-JI)=SRSIG0(18-JI)+CORSR(JI)
C   CORSR IS THE CORRECTION TO THE CROSS SECTION ADVISED BY
C   JOHN FREDERICK 5/26/83.

  260 SUM(JI)=SR2(1,JI)


C     CROSS-SECTIONS FOR CURRENT SOLAR ZENITH ANGLE CHID.
      NCO2=1.0
      DO 270 I=2,6
      NCO2=NCO2*NCL
      DO 270 JI=1,17
  270 SUM(JI)=SUM(JI)+SR2(I,JI)*NCO2

C     STORE ON ACKERMAN INTERVAL NUMBER.
      DO 280 JI=1,17
      SR3=10.**SUM(JI)
      IF(SR3 .LT. 1.e-15)SR3=0.00
  280 SRSIG(18-JI)=SRSIG0(18-JI)*SFAC(IJC)**(-SR3)


C     STORE CROSS-SECTIONS IN XSCHRUN AND XO2 (FOR COMPUTING OPTICAL DEPTH).
	DO 290 JL=5,18
	XSCHRUN(JL,IJC,IKC)=SRSIG(JL-5+1)
c        if(ikc.eq.20 .and. ijc.eq.9 .and. jl.eq.10)then
c           write(14,4101)jl,ijc,ikc,xschrun(jl,ijc,ikc),ncl,
c     *  prp,ts
c 4101      format(' jl=',i3,' ijc=',i3,' ikc=',i3,' xsch=',
c     *  1pe14.6,' ncl=',1pe14.6,' prp=',1pe14.6,' ts=',
c     *  1pe14.6)
c           endif
290	XO2(JL)=SRSIG(JL-5+1)



C     NO PHOTODISSOCIATION CROSS-SECTIONS.
C     CROSS SECTIONS FOR SOLAR ZENITH ANGLE = 0 DEGREES.
      SUMNO1=SRNO1(1,1)
      SUMNO2=SRNO1(1,2)
C     PRESSURE EXPANSION: UPPER ALTITUDE CUT OFF IS 90 KM; FOR
C     LOWER ALTITUDES THE CROSS-SECTION LIMIT IS 1.E-15 CM**2.
      IF(PRIN.LT.PRA(4))THEN
        P0=PRAL(4)
      ELSE
        P0=PRINL
      END IF
      PRI=1.0
      DO 300 I=2,9
      PRI=PRI*P0
      SUMNO1=SUMNO1+SRNO1(I,1)*PRI
      SUMNO2=SUMNO2+SRNO1(I,2)*PRI
  300 CONTINUE
C     STORE IN ACKERMAN INTERVAL ORDER.
      IF(SUMNO1.GT.ESIGLM)THEN
        NOSIG0(2)=SIGLM
      ELSE
        NOSIG0(2)=10.**SUMNO1
      END IF
      IF(SUMNO2.GT.ESIGLM)THEN
        NOSIG0(1)=SIGLM
      ELSE
        NOSIG0(1)=10.**SUMNO2
      END IF
C     CROSS-SECTIONS FOR CURRENT ZENITH ANGLE. RESTRICT
C     CROSS-SECTIONS TO A MAXIMUM VALUE OF SIGLM. THIS
C     APPROXIMATION IS VALID FOR CHID<85 DEGREES.
      IF(CHID(IJC).GT.85.0)THEN
        CHX=10.260
      ELSE
        CHX=SFAC(IJC)
      END IF
      SUMNO1=SRNO2(1,1)
      SUMNO2=SRNO2(1,2)
      NCO2=1.0
      DO 320 I=2,5
      NCO2=NCO2*NCL
      SUMNO1=SUMNO1+SRNO2(I,1)*NCO2
      SUMNO2=SUMNO2+SRNO2(I,2)*NCO2
  320 CONTINUE
C     STORE ON ACKERMAN INTERVAL.
      NOSIG(1,IJC,IKC)=NOSIG0(1)*CHX**SUMNO2
      IF(NOSIG(1,IJC,IKC).GT.SIGLM)NOSIG(1,IJC,IKC)=SIGLM
      NOSIG(2,IJC,IKC)=NOSIG0(2)*CHX**SUMNO1
      IF(NOSIG(2,IJC,IKC).GT.SIGLM)NOSIG(2,IJC,IKC)=SIGLM

      RETURN
      END