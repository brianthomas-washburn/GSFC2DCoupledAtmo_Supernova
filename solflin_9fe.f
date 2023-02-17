	SUBROUTINE SOLFLIN

        include "com2d.h"

      DATA D1/1.060693/,D2/0.55643831/,D3/1.0619896/,
     .  D4/1.7245609/,D5/0.56498823/,D6/0.06651874/,
     .  DEC1/279.575/,DEC2/0.985647/,DEC3/356.967/,
     .  DEC4/0.985600/,DEC5/1.916/,DEC6/0.020/,
     .  DEC7/0.39782/,DEC8/0.0167/
        DATA RSC/1.12,1.10,1.04,1.03,1.01/
        DATA R27/1.07,1.06,1.03,1.01/
        DATA RLAMB/1750.,1900.,2100.,2400.,2600.,3000./

	SAVE


C       INITIALIZE CONSTANTS.
        DEC1=DEC1*DTR
        DEC2=DEC2*DTR
        DEC3=DEC3*DTR
        DEC4=DEC4*DTR
        DEC5=DEC5*DTR
        DEC6=DEC6*DTR

C INPUT WAVELENGTHS AND FLUXES

C UNIT=92,FILE='WFLUX85MA.DAT'
	WRITE(6,3102)
3102   FORMAT(' READ WFLUX85MA.DAT')

	READ(92,18)
C INPUT NUMBER OF WAVELENGTH INTERVALS; ACTUALLY IL$ IN PARAMETER STATEMENT

	READ(92,11)(WVL(JL),JL=1,IL$)
C INPUT WAVELENGTHS

	READ(92,18)
C DUMMY READ TO SKIP HEADER CARD.

	READ(92,13)(FLUX0(JL),JL=1,IL$)
C INPUT SOLAR FLUXES


C     TABLES USED FOR DIURNAL AVERAGING OF SOLAR FLUXES. REFERENCE
C     LATITUDES ARE ALL POSITIVE.
C      UNIT=94,FILE='DAVG.DAT'
       WRITE(6,3104)
3104   FORMAT(' AFTER OPEN DAVG.DAT')
      READ(94,6)ZREF,ILATR,IDA,ITA
      READ(94,11)(LATR(I),I=1,ILATR)
      READ(94,11)(DECA(I),I=1,IDA)
      READ(94,11)(TAUA(I),I=1,ITA)
      READ(94,19)ITAUA
      READ(94,11)(((FAVG(I,JI,JL),I=1,ITA),JI=1,IDA),JL=1,ILATR)
      WRITE(6,8)(LATR(I),I=1,ILATR)

C     CONVERT OPTICAL DEPTHS AND SOLAR FLUX RATIOS TO LOG FORM FOR
C     INTERPOLATION.
      R=aLOG(1.e-20)
      DO 500 I=1,ITA
      RT=TAUA(I)
      RN=R
      IF(RT.NE.0.e0)RN=aLOG(RT)
      TAUA(I)=RN
      DO 500 JI=1,IDA
      DO 500 JL=1,ILATR
      RT=FAVG(I,JI,JL)
      RN=R
      IF(RT.NE.0.e0)RN=aLOG(RT)
  500 FAVG(I,JI,JL)=RN


	RETURN
    6 FORMAT(44X,F4.0,17X,3I5)
    8 FORMAT('0*** REFERENCE LATITUDES FOR DIURNAL AVERAGING OF SOLAR
     . FLUXES ***'/(5X,20F5.0))
   11 FORMAT(10F8.2)
   13 FORMAT(1P8E10.4)
   14 FORMAT(1P4E14.6/1P4E14.6/1PE14.6)
   15 FORMAT(1P4E14.6/1P2E14.6)
   16 FORMAT(1P4E14.6/1PE14.6)
   18 FORMAT(I10)
   19 FORMAT(20I4)
      END
