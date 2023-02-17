	SUBROUTINE SPDR(IJC,IKC,ICLOCK)

        include "com2d.h"
        include "comphot.h"

	REAL*8 ratjo2, JX, XSRFU, JX39, JZ39, J39
	COMMON/JALLTIME/JX(PH$,L$,Z$)

C  wavelength dependent J-arrays (for diagnostics ONLY)

        COMMON/CJX39/JX39(PH$,5,L$,Z$), JZ39(18,PH$,5,L$,Z$),
     >                J39(PH$,5,L$,Z$)


C  commons for O3, O2, NO2  cross sections*RFLUX for consistent photolysis/heating rates: 
C  also, use xsectno2(IL$) which has just been loaded in CROSEC for the temperture at the current grid pt

        COMMON/CXSRFU/XSRFU(4,L$,Z$)
        COMMON/CXNO2/xsectno2t(IL$,201), xsectno2(IL$)


c   solar flux data from Judith Lean:
C   index 1 = time;  2-40 are the 39 model wavelength bins (in ph/cm^2/sec); 41=TOA flux (W/m2)
C    SJFLUXDAY(41) are TD values interpolated to the current day ;  FLXA(41) are the long term avgs
C
C     if ISOLCYC = 0, uses constant avg solar flux
C     if ISOLCYC = 1, uses time dependent flux data for current day

        REAL*8 sjfluxday, hjfluxday, flxa, scfac
        COMMON/CSJFLUXD/ SJFLUXDAY(41), HJFLUXDAY(12), FLXA(41), ISOLCYC


        SAVE


C     SOLAR PHOTODISSOCIATION 
C     RATES FOR CHEMISTRY MODULE IN POINT CALL FORMAT.

C     NO CALCULATION OF J COEFFICIENTS OR REDUCED FLUXES IF THE SUN IS
C     NOT UP AT THE CURRENT POINT.
      IF(ZGRZ(IJC).LE.0.0)THEN
        DO 200 JL=1,IL$
200	RFLUX(JL,IJC,IKC) = 0.0
        DO 220 JP=1,PH$
220	JX(JP,IJC,IKC) = 0.0
        RETURN
      END IF

C     SOLAR FLUX AND ABSORPTION CROSS-SECTIONS (PHOTODISSOCIATION DATA)
C     FOR SPECIES IN OPTICAL DEPTH CALCULATION ARE DEFINED AT IL$
C     WAVELENGTHS. 

C     INITIALIZE ARRAYS FOR J-COEFFICIENTS AND JX39(PH$,5,L$,Z$)
	DO 320 JP=1,PH$
           JX(JP,IJC,IKC) = 0.0

           do 321 ijj0=1,5
 321	     JX39(JP,ijj0,ijc,ikc) = 0.0
320 	CONTINUE


C     J COEFFICIENTS FOR EACH PROCESS WHEN GRAZING HEIGHT > 0KM.
C     COMPUTE EACH J-COEFFICIENT OVER ITS DEFINED WAVELENGTH RANGE.
C     I IS THE WAVELENGTH INDEX; J, THE CROSS-SECTION INDEX. DO NOT
C     COMPUTE J-COEFFICIENTS CORRESPONDING TO OMITTED PROCESSES.
	DO 400 JL=1,IL$
	DO 400 JP=1,PH$
        JX(JP,IJC,IKC) = JX(JP,IJC,IKC) + XSECT(JP,JL)*RFLUX(JL,IJC,IKC)
c	if(ijc.eq.9 .and. jp.eq.30 .and. ikc.eq.4)then
c	print *,' jl=',jl,' j=',jx(jp,ijc,ikc),
c     *  ' xsect(1,jl)=',xsect(jp,jl),' rflux=',rflux(jl,ijc,ikc)
c	endif
400	CONTINUE


C  LIFE2 - save J's for 5 wavelength bins: JX39(PH$,5,L$,Z$)
C   1=Lyman alpha;   2=1695-1905;   3=1905-2299;  4=2299-2857;   5=2857-7350

        DO 500 JP=1,PH$
              JL=1
              JX39(JP,1,IJC,IKC) = XSECT(JP,JL)*RFLUX(JL,IJC,IKC)


           do 510 JL=2,13
              JX39(JP,2,IJC,IKC) = JX39(JP,2,IJC,IKC)
     >                           + XSECT(JP,JL)*RFLUX(JL,IJC,IKC)
 510	   CONTINUE


           do 520 JL=14,22
              JX39(JP,3,IJC,IKC) = JX39(JP,3,IJC,IKC)
     >                           + XSECT(JP,JL)*RFLUX(JL,IJC,IKC)
 520	   CONTINUE


           do 530 JL=23,27
              JX39(JP,4,IJC,IKC) = JX39(JP,4,IJC,IKC)
     >                           + XSECT(JP,JL)*RFLUX(JL,IJC,IKC)
 530	   CONTINUE


           do 540 JL=28,39
              JX39(JP,5,IJC,IKC) = JX39(JP,5,IJC,IKC)
     >                           + XSECT(JP,JL)*RFLUX(JL,IJC,IKC)
 540	   CONTINUE
 500	   CONTINUE



	CALL NOPHOT(ENOP,IJC,IKC,ICLOCK)
	jx(16,ijc,ikc)=enop


C
C  9HB: (Feb 2012) - this works with and without solar cycle (ISOLCYC = 1/0)
C
C   w/ solar cycle, do JO2 properly (as in Charley's spdr_ce.f)
C
C   with Judith's solar flux data, use scaling factor for current day (SOLCYCR are scaling factors)
C      ie, divide by long term average:  SJFLUXDAY(41)/FLXA(41), 2-40 are the 39 model wavelength bins
C                                        this is the ONLY change needed here (replace SOLCYCR)
	O2JCOMP=0.0E0
	DO 410 JL=1,37
           scfac = SJFLUXDAY(JL+1)/FLXA(JL+1)

	   O2JCOMP = O2JCOMP + (XSECT(1,JL)*
     *         RFLUX(JL,IJC,IKC)/scfac)                     !  /SOLCYCR(JL,iday360,iysc))
c	if(ijc.eq.1)then
c           write(6,1105)jl,ikc,o2jcomp,xsect(1,jl),rflux(jl,ijc,ikc),
c     *  fluxratmult(jl)
c 1105	   format(' jl=',i5,' ikc=',i5,' o2jc=',1pe13.5,' xs=',
c     *  1pe13.5,' rfl=',1pe13.5,' frm=',1pe13.5)
c	endif

 410	   CONTINUE


C  NOTE: To do solar cycle, use Jx(1,IJC,IKC) from 400 loop above, which contains the solar cycle via RFLUX.
C        O2JCOMP is J[O2] with the solar cycle REMOVED IN RFLUX. 
C
C        So RATO2J is the ratio of J[O2] with solar cycle/J[O2] without solar cycle.
C           if NO solar cycle (ISOLCYC = 0) then RATO2J = 1.

C        and O2JINT is J[O2] from the look-up table.

	   RATO2J=1.0
c	if(ijc.eq.1)then
c           write(6,1103)j(1,ijc,ikc),rato2j,fluxratmult(19)
c	endif

	   IF(Jx(1,IJC,IKC) .GT. 1.E-32) RATO2J = Jx(1,IJC,IKC)/O2JCOMP
           if (ISOLCYC .eq. 0)  RATO2J = 1.

        jx(1,ijc,ikc)=o2jint(ijc,ikc)*RATO2J
c J(O2) from look-up table - 4/5/95


C   RATJO2 = ratio of J[O2] from look-up table (O2JINT) to that of JX(1) WITHOUT the solar cycle (O2JCOMP)

        ratjo2 = 1.d0
        IF (O2JCOMP .GT. 1.E-32) ratjo2 = o2jint(ijc,ikc)/O2JCOMP

        if (ISOLCYC .eq. 0  .and.  jx(1,ijc,ikc) .GT. 1.E-32)
     >      ratjo2 = o2jint(ijc,ikc)/jx(1,ijc,ikc)



c        if(ijc.eq.9 .and. ikc.eq.46)then
c           write(6,1101)(jx(30,ijc,ikout),ikout=1,z$)
c 1101      format(' j1301=',/,1p6e12.6,/,1p6e12.6,/,1p6e12.6,/,
c     *  1p6e12.6,/,1p6e12.6,/,1p6e12.6,/,1p6e12.6)
c           endif

C
C   load in XSECT*RFLUX here for O3, O2, NO2 xsections (for heating rate calculation):   XSRFU is in photons/second/nm
C   XSRFU(4,L$,Z$), 1=Ozone UV; 2=O3-vis;   3=O2;  4=NO2  ;   XSECT(PH$,IL$), xsectno2(IL$), RFLUX(IL$,L$,Z$)
C   to account for lookup table for J[O2], apply ratio RATJO2 to XSECT(1), XSECT(46) is NOT computed from lookup table
C                           note that RATJO2 DOES NOT have the solar cycle effect, the solar cycle is contained in RFLUX
C
C   also, integrate over wavelength to include 1/wvmid(IL$) factor,  start w/ Lyman alpha - 121.567nm 
C       wvmid(IL$) (in nm, in COMMON) which is the average of 1/lambda for each bin (defined in INPUT)
C
C   XSRFU are UNCORRECTED values - insolation correction is applied in PDFINSOL to be consistent w/ the J correction
C                                           first initialize
C
        do 557 ii=1,4
 557	   xsrfu(ii,IJC,IKC) = 0.d0

C                                                      ! start w/ Lyman alpha - 121.567nm
        DO 750 JL=1,37
           xsrfu(1,IJC,IKC) = xsrfu(1,IJC,IKC) + 
     >           (XSECT(2,JL) + XSECT(3,JL))*RFLUX(JL,IJC,IKC)/wvmid(JL)
 750	CONTINUE


        DO 760 JL=38,39
           xsrfu(2,IJC,IKC) = xsrfu(2,IJC,IKC) + 
     >           (XSECT(2,JL) + XSECT(3,JL))*RFLUX(JL,IJC,IKC)/wvmid(JL)
 760	CONTINUE


        DO 770 JL=1,IL$
           xsrfu(3,IJC,IKC) = xsrfu(3,IJC,IKC) + 
     >   (XSECT(1,JL)*ratjo2 + XSECT(46,JL))*RFLUX(JL,IJC,IKC)/wvmid(JL)

           xsrfu(4,IJC,IKC) = xsrfu(4,IJC,IKC) + 
     >                      XSECTNO2(JL)*RFLUX(JL,IJC,IKC)/wvmid(JL)
 770	CONTINUE


	RETURN
	END
