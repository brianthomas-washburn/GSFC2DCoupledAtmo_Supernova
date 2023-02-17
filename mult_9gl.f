	SUBROUTINE MULTSC(ICLOCK)
C  this mult.f uses variable beta1 instead of beta (as in multme.f)

	include "com2d.h"
	include "comphot.h"
	
	
	DIMENSION D(Z$),TAU(Z$),TA(Z$),RA(Z$),LL(3,Z$)
	DIMENSION U(Z$),WW(Z$),XSS(3),XSA(3),XST(3),POL(3)
	DIMENSION WVLB(20),XSO2(IL$,Z$),XUSO2(Z$),DELTASR(Z$)

	REAL*8 REDUCE(IL$,Z$),QSO(IL$,Z$1),ylam(Z$1),qp(Z$1),tau0(Z$1)
	REAL*8 xk(Z$1),qpp(Z$1),alp(Z$1),bet(Z$1),a(Z$1),b(Z$1),xid(Z$1)
	REAL*8 xip(Z$1),cxxx(94,94),qvec(94),ansa(Z$1,94),ansb(Z$1,94)
	REAL*8 DTAU(Z$),RSCAT(IL$,Z$),REDFLUX(IL$,Z$1)
	REAL*8 TAUUSE,SCATTERUSE,TAUCALC,ABSORBUSE,XMU,TAU1,ALBUSE
	
	DATA BETA1/1.32293e-13/,POL/1.74,0.0,1.57/
	DATA WVLB/7300.,19*0./,ICNT/0/,XMU/.577/
C       BETA1 = CONSTANT FOR OBTAINING SCATTERING CROSS-SECTIONS =
C       128*PI**3/3*1.E-16 (POLARIZABILITIES IN 1.E24*CM**3,
C       WAVELENGTH IN ANGSTROMS).
C       xmu=.577 From Gaussian quadrature, Chandrasekhar, Radiative Transfer, p.62

        SAVE
	
	
	IZTWO=2*Z$1
	DO 110 I=1,Z$
110       ALBB(I)=0.0

	DO 1600 IK=1,Z$
	DO 1600 JL=1,IL$
	REDUCE(JL,IK)=0.0D0
	RSCAT(JL,IK)=0.0D0
1600	CONTINUE

c	print *,' press in mult before 1700 =',press

	DO 1700 IK=1,Z$
	DO 1700 IJ=1,L$
	DO 1700 JL=1,IL$
	FLUXMULT(JL,IJ,IK)=0.0E0
	TAUMULT(JL,IJ,IK)=1.E5
1700	CONTINUE
C INITIALIZE RSCAT AND REDUCE


	DO 1751 IJ=1,L$
C SET UP SOLAR ZENITH ANGLES FOR TABLE LOOK UP - 1/26/95
	   SZA2D(IJ)=CHID(IJ)
 1751	   CONTINUE

c	print *,' press in mult after 1700 =',press

c        print *,' il$=',il$
c        print *,' in multsc  flux=',(flux(ii),ii=1,il$)

	ALBB(1)=0.3
C SET ALBEDO AT GROUND EQUAL TO 0.3

	IU=Z$-1

C     CALCULATION WITH RADIATIVE TRANSFER.
	DO 1753 IJ=1,L$
	DO 1753 IK=1,Z$
C SET UP OVERHEAD COLUMN OZONE FOR TABLE LOOK UP - 1/26/95 - use diurnally varying O3 column
C                                                    NCOLD(2,18+3,L$,Z$) (COMMON); 1=Ozone;   2=NO

           OVHO3(IJ,IK) = NCOLD(1,iclock,IJ,IK)    ! NCOLGD(2,IJ,IK)
 1753	CONTINUE

C DO TABLE LOOK UP - 1/26/95
	CALL RFLUXINTERP

	XSA(1)=0.
C  XSA(1) IS FOR N2

       DO 2000 IJ=1,L$
       DO 2005 IK=1,Z$
       SNDZ(1,IK)=NCOLGD(1,IJ,IK)*(0.781/0.209)
C  SNDZ(1) IS FOR N2

       SNDZ(2,IK) = NCOLD(1,iclock,IJ,IK)    ! NCOLGD(2,IJ,IK)
C  SNDZ(2) IS FOR O3

       SNDZ(3,IK)=NCOLGD(1,IJ,IK)
C  SNDZ(3) IS FOR O2

       XNTAU(1,IK)=N2(IJ,IK)
C  XNTAU(1) IS FOR N2

       XNTAU(2,IK)=C(4,IJ,IK)
C  XNTAU(2) IS FOR O3

       XNTAU(3,IK)=C(3,IJ,IK)
C  XNTAU(3) IS FOR O2
2005       CONTINUE


c           print *,' sndz(1,*)=',(sndz(1,i),i=1,z$),' sndz(2,*)=',
c     *  (sndz(2,i),i=1,z$),' sndz(3,*)=',(sndz(3,i),i=1,z$)
c           print *,' xntau(1,*)=',(xntau(1,i),i=1,z$),
c     *  ' xntau(2,*)=',(xntau(2,i),i=1,z$),' xntau(3,*)=',
c     *  (xntau(3,i),i=1,z$)

	SFACMULT=0.0

       DO 2015 I=1,Z$
       IK=Z$-I+1

       IF(IK.NE.Z$ .AND. SFAC(IJ).LE.0.0)SFAC(IJ)=SFACMULT

       IF(SFAC(IJ) .LE. 0.0)GO TO 2025

       CALL SRBAND(IJ,IK)

       DO 2020 JL=1,IL$
c	if(ik .eq. z$) print *,'xo2(jl),jl',xo2(jl),jl
2020       XSO2(JL,IK)=XO2(JL)

2025	SFACMULT=SFAC(IJ)
	ZGRZMULT=ZGRZ(IJ)

c           print *,' ij=',ij,' ik=',ik,' sfacmult=',sfacmult,
c     *  ' zgrzmult=',zgrzmult

2015       CONTINUE

c           print *,' ij=',ij,' sfacmult=',sfacmult,
c     *  ' zgrzmult=',zgrzmult

C                              !! not sure O3 xsection (XO3) is used here since table look up is used below, 
C                              !! but load in here just in case - just use 240K value
       DO 2030 JL=1,IL$
       PHLUX=FLUX(JL)
       LAMMS=WVL(JL)
       XSA(2)= xsecttd(2,JL,121) + xsecttd(3,JL,121)         !!   XO3(I)
C  XSA(2) IS FOR O3


c	if(jl.eq.36 .or. jl.eq.4 .and. ij.eq.9)then
c	print *,' jl=',jl,' phlux=',phlux,' lamms=',lamms,
c    *  ' xsa(2)=',xsa(2)
c	endif

      DO 1200 JI=1,3
	WVLMICRON=WVL(JL)*1.E-4
C WVLMICRON is the wavelength in microns.
	TERMSS=3.916+(0.074*WVLMICRON)+(0.05/WVLMICRON)
	XSS(JI)=4.E-28/(WVLMICRON**(TERMSS))
c	if(jl.eq.36 .and. ij.eq.9)then
c	print *,' wvl=',wvl(jl),' wvlmicron=',wvlmicron,
c     *  ' ji=',ji,' xss(ji)=',xss(ji),' termss=',termss
c	endif
1200	CONTINUE
C Rayleigh scattering cross section from Brasseur and Solomon, Aeronomy
C of the Middle Atmosphere, p. 110


       XST(1)=XSA(1) + XSS(1)
C  XST(1) IS FOR N2

       XST(2)=XSA(2) + XSS(2)
C  XST(2) IS FOR O3
c	if(jl.eq.36 .and. ij.eq.9)then
c	print *,' xst2=',xst(2),' xsa2=',xsa(2),
c     *  ' xss2=',xss(2)
c	endif

        DO 1250 IK=1,Z$
        XUSO2(IK)=XSO2(JL,IK) + XSS(3)
c	if(ij .eq. 9 .and. ik .eq. 20)
c     *       print *,' ik=',ik,'jl=',jl,' wvl=',wvl(jl),
c     *  ' xuso2=',xuso2(ik),' xso2(jl,ik)=',xso2(jl,ik),' xss(3)=',
c     *	xss(3)
1250	CONTINUE
C  XUSO2 IS FOR O2

      ALB0=ALBB(1)
C     LOOP DOWN OVER ALTITUDE TO DETERMINE INTERMEDIATE QUANTITIES,
C     INCLUDING OPTICAL DEPTH TAU.
      DO 1300 I=1,Z$
      IK=Z$-I+1
      S=0.
       SFACMS=SFACMULT
       ZGRZMS=ZGRZMULT
       IF(ZGRZMS .LE. 0.0)SFACMS=1.e6

      DO 1320 JI=1,2 
	S=S+XST(JI)*SNDZ(JI,IK)
c	if(ij.eq.9 .and. jl.eq.4)then
c	print *,' s=',s,' xst(ji)=',xst(ji),' sndz(ji,ik)=',
c     *  sndz(ji,ik)
c	endif
1320	CONTINUE

	IF(JL.LT.5 .OR. JL.GT.18)THEN
       S=S + XUSO2(IK)*SNDZ(3,IK)
c	if(ij.eq.9 .and. jl.eq.4)then
c	print *,' s=',s,' xuso2(ik)=',xuso2(ik),' sndz(3,ik)=',
c     *  sndz(3,ik)
c	endif
	ENDIF


C  FOR SCHUMANN-RUNGE BANDS - NEED TO DIFFERENTIAL INTEGRATION
C  12/9/88
       IF(JL.GE.5 .AND. JL.LE.18)THEN
	S=S + XSS(3)*SNDZ(3,IK)
C S CONTAINS ALL OF OVERHEAD OPTICAL DEPTH EXCEPT SCHUMANN-RUNGE BANDS

       IF(IK.EQ.Z$)DELNCOL=.21*(PRESS(Z$)*1000./980.)*
     * (6.02E23/28.964)
       IF(IK.LT.Z$)THEN
       DELTAA=.21*(PRESS(IK)*1000./980.)*(6.02E23/28.964)
       DELTAB=.21*(PRESS(IK+1)*1000./980.)*(6.02E23/28.964)
       DELNCOL=DELTAA-DELTAB

c	if(ij.eq.2 .and. jl.eq.6)then
c	print *,' ik=',ik,' deltaa=',deltaa,' deltab=',
c     *  deltab,' press(ik)',press(ik),' press(ik+1)=',
c     *  press(ik+1)
c	endif
       ENDIF
       DELTASR(IK)=XSCHRUN(JL,IJ,IK)*DELNCOL

c	IF(IJ.EQ.9 .AND. IK.GE.20)THEN

c	WRITE(6,8182)DELTASR(IK),XSCHRUN(JL,IJ,IK),
c     *  DELNCOL,JL,IJ,IK
c8182	FORMAT(' DELTASR=',1PE11.3,' XSCHRUN=',1PE11.3,
c     *  ' DELNCOL=',
c     *  1PE11.3,' JL=',I5,' IJ=',I5,' IK=',I5)

c	WRITE(14,8182)DELTASR(IK),XSCHRUN(JL,IJ,IK),
c     *  DELNCOL,JL,IJ,IK

c	ENDIF

c	if(ij.eq.2 .and. jl.eq.6)then
c	print *,' ik=',ik,' delncol=',delncol,' deltasr(ik)=',
c     *  deltasr(ik),' xschrun(jl,ij,ik)',xschrun(jl,ij,ik)
c	endif
C COMPUTE SCHUMANN-RUNGE CONTRIBUTION AT EACH ALTITUDE

	ADDSR=0.0E0
	DO 7700 IKK=IK,Z$
	ADDSR=ADDSR+DELTASR(IKK)

c	if(ij.eq.9 .and. jl.eq.4)then
c	print *,' ikk=',ikk,' addsr=',addsr,' deltasr(ikk)=',
c     *  deltasr(ikk)
c	endif

c	IF(IJ.EQ.9 .AND. IK.EQ.20)THEN

c	WRITE(6,8184)IKK,DELTASR(IKK),ADDSR,
c     *  JL,IJ,IK
c 8184   FORMAT(' IKK=',I5,' DELTASR=',1PE11.3,' ADDSR=',1PE11.3,
c     *  ' JL=',I5,' IJ=',I5,' IK=',I5)

c	WRITE(14,8184)IKK,DELTASR(IKK),ADDSR,
c     *  JL,IJ,IK

c	ENDIF


7700	CONTINUE

c	IF(IJ.EQ.9 .AND. IK.EQ.20)THEN

c	WRITE(6,8188)S,JL,IJ,IK
c 8188   FORMAT(' S=',1PE11.3,
c     *  ' JL=',I5,' IJ=',I5,' IK=',I5)

c	WRITE(14,8188)S,JL,IJ,IK

c	ENDIF

	S=S+ADDSR
C ADD SCHUMANN-RUNGE CONTRIBUTION TO THE OTHER OPTICAL DEPTH COMPONENTS

c	IF(IJ.EQ.9 .AND. IK.EQ.20)THEN

c	WRITE(6,8186)S,ADDSR,
c     *  JL,IJ,IK
c 8186   FORMAT(' S=',1PE11.3,' ADDSR=',1PE11.3,
c     *  ' JL=',I5,' IJ=',I5,' IK=',I5)

c	WRITE(14,8186)S,ADDSR,
c     *  JL,IJ,IK

c	ENDIF

c	if(ij.eq.9 .and. jl.eq.4)then
c	print *,' ik=',ik,' addsr=',addsr,' deltasr(ik)=',
c     *  deltasr(ik),' xss3=',xss(3),' sndz(3,ik)=',sndz(3,ik)
c	endif

	ENDIF

	TAU(IK)=S*SFACMS

c	if(jl.eq.4 .and. ij.eq.9)then
c	if(jl.eq.36 .and. ij.eq.9 .and. ik.eq.20)then
c	print *,' ik=',ik,' s=',s,' sfacms=',sfacms,' tau=',tau(ik),
c     *  ' day=',day
c	endif

c	if(ij.eq.2 .and. jl.eq.6)then
c	print *,' ik=',ik,' s=',s,' sfacms=',sfacms
c	endif

C This formulation taken from R. S. Stolarski  -  2/19/92
c Calculate fraction by which solar flux is reduced due to ozone
c and O2 absorption and Rayleigh scattering at each wavelength

ccc        TAUUSE=DFLOAT(TAU(IK))
	TAUUSE = DBLE(TAU(IK))
	REDUCE(JL,IK)=DEXP(-TAUUSE)

c        if(jl.eq.10 .and. ik.eq.20 .and. ij.eq.9)then
c           print *,' jl ik ij tau(ik) reduce(jl,ik) s sfacms',
c     *  jl,ik,ij,tau(ik),reduce(jl,ik),s,sfacms
c           write(14,5101)jl,ik,ij,tau(ik),reduce(jl,ik),s,
c     *  sfacms
 5101      format(' jl=',i3,' ik=',i3,' ij=',i3,' tau=',
     *  1pe11.3,' reduce=',1pe11.3,' s=',1pe11.3,
     *  ' sfacms=',1pe11.3)
c           endif

C This formulation taken from R. S. Stolarski  -  2/19/92
C Calculate the source function which contributes to the scattering
C component
	SCATTER=XNTAU(1,IK)*XSS(1) + XNTAU(3,IK)*XSS(3)
	SCATTERUSE=DBLE(SCATTER)
	QSO(JL,IK)=REDUCE(JL,IK)*SCATTERUSE

C Set up optical depth in layer
	IF(IK.NE.Z$)TAUINT=ABS(TAU(IK)-TAU(IK+1))/SFACMS
	IF(IK.EQ.Z$)TAUINT=TAU(IK)/SFACMS
	DTAU(IK)=DBLE(TAUINT)

c	if(jl.eq.4 .and. ij.eq.9)then
c	print *,' ik=',ik,' dtau=',dtau(ik),' qso=',qso(jl,ik),
c     *  ' scatter=',scatter,' xntau1=',xntau(1,ik),' xss1=',xss(1),
c     *  ' xntau3=',xntau(3,ik),' xss3=',xss(3),' tau(ik)=',
c     *  tau(ik),' tau(ik+1)=',tau(ik+1),' sfacms=',sfacms
c	endif



1300	CONTINUE
      S=0.
       SFACMS=SFACMULT
       ZGRZMS=ZGRZMULT
       IF(ZGRZMS .LE. 0.0)SFACMS=1.e6
C
	S=XST(1)*(0.78*1013.*(1000./980.)*(6.02E23/28.964)) + S
	S=XST(2)*SNDZ(2,1) + S
	S=XUSO2(1)*(0.21*1013.*(1000./980.)*(6.02E23/28.964)) + S
C
	TAUGROUND=S*SFACMS

C This formulation taken from R. S. Stolarski  -  2/19/92
c Calculate fraction by which solar flux is reduced due to ozone
c and O2 absorption and Rayleigh scattering at each wavelength
	TAUUSE=DBLE(TAUGROUND)
	REDFLUX(JL,1)=DEXP(-TAUUSE)

C     DETERMINE REDUCED FLUXES AT EACH ALTITUDE FOR THIS WAVELENGTH.
      DO 1150 IK=1,Z$
      FLUXMULT(JL,IJ,IK)=0.0
       IF(ZGRZMULT .LE. 0.0E0)THEN
         TAUMULT(JL,IJ,IK)=1.e5
         GO TO 1150
       END IF
       TAUMULT(JL,IJ,IK)=TAU(IK)/SFACMS
C   TAU(IK) IS REALLY THE OPTICAL DEPTH TIMES THE SFAC.  THUS IN ORDER
C   TO GET THE OPTICAL DEPTH, TAUMS, USED IN CALCULATION OF THE DIURNAL
C   AVERAGE SOLAR FLUX, THIS IS THE CORRECT EXPRESSION.

	DIR=SNGL(REDUCE(JL,IK))
      IF(JL.LE.4)THEN
	 FLUXMULT(JL,IJ,IK)=PHLUX*(DIR)
      ENDIF
      IF(JL.GE.5 .AND. JL.LE.19)THEN
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFLUX(IJ,IK,JL-4)
      ENDIF
      IF(JL.EQ.20)THEN
	 J1=16
	 J2=17
	 J3=18
	 J4=19
	 SFUSE=(SFLUX(IJ,IK,J1)+SFLUX(IJ,IK,J2)+SFLUX(IJ,IK,J3)+
     *    SFLUX(IJ,IK,J4))/4.
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFUSE
      ENDIF
      IF(JL.EQ.21)THEN
	 J1=20
	 J2=21
	 J3=22
	 J4=23
	 SFUSE=(SFLUX(IJ,IK,J1)+SFLUX(IJ,IK,J2)+SFLUX(IJ,IK,J3)+
     *    SFLUX(IJ,IK,J4))/4.
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFUSE
      ENDIF
      IF(JL.EQ.22)THEN
	 J1=24
	 J2=25
	 J3=26
	 J4=27
	 SFUSE=(SFLUX(IJ,IK,J1)+SFLUX(IJ,IK,J2)+SFLUX(IJ,IK,J3)+
     *    SFLUX(IJ,IK,J4))/4.
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFUSE
      ENDIF
      IF(JL.EQ.23)THEN
	 J1=28
	 J2=29
	 J3=30
	 J4=31
	 SFUSE=(SFLUX(IJ,IK,J1)+SFLUX(IJ,IK,J2)+SFLUX(IJ,IK,J3)+
     *    SFLUX(IJ,IK,J4))/4.
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFUSE
      ENDIF
      IF(JL.EQ.24)THEN
	 J1=32
	 J2=33
	 J3=34
	 J4=35
	 SFUSE=(SFLUX(IJ,IK,J1)+SFLUX(IJ,IK,J2)+SFLUX(IJ,IK,J3)+
     *    SFLUX(IJ,IK,J4))/4.
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFUSE
      ENDIF
      IF(JL.EQ.25)THEN
	 J1=36
	 J2=37
	 J3=38
	 J4=39
	 SFUSE=(SFLUX(IJ,IK,J1)+SFLUX(IJ,IK,J2)+SFLUX(IJ,IK,J3)+
     *    SFLUX(IJ,IK,J4))/4.
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFUSE
      ENDIF
      IF(JL.EQ.26)THEN
	 J1=40
	 J2=41
	 J3=42
	 J4=43
	 SFUSE=(SFLUX(IJ,IK,J1)+SFLUX(IJ,IK,J2)+SFLUX(IJ,IK,J3)+
     *    SFLUX(IJ,IK,J4))/4.
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFUSE
      ENDIF
      IF(JL.EQ.27)THEN
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFLUX(IJ,IK,44)
      ENDIF
      IF(JL.EQ.28)THEN
	 J1=45
	 J2=46
	 J3=47
	 SFUSE=(SFLUX(IJ,IK,J1)+SFLUX(IJ,IK,J2)+SFLUX(IJ,IK,J3))/3.
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFUSE
      ENDIF
      IF(JL.EQ.29)THEN
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFLUX(IJ,IK,48)
      ENDIF
      IF(JL.EQ.30)THEN
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFLUX(IJ,IK,49)
      ENDIF
      IF(JL.EQ.31)THEN
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFLUX(IJ,IK,50)
      ENDIF
      IF(JL.EQ.32)THEN
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFLUX(IJ,IK,51)
      ENDIF
      IF(JL.EQ.33)THEN
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFLUX(IJ,IK,52)
      ENDIF
      IF(JL.EQ.34)THEN
	 J1=53
	 J2=54
	 J3=55
	 SFUSE=(SFLUX(IJ,IK,J1)+SFLUX(IJ,IK,J2)+SFLUX(IJ,IK,J3))/3.
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFUSE
      ENDIF
      IF(JL.EQ.35)THEN
	 J1=56
	 J2=57
	 J3=58
	 J4=59
	 SFUSE=(SFLUX(IJ,IK,J1)+SFLUX(IJ,IK,J2)+SFLUX(IJ,IK,J3)+
     *    SFLUX(IJ,IK,J4))/4.
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFUSE
      ENDIF
      IF(JL.EQ.36)THEN
	 J1=60
	 J2=61
	 J3=62
	 J4=63
	 SFUSE=(SFLUX(IJ,IK,J1)+SFLUX(IJ,IK,J2)+SFLUX(IJ,IK,J3)+
     *    SFLUX(IJ,IK,J4))/4.
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFUSE
      ENDIF
      IF(JL.EQ.37)THEN
	 J1=64
	 J2=65
	 J3=66
	 J4=67
	 SFUSE=(SFLUX(IJ,IK,J1)+SFLUX(IJ,IK,J2)+SFLUX(IJ,IK,J3)+
     *    SFLUX(IJ,IK,J4))/4.
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFUSE
      ENDIF
      IF(JL.EQ.38)THEN
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFLUX(IJ,IK,68)
      ENDIF
      IF(JL.EQ.39)THEN
	 FLUXMULT(JL,IJ,IK)=PHLUX*SFLUX(IJ,IK,69)
      ENDIF

      IF(FLUXMULT(JL,IJ,IK).LT.1.e-20.AND.
     .  FLUXMULT(JL,IJ,IK).NE.0.)FLUXMULT(JL,IJ,IK)=1.e-20

c	if(jl.eq.38 .and. ij.eq.9)then
c	print *,' ij=',ij,' ik=',ik,' jl=',jl
c	print *,' dir=',dir,' scat=',scat,
c     *  ' phlux=',phlux,' tau=',tau(ik),' sfacms=',
c     *  sfacms
c	print *,' fluxmult=',fluxmult(jl,ij,ik),' taumult=',
c     *  taumult(jl,ij,ik)
c	endif

c	if(ij.eq.2 .and. jl.eq.6)then
c	print *,' ij=',ij,' ik=',ik,' jl=',jl
c	print *,' fluxmult=',fluxmult(jl,ij,ik),' taumult=',
c     *  taumult(jl,ij,ik)
c	endif

 1150 CONTINUE
c	if(jl.eq.4 .and. ij.eq.9)then
c	if(1.gt.0)stop
c	endif
2030       CONTINUE
2000       CONTINUE

c           print *,' il$=',il$
c	print *,' fluxmult(*,9,20)=',
c     *  (fluxmult(ii,9,20),ii=1,il$),il$
c	print *,' taumult=',taumult

      RETURN
      END

