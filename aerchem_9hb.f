C
C
        SUBROUTINE AERCHEM

C
C  routine to do the AER FAST CHEMISTRY CALCULATIONS    - EF 12/8/03
C     code is just been taken out of MAIN and put into this separate routine
C


         include "com2d.h"
         include "com_aerd.h"
         include "com_aerg.h"
         include "com_aers.h"



C  Common for mapping into GSFC code - NOTE: CAN't have any L$, Z$ here in these COMMONS
C  

         COMMON/CAER1/ISPMAP(NFSP)
         COMMON/CAER2/DTIMEA(NTIME), DTIMEB(L$,NTIME), 
     >                PMTIME(NTIME), DLITE(NHT)
         COMMON /AERHET1/KHETA(NKR,18,46), GAMAER(NKR,18,46)

         REAL KHETA, GAMAER, PMTIME, DLITE
         DOUBLE PRECISION DTIMEA, DTIMEB


C  JAVAER are the diurnally avg J-coeffs,      JAER are the diurnally varying Js
C  CNDCA = CN array for 18+3 diurnal time steps for up to 1 daily iterations, CNDCA is in Num dens., NFSP=53
C
         REAL*8 JAVAER(NJR+3,L$,Z$), CNDCA(S$,ntime+3,L$,Z$), JX39, JZ39
         REAL*8 JAER(NJR+3,NTIME,L$,Z$), OXG(L$,Z$), ch4out(L$,Z$), J39


Cj  common for diurnally avged J-coeffs - ONLY USE if dayime avg J's are needed
cj      COMMON/JAER1/JAVAER(NJR,18,46), AJOUT(2,ntime,18), 
cj     >     AJAVDAY(3,NJR,18,46)
ccccccc   REAL AJOUT(2,ntime,18), AJAVDAY(3,NJR,18,46)
cj
cj      COMMON/JAER1/AJOUT(2,ntime,18)


         COMMON/CHIDSFAC/CHIDA(18,L$),SFACA(18,L$)


         COMMON/OXCC/OXTR(L$,Z$)
         REAL OXB(L$,Z$)


C  wavelength dependent J-arrays (for diagnostics ONLY), from SPDR and PDFINSOL, TIMEAVJ

         COMMON/CJX39/JX39(PH$,5,L$,Z$), JZ39(18,PH$,5,L$,Z$),
     >                 J39(PH$,5,L$,Z$)



C  indicies for mapping AER FAST species indicies into GSFC indicies, 
C    ISPMAP(NFSP=54) in COMMON,  for the AER hydrocarbons and Iodine family 
C    that GSFC doesn't have, just set index to 0
C
      DATA ISPMAP/2,  1, 41,  4, 12, 13, 14, 16, 22, 24, 23,  0,  0,  0, 
     >            0,  0,  0,  0,  0,  0, 27, 28, 29, 25, 30, 26, 62, 65,
     >           61, 64, 63,  9,  5,  6,  7,  8, 38, 42, 10, 89, 45, 44,
     >           46, 68, 47, 79, 43, 60,  0,  0,  0,  0,  0,  0/


C
C  *******************   DO  CHEMISTRY w/ NEW AER DIURNAL CYCLE  **********************
C
C
c  CALL CLOCKAER(DTIMEA,PMTIME,DLITE,DL90,JHT)

         CALL SETUPA(DTIMEB)
 


C  before photolysis calc, initialize final INSOLATION CORRECTED photolysis array JZ
C  and the corrected XSRF array here (for coupled model heating rate calc consistent w/ photolysis)

         IF (icoup .EQ. 1) CALL PHOTHEATIN

C
C  ICLOCK/JCLOCK ARE the diurnal time step counters, which go from 1,NTIME (18)
C
C  NOTE:  SUNSET=7 and SUNRIS=12 ALWAYS (at ALL latitudes, days of the year)
C                                      
        
       DO 500 ICLOCK=1,18   ! SUNSET  - skip over nightime

          if (iclock .ge. 8  .and.  iclock .le. 11) goto 500

C                              DIURNALLY varying O3, NO columns in NCOLD(2,18+3,L$,Z$) (COMMON)
          CALL COLDIUR(ICLOCK)

          DO 510 IJ=1,L$
          CHID(IJ)=CHIDA(ICLOCK,IJ)
          SFAC(IJ)=SFACA(ICLOCK,IJ)
 510      CONTINUE

              CALL MULTSC(ICLOCK)
ccc                type *,'multsc'
              CALL RADT
ccc                type *,'radt'

		DO 1000 IJ=1,L$
		DO 1000 IK=1,Z$
			CALL FILL_COL(IJ,IK,ICLOCK)
			CALL CROSEC(IJ,IK)
			CALL SPDR(IJ,IK,ICLOCK)	
1000	CONTINUE
cccccc        print *,' ICLOCK(MAIN)=',ICLOCK

C
C  Now do INSOLATION CORRECTION on diurnally varying J-values: 
C   JX(PH$,L$,Z$) - get new array JZ(18,PH$,L$,Z$), which gets mapped for the AER chemistry in JMAP 
C                                      also do for heating rates: XSRFU(4,L$,Z$)-> XSRF(18,4,L$,Z$)
              CALL PDFINSOL(ICLOCK)

 500    CONTINUE


C
C  Call NATICE BEFORE AER FAST CHEMISTRY to set up NAT and ICE PSCs
C      HETCHEM reactions now computed in HETCHEM called within the DIURNAL CHEMISTRY
C
C         (OLD: correct problem w/ defining KHs after SOLVER call (EF,3/99))
C

         CALL NATICE(idor)


cccccc       CALL HETCHEM



C   load  CN array into AER's ALLSSP and ALLFSP arrays for FAST CHEMISTRY - use NOONRAT array

          CALL FILLA
cccc           type *,'FILLA'


cc   load Ox into OXG(L$,Z$) array for use in FAST CHEM

cc             do 811 ik=1,Z$
cc             do 811 ij=1,L$
cc 811            oxg(ij,ik) = cn(39,ij,ik)


C  driver for AER fast chemistry 
c                                                            !  CNDCA(S$,ntime+3,L$,Z$) is NOT in COMMON

          CALL FASTCHAER(DTIMEB,TFD,L$,Z$,S$,ikmes,JAVAER,JAER,CNDCA,
     >       iday360, iyr, aerosol, nataer, aerice, PH$, IL$, JZ39, J39)

ccccc          CALL FASTCHAER(DTIMEB,TFD,L$,Z$,S$,ikmes,JAVAER,JAER,CNDCA) 
ccccc             type *,'FASTCHEM'



C  store noontime ratios into array NOONRAT(S$,L$,Z$) for FAST CHEM initialization on next day
C  load updated diurnally avged fast species into C and CN arrays
C      and sum up CHx and HOx families for transport (NO additional chemistry) 
C      also load CNDCA into CNDC which is in COMMON

          CALL UNFILLA(CNDCA)

ccccccc      type *,'UNFILLA'



C  for coupled model, compute diurnally avged heating rates consistent w/ photolysis,
C  use xsect*rflux arrays from SPDR with INSOLATION correction applied in PDFINSOL
C  and updated diurnal variations of ozone, NO2, adjusted in UNFILLA:
C  these will be used in the call to the COUPLED model on the following day:
C
          IF (icoup .EQ. 1) CALL PHOTHEAT(DTIMEB)




c  load diurnally averaged Js from the AER array, JAVAER(NJR+3,L$,Z$) into the GSFC array J(PH$,L$,Z$) array
C    and the diurnally varying Js - JAER(NJR+3,NTIME,L$,Z$) - into the GSFC array JDC(PH$,NTIME,L$,Z$) array

          CALL UNJMAP(JAVAER, JAER, NJR, NTIME)
cccccc             type *,'UNJMAP' 



cccc   JAVAER(NJR+3,L$,Z$), AJOUT(2,ntime,L$),  AJAVDAY(3,NJR,L$,Z$)
cccj         write (17) javaer, ajout
CCCC         write (17) javaer, ajout, ajavday



C  ALL TRANSPORT (ADVECTION + DIFFUSION) Now done from separate routine XTRANSPORT (MAIN is less complicated)
C      also, store Ox values BEFORE transport, and then take the difference after transport to get 
C      the transport term, in num den/sec;    OXB(L$,Z$),   OXTR(L$,Z$)
C
             do 954 ik=1,Z$
             do 954 ij=1,L$
 954            oxb(IJ,IK) = c(39,IJ,IK)


          CALL XTRANSPORT



             do 964 ik=1,Z$
             do 964 ij=1,L$
 964            oxtr(IJ,IK) = (c(39,IJ,IK) - oxb(IJ,IK))/DT8


ccth  -  UPDATE CH4 and H2O AFTER TRANSPORT - NO CHEMISTRY (NO SOLVER)
cc
cc        do 904 ik=1,Z$
cc        do 904 ij=1,L$
cc           cn(18,IJ,IK) = c(18,IJ,IK)
cc           cn(15,IJ,IK) = c(15,IJ,IK)
cc 904    CONTINUE
cc
cc       	do 924 ik=1,Z$X
cc     	do 924 ij=1,L$
cc          cn58(1,ij,ik) = c58(1,ij,ik) 
cc          cn58(2,ij,ik) = c58(2,ij,ik) 
cc 924    CONTINUE
cc
C
C  use diurnally averaged J's and radicals for SOLVER (long-lived species only)
C    also need diurnally VARYING Js and radicals where necessary  -  the IJ,IK loop is now in the SUBROUTINE
C  
ccccccccc   do 900 ik=1,Z$
ccccccccc   do 900 ij=1,L$
ccccccccc		CALL SOLVER(IJ,IK)
cccccccccc            CALL OX(IJ,IK)
ccccccccccc 900	continue


C    Ox-39 now updated in SOLVER w/ ozone as a fast species - May 2003
C         routine OX just computes Ox Prod and loss diagnostically using counting scheme
C
            CALL SOLVER(DTIMEB, NTIME)
            CALL OX(DTIMEB, NTIME)


	RETURN
         END

