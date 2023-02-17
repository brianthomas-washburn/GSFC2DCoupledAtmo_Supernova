
	SUBROUTINE SETDAILY


C  THIS routine is for the new circulation (BASE8A), and is a reduced version 
C       of the old CONTIN routine. It only has the set up for basic model
C       parameters, eg, pressure, zkm arrays, major constituent set up
C       and is CALLED EVERY DAY to INTERPOLATE TEMPS AND TRANSPORT ARRAYS TO DAILY VALUES
C       *** ALSO, NEED TO RE-COMPUTE V* FROM INTERPOLATED W* FIELD EVERYDAY .........
C

        include "com2d.h"


        REAL IDAYX(74), idaypr(72), tempcg(45,76), TPROBST(211)
        REAL TPROB1(211), TPROB2(211), LPROB1(181), LPROB2(181)
        REAL TPROB_PSCD(211,45,76), TPROBHYB(211,45,76), TPROBFR(211)
        REAL TPROBC(211,45,76)
        REAL KZZ0(45,89), KZZDAY(45,89), KZZTH(45,89), EKZZA(L$,Z$X1)

        DIMENSION tzonalw(45,76),iadjw(45,76),ffow(45,76),tadjzmw(45,76)
        DIMENSION ic1w(45,76), ic2w(45,76), id1w(45,76), id2w(45,76)
ccccccc        DIMENSION tprobfrw(211,45,76)


C      common block for the coupled model dynamics,  from MAIN

        COMMON/TRCOUP/tcoup(L$,Z$X), ubard(L$,Z$X), kyycoup(L$+1,Z$X), 
     >              kzzcoup(L$,Z$X+1), wcoup(L$,Z$X+1), DRAGOUT(7,L$,Z$)


C   common for tropospheric Kzz adjustment (45x89 grid)

        COMMON/CFKZZ/ylatf(45), zaltf(89), timefz(14)


c  common for coupled model temperature - NCEP offset (lat-hgt-season) from 1985-2005:  (from TEMPIN)

        COMMON/CTOFFS/ xlatoffs(45), zzoffs(76)
ccccccc        COMMON/CTOFFS/ TEMPOFFS(360,45,76), xlatoffs(45), zzoffs(76)


C   common of Coupled model longitudinal T', from PWave geop hgts; interpolated to chemistry grid
C
        COMMON/CTPDF/tcc(72,L$,Z$)


	SAVE


C   this is for the FIXED MODEL:

      IF (ICOUP .NE. 1) THEN

C    on 1st day of each year, read in appropriate transport fields
C       also get UBAR, DELF, QBARY fields for output - UBARCG, EPCG, QYCG(L$,Z$X,74)

        IF (iday360 .eq. 1) then
             CALL GETDYN
             CALL UBARDIAG
        ENDIF
C
c
C   Interpolate temps, w*, Kyy, Kzz to daily values (from the 5-day means) for the model grid
c       IDAY360 goes from 1-360
C
C   NOTE: the interpolation here is INDEPENDENT OF INTERANNUAL/CLIM specification, so it's MUCH SIMPLER!!
C         SPECIFIC years to use are NOW defined in GETDYN/GETPROB
C
C
C    Note: Gwave and tropospheric Kzzs are now combined in KZZTOT(L$,Z$X1,74) ,  EKZZ(L$,Z$X1)
C
C    always use WALL(L$,Z$X1,74), KYYALL(L$+1,Z$X,74), KZZTOT(L$,Z$X1,74), TEMPALL(L$,Z$X,74) 
C      read in for CURRENT YEAR, interpolate to: 
C
C   TEMP(L$,Z$X) is REAL*8,  EKYY(L$+1,Z$X),  W1(L$,Z$X1), EKZZ(L$,Z$X1),   UBARD(L$,Z$X), DRAGOUT(7,L$,Z$)
C               

         do 9772 ijk=1,74
9772     IDAYX(ijk) = (ijk-1)*5.- 2.

         ILM = INT((IDAY360 + 2)/5) + 1


         DO 8505 IK=1,Z$X
         DO 8505 IJ=1,L$
           TEMP(IJ,IK) = (TEMPALL(IJ,IK,ILM+1) - TEMPALL(IJ,IK,ILM))/5.*
     >       (IDAY360 - IDAYX(ILM)) + TEMPALL(IJ,IK,ILM)

           UBARD(IJ,IK) = (UBARCG(IJ,IK,ILM+1) - UBARCG(IJ,IK,ILM))/5.*
     >       (IDAY360 - IDAYX(ILM)) + UBARCG(IJ,IK,ILM)
 8505   CONTINUE


       DO 8507 IK=1,Z$X
       DO 8507 IJ=1,L$+1
           EKYY(IJ,IK) = (KYYALL(IJ,IK,ILM+1) - KYYALL(IJ,IK,ILM))/5.*
     >       (IDAY360 - IDAYX(ILM)) + KYYALL(IJ,IK,ILM)
 8507   CONTINUE


       DO 8506 IK=1,Z$X1
       DO 8506 IJ=1,L$
         W1(IJ,IK) = (WALL(IJ,IK,ILM+1) - WALL(IJ,IK,ILM))/5.*  
     >       (IDAY360 - IDAYX(ILM)) + WALL(IJ,IK,ILM)

         EKZZ(IJ,IK) = (KZZTOT(IJ,IK,ILM+1) - KZZTOT(IJ,IK,ILM))/5.*
     >       (IDAY360 - IDAYX(ILM)) + KZZTOT(IJ,IK,ILM)
 8506  CONTINUE


         DO 9505 IK=1,Z$
         DO 9505 IJ=1,L$
           DRAGOUT(2,IJ,IK) = (EPCG(IJ,IK,ILM+1) - EPCG(IJ,IK,ILM))/5.*
     >       (IDAY360 - IDAYX(ILM)) + EPCG(IJ,IK,ILM)

           DRAGOUT(1,IJ,IK) = (QYCG(IJ,IK,ILM+1) - QYCG(IJ,IK,ILM))/5.*
     >       (IDAY360 - IDAYX(ILM)) + QYCG(IJ,IK,ILM)
 9505   CONTINUE



C  do here for just using the first time period dynamics                               
C            
C  KYYALL(L$+1,Z$X,74), TEMPALL(L$,Z$X,74), TEMP(L$,Z$X) is REAL*8    ;   UBARD(L$,Z$X)
c    EKYY(L$+1,Z$X),    WALL(L$,Z$X1,74), W1(L$,Z$X1)
C
C  Note: Gwave and tropospheric Kzzs are now combined in KZZTOT(L$,Z$X1,74) ,  EKZZ(L$,Z$X1)
C

      IF (ICLIM .EQ. 3) THEN 

       DO 7005 IK=1,Z$X
       DO 7005 IJ=1,L$
         TEMP(IJ,IK) = TEMPALL(IJ,IK,1)
         UBARD(IJ,IK) = UBARCG(IJ,IK,1)
 7005   CONTINUE


       DO 7017 IK=1,Z$
       DO 7017 IJ=1,L$
         DRAGOUT(2,IJ,IK) = EPCG(IJ,IK,1)
 7017	 DRAGOUT(1,IJ,IK) = QYCG(IJ,IK,1)


       DO 7007 IK=1,Z$X
       DO 7007 IJ=1,L$+1
         EKYY(IJ,IK) = KYYALL(IJ,IK,1)  
cef          EKYY(IJ,IK) = 5.e10
 7007   CONTINUE

       DO 7006 IK=1,Z$X1
       DO 7006 IJ=1,L$
         W1(IJ,IK) = WALL(IJ,IK,1) 
         EKZZ(IJ,IK) = KZZTOT(IJ,IK,1) 
cef          EKZZ(IJ,IK) = 7.e4
 7006  CONTINUE

      ENDIF


C   do tropospheric Kzz adjustment identical to coupled model for chemistry ONLY
C   first interpolate to 45x89 grid (do log interpolation in vertical);   timefz(14)
C   EKZZ(L$,Z$X1) ;  LAT4(L$), ZALTE(Z$X1) ->  KZZ0(45,89),  ylatf(45), zaltf(89)
C
C   interpolate adjusted KZZDAY(45,89) back to EKZZA(L$,Z$X1); KZZTH(45,89) not used here
C
C
CCWCONV276 - DON'T DO this, not as good as original Kzz (too much ozone, age too old)
CCWCONV276                 probably can adjust Kzz here somehow though, maybe use CCM OTs?
CCWCONV276
CCWCONV276
CCWCONV276        CALL BINTERP(lat4, L$, zalte, Z$X1, EKZZ,
CCWCONV276     >               ylatf, 45, zaltf, 89, 0, 1, KZZ0)
CCWCONV276
CCWCONV276
CCWCONV276        CALL KZZADJ(IDAY360, TIMEFZ, ylatf, zaltf, KZZ0, KZZTH, KZZDAY)
CCWCONV276
CCWCONV276
CCWCONV276        CALL BINTERP(ylatf, 45, zaltf, 89, KZZDAY,
CCWCONV276     >               lat4,  L$, zalte, Z$X1, 0, 1, EKZZA)
CCWCONV276
CCWCONV276
CCWCONV276        DO 717 IK=1,Z$X1
CCWCONV276        DO 717 IJ=1,L$
CCWCONV276 717	   EKZZ(IJ,IK) = EKZZA(IJ,IK)
CCWCONV276


C
C  Now write out every 5 days on days 5,10,15,20,25, etc for years 56-65 for checking
C
cc       if (iyr .ge. 56) then 
cc          if (mod(iday360,5) .eq. 0.0) then
cc             write(36) iyr, iday360, idaytot
cc             write(36) temp
cc             write(36) w1
cc             write(36) ekyy
cc          endif
cc       endif


      ENDIF

C
C   *************************  END  FIXED MODEL  INTERPOLATION  BLOCK  **********************
C
C
C
C
C   *************************  LOAD  IN  COUPLED  MODEL  DYNAMICS  HERE *********************
C
C   from COMMON block TRCOUP: passed from MAIN via subroutine XCOUPLED
C              tcoup(L$,Z$X),  kyycoup(L$+1,Z$X), kzzcoup(L$,Z$X+1), wcoup(L$,Z$X+1)
C
C
       IF (ICOUP .EQ. 1) THEN

           DO 7505 IK=1,Z$X
           DO 7505 IJ=1,L$
 7505	      TEMP(IJ,IK) = TCOUP(IJ,IK)

           DO 7506 IK=1,Z$X
           DO 7506 IJ=1,L$+1
 7506         EKYY(IJ,IK) = KYYCOUP(IJ,IK)

           DO 7507 IK=1,Z$X1
           DO 7507 IJ=1,L$
              W1(IJ,IK)   = WCOUP(IJ,IK) 
              EKZZ(IJ,IK) = KZZCOUP(IJ,IK)
 7507  CONTINUE

C   also interpolate current temps to 45x76 grid for PDFs below:    xlatoffs(45), zzoffs(76). TEMPCG(45,76)
C                                                  use LAT4(L$), ZALT(Z$X), TCOUP(L$,Z$X) which are REAL*4
        CALL BINTERP(lat4, L$, zalt, Z$X, TCOUP, 
     >               xlatoffs, 45, zzoffs, 76, 0, 0, TEMPCG)

       ENDIF


C  *************************  END  COUPLED  MODEL  LOADING  *****************************
C
c
C  
c   Now compute EKYZ from the interpolated TEMP, EKYY fields, 
C     and adjust the slope, EKYY and/or EKZZ as necessary to ensure the diffusion matrix is positive
C
C                       NOTE: this is used for the COUPLED model similarly to that done for the FIXED model

         CALL GETKYZ


c   
c   Now re-adjust the interpolated w* field, and compute v*, and load into w and v arrays for transport
C     Note: WALL array is already in cm/sec, so W1,  W and V  are in cm/sec here  *****


         CALL VWMASS

c
c  save mixing ratios; multiply by the new M after updating for new month, 
C    Still NEED THIS  WITH  DAILY UPDATES OF M, so that transport (which keys on 
C        mixing ratio) does not get out of whack
C  ALSO do for CN array, although this shouldn't make a difference, but do it for completeness  - 1/7/97
c           C(S$,L$,Z$), CN(S$,L$,Z$), C58(T$,L$,Z$X), CN58(T$,L$,Z$X), M(L$,Z$X)
C
	DO 210 IK=1,Z$
	DO 210 IJ=1,L$
	DO 210 IT1=1,S$
		CN(IT1,IJ,IK) = CN(IT1,IJ,IK)/M(IJ,IK)
		C(IT1,IJ,IK) = C(IT1,IJ,IK)/M(IJ,IK)
210	CONTINUE

        do 212 ik=1,Z$X
        do 212 ij=1,L$
        do 212 ist=1,ITRANS
    	   cn58(ist,ij,ik) = cn58(ist,ij,ik)/m(ij,ik)
 212	   c58(ist,ij,ik) = c58(ist,ij,ik)/m(ij,ik)
 
C
C  Now compute the following arrays from DAILY TEMPS which  are in COMMON:
C      M(L$,Z$X), ZKM(L$,Z$X), Q(L$,Z$X), DELTAZ(L$,Z$X), PRESS8(Z$X) and TEMP(L$,Z$X) are REAL*8
c

        DO 601 IK=1,Z$X
        DO 601 IJ=1,L$
 601	   m(ij,ik) = PRESS8(IK)*1.d3/(WTA*RD*TEMP(IJ,IK)*1.66d-24)


	DO 30 IJ=1,L$
        DO 5620 IK=1,Z$X

C   SCALE HEIGHT SHOULD PROPERLY DEPEND ON VIRTUAL TEMP.,  TEMP(L$,Z$X) is REAL*8 in COMMON 
C   THE TERM 29.3*1.e-3 IS K/MG IN UNITS OF KM*DEG**(-1),  SCALE HEIGHT IN (KM)

         H = 29.3*TEMP(IJ,IK)*1.D-3

C  integrating up to compute ZGP (geopotential height), from Wallace and Hobbs - EF (Aug. 2001)
C    this is very slightly different than the previous method of integrating the scale height
C           (~0.5% or less difference - not a problem)
C    for ik=1, just assume avg temp in layer is TEMP(1) , ie, just use H computed above
C    for ik=2 and above, compute avg temp in layer
C    PRESS(Z$X) array are the pressures of the grid points (box centers),   ZGP(L$,Z$X)

       if (ik .eq. 1) ZGP(IJ,IK) = H*alog(1013./press(1))

       if (ik .ge. 2) ZGP(IJ,IK) = 29.3*(TEMP(ij,ik-1)+TEMP(ij,ik))/2.
     >                *1.D-3*alog(press(ik-1)/press(ik)) + ZGP(IJ,IK-1)

       ZKM(IJ,IK) = RE*(1./(1.-ZGP(IJ,IK)/RE)-1.)
C                                                    CHANGE FROM GEOPOTENTIAL TO REAL HEIGHT, ZKM(L$,Z$X)

C   Q(IJ,IK)=-1./(H*1.e5) ;;    Q MUST BE IN CGS UNITS
C  --  NOTE:  Q(L$,Z$X) is in COMMON,  used  in  NEWDIF  *********

       Q(IJ,IK) = -1.d0/(H*1.d5)

C  DEFINE DELTA Z FOR GRID BOXES   11/9/88   DELTAZ(L$,Z$X) (in KM) is in COMMON 
C  - this is correct, DELTAZ is defined at the model grid points, and is the difference between the 
C    geopotential height of the top edge of the box minus the geop. hgt. of the bottom edge of the box
C    which is just the scale height in the box (ie, using the temp. at the grid point as the 
C    avg temp. in the box) times alog(presse(bottom edge)/presse(top edge)), which is just -DP(Z$X) array
C    computed from PRESSE(Z$X1), the pressures of the box edges, need to include neg. sign here since the 
C    DP(Z$X) array is always negative, alog(presse(top edge)/presse(bottom edge))   - all checked, 8/01 - EF

       DELTAZ(IJ,IK) = -H*DP(IK)

5620	CONTINUE
30	CONTINUE


C   convert back to number density with the new M, do for both C and CN arrays
C
	DO 230 IK=1,Z$
	DO 230 IJ=1,L$
	DO 230 IT1=1,S$
     	   CN(IT1,IJ,IK) = CN(IT1,IJ,IK)*M(IJ,IK)
 230	   C(IT1,IJ,IK) = C(IT1,IJ,IK)*M(IJ,IK)

        do 232 ik=1,Z$X
        do 232 ij=1,L$
        do 232 ist=1,ITRANS
     	   cn58(ist,ij,ik) = cn58(ist,ij,ik)*m(ij,ik)
 232	   c58(ist,ij,ik) = c58(ist,ij,ik)*m(ij,ik)
 
c 
c  Update  major  constituents - now for N2 use N2(L$,Z$) - REAL*8 in COMMON (9fk run)

	do 500 IK=1,Z$
	do 500 IJ=1,L$
c  N2
		n2(IJ,IK) = 0.79d0*m(IJ,IK)

c  O2
		c(3,IJ,IK) = 0.21d0*m(IJ,IK)
                cn(3,ij,ik) = c(3,ij,ik)

c  CO2  (set to 350 ppmv for 1990 from UNEP-94), NOT done now w/ CO2 calculated
c
c		c(20,IJ,IK) = 350.*1.e-6*m(IJ,IK)
c                cn(20,ij,ik) = c(20,ij,ik)

500      CONTINUE


C
c  if ozone is fixed to data, insert data here
c	IF(YSO3FIX)THEN
c		DO 220 IJ=1,L$
c		DO 220 IK=1,Z$
c			c(4,IJ,IK)=SBUVMRO3(I,IJ,IK)*M(IJ,IK)*1.E-6
c220		CONTINUE
c	ENDIF
c



C    **********************    BEGIN    PDF   LOADING   *********************************
C
C
C  IF USING PARCEL PDFs of TEMPERATURE and/or LATITUDE, THEN DO THIS LOOP, OTHERWISE, SET TO ZONAL MEANS
C
c  get proper TPROB/LPROB arrays for current YEAR and MODEL GRID, called ON FIRST DAY OF THE YEAR ONLY
c     and interpolate HERE to current daily value (this works OK, always get total of 1.0),  idaypr(72)
C     first expand input TPROB45/LPROB45 to 211/181 for 1st and 2nd time periods: 
C         TPROB1(211), TPROB2(211), LPROB1(181), LPROB2(181)
C
C        TPROB45(125,45,76,72) --> TPROB1(211)/TPROB2(211) --> TPROBD(211,45,76)
C        LPROB45(70,45,76,72) --> LPROB1(211)/LPROB2(211) --> LPROBD(181,45,76)
C
C                        ! if using coupled model longitudinal temp PDF, then ICOUPDF=1
      ICOUPDF = 0


      IF (ITPDF .EQ. 1   .OR.   ILPDF .EQ. 1) THEN


       IF (iday360 .eq. 1)  CALL GETPROB


C    now interpolate to current daily value - for days 1-2, and 358-360, just use 1st and last periods
C    TPROB1(211), TPROB2(211), LPROB1(211), LPROB2(211), TPROB45(125,45,76,72);  LPROB45(70,45,76,72)
C       itr1(45,76,72), itr2(45,76,72), ilr1(45,76,72), ilr2(45,76,72) are in COMMON
c       TPROBD(211,45,76); LPROBD(181,45,76)

        do 4501 ijk=1,72
 4501	   idaypr(ijk) = ijk*5.- 2.

  ! just use  1st and last periods for days 1-2, 358-360, ie, TPROB1/TPROB2 are identical

           ilmp1 = INT((iday360 + 2)/5)
           ilmp2 = ilmp1 + 1

           if (iday360 .le. 2) then 
              ilmp1 = 1
              ilmp2 = 1
           endif

           if (iday360 .ge. 358) then 
              ilmp1 = 72
              ilmp2 = 72
           endif



           do 4520 ik=1,76
           do 4520 ij=1,45
                                        ! regenerate full PDF(211)/(181) for 2 time periods, first initialize
              do 4522 itt=1,211
                 tprobd(itt,ij,ik) = 0.0
                 tprob1(itt) = 0.0
                 tprob2(itt) = 0.0
 4522         CONTINUE

            do 4524 itt=itr1(ij,ik,ilmp1), itr2(ij,ik,ilmp1)
 4524         TPROB1(itt) = TPROB45(itt-itr1(ij,ik,ilmp1)+1,ij,ik,ilmp1)

            do 4525 itt=itr1(ij,ik,ilmp2), itr2(ij,ik,ilmp2)
 4525	      TPROB2(itt) = TPROB45(itt-itr1(ij,ik,ilmp2)+1,ij,ik,ilmp2)

              do 4527 itt=1,211
                 tprobd(itt,ij,ik) = (tprob2(itt) - tprob1(itt))/5.
     >                          *(iday360 - idaypr(ilmp1)) + tprob1(itt)

                 tprob_pscd(itt,ij,ik) = tprobd(itt,ij,ik)
 4527         CONTINUE

c  TPROB_PSCD(211,45,76) is separate array for PSCs for coupled model (when using NCEP PDFs for PSCs)



cccoup  **************************************************************************************************
cccoup
cccoup - THIS USES COUPLED MODEL ZONAL MEAN TEMPS + NCEP T' (ORIGINAL METHOD - uses integer offset ONLY)
cccoup
cccoupC
cccoupC   for coupled model do ADJUSTED PDFs:
cccoupC         1) get zonal mean and perturbation temps from TPROB  -  TPROBD(211,45,76)
cccoupC         2) remove NCEP zonal mean,  add in coupled model zonal mean - TEMPCG(45,76)
cccoupC         3) apply prespecified MODEL - NCEP adjustment  -  TEMPOFFS(360,45,76), xlatoffs(45), zzoffs(76)
cccoupC         
cccoupC                                     also store current PDF in tprobst(211), then TPROBHYB(211,45,76)
cccoup       IF (ICOUP .EQ. 1) THEN
cccoup
cccoup             tzonal = 0.0
cccoup             do 6527 itt=1,211
cccoup                 tprobst(itt) = tprobd(itt,ij,ik)
cccoup 6527	         tzonal = tzonal + (itt + 119.)*tprobst(itt)
cccoup
cccoupC                                 determine zonal mean temp adjustment, but DON'T use TEMPOFFSET here
cccoupCC                                                                           
cccoupcctoffs                                      tadjzm = TEMPCG(ij,ik) - tzonal - TEMPOFFS(iday360,ij,ik)
cccoup             tadjzm = TEMPCG(ij,ik) - tzonal
cccoup
cccoup             iadj = 0
cccoup             if (tadjzm .lt. 0.) iadj = INT(tadjzm - .5)
cccoup             if (tadjzm .ge. 0.) iadj = INT(tadjzm + .5)
cccoup
cccoupC  adjust PDFs - move temperatures by amount IADJ, ensure it's 1-211, extreme temps shouldn't be a problem
cccoupC      find 1st, last temp in PDF  ; initialize tprobhyb
cccoup
cccoup             ic1 = 1
cccoup             ic2 = 211
cccoup
cccoup             do 6537 itt=1,211
cccoup                if (tprobst(itt) .ne. 0.) ic2 = itt
cccoup 6537        CONTINUE
cccoup
cccoup             do 6538 itt=211,1,-1
cccoup                if (tprobst(itt) .ne. 0.) ic1 = itt
cccoup 6538        CONTINUE
cccoup
cccoup
cccoup             do 6547 itt=1,211
cccoup 6547		tprobhyb(itt,ij,ik) = 0.0
cccoup
cccoupC  ensure that adjustment doesn't go out of bounds for coding purposes, but this should NOT be a problem
cccoupC                                   in case it does go out of bounds, PDF is NOT CHANGED (reload original)
cccoup             idiff1 = ic1 + iadj
cccoup             idiff2 = ic2 + iadj
cccoup
cccoup             if (idiff1 .ge. 1  .and.  idiff2 .le. 211) then
cccoup                do 6577 itt=ic1,ic2
cccoup                    inew = itt + iadj
cccoup                    tprobhyb(inew,ij,ik) = tprobst(itt)
cccoup 6577            CONTINUE
cccoup             else
cccoup                 do 6587 itt=1,211
cccoup                    tprobhyb(itt,ij,ik) = tprobst(itt)
cccoup 6587            CONTINUE
cccoup             endif
cccoupc                                        end coupled model temperature adjustment loop
cccoup       ENDIF
cccoup
cccoup  **************************************************************************************************


ccc
ccc - THIS USES COUPLED MODEL ZONAL MEAN TEMPS + NCEP T' - HYBRID PDF
ccc
ccc    WCONV260 - does integer temperature offset, then accounts for fractional degree offset
ccc
C    for coupled model do ADJUSTED PDFs:
C      1) get zonal mean and perturbation temps from TPROB  -  TPROBD(211,45,76)
C      2) remove NCEP zonal mean,  add in coupled model zonal mean - TEMPCG(45,76)
C      3) adjust PDF: 1) integer offset, then 2) fractional degree offset
C
C                                store current PDF in tprobst(211), then TPROBHYB(211,45,76)
       IF (ICOUP .EQ. 1) THEN

          tzonal = 0.0
          do 6527 itt=1,211
             tprobst(itt) = tprobd(itt,ij,ik)
 6527        tzonal = tzonal + (itt + 119.)*tprobst(itt)
C                                                   determine zonal mean temp adjustment
          tadjzm = TEMPCG(ij,ik) - tzonal


C  adjust PDFs - move temperatures by amount IADJ
C     then do fractional degree offsett, do fractional offset separately for +/- offset
C     ensure it's 1-211, extreme temps shouldn't be a problem
C       find 1st, last temp in PDF  ; initialize TPROBHYB(211,45,76), TPROBFR(211) 

          do 6547 itt=1,211
             tprobfr(itt) = 0.0
     	     tprobhyb(itt,ij,ik) = 0.0
ccccccc             tprobfrw(itt,ij,ik) = 0.0
 6547	  CONTINUE

             iadj = INT(tadjzm)
             ffo  = ABS(tadjzm - iadj)

C                                        ! tzonalw, iadjw, ffow, tadjzmw(45,76) are for output checking only
             tzonalw(ij,ik) = tzonal
             iadjw(ij,ik) = iadj
             ffow(ij,ik) = ffo
             tadjzmw(ij,ik) = tadjzm


             ic1 = 1
             ic2 = 211

             do 6537 itt=1,211
                if (tprobst(itt) .ne. 0.) ic2 = itt
 6537        CONTINUE

             do 6538 itt=211,1,-1
                if (tprobst(itt) .ne. 0.) ic1 = itt
 6538        CONTINUE


             ic1w(ij,ik) = ic1
             ic2w(ij,ik) = ic2


C  ensure that adjustment doesn't go out of bounds for coding purposes, but this should NOT be a problem
C    in case it does go out of bounds, just use COUPLED model Zonal mean temps
C    NOTE: idiff1, idiff2 have +/- 1 added here for the fractional PDF change (extra bin is used)
C
C  ic1w(45,76), ic2w(45,76), id1w(45,76), id2w(45,76), tprobfrw(211,45,76) are for output checking only

             idiff1 = ic1 + iadj - 1
             idiff2 = ic2 + iadj + 1

             if (idiff1 .ge. 1  .and.  idiff2 .le. 211) then

                do 6577 itt=ic1,ic2
                    inew = itt + iadj
                    tprobhyb(inew,ij,ik) = tprobst(itt)
                    tprobfr(inew) = tprobhyb(inew,ij,ik)
ccccccc                    tprobfrw(inew,ij,ik) = tprobfr(inew)
 6577           CONTINUE

C  fractional adjustment:  TPROBFR(211), get starting, ending indicies of NEW PDF
C
                id1 = 1
                id2 = 211

                do 6737 itt=1,211
                   if (tprobfr(itt) .ne. 0.) id2 = itt
 6737           CONTINUE

                do 6738 itt=211,1,-1
                   if (tprobfr(itt) .ne. 0.) id1 = itt
 6738           CONTINUE


                id1w(ij,ik) = id1
                id2w(ij,ik) = id2


                if (tadjzm .gt. 0.) then
                   tprobhyb(id1,ij,ik) = tprobfr(id1)*(1.-ffo)
                   tprobhyb(id2+1,ij,ik) = tprobfr(id2)*ffo
        
                   do 9577 itt = id1+1,id2
 9577   tprobhyb(itt,ij,ik) = tprobfr(itt-1)*ffo + tprobfr(itt)*(1.-ffo)
                endif


                if (tadjzm .lt. 0.) then
                   tprobhyb(id1-1,ij,ik) = tprobfr(id1)*ffo
                   tprobhyb(id2,ij,ik) = tprobfr(id2)*(1.-ffo)

                   do 9677 itt = id1,id2-1
 9677   tprobhyb(itt,ij,ik) = tprobfr(itt)*(1.-ffo) + tprobfr(itt+1)*ffo
                endif

             endif

C  just in case PDF goes out of bounds, load in coupled model zonal mean temps here

             if (idiff1 .lt. 1  .or.  idiff2 .gt. 211) then
                iit = NINT(TEMPCG(ij,ik)) - 119
                tprobhyb(iit,ij,ik) = 1.
             endif


C  reload TPROBD(211,45,76) with TPROBHYB(211,45,76) for use below

          do 7576 itt=1,211
 7576	     tprobd(itt,ij,ik) = tprobhyb(itt,ij,ik)

c                                                    end coupled model temperature adjustment loop
       ENDIF


ccc  **************************************************************************************************


C
C   If using coupled model PWAVE PERTURBATION longitudinal temps - tcc(72,L$,Z$) (T') - ICOUPDF = 1
C      load into TPROBC(211,45,76) is for 120-330K
C          TEMPCG(45,76) are the coupled model zonal mean temps            
C
C          add in coupled model T', re-generate PDFs
C
C
C   NOTE: load TPROBC below for use in gas and het reactions
C         TPROB_PSC will still be from NCEP PDFs for coupled model, if requested
C

       IF (ICOUP .EQ. 1) THEN   !   .and.  ICOUPDF .eq. 1) THEN
C                                                     initialize TPROBC(211,45,76)
             do 7775 itt=1,211
 7775		 tprobc(itt,ij,ik) = 0.0

C  first loop over all longitudes to check if Temps go <120K or >330K, if they do, just use ZONAL MEAN
C
             iich = 0
             do 7776 ix=1,72
                 ttc = TEMPCG(ij,ik) + tcc(ix,ij,ik)
                 iit = NINT(ttc) - 119
                 if (iit .lt. 1  .or.  iit .gt. 211) iich = 1
 7776        CONTINUE


             if (iich .eq. 0) then
                do 7777 ix=1,72
                   ttc = TEMPCG(ij,ik) + tcc(ix,ij,ik)
                   iit = NINT(ttc) - 119
                   tprobc(iit,ij,ik) = tprobc(iit,ij,ik) + 1./72.
 7777           CONTINUE
             endif

             if (iich .eq. 1) then
                ttc = TEMPCG(ij,ik)
                iit = NINT(ttc) - 119
                tprobc(iit,ij,ik) = 1.
             endif

C                                            if using Coupled model T' for PSCs, load here
ccoupsc             do 7774 itt=1,211
ccoupsc 7774		 tprob_pscd(itt,ij,ik) = tprobc(itt,ij,ik)

       ENDIF
C                      end coupled model T' PDF loop


              do 4521 itt=1,181
                 lprobd(itt,ij,ik) = 0.0
                 lprob1(itt) = 0.0
                 lprob2(itt) = 0.0
 4521         CONTINUE

            do 4544 itt=ilr1(ij,ik,ilmp1), ilr2(ij,ik,ilmp1)
 4544         LPROB1(itt) = LPROB45(itt-ilr1(ij,ik,ilmp1)+1,ij,ik,ilmp1)

            do 4545 itt=ilr1(ij,ik,ilmp2), ilr2(ij,ik,ilmp2)
 4545         LPROB2(itt) = LPROB45(itt-ilr1(ij,ik,ilmp2)+1,ij,ik,ilmp2)

            do 4547 itt=1,181
                 lprobd(itt,ij,ik) = (lprob2(itt) - lprob1(itt))/5.
     >                          *(iday360 - idaypr(ilmp1)) + lprob1(itt)
 4547         CONTINUE

 4520   CONTINUE


C
C   NOW convert to proper grid if necessary, using IJL(L$), IKZ(Z$) in COMMON already defined in INPUT
C                 TPROBD(211,45,76)     --> TPROB(211,L$,Z$)
C                 TPROB_PSCD(211,45,76) --> TPROB_PSC(211,L$,Z$)
C                 LPROBD(181,45,76)     --> LPROB(181,L$,Z$)
C
C  WCONV260 - use TPROBHYB(211,45,76) - loaded into TPROBD(211,45,76) 
C            - MODEL ZONAL MEAN TEMPS + NCEP T' (Hybrid)
C
      IF (L$ .eq. 45  .and.  Z$ .eq. 76) THEN 

           do 755 ik=1,Z$
           do 755 ij=1,L$
                 do 751 itt=1,211
    		    tprob(itt,ij,ik)  = tprobd(itt,ij,ik)
    		    tprob_psc(itt,ij,ik) = tprob_pscd(itt,ij,ik)
 751		 CONTINUE

                 do 752 itt=1,181
 752		    lprob(itt,ij,ik) = lprobd(itt,ij,ik)
 755    continue  


        ELSE


C  if model grid is NOT 45x76, find nearest model grid points to the 45x76 grid
C      using IJL(L$), IKZ(Z$) in COMMON already defined in INPUT;  
C
C      NOTE: the limit on ijl(L$) is 45,   the limit on ikz(Z$) is 76
C

          do 1991 ik=1,Z$
          do 1991 ij=1,L$
             do 112 itt=1,211
    		tprob(itt,ij,ik) = tprobd(itt,ijl(ij),ikz(ik))
    		tprob_psc(itt,ij,ik) = tprob_pscd(itt,ijl(ij),ikz(ik))
 112         CONTINUE


             do 113 itt=1,181
 113		lprob(itt,ij,ik) = lprobd(itt,ijl(ij),ikz(ik))
 1991     CONTINUE
c                                        end loop converting tprob/lprob to proper grid
      ENDIF


c
C     and define starting and ending temperatures in ITP1(L$,Z$), ITP2(L$,Z$), tprob(211,L$,Z$)
C                         and latitudes in  ILP1(L$,Z$), ILP2(L$,Z$),  lprob(181,L$,Z$)  - first initialize
C   
C   also do for ITPSC1(L$,Z$), ITPSC2(L$,Z$) in TPROB_PSC(211,L$,Z$) (in COMMON) for use in PSC calcs in NATICE
C                                    although itp1, itp2 are NOT needed in NATICE
                                   
         do 142 ik=1,Z$
         do 142 ij=1,L$
             itp1(ij,ik) = 1
             itp2(ij,ik) = 211

             itpsc1(ij,ik) = 1
             itpsc2(ij,ik) = 211

             ilp1(ij,ik) = 1
             ilp2(ij,ik) = 181
 142	 CONTINUE


         do 145 ik=1,Z$
         do 145 ij=1,L$
c                                  count forward to find last temperture in PDF
            do 146 itt=1,211
              if (tprob(itt,ij,ik) .ne. 0.) itp2(ij,ik) = itt
 146	    CONTINUE
c                                  count backward to find first temperture in PDF
            do 147 itt=211,1,-1
              if (tprob(itt,ij,ik) .ne. 0.) itp1(ij,ik) = itt
 147	    CONTINUE


c                                  count forward to find last temperture in PDF
            do 346 itt=1,211
              if (tprob_psc(itt,ij,ik) .ne. 0.) itpsc2(ij,ik) = itt
 346	    CONTINUE
c                                  count backward to find first temperture in PDF
            do 347 itt=211,1,-1
              if (tprob_psc(itt,ij,ik) .ne. 0.) itpsc1(ij,ik) = itt
 347	    CONTINUE


c                                  count forward to find last latitude in PDF
            do 246 itt=1,181
              if (lprob(itt,ij,ik) .ne. 0.) ilp2(ij,ik) = itt
 246	    CONTINUE
c                                  count backward to find first latitude in PDF
            do 247 itt=181,1,-1
              if (lprob(itt,ij,ik) .ne. 0.) ilp1(ij,ik) = itt
 247	    CONTINUE

 145	 CONTINUE

c                          end PDF loading loop       
      ENDIF



C  if using ZONAL MEAN TEMPERATURES for CHEMISTRY,  load in here:   TEMP(L$,Z$X)
C     NOTE: NATICE ALWAYS USES TPROB_PSC (NEVER Zonal means)       TPROB(211,L$,Z$), ITP1(L$,Z$), ITP2(L$,Z$)

ccccccc       IF (ITPDF .EQ. 0   .OR.   ICOUP .EQ. 1) THEN - insert this if using COUPLED MODEL zonal mean temps 
c                                                         this will overwrite above, but also need to set TPROB_PSC
      IF (ITPDF .EQ. 0) THEN 

         IZONAVGT = 1

         do 170 ik=1,Z$
         do 170 ij=1,L$
                                                    ! first initialize entire PDF to zero
           do 166 itt=1,211
 166           tprob(itt,ij,ik) = 0.
                                                    ! use first index as zonal mean temp w/ probability of 1.
                                                    !   and use TEMP(L$,Z$X) array wherever TPROB is used
           tprob(1,ij,ik) = 1.

           itp1(ij,ik) = 1
           itp2(ij,ik) = 1
 170	 CONTINUE

      ENDIF


C  if using ZONAL MEAN INSOLATION for CHEMISTRY,  load in here:   first initialize to zero;   LAT4(L$)
C                                                                 LPROB(181,L$,Z$), ILP1(L$,Z$), ILP2(L$,Z$)
C
      IF (ILPDF .EQ. 0) THEN 

         do 7176 ik=1,Z$
         do 7176 ij=1,L$
         do 7176 itt=1,181
 7176	   lprob(itt,ij,ik) = 0.

         do 370 ik=1,Z$
         do 370 ij=1,L$
           ilat0 = INT(LAT4(ij) + 91)

           lprob(ilat0,ij,ik) = 1.0
           ilp1(ij,ik) = ilat0
           ilp2(ij,ik) = ilat0
 370	 CONTINUE

      ENDIF


C   wconv260 - write out TPROBs from coupled model/NCEP hybrid PDFs to fort.457
C        tzonalw(45,76), iadjw(45,76), ffow(45,76), tprobfrw(211,45,76)
C        ic1w(45,76), ic2w(45,76), id1w(45,76), id2w(45,76)
c

ccccccc        if (mod(iday360-15, 30) .eq. 0.0) then
ccccccc          write (457) iday360, TEMPCG, tzonalw, iadjw, ffow, tadjzmw
ccccccc          write (457) tprob_psc, tprobc, tprobhyb, tprobfrw
ccccccc          write (457) ic1w, ic2w, id1w, id2w
ccccccc        endif

cccc            write (555) TEMPCG, tcc, tprobd, tprob


C    ENHANCED HNO3 washout in troposphere:
C 
C  load in HNO3 washout array HNO3W(L$,Z$), LAT4(L$), ZALT90(Z$) are in COMMON
C    use current daily lightning array - DLIGHT(L$,Z$) in COMMON, scaled by 5.e-10
C       which gives time scale of ~5 days for production rate of 5000. #/sec
C       larger scale washout of ~10 days in tropical troposphere per Ken P.
C       set minimum time scale to 30 days (4.e-7 sec^-1) everywhere
C                     NOTE, in RAIN:  if (ZALT90(ik) .gt. 16.) cn(59,ij,ik) = 0.d0 
      do 7756 ik=1,Z$
      do 7756 ij=1,L$

         HNO3W(ij,ik) = cn(59,ij,ik)

         rrlat = ABS(LAT4(ij))
         rrz   = ZALT90(ik)
         rrmin = DLIGHT(ij,ik)*5.e-10           ! 1.2e-6
         if (rrmin .le. 4.e-7) rrmin = 4.e-7

         if (rrlat .ge. 35. .and. rrz .le. 8.) then
            if (hno3w(ij,ik) .le. rrmin) hno3w(ij,ik) = rrmin
         endif

         if (rrlat .ge. 30. .and. rrlat .lt. 35. .and. rrz .le. 9.) then
            if (hno3w(ij,ik) .le. rrmin) hno3w(ij,ik) = rrmin
         endif

         if (rrlat .ge. 25. .and. rrlat .lt. 30. .and. rrz .le. 11.)then
            if (hno3w(ij,ik) .le. rrmin) hno3w(ij,ik) = rrmin
         endif

         if (rrlat .ge. 20. .and. rrlat .lt. 25. .and. rrz .le. 13.)then
            if (hno3w(ij,ik) .le. rrmin) hno3w(ij,ik) = rrmin
         endif

         if (rrlat .lt. 20. .and. rrz .le. 15.) then
            if (hno3w(ij,ik) .le. rrmin) hno3w(ij,ik) = rrmin
         endif
C                            ramp down to min 100 day time scale at 15.5 km in tropics
         if (rrlat .lt. 20. .and. rrz .le. 16.) then
            if (hno3w(ij,ik) .le. 1.2e-7) hno3w(ij,ik) = 1.2e-7
         endif

 7756   CONTINUE


      RETURN
      END
