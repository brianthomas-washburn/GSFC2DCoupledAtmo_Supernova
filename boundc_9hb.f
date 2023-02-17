        SUBROUTINE BOUNDC

        include "com2d.h"


C   COMMON for SPARC OH field in (#/cm3)
C
       COMMON/CCOH/ohm(14,L$,Z$), ohday(L$,Z$)



C  COMMON for GMI sfc deposition, emissions, interpolated to L$ latitude grid (from TEMPIN):
C    DEP0 in cm/sec;   EMIS0 in #/cm2/sec - need to divide by DELTAZ(L$,Z$X) (is in KM) in COMMON 
C
       COMMON/CDEP/tdep(14), DEP0(14,L$,7), EMIS0(14,L$,3), 
     >             DEPDAY(L$,7), EMISDAY(L$,3)


C   common for CGCM-COMBO T1R8 sfc deposition and wet scavenging, interpolated to L$,Z$ grid in TEMPIN
C     DEPR8 is in cm/sec;  - need to divide by DELTAZ(L$,Z$X) to get 1/sec
C
C     SCAVR8 is in 1/sec ;  SCAVDAY is the combined surface deposition + scavenging  in  1/sec
C
       COMMON/CDEPR8/DEPR8(10,14,L$), SCAVR8(20,14,L$,Z$),
     >               SCAVDAY(20,L$,Z$)



c  COMMON for GMI lightning, sfc constits: FORTRAN indicies are:
c
C    1=HNO2  2=HNO3   3=HNO4   4=NO    5=NO2   6=NO3   7=N2O5   8=CH2O   9=CO   10=H 
C   11=OH   12=HO2   13=H2O2   14=ClONO2  15=BrONO2   16=N   17=O3   18=DU/km   19=NOy
C
        COMMON/CGMI/ gmi2dbc(14,19,L$,10), timesb(14),
     >               gmilt(14,45,20), xlatlt(45), zzlt(20)


        REAL ttime(1), ttout(1), bb14(14), glday(45,20)
        REAL DEPDAYR8(10,L$), SCAVDAYR8(20,L$,Z$)


        DATA IGSPEC/42, 10, 38,  5, 6, 7, 8, 23, 19, 12, 13, 14, 16,
     >              30, 47,  9,  4, 80, 31/


	SAVE

C SET UP BOUNDARY CONDITIONS FOR TIME-DEPENDENT SPECIES
C   Note: IYR=0 corresponds to 1935, IYR=15 => 1950, IYR=65 => 2000, IYR=115 => 2050, IYR=165 => 2100 
C   new TD BCs from WMO-2006 now corresponds to the BEGINNING of the year when interpolating
C                               GHGs also corresponds to the 1st of the year 

        DO 200 II=1,165
        IUSE=II
        IF(IYR .GE. IYRBC(II) .AND. IYR .LT. IYRBC(II+1))GO TO 250
200     CONTINUE


250     DO 300 JJ=1,30

	IS=ISPBC(JJ)

        MIXRATBC = (BCTDINPUT(IS,IUSE+1)-BCTDINPUT(IS,IUSE))/360.
     >                   *(iday360-1.) + BCTDINPUT(IS,IUSE)


        if (iyr .le. 14)  MIXRATBC = BCTDINPUT(IS,16)
        if (iyr .ge. 165) MIXRATBC = BCTDINPUT(IS,166)


C   if steady-state run, overwrite MIXRATBC with proper year from TD TABLE
C     choose any year in the table (1950-2100); this is repeated at each time step 

        IF (LSTSTATE) MIXRATBC = BCTDINPUT(IS,iyearss-1934)



cef        DELTABC=(BCTDINPUT(IS,IUSE+1)-BCTDINPUT(IS,IUSE))/
cef     c          (1.*iyearbctd(iuse+1)-1.*iyearbctd(iuse))

c	if(iyr .gt. 90)then
c	if(is.eq.34)print *,' is=',is,' iuse=',iuse,' deltabc=',
c     *  deltabc
c	if(is.eq.34)print *,' bciuse+1=',bctdinput(is,iuse+1),
c     *  ' bciuse=',bctdinput(is,iuse)
c	endif

cef        TEMP1USE=IYR + (DAY360/360.)
cef        TEMP2USE=IYRBC(IUSE)

	IF(LBCMRTD(IS))THEN
C THIS SECTION IS FOR MIXING RATIO BOUNDARY CONDITIONS
c                                                           BVAL(S$,L$), TDLAT(S$,L$)
cef	MIXRATBC=(TEMP1USE-TEMP2USE)*DELTABC + 
cef     *  BCTDINPUT(IS,IUSE)
c	if(iyr .gt. 90)then
c	if(is.eq.34)print *,' mixrat=',mixratbc,' iyr=',iyr,
c     *  ' day360=',day360,' temp1use=',temp1use,' temp2use=',
c     *  temp2use,' iuse=',iuse
c	endif
	DO 500 IL=1,L$
	BVAL(IS,IL)=TDLAT(IS,IL)*MIXRATBC*1.E-9
            if (IS .eq. 20) BVAL(IS,IL) = TDLAT(IS,IL)*MIXRATBC*1.E-6             ! CO2 BC's are in PPMV
c	if(iyr.gt.90 .and. il.eq.9)then
c	if(is.eq.34)print *,' bval=',bval(is,il),' il=',il,
c     *  ' tdlat=',tdlat(is,il),' mixrat=',mixratbc
c	endif
500	CONTINUE

	ENDIF


	IF(.NOT.LBCMRTD(IS))THEN
C THIS SECTION IS FOR FLUX BOUNDARY CONDITIONS

	FLUXBC=(TEMP1USE-TEMP2USE)*DELTABC + 
     *  BCTDINPUT(IS,IUSE)

C NORMALIZE BY MOLECULAR WEIGHT OF MOLECULE
	FLUXBC=FLUXBC*(FRMOLWT(1)/FRMOLWT(JJ))

	DO 1500 IL=1,L$
C USE LATITUDE DEPENDENCE READ IN FOR CFC11 IN ALL CFCs
	BVAL(IS,IL)=TDLAT(ISPBC(1),IL)*FLUXBC
1500	CONTINUE

	ENDIF

300	CONTINUE




C  now set up boundary conditions for NOx species, CO, and CH2O from GMI simulations
C  - gmi2dbc(14,19,L$,10) has seasonal variations, but no long term TD 
C  - interpolate to current day:    bb14(14), timesb(14) 
C    BVALG(S$,L$,10) is in COMMON,  unit are all ppp
C    load in using igspec(19) (in COMMON) for 2D model species,  ik = 0.5 - 9.5 km
C
         ttime(1) = IDAY360*1.

         DO 588 ik=1,10
         DO 588 ij=1,L$
         DO 588 iis=1,19

            do 589 iii=1,14
 589            bb14(iii) = gmi2dbc(iii,iis,ij,ik)
 
            CALL LINTERP(timesb, 14, bb14, ttime, 1, 0, ttout)

            bvalg(igspec(iis),ij,ik) = ttout(1)
 588   CONTINUE


C   interpolate GMI lightning to current day:  gmilt(14,45,20) -> glday(45,20)
C         units are  N production in #/cm^3/sec
C                             interpolate to DLIGHT(L$,Z$) (in COMMON) - xlatlt(45), zzlt(20)
         DO 688 ik=1,20
         DO 688 ij=1,45

            do 689 iii=1,14
 689            bb14(iii) = gmilt(iii,ij,ik)
 
            CALL LINTERP(timesb, 14, bb14, ttime, 1, 0, ttout)

            glday(ij,ik) = ttout(1)
 688   CONTINUE


        CALL BINTERP(xlatlt, 45, zzlt, 20, GLDAY,
     >               lat4, L$, zalt90, Z$, 0, 0, DLIGHT)



C  interpolate GMI surface dep, emissions to current day: divide by DELTAZ(L$,Z$X) at sfc (in KM)
C      tdep(14), DEP0(14,L$,7) is in cm/sec      > DEPDAY(L$,7) is in 1/sec
C                EMIS0(14,L$,3) is in #/cm^2/sec > EMISDAY(L$,3) is in #/cm^3/sec
C
         DO 788 iis=1,7
         DO 788 ij=1,L$

            do 789 iii=1,14
 789            bb14(iii) = dep0(iii,ij,iis)
 
            CALL LINTERP(tdep, 14, bb14, ttime, 1, 0, ttout)

            depday(ij,iis) = ttout(1)/(DELTAZ(ij,1)*1.e5)
 788   CONTINUE


CCWCONV200 
C  DEPOSITIONS:  1=CH2O  2=HNO3  3=H2O2  4=CH3OOH (MP)  5=N2O5  6=NO2   7=O3
C  modify Ozone surface dep to get better trop Ozone: DEPDAY(IJC,7) is ozone; LAT4(L$) in COMMON
C    NH, ramp up to minimum of .12-.18 1/days

         DO 577 ij=1,L$
            dep2 = DEPDAY(ij,7)*86400.

            if (LAT4(ij) .ge.  0. .and. dep2 .le. .12) dep2 = .12
            if (LAT4(ij) .ge.  5. .and. dep2 .le. .13) dep2 = .13
            if (LAT4(ij) .ge. 10. .and. dep2 .le. .14) dep2 = .14
            if (LAT4(ij) .ge. 15. .and. dep2 .le. .15) dep2 = .15
            if (LAT4(ij) .ge. 20. .and. dep2 .le. .16) dep2 = .16
            if (LAT4(ij) .ge. 25. .and. dep2 .le. .17) dep2 = .17
            if (LAT4(ij) .ge. 30. .and. dep2 .le. .18) dep2 = .18

ccc            if (dep2 .ge. .15) dep2 = .15


ccc             if (LAT4(ij) .lt. 15. .and. dep2 .le. .09) dep2 = .09
ccc 
ccc             if (LAT4(ij) .lt. 10. .and. dep2 .ge. .08) dep2 = .08
ccc             if (LAT4(ij) .lt.  5. .and. dep2 .ge. .07) dep2 = .07
ccc             if (LAT4(ij) .lt.  0. .and. dep2 .ge. .06) dep2 = .06
ccc             if (LAT4(ij) .lt. -5. .and. dep2 .ge. .05) dep2 = .05

            DEPDAY(ij,7) = dep2/86400.

            DEPDAY(ij,1) = DEPDAY(ij,1)/5.   ! reduce CH2O dep to better match GMI
 577     CONTINUE



         DO 888 iis=1,3
         DO 888 ij=1,L$

            do 889 iii=1,14
 889            bb14(iii) = emis0(iii,ij,iis)
 
            CALL LINTERP(tdep, 14, bb14, ttime, 1, 0, ttout)

            emisday(ij,iis) = ttout(1)/(DELTAZ(ij,1)*1.e5)
 888   CONTINUE


C
CCWCONV250
C
C  interpolate CGCM/COMBO T1R8 surface depositions to current day: divide by DELTAZ(L$,Z$X) at sfc (in KM)
C      tdep(14), DEPR8(10,14,L$) is in cm/sec ->  DEPDAYR8(10,L$) is in 1/sec;   bb14(14)
C
C  DEPDAYR8:  1=brono2  2=ch2o   3=clo  4=clono2   5=h2o2   6=hbr   7=hcl   8=hocl   9=mp    10=o3
C
         DO 988 ij=1,L$
         DO 988 iis=1,10

            do 989 iii=1,14
 989            bb14(iii) = depr8(iis,iii,ij)
 
            CALL LINTERP(tdep, 14, bb14, ttime, 1, 0, ttout)

            depdayr8(iis,ij) = ttout(1)/(DELTAZ(ij,1)*1.e5)
 988   CONTINUE


C
C  interpolate CGCM/COMBO T1R8 wet scavenging to current day;  tdep(14), bb14(14)
C     SCAVR8(20,14,L$,Z$) -> SCAVDAYR8, SCAVDAY(20,L$,Z$) are in 1/sec, ensure 0. above 16 km
C
C  SCAVDAY:  1=ch2o  2=hno2     3=hno3   4=ho2no2  5=ho2   6=h2o2  7=mp    8=no2  9=no3  10=n2o5
C           11=o3   12=brono2  13=hbr   14=br     15=brcl 16=hobr 17=clo  18=clono2  19=hcl  20=hocl
C
         DO 1088 ik=1,Z$
         DO 1088 ij=1,L$
         DO 1088 iis=1,20

            do 1089 iii=1,14
 1089            bb14(iii) = scavr8(iis,iii,ij,ik)
 
            CALL LINTERP(tdep, 14, bb14, ttime, 1, 0, ttout)

            scavdayr8(iis,ij,ik) = ttout(1)
            if (zalt90(ik) .ge. 16.) scavdayr8(iis,ij,ik) = 0.
 1088   CONTINUE


C  combine surface deposition and wet scavenging into one array here: 
C           SCAVDAY(20,L$,Z$) in 1/sec;   DEPDAYR8(10,L$) is in 1/sec
C
C  DEPDAYR8:  1=brono2  2=ch2o   3=clo  4=clono2   5=h2o2   6=hbr   7=hcl   8=hocl   9=mp    10=o3
C
         DO 1988 ik=1,Z$
         DO 1988 ij=1,L$
         DO 1988 iis=1,20
 1988       scavday(iis,ij,ik) = scavdayr8(iis,ij,ik) 


C  9HA - need to increase O3 sfc deposition, otherwise winter surface ozone too large

         DO 1991 ij=1,L$
            scavday(1,ij,1) = scavdayr8(1,ij,1) + depdayr8(2,ij)
            scavday(6,ij,1) = scavdayr8(6,ij,1) + depdayr8(5,ij)
            scavday(7,ij,1) = scavdayr8(7,ij,1) + depdayr8(9,ij)
           scavday(11,ij,1) = scavdayr8(11,ij,1) + depdayr8(10,ij)*10.

           scavday(12,ij,1) = scavdayr8(12,ij,1) + depdayr8(1,ij)
           scavday(13,ij,1) = scavdayr8(13,ij,1) + depdayr8(6,ij)
           scavday(17,ij,1) = scavdayr8(17,ij,1) + depdayr8(3,ij)
           scavday(18,ij,1) = scavdayr8(18,ij,1) + depdayr8(4,ij)
           scavday(19,ij,1) = scavdayr8(19,ij,1) + depdayr8(7,ij)
           scavday(20,ij,1) = scavdayr8(20,ij,1) + depdayr8(8,ij)
 1991    CONTINUE

C
C  SCAVDAY:  1=ch2o  2=hno2     3=hno3   4=ho2no2  5=ho2   6=h2o2  7=mp    8=no2  9=no3  10=n2o5
C           11=o3   12=brono2  13=hbr   14=br     15=brcl 16=hobr 17=clo  18=clono2  19=hcl  20=hocl
C
C  9HA - limit Bry washout to get better agreement w/ COMBO (HCl 1 day limit not needed)
C              BrONO2 - 5 days;  HOBr - 5 days
C
C   rbry - for testing, DONT limit BrONO2, HOBr washout

ccrbry         DO 1998 ik=1,Z$
ccrbry         DO 1998 ij=1,L$

ccccccc           sclim = 1./(1.*86400.)
ccccccc           if (scavday(19,ij,ik) .ge. sclim) scavday(19,ij,ik) = sclim

ccrbry            sclim = 1./(5.*86400.)
ccrbry            if (scavday(12,ij,ik) .ge. sclim) scavday(12,ij,ik) = sclim

ccrbry            sclim = 1./(5.*86400.)
ccrbry            if (scavday(16,ij,ik) .ge. sclim) scavday(16,ij,ik) = sclim
ccrbry 1998    CONTINUE


C  interpolate SPARC OH field to current day, tdep(14), ttime(1), ttout(1)
C     ohm(14,L$,Z$) -> ohday(L$,Z$) (in COMMON) in #/cm3,  ensure OHDAY > 0.
C
         DO 791 ik=1,Z$
         DO 791 ij=1,L$

            do 795 iii=1,14
 795            bb14(iii) = ohm(iii,ij,ik)
 
            CALL LINTERP(tdep, 14, bb14, ttime, 1, 0, ttout)

            ohday(ij,ik) = ttout(1)
            if (ohday(ij,ik) .lt. .001) ohday(ij,ik) = .001
 791   CONTINUE


C
C  modify Ozone surface dep to get better trop Ozone: DEPDAYR2(10,ij) is ozone; LAT4(L$) in COMMON
C  NH, ramp up to minimum of .12-.18 1/days, then increase by 50% everywhere (WCONV202)
ccc
ccc         DO 677 ij=1,L$
ccc            dep3 = DEPR2DAY(10,ij)*86400.
ccc
ccc            if (LAT4(ij) .ge.  0. .and. dep3 .le. .12) dep3 = .12
ccc            if (LAT4(ij) .ge.  5. .and. dep3 .le. .13) dep3 = .13
ccc            if (LAT4(ij) .ge. 10. .and. dep3 .le. .14) dep3 = .14
ccc            if (LAT4(ij) .ge. 15. .and. dep3 .le. .15) dep3 = .15
ccc            if (LAT4(ij) .ge. 20. .and. dep3 .le. .16) dep3 = .16
ccc            if (LAT4(ij) .ge. 25. .and. dep3 .le. .17) dep3 = .17
ccc            if (LAT4(ij) .ge. 30. .and. dep3 .le. .18) dep3 = .18
ccc
ccc            DEPR2DAY(10,ij) = dep3/86400.*1.5
ccc
ccc            DEPR2DAY(2,ij) = DEPR2DAY(2,ij)/5.   ! reduce CH2O dep to better match GMI
ccc 677     CONTINUE

	RETURN
	END
