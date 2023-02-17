
        SUBROUTINE DAILYINT


C    THIS routine interpolates various arrays to the current day's value
C         for things that used by both the fixed and coupled models, so it's 
C         the first routine that is called each time step (before the coupled model)
C         eg, aerosols, solar cycle, etc. - Feb. 2012
C         (as opposed to XCOUPLED which is just for the coupled model inputs)


        include "com2d.h"


	common/aeros/aerosol_all(2,L$,Z$)
	common/aer/aer_dbc(264,L$,Z$)


C   common for CCMVal SAD data set, 1950-2010, plus background

        COMMON/CSAD/sadbk7(L$,Z$,14), sad7(L$,Z$,492), times7(492)


c  solar flux data from Judith Lean, ALL in REAL*8 ; SJFLUXA = FLXA are the 1954-2008 averages
C  index 1 = time;  2-40 are the 39 model wavelength bins (in ph/cm^2/sec); 41 = TOA flux (W/m2)
C
C  HJFLUX(12,1524) are for the heating rate intervals, HJFLUXA(12) are 1954-2008 avgs (in mW/m2)
C   SJFLUXDAY(41), HJFLUXDAY(12) are values interpolated to the current day or long term avgs

        COMMON/CSJFLUX/ SJFLUX(41,1524), SJFLUXA(41),
     >                  HJFLUX(12,1524), HJFLUXA(12)

        COMMON/CSJFLUXD/ SJFLUXDAY(41), HJFLUXDAY(12), FLXA(41), ISOLCYC


        REAL*8 sjflux, sjfluxa, hjflux, hjfluxa, sjfluxday, hjfluxday
        REAL*8 ttimes(1), tsol(1524), gsol(1524), ggout(1), flxa

        REAL ttime(1), ttout(1), g4in(492), g4ss(14)


	SAVE


C
C  **********************   BEGIN   AEROSOL   LOADING   *************************************
C
C
c set up aerosol array  -  divide by 4 as per M. Prather for lower limit to aerosol affect
c  select appropriate distribution based on MONTH,  
C  AEROSOL(L$,Z$) in COMMON,  AEROSOL_ALL(2,L$,Z$) above

	if (iday360 .le. 180) then
		do 7127 ik=1,Z$
		do 7127 ij=1,L$
			aerosol(ij,ik)=aerosol_all(1,ij,ik)*1.e-8/4.
7127		continue
	else
		do 7128 ik=1,Z$
		do 7128 ij=1,L$
			aerosol(ij,ik)=aerosol_all(2,ij,ik)*1.e-8/4.
7128		continue
	endif

C note: the ik limit here for loops 7127 and 7128 was previously 30, ie, 58 km
C   but we have now set the AEROSOL_ALL array to 0.0 above 58 km (in the AEROSOLS routine)
C   so now just loop through all the way to Z$
C
C   but if we define an ik58 index for the PSC routine then we should put it in 
C   COMMON and use it here, or maybe make it lower in the atmosphere?? 
C   so things won't look strange in upper strat.?????
C
C
C   plug in SADs for 1979-2000 using David's new SAD data set, use iimon, aer_dbc(264,L$,Z$)
C         IYR=44 => 1979;    IYR=60 => 1995,   IYR=65 => 2000
C
C    for all years prior to 1979, use 1979 values
C    for 1979-1999 use supplied values      
C    for future years, use 1997 values for 2000-2050
C
C
        if (iyr .le. 43) iimon = int((iday360 + 29)/30)

        if (iyr .ge. 44  .and.  iyr .le. 64) 
     >            iimon = int((iyr - 44)*12 + (iday360 + 29)/30)

        if (iyr .ge. 65) iimon = int((iday360 + 29)/30) + 216

C
C   for Background (Clean) SAD, use 1979 values for ALL YEARS here:  (iimon = 1->12)
C       set ISADBK = 1 (in COMMON)
C

        ISADBK = 1

        if (ISADBK .eq. 1) iimon = INT((iday360 + 29)/30)


        do 12 ik=1,Z$
        do 12 ij=1,L$
              aerosol(ij,ik) = aer_dbc(iimon,ij,ik)
 12     continue


C
C  9HB:
C     IF ISADCCMV = 1, use CCMVal SAD data set for 1950-2010
C     interpolate to current day, OVERWRITE AEROSOL(L$,Z$) (in COMMON)

C     SAD7(492) is for 1963-2003 - use this for time dependent runs, and 
C     for year 2000 (ONLY) steady state (set ISADBK = 0)
C     OTHERWISE use background values (set ISADBK = 1)
C
C     also define year for SAD data set (YRSAD)
C
C   sadbk7(L$,Z$,14), sad7(L$,Z$,492), times7(492), timefz(14) are all in COMMON
C     load into AEROSOL(L$,Z$) in COMMON,   g4ss(14), g4in(492), times7(492)

        ISADCCMV = 1

        ISADBK = 1
        YRSAD = (IDAY360*1. - .5)/360. + iyr + 1935.

       if (.NOT. LSTSTATE .and. YRSAD .ge. 1963. .and. YRSAD .lt. 2004.) 
     >        ISADBK = 0

                                                !   use 2000 SAD for 2000 and 2100 TSlice

        if (LSTSTATE .and. IYEARSS .eq. 2000 .or. IYEARSS .eq. 2100)then
            ISADBK = 0
            YRSAD = 2000. + (IDAY360*1. - .5)/360.
        endif


        IF (isadccmv .eq. 1) THEN

          IF (ISADBK .eq. 1) then
              ttime(1) = IDAY360*1.

              DO 588 ik=1,Z$
              DO 588 ij=1,L$

                 do 589 iii=1,14
 589                g4ss(iii) = sadbk7(ij,ik,iii)
 
                 CALL LINTERP(timefz, 14, g4ss, ttime, 1, 0, ttout)

              aerosol(ij,ik) = ttout(1)
 588   CONTINUE

          ELSE

              ttime(1) = YRSAD

              DO 788 ik=1,Z$
              DO 788 ij=1,L$

                 do 789 iii=1,492
 789                g4in(iii) = sad7(ij,ik,iii)
 
                 CALL LINTERP(times7, 492, g4in, ttime, 1, 0, ttout)

              aerosol(ij,ik) = ttout(1)
 788   CONTINUE

          ENDIF

        ENDIF

C
C    **********************    END   AEROSOL   LOADING   **********************
C
C


C
C    **********************   BEGIN  SOLAR  CYCLE  LOADING   **********************
C
C  9HB - use Judith Lean's solar flux data:
C
C  IF ISOLCYC = 1, use time dependent data for 1882-2008, interpolate to current day
C  LINTERP8 works fine - OK to interpolate before/after first/last time index keeps constant value
c     ttimes(1), ggout(1), gsol(1524), tsol(1524)
C
C  index 1 = time (1882.04-2008.96)
C    2-40 are the 39 model wavelength bins (in ph/cm^2/sec);  41 = TOA flux (W/m2)
C
C    SJFLUX(41,1524) => SJFLUXDAY(41) ;   HJFLUX(12,1524) => HJFLUXDAY(12)
C
C
C  If NOT DOING SOLAR CYCLE (default, ISOLCYC = 0), eg, steady state runs, etc.,
C     then use SJFLUXA(41) and HJFLUXA(12) (1954-2008 avgs), load into daily arrays
C
C   except for 2000 time slice, use year 2000 solar flux conditions repeated


        do 5591 ij=1,41
           sjfluxday(ij) = sjfluxa(ij)
 5591      flxa(ij) = sjfluxa(ij)

        do 7791 ij=1,12
 7791      hjfluxday(ij) = hjfluxa(ij)


        ISOLCYC = 0

        ttimes(1) = (IDAY360*1. - .5)/360. + iyr + 1935.

        IF (LSTSTATE  .and.  IYEARSS .eq. 2000) then
           ISOLCYC = 1
           ttimes(1) = (IDAY360*1. - .5)/360. + 2000.
        ENDIF


        IF (ISOLCYC .eq. 1) then

             do 4491 iii=1,1524
 4491           tsol(iii) = sjflux(1,iii)

             do 4488 ij=2,41

                do 4489 iii=1,1524
 4489               gsol(iii) = sjflux(ij,iii)
 
                CALL LINTERP8(tsol, 1524, gsol, ttimes, 1, 0, ggout)

             sjfluxday(ij) = ggout(1)
 4488     CONTINUE

C  
C  and interpolate flux data for the heating rate grid to current day:   
C     HJFLUX(12,1524) => HJFLUXDAY(12) (mW/m2)
C  index 1=time (identical to above);  2=TSI (W/m2) (identical to SJFLUX(41) above); 
C  indicies 3-11 = are the heating rate grid;  12 is blank;  
C   use tsol(1524) and ttimes(1) from above
C
             do 5488 ij=2,11

             do 5489 iii=1,1524
 5489            gsol(iii) = hjflux(ij,iii)
 
                 CALL LINTERP8(tsol, 1524, gsol, ttimes, 1, 0, ggout)

              hjfluxday(ij) = ggout(1)
 5488     CONTINUE

        ENDIF


C    **********************   END  SOLAR  CYCLE  LOADING   **********************
C

        RETURN
        END
