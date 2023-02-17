C
C  NEW 2D MODEL Routine for computing  RAINOUT, sets trop values to Oort, 83  --  EF  4/18/94
C               NOW set to ERA-40/HALOE clim. -- EF, 11/2008          
C
C   NO MORE RELATIVE HUMIDITY FIX at TROPICAL TROPOPAUSE NEEDED WITH PPM
C   Also, Oort values at mid-high latitudes seem too wet much of the year compared to saturation
C   vapor pressure over ice based on the NMC temps, so rainout also at levels 3-4, 
C
C    NOTE: Oort values for NOW USED FOR ONLY LEVELS 1-2, although this is very similar to setting 
C          levels 1-4, since its determined by the RH Rainout limit......
C
C
        SUBROUTINE RAINOUT


        include "com2d.h"

C
C   ITRAIN IS THE INDEX OF THE TOP MODEL LEVEL TO COMPUTE RAINOUT FOR EACH LATITUDE, -- NO (OLD)
C
C   For BASE8SB, now have A RAINOUT PROFILE FOR the TROPICS (15S-15N - needs to be somewhat wet),
C   and outside the tropics (dry), and interpolate at 25-35N,S,
C   so relative humidity limit (RHRAIN) is a function of latitude and height
C   ****  NOW specified based on the daily climalogical trop hgts for up to 9 levels in the 
C         tropics, lower levels at mid-high latitudes, uses ITROP360(L$,360) in COMMON *******
C
C   Note: need to make stratosphere/mesosphere wetter with new PPM transport, so raise RH level, 
C   keep constant w/ height at 45% (60% too wet), note levels 1-2 are not used
C

         DIMENSION ITRAIN(L$), RHRAIN(L$,9), RHTR(9), RHML(9)

         REAL*4 ZZHH(45,92), ZZHHM(L$,Z$), H2OFIX(L$,Z$), co2sens(L$,Z$)
         REAL*4 RAINSTR(L$S,Z$S), RAINCHEM(L$,Z$), g4hday(L$,Z$)
         REAL*4 g4h2o1(91,76), ttime(1), ttout(1), g4ss(14)

         REAL*8 PRECIP, XSAT, VPS, WVLIM, krain, rhum, tvap


C   latent heating COMMON block from TEMPIN:   WLH, WHACK are from WACCM, in K/sec 
C
       COMMON/CCLH/WLH(91,117,360)
ccwconv175 - OLD       COMMON/CCLH/LHC(91,117,360), WLH(91,117,360), WHACK(91,117,360)


c  common for TD GEOS 4 tropospheric H2O 1950-2100 in ppmv, includes sensitivity to CO2, w/ seas cycle
                                                                              from TEMPIN
       COMMON/CG4H2O/ g4h2o(16,91,76), latg4h(91), zh76(76), timesfh(14)


C
CC         DATA ITRAIN/7,7,7,7, 8,8, 9,9,9, 9,9,9, 8,8, 7,7,7,7/
CC  OLD       DATA ITRAIN/7,7,8,9,10,11,11,11,11,11,11,11,11,10,9,8,7,7/

         DATA RHTR/.30, .30, .30, .30, .30, .30, .30, .30, .30/
         DATA RHML/.20, .20, .20, .20, .20, .20, .20, .20, .20/

         SAVE


C  first interpolate daily ERA-40/HALOE High res H2O climatology - used for tropical tropopause
C   HALH2OR(45,92,360), LATHAL(45), ZZHAL(92) are in COMMON to PROPER GRID FOR THE CURRENT DAY 
C   do LINEAR/LOG INTERPOLATION in Latitude/pressure-alt
C   zzhh(45,92), lat4(L$), zalt90(Z$), zzhhm(l$,z$)  -->  HALH2O(L$,Z$) (ppmv) is in COMMON
C

          do 475 ik=1,92 
          do 475 ij=1,45
 475	       zzhh(ij,ik) = HALH2OR(ij,ik,iday360)

            CALL BINTERP(LATHAL, 45, ZZHAL, 92, ZZHH, 
     >                   LAT4, L$, ZALT90, Z$, 0, 1, ZZHHM)

            DO 476 ik=1,Z$
            DO 476 ij=1,L$
 476	      HALH2O(ij,ik) = zzhhm(ij,ik)


C  ****************************************************************
C  ****************************************************************
C
C  similarly, interpolate GEOS 4 H2O to current day, then interpolate to PROPER GRID as above
C    g4h2o(16,91,76), 1-14 are months, 15-16 are the sensitivity to CO2 in ppmv H2O/ppmv CO2
C                                      15 is the GEOS 4 sensitivity scaled to the ERA-40 values
C                                      16 is the GEOS 4 sensitivity UNSCALED
C    latg4h(91), zh76(76), timesfh(14)   --> g4h2o1(91,76) --> g4hday(L$,Z$)
C
C    use CO2 BC for current day - CN(20,ij,ik), CONVERT TO PPMV
C       use level 1 for any latitude (use 9 so it won't blow up w low res)
C       and get change from 1950-2100 avg (468.35 ppmv)
C
C
C   NOTE:  GEOS 4 tropospheric H2O is WAY too WET - DON'T USE
ccg4h
ccg4h         ddco2 = CN(20,9,1)/M(9,1)*1.e6 - 468.35
ccg4h         ddco2 = 0.0    ! try first w/o long term trend
ccg4h
ccg4h         ttime(1) = IDAY360*1.
ccg4h
ccg4h         DO 588 ik=1,76
ccg4h         DO 588 ij=1,91
ccg4h
ccg4h            do 589 iii=1,14
ccg4h 589            g4ss(iii) = g4h2o(iii,ij,ik)
ccg4h 
ccg4h            CALL LINTERP(timesfh, 14, g4ss, ttime, 1, 0, ttout)
ccg4h
ccg4h            g4h2o1(ij,ik) = ttout(1) + ddco2*g4h2o(15,ij,ik)
ccg4h 588   CONTINUE
ccg4h
ccg4h          CALL BINTERP(LATG4H, 91, ZH76, 76, g4h2o1,
ccg4h     >                 LAT4, L$, ZALT90, Z$, 0, 0, G4HDAY)
ccg4h  ****************************************************************





C
C  ****************************************************************
C
C   USE ERA40/HALOE seasonal cycle H2O + CO2 sensitivity - use average of SCALED/UNSCALED 
C    g4h2o(16,91,76), 1-14 are months, 15-16 are the sensitivity to CO2 in ppmv H2O/ppmv CO2
C                                      15 is the GEOS 4 sensitivity scaled to the ERA-40 values
C                                      16 is the GEOS 4 sensitivity UNSCALED
C
C    interpolate CO2 sensitivity factors to proper grid (do LOG INTERP in ALTITUDE for H2O sens)
C    latg4h(91), zh76(76),  g4h2o1(91,76) --> CO2SENS(L$,Z$)
C
C    to mimic CCM lower strat H2O trends (~.5 ppmv for 2000-2100):
C    set CO2 sensitivity factor = .0015 ppm/ppm for the strat ENTRY VALUE
C       ie, the CCM trend just above tropical TTL (~17-18 km)
C

         DO 788 ik=1,76
         DO 788 ij=1,91
 788        g4h2o1(ij,ik) = g4h2o(15,ij,ik)        ! + g4h2o(16,ij,ik))/2.


          CALL BINTERP(LATG4H, 91, ZH76, 76, g4h2o1,
     >                 LAT4, L$, ZALT90, Z$, 0, 1, CO2SENS)


         DO 789 ik=1,Z$
         DO 789 ij=1,L$
            if (ZALT90(ik) .ge. 16.8  .and.  ZALT90(ik) .le. 18.) then
               if (ABS(LAT4(ij)) .lt. 30.) co2sens(ij,ik) = .0015
            endif
 789     CONTINUE


C  Use CO2 BC for current day - CN(20,ij,ik) -> load into G4HDAY(L$,Z$) is in PPMV
C     use level 1 for any latitude (use 9 so it won't blow up w low res)
C     ERA-40 data is 1981-2001 avg, so get change from 1991 (353.47 ppmv)
C
         ddco2 = CN(20,9,1)/M(9,1)*1.e6 - 353.47

         DO 791 ik=1,Z$
         DO 791 ij=1,L$
 791        g4hday(ij,ik) = ddco2*co2sens(ij,ik)


C  ****************************************************************



         DO 1001 ik=1,9
            DO 1005 ij=1,5
 1005          RHRAIN(ij,ik) = RHML(ik)

            DO 1010 ij=6,7
 1010          RHRAIN(ij,ik) = (RHTR(ik)-RHML(ik))/3.*(ij-5) + RHML(ik)

            DO 1015 ij=8,11
 1015          RHRAIN(ij,ik) = RHTR(ik)

            DO 1020 ij=12,13
 1020          RHRAIN(ij,ik) = (RHML(ik)-RHTR(ik))/3.*(ij-14) + RHML(ik)
   
            DO 1025 ij=14,18
 1025          RHRAIN(ij,ik) = RHML(ik)
 1001  CONTINUE

C
C   Area of the globe is in  CM-2  (in COMMON)
C
cccccccccccc        GLOBE = 5.1E18 
        IDAY = IDAY360
	DRPLFT = DT       
C

C  Now compute rainout, SET INDEX OF LOWEST MODEL LEVEL TO COMPUTE RAINOUT  EG, 1 = 879 mb 
C   and rainout up to the NMC tropopause hgt (ITROP360) + 1 since with only a 2km grid spacing
C   we'll miss some of the troposphere, especially in the tropics, 
C   NOTE:  levels 1 and 2 are reset to OORT climatology below, set IBW=1 here for diagnostics

        IBW = 1

        DO 444 IJ = 1,L$
        DO 444 IK = IBW, Z$
                                                      !! only do this loop for below 34 km
             if (ZALT90(ik) .ge. 34.) GO TO 444
ccccccc             if (ZALT90(ik) .ge. 10.) GO TO 444


ccccc        DO 444 IK = IBW, ITROP360(IJ,IDAY360)+1
CC        DO 444 IK = IBW,ITRAIN(IJ)


C   loop through temperature PDF here to get ZONAL MEAN XSAT, then apply RH to get WVLIM
C       use ITP1(L$,Z$), ITP2(L$,Z$), TPROB(211,L$,Z$) as defined for current day in SETDAILY

           XSAT = 0.d0
           
           do 4000 itemp = itp1(ij,ik)+119, itp2(ij,ik)+119
               tvap = DBLE(itemp)                                ! get real*8 value
                                                        !if using zonal mean temps, use TEMP(L$,Z$X) - REAL*8
               if (IZONAVGT .eq. 1) tvap = TEMP(ij,ik)

               VPS = 1013./760.*exp(-88.44375824 + 0.7641223392*tvap -
     >             0.002285011125*(tvap**2) + 2.539292763e-6*(tvap**3))

c
c   calculate LOCAL saturation vapor pressure over ice from analytic form since temps in region
C      of rainout are 180-239K, and convert from Torr (mm of Hg) to mb => 1013./760.
C                                                 TEMP(L$,Z$X), PRESS(Z$X), ZALT90(Z$)
C
C   PRECIP IS THE NUMBER OF H2O MOLECULES WHICH ARE MORE THAN the Rel. hum. limit OF THE
C      H2O MOLECULES IN SATURATED AIR - then subtract from the total density,
C
C      *** NOTE: there should be NO .622 factor here (ratio of molecular weight of water to dry air), 
C          this is now corrected in base8zo, and the relative humidity limit in the 
C          tropics should be .622*.48= ~.3  to give the same mixing ratios in the strat/meso. ************
C     -- need to adjust RHRAIN accordingly - (Jan. 2000)  

             XSAT = XSAT + VPS/PRESS(IK)*TPROB(itemp-119,ij,ik)
 4000     CONTINUE


cef        WVLIM = RHRAIN(ij,ik)*XSAT*M(IJ,IK)

c  to get annual mean tropical tropopause mixing ratio of 3.84 ppmv, need RH=41.41% based on 1992 UKMO temps
C  and RH=44.77% w/ 1992-2000 UKMO temps, so for now just set RH=42% everywhere, up from previous 30% value
C     SET TO RH=35% W/ HALOE RESET BELOW - 11/04
C
C
C  New EXTRATROPICAL RAINOUT  -   (EF, AUg. 2005):
C
C     for tropics, do rainout as usual, outside of tropics, also do rainout of 10% RH
C     based on RH computed from UARS water vapor climatology and met temperatures 
C     but with this convective process, it only cares about the lowest temperature 
C     in the profile to determine the amount of water entering the stratosphere. 
C     Don't use the tropopause climatology since this may not be the lowest
C     temperature (eg, could be a little colder above, especially in polar night). So just do
C     rainout at 10% everywhere below 34 km - this will ensure getting the minimum temperature, 
C     but then also need to "rainout" in polar vortex, but need to ramp up to 100% RH above ~21 km
C     since this is NOT a convective process, but rather the cold temperatures are caused by diabatic 
C     cooling, and PSCs form. 
C
C     also for latitudes > 40, put excess water vapor into solid H2O above 13 km
C
C
C     Problem with using analyzed temps for rainout: SAT Mixing rat is non-linear in temperature, and
C     will be very sensitive to errors, also rising motion in over-shooting convective tops, where
C     the parcel is now COLDER than the environment will give lower mixing ratio at the tropopause
C     compared to analyzed temps, which are averages of many profiles at any one grid point (warmer 
C     upper trop profiles will preferentially occur in regions of sinking, so I can see where the 
C     water vapor mix ratio will be OVER-estimated using analyzed temperatures
C     e.g., at 227K, sat mix ratio = 280 ppmv (~10 km), while satur. at 212K = 42 ppmv, so that gives
C     a variation of 15% RH vs. 100% RH (42/280 = 15%) over a 15K temperature range 
C     ie, a bias/error of 15K at any one grid point corresponds to a error of 85% in RH.
C
C
          rhum = 1.0d0       !  0.35d0
                                                       ! set rhum = 20% for troposphere for washout rate - NO
ccwash          if (ZALT90(ik) .le. 9.) rhum = 0.20d0


          WVLIM = rhum*XSAT*M(IJ,IK)

          PRECIP = CN(15,IJ,IK) - WVLIM
          PRECIP = DMAX1(PRECIP,0.0)


CC  for 9FN run, don't remove any rain from CN15 since it gets reset to CLIM below
CC     , just compute PRECIP here for WASHOUT rate below
        
cccccccccc          CN(15,IJ,IK) = CN(15,IJ,IK) - PRECIP


cccc          if (ABS(LAT4(ij)) .ge. 40.  .and.  ZALT90(ik) .ge. 10.)
cccc     >           CN(67,IJ,IK) = CN(67,IJ,IK) + PRECIP


C
C   Also, load the precip into two constituents: CN(58)  (or CN(42) for the old GSFC chemistry)
C       and also CN(59), which loss frequency for rainout of certain species (the number density 
C       multiplied by some 'reaction rate' so that the wash-out time scale is ~1-3 days in 
C       the lower troposphere and decreases above. Also, we specify CN(59) at levels 1-3 at 
C       low latitudes (specify as a mixing ratio in parts per part: .0007, .0003, .0001) 
C       where PRECIP=0 since the RH limit we are using (~48%) is not low enough. 
C       We can fix this up eventually if and when we work on the troposphere more carefully.
C       e.g., use cloud/latent heat/precip climatologies, etc....
C       also, set rainout loss frequency = 0 above tropopause
C
C  for 9GL, washout rate based on WACCM latent heating, set below
C

cctfam  DON'T SET CN(58) here - it's now used for Cly/Bry family transport
cctfam
cctfam       CN(58,IJ,IK) = PRECIP
cctfam       C(58,IJ,IK) = CN(58,IJ,IK)


       krain = 5.D-22
ccwash       CN(59,IJ,IK) = PRECIP*krain
ccwash           if (ik .eq. 1  .and.  PRECIP .eq. 0.0) 
ccwash     >                             CN(59,ij,ik)=.0007*m(ij,ik)*krain
ccwash
ccwash           if (ik .eq. 2  .and.  PRECIP .eq. 0.0) 
ccwash     >                             CN(59,ij,ik)=.0003*m(ij,ik)*krain
ccwash
ccwash           if (ik .eq. 3  .and.  PRECIP .eq. 0.0) 
ccwash     >                             CN(59,ij,ik)=.0001*m(ij,ik)*krain
ccwash
ccwash       if (ik .gt. ITROP360(IJ,IDAY360)+1) CN(59,IJ,IK) = 0.0


       RAIN=PRECIP/DRPLFT

ccrain       PRECIPYR(IJ,IK,IDAY) = PRECIP/1.E20

C  Load in rainout loss in Mixing Ratio per second for output -  RAIN1(L$,Z$) is in COMMON

        RAIN1(ij,ik) = -PRECIP/M(IJ,IK)/DT

444     CONTINUE

C
C   RESET TROPOSPHERIC VALUES (below 3 km) TO THOSE OF OORT, 1983, from DAILY input values
C          TROPWV(L$,Z$X,360) has been interpolated in TEMPIN from original 18x11 grid , ZALT8(Z$X)
C
C   NOTE: this is now set using the HALOE clim. below (which has ERA-40 in troposphere)
C                                   don't use TROPWV and precipyr to save memory
Crain        DO 8805 IJ=1,L$
Crain            DO 8805 IK=1,2
cccccc            DO 8805 IK=1,ITROP360(IJ,IDAY360)-1
Crain 8805       CN(15,IJ,IK) = TROPWV(IJ,IK,IDAY)*M(IJ,IK)/1.D06
c

C
C  Reset tropospheric H2O to ERA40/HALOE high resolution DAILY climatology up to 18 km in tropics
C    HALH2O(L$,Z$) in COMMON (ppmv) - LAT4(L$), ZALT90(Z$)
C    include contribution to CO2 sensitivity - G4HDAY(L$,Z$) is ppmv, load into H2OFIX(L$,Z$)
C
       do 8891 ik=1,Z$
       do 8891 ij=1,L$
          h2ofix(ij,ik) = (HALH2O(ij,ik) + G4HDAY(ij,ik))*1.d-6*m(ij,ik)

          if (ABS(LAT4(ij)) .lt. 25.) then
              if (ZALT90(ik) .le. 18.)  cn(15,ij,ik) = H2OFIX(ij,ik)
          endif

          if (ABS(LAT4(ij)) .ge. 25.  .and.  ABS(LAT4(ij)) .le. 30.)then
             if (ZALT90(ik) .le. 17.)  cn(15,ij,ik) = H2OFIX(ij,ik)
          endif

          if (ABS(LAT4(ij)) .gt. 30.  .and.  ABS(LAT4(ij)) .le. 50.)then
             if (ZALT90(ik) .le. 14.)  cn(15,ij,ik) = H2OFIX(ij,ik)
          endif

          if (ABS(LAT4(ij)) .gt. 50.  .and.  ABS(LAT4(ij)) .le. 60.)then
             if (ZALT90(ik) .le. 12.)  cn(15,ij,ik) = H2OFIX(ij,ik)
          endif

          if (ABS(LAT4(ij)) .gt. 60.) then
             if (ZALT90(ik) .le. 11.)  cn(15,ij,ik) = H2OFIX(ij,ik)
          endif
 8891   CONTINUE



ccccccc        do 8889 ik=1,Z$
ccccccc       do 8889 ij=1,L$
ccccccc          if (ZALT90(ik) .le. 7.)
ccccccc     >        cn(15,ij,ik) = HALH2O(ij,ik)*1.d-6*m(ij,ik)
ccccccc 8889   CONTINUE


C  for 9GL run, reset washout rate based on WACCM latent heating distribution 
C         and it's the SAME for the FIXED and COUPLED MODELS:
C
C   WLH(91,117,360), WHACK(91,117,360) are in K/sec,  -> rainstr(L$S=91,Z$S=117), RAINCHEM(L$,Z$) is in 1/sec
C   
C  to get rainout rate for current day:
C    take absolute value of sum, divide by 50, interpolate to chemistry grid, set to 0 ABOVE 15 km
C                                   LATST4(L$S), ZSTR4(Z$S) -> LAT4(L$), ZALT90(Z$) all are REAL*4 in COMMON

      do 5677 ik=1,117
      do 5677 ij=1,91
 5677  rainstr(ij,ik)=ABS(WLH(ij,ik,iday360))/50.    !  + WHACK(ij,ik,iday360))/50.


            CALL BINTERP(LATST4, L$S, ZSTR4, Z$S, RAINSTR,
     >                   LAT4,   L$, ZALT90, Z$, 0, 0, RAINCHEM)


      do 7756 ik=1,Z$
      do 7756 ij=1,L$
           cn(59,ij,ik) = RAINCHEM(ij,ik)
           if (ZALT90(ik) .gt. 16.) cn(59,ij,ik) = 0.d0
 7756   CONTINUE




C
C   Now compute rainfall for a year for each latitude, and total # of H2O 
C       molecules rained out in a year for the entire globe for LAST 2 YEARS OF RUN
C          
cegf       IF (ITDAY .GE. INDAYS-360-359) THEN
cegf 
cegf            DO 201 ILAT=1,L$
cegf            SUMZ = 0.0
cegf            DO 200 IZ = IBW,11
cegf 200        SUMZ = SUMZ + (PRECIPYR(ILAT,IZ,IDAY)*DELTAZ(ILAT,IZ)*1.E5)
cegf 201        RAINDAY(ILAT,IDAY) = SUMZ
C
C        SUML = 0.0
C        DO 203 ILAT=1,18
C203     SUML = SUML + (RAINDAY(ILAT,IDAY)*AREA(ILAT)*GLOBE)
C        TOTRAINDAY = SUML
C        WRITE(52,709)
C709    FORMAT(/,5X,'LATITUDE  RAINOUT  1.E+20 # mole / cm2')
C       WRITE(52,707) (RAINDAY(IL,IDAY),IL=1,18)
C        ENDIF
C
cegf 
cegf       IF (IDAY .EQ. 360) THEN
cegf         DO 668 JLAT=1,18
C 
cegf         SUM = 0.0
cegf         DO 667 JDAY=1,360
cegf 667     SUM = SUM + RAINDAY(JLAT,JDAY)
C
cegf 668     RAINYR(JLAT) = SUM
C
cegf         TOTH2O = 0.0
cegf         DO 670 ILAT=1,18
cegf 670     TOTH2O = TOTH2O + (RAINYR(ILAT)*AREA(ILAT)*GLOBE)
C
cegf        WRITE(52,810) IYR, DAY360, IDAY, ITDAY
cegf 810    FORMAT(//,10X,'IYR =',I4,7X,'DAY360 =',F6.1,5X,
cegf      c      'IDAY=',I4,5X,'ITDAY =',I6)
C
cegf         WRITE(52,811) 
cegf 811     FORMAT(/,3X,'TOTAL RAINOUT for the year at each Latitude is',
cegf      *  '  --  E+20 # molecules H2O/cm-2 :')
cegf         WRITE(52,707) (RAINYR(IJ),IJ=1,18)
cegf 707     FORMAT(5X,1P6E14.4)
C
cegf         WRITE(52,813) TOTH2O 
cegf 813   FORMAT(/,1X,'Total  RAINOUT for the globe',
cegf      *' for 1 year  -- E+20 # molecules H2O = ',1PE14.4,/)
C
cegf         DO 901 JDAY=1,360
cegf         DO 901 ILAT=1,18
cegf  901    RAINDAY(ILAT,JDAY) = 0.0
cegf         
cegf         DO 902 ILAT=1,18
cegf  902    RAINYR(ILAT) = 0.0
C
cegf 
cegf       ENDIF
cegf       ENDIF
cegf 
cegf 

ccrain      IF (IDAY .EQ. 360) THEN
ccrain        DO 910 IKD = 1,360 
ccrain        DO 910 IJD = 1,Z$
ccrain        DO 910 ID = 1,L$
ccrain 910    PRECIPYR(ID,IJD,IKD) = 0.0
ccrain      ENDIF
C

	RETURN
	END
