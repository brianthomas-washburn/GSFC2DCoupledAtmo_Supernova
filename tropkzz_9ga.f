C
          SUBROUTINE TROPKZZ(im10, iiyr, im324) 
c
C
c  file TROPKZZ_9dg.f   FORTRAN routine TO  PRESCRIBE  Tropospheric/lower stratospheric Kzz's, 
C      and STRATOSPHERIC BGD Kzz,  BASED ON  temperature lapse rate and Kyy/Kyz  
C
C    Kyz is first computed at the grid box corners for use in AER adapted diffusive flux routiune (NEWDIFYZ)
C
C  This is updated (BASE8SA) to make sure that the diffusion matrix is always positive (diffusive)
C    such that, Kyy*Kzz > Kyz*Kyz    
C    
C
C    updated Sept. 2005 to better simulate 10-20 km (IMPORTANT for TOTAL OZONE)
C
C    updated again, May 2006 to better simulate subtropics/midlat (12-50N,S), at 15-30 km
C                   and eliminate separate 35-55N,S loop


        include "com2d.h"


 	SAVE


        INTEGER TROPIND(L$), P320IND(L$)
        REAL PT(L$,Z$X), kfac, kyysh(72), kyynh(72), kyyun(L$+1,Z$X)

        REAL*8 yylf, ubark(L$+1,Z$X), xn2k(L$+1,Z$X), ubarkz8(L$,Z$X1)
        REAL*8 vps0, vps1, wsat0, wsat1, ww0, ww1, wz

        REAL*4 ubarkz(L$,Z$X1), tempkz(L$,Z$X1)


C  define minimum Kyy in polar regions for 10-20 km: (function of season)

        DATA kyysh/18*20., 30., 40., 50., 47*50., 40., 30., 2*20./

        DATA kyynh/21*50., 40., 30., 20., 33*20., 30., 40., 50., 12*50./



cccc        REAL dptdz(L$+1,Z$X1), dptdy(L$+1,Z$X1), tslope(L$+1,Z$X1)
cccc        REAL kzzmin(L$+1,Z$X1), TEMP4(L$,Z$X1), dtdz2(L$,Z$X1)



C  first compute potential temperature, TEMPALL(L$,Z$X,74), PRESS(Z$X)

      do 33 ikt=1,Z$X
      do 33 ij=1,L$
 33      pt(ij,ikt) = TEMPALL(IJ,IKT,im10+1)*(1013./press(ikt))**.286


C  **************   FIND  TROPOPAUSE   **************************************************
c
C find tropopause  and  320 K pot temp sfc (not currently used) ;;;tropind = intarr(18), p320ind = intarr(18)
C    For re-specified tropospheric Kyy's based on 320 theta sfc and trop hgt, see  KYYBASE6C.F in ../base6
C  ALSO, if TROPHTF is between 2D model pressure sfcs, should define the TROPIND to be the LOWER (altitude)
C   of the two surfaces, as we have done here, to avoid spurrious transport across tropopause
C
C First find tropopause height index for all 360 days, TROPHTF, ITROP360(L$,360) - do this only once (im10=1)
C     THIS IS FOR DIFFUSION (KYY, KZZ), SO DO ON CONSTITUENT GRID, L$, Z$X, LAT(L$), ZALT(Z$X), PRESS(Z$X)
C
C   TROPHTF range is 333.128 - 98.8221 mbar
C

       IF (im10 .eq. 1) THEN
          DO 700 idy=1,360
          DO 700 ij=1,L$ 
                                                      ! initialize to level 5 - DON'T NEED TO DO THIS
cef         ITROP360(ij,idy) = 5

            DO 701 IK=1,Z$X-1
             IF (TROPHTF(ij,idy) .le. press(ik) .and.
     >           TROPHTF(ij,idy) .ge. press(ik+1)) ITROP360(ij,idy) = ik
 701	    continue
 700      CONTINUE
       ENDIF

C                                      need to define tropind(L$) for the current 5-day period
        ij360 = int((im10-1)*5 + 3)   

        DO 400 ij=1,L$
          tropind(ij) = ITROP360(ij,ij360) 

          DO 402 IK=1,16
 402  IF(PT(ij,ik) .le. 320. .and. PT(ij,ik+1) .ge. 320.) P320IND(ij)=ik
 400	CONTINUE



C
C  *************  KYY  upgrade  CALCULATION   ***********************************
C    
C    use UBAR:   UBAR1(L$S,Z$S) in COMMON
C  
C     first interpolate UBAR, XN2 (REAL*8) to  (L$+1, Z$X) -->  ubark(L$+1,Z$X), xn2k(L$+1,Z$X) (REAL*8) 
C     LATST(L$S), ZSTR(Z$S), LATEG(L$+1), ZALT8(Z$X)  are all REAL*8 in COMMON
C
cccc
cccc            CALL BINTERP8(LATST, L$S, ZSTR, Z$S, UBAR1,
cccc     >                    LATEG, L$+1, ZALT8, Z$X, 0, 0, UBARK)
cccc
CUN
CUN            CALL BINTERP8(LATST, L$S, ZSTR, Z$S, XN2,
CUN     >                    LATEG, L$+1, ZALT8, Z$X, 0, 0, XN2K)
CUN
C
C  compute KYYUN(L$+1,Z$X) REAL*4, just go up to Z$, convert to cm2/sec
C
CUN         do 791 ik=1,Z$
CUN         do 791 ij=1,L$+1
CUN             if (ubark(ij,ik) .ge. 0.) kyyun(ij,ik) = 1.e8*
CUN     >          (250. - ubark(ij,ik)*4. - (xn2k(ij,ik)*86400.- 7.)*5.)
CUN
CUN             if (ubark(ij,ik) .lt. 0.) kyyun(ij,ik) = 1.e8*
CUN     >      (250. - ABS(ubark(ij,ik))*20. - (xn2k(ij,ik)*86400.- 7.)*5.)
CUN
CUN             if (kyyun(ij,ik) .lt. 1.e8) kyyun(ij,ik) = 1.e8
CUN 791     CONTINUE

C
C
C  **************   TROPOSPHERIC  Kyy  ADJUSTMENT   ******************************

C   adjust Kyy in troposphere to get correct interhemispheric transport for SF6 
C   and CO2 in tropical upper troposphere, set Kyy to minimum defined above blend over 1 level 
C   to get a smooth transistion. Don't worry about consistency with QBARY, 
C   since this is just the TROPOSPHERE!        - ONLY DO UP TO 8 km, don't adjust upper tropical troposphere
C
C   TROPIND(L$), PRESS(Z$X), ZALT(Z$X), KYYALL(L$+1,Z$X,74)
C
C - OLD:
C      DO 111 ik=1,Z$X
C        IF (zalt(ik) .le. 8.) then 
C
C          if (zalt(ik) .le. 4.) kyyad=200.e8
C          if (zalt(ik) .gt. 4. .and. zalt(ik) .le. 8.) kyyad=100.e8
C
C          do 112 ij=1,L$+1
C           if (kyyall(ij,ik,im10+1) .lt. kyyad) kyyall(ij,ik,im10+1)=kyyad
C 112      continue
C
C       ENDIF
C 111   CONTINUE
C
C
C   increase winter-spring Kyy at 30-50S,N, do from 12-25 km
C      assume that QBARY min is smaller (instead of 0.3e-11, reasonable assumption), so NO change in DELF
C
C       KYYALL(L$+1,Z$X,74), in cm2/sec ;  LATEG4(L$+1), ZALT(Z$X)    - key on UBARK(L$+1,Z$X)
C                            also blend with latitude, altitude transition
cccc
cccc      do 7730 ik=1,Z$
cccc       IF (zalt(ik) .ge. 12.  .and.  zalt(ik) .le. 25.) then
cccc          zfac = 1.
cccc         if (zalt(ik) .gt. 22.)zfac=(.75-1.)/(25.-22.)*(zalt(ik)-22.)+1.
cccc         if (zalt(ik) .lt. 15.)zfac=(.75-1.)/(12.-15.)*(zalt(ik)-15.)+1.
cccc
cccc
cccc       do 7732 ij=1,L$+1
cccc        IF (ABS(lateg4(ij)) .le. 55. .and. ABS(lateg4(ij)) .ge. 25.)then
cccc            lfac = 1.
cccc
cccc          IF (ABS(lateg4(ij)) .gt. 50.)
cccc     >      lfac = (.75-1.)/(55.-50.)*(ABS(lateg4(ij)) - 50.) + 1.
cccc
cccc          IF (ABS(lateg4(ij)) .lt. 30.)
cccc     >      lfac = (.75-1.)/(25.-30.)*(ABS(lateg4(ij)) - 30.) + 1.
cccc
cccc
ccccC                                ! set Kyy factor of 1.3 based on UBAR, ramp down for < -3,   > 30.
cccc         fmax = 1.3
cccc         kfac0 = 1.
cccc         if (ubark(ij,ik) .gt. -3.  .and.  ubark(ij,ik) .lt. 1.)
cccc     >       kfac0 = (fmax- 1.)/(1. - (-3.))*(ubark(ij,ik) - (-3.)) + 1.
cccc
cccc         if (ubark(ij,ik) .ge. 1. .and. ubark(ij,ik) .le. 25.)kfac0=fmax
cccc
cccc         if (ubark(ij,ik) .gt. 25.  .and.  ubark(ij,ik) .lt. 30.)
cccc     >        kfac0 = (1. - fmax)/(30. - 25.)*(ubark(ij,ik) - 25.) +fmax
cccc
cccc         kfac = kfac0*lfac*zfac
cccc         if (kfac .lt. 1.) kfac = 1.
cccc         if (kfac .gt. fmax) kfac = fmax
cccc
cccc         kfac = 1.
cccc         kyyall(ij,ik,im10+1) = kyyall(ij,ik,im10+1)*kfac
cccc
cccc       ENDIF
cccc 7732  CONTINUE
cccc
cccc       ENDIF
cccc 7730  CONTINUE
cccc



C
C  ***********************************************************************************
C
C
C  NEW KYY ADJUSTMENTS -  March 2007  (EF);     KYYALL(L$+1,Z$X,74) is in cm2/sec 
C
C     set MINIMUM KYY to 1.e8 EVERYWHERE, reduce Kyy 7-15 km, 34N-90N, 
C                              and key on  LATEG4(L$+1),  ZALT(Z$X)
C
        do 930 ik=1,Z$X
        do 930 ij=1,L$+1

          kyf = 1.
          if (lateg4(ij) .ge. 33.) then
              if (zalt(ik) .le. 17.) kyf = .95
              if (zalt(ik) .le. 16.) kyf = .9
              if (zalt(ik) .le. 15.) kyf = .7
              if (zalt(ik) .le. 14. .and. lateg4(ij) .gt. 33.) kyf = .6
              if (zalt(ik) .le. 14. .and. lateg4(ij) .ge. 35.) kyf = .5

              if (zalt(ik) .le. 8.) kyf = .75
              if (zalt(ik) .le. 7.) kyf = 1.
          endif

          kyyall(ij,ik,im10+1) = kyyall(ij,ik,im10+1)*kyf


          kyyl = 1.e8
          if (kyyall(ij,ik,im10+1) .lt. kyyl) kyyall(ij,ik,im10+1)= kyyl

ccccccc   OLD:          kyf = (COS(lateg4(ij)*DTR))**3  -  cos^3 is < .5 poleward of 35N, so do simple reduc.
ccccccc                 if (kyf .le. .5) kyf = .5
 930    CONTINUE


C
C  March 2007 - Kyy INCREASE at 21S-55S for ~13.5-27.5 km - ONLY DO FOR SH - see AURA/MLS data
C    do for May 17 - Nov 2, ramp up/down - use factor of 2 max  (No change for Nov 7 - May 12)
C    ramp up/down w/ latitude/height as usual, KYYALL(L$+1,Z$X,74) is in cm2/sec ;   LATEG4(L$+1),  ZALT(Z$X)
C
C             
C         FOR 9FQ run, turn this OFF.....
C 
cckyy
cckyy      IF (im10 .ge. 28   .and.  im10 .le. 61) then
cckyy
cckyy       tfac = 1.
cckyy       if (im10 .eq. 28) tfac = .6
cckyy       if (im10 .eq. 29) tfac = .7
cckyy       if (im10 .eq. 30) tfac = .8
cckyy       if (im10 .eq. 31) tfac = .9
cckyy
cckyy       if (im10 .eq. 59) tfac = .9
cckyy       if (im10 .eq. 60) tfac = .75
cckyy       if (im10 .eq. 61) tfac = .6
cckyy
cckyy       do 7717 ik=1,Z$X
cckyy       do 7717 ij=1,L$+1
cckyy
cckyy          IF (zalt(ik) .ge. 13.  .and.  zalt(ik) .le. 28.) then
cckyy          IF (lateg4(ij) .ge. -55.  .and.  lateg4(ij) .le. -21.) then
cckyy
cckyy              kfac = 1.1
cckyy              if (zalt(ik) .le. 27.) kfac = 1.4
cckyy              if (zalt(ik) .le. 26.) kfac = 1.7
cckyy              if (zalt(ik) .le. 25.) kfac = 2.
cckyy
cckyy              lfac = 1.
cckyy              if (lateg4(ij) .le. -51.) lfac = .75
cckyy              if (lateg4(ij) .ge. -29.) lfac = .75
cckyy              if (lateg4(ij) .ge. -25.) lfac = .55
cckyy
cckyy              if (zalt(ik) .le. 15.) then 
cckyy                 kfac = 1.5
cckyy                 lfac = 1.
cckyy                 if (lateg4(ij) .ge. -27.) kfac = 1.
cckyy              endif
cckyy
cckyy              if (zalt(ik) .le. 14.) then
cckyy                 kfac = 1.2
cckyy                 lfac = 1.
cckyy                 if (lateg4(ij) .ge. -31.) kfac = 1.
cckyy              endif
cckyy
cckyy              totfac = kfac*tfac*lfac
cckyy              if (totfac .le. 1.) totfac = 1.
cckyy
cckyy              kyyall(ij,ik,im10+1) = kyyall(ij,ik,im10+1)*totfac
cckyy          ENDIF
cckyy          ENDIF
cckyy
cckyy 7717  CONTINUE
cckyy
cckyy      ENDIF


C
C    March 2007 - bring back OT KYYs!!!!!  -  KYYOT(L$+1,Z$,72) is in cm2/sec in COMMON - use as minimum Kyy
C       MLS N2O is VERY flat across tropics below 19 km throughout the year - 
C       NCEP data doesn't get large enough Kyys across the equator;    KYYALL(L$+1,Z$X,74)
C       KYYOT is good for 7-23 km, use for 7-16, blend over 17.5-18.5 km
C                              and use for 30S-30N, blend at 22-30S,N (>32N reduced above)
C
C                              set these here in a separate loop;  LATEG4(L$+1);  ZALT(Z$X)
C                                             minimums are checked below

       do 7756 ik=1,Z$
       do 7756 ij=1,L$+1

         IF (zalt(ik) .ge. 7.  .and.  zalt(ik) .le. 19.) then
           IF (ABS(lateg4(ij)) .lt. 32.) then

             if (zalt(ik) .gt. 18.) then
                kyyav = .75*kyyall(ij,ik,im10+1) + .25*kyyot(ij,ik,im10)

               if (ABS(lateg4(ij)) .lt. 20.) 
     >           kyyav = .9*kyyall(ij,ik,im10+1) + .1*kyyot(ij,ik,im10)

               if (ABS(lateg4(ij)) .lt. 15.) kyyav= kyyall(ij,ik,im10+1)
             endif


             if (zalt(ik) .gt. 17.  .and.  zalt(ik) .le. 18.) then
                 kyyav = .5*kyyall(ij,ik,im10+1) + .5*kyyot(ij,ik,im10)

               if (ABS(lateg4(ij)) .lt. 15.) 
     >           kyyav = .8*kyyall(ij,ik,im10+1) + .2*kyyot(ij,ik,im10)
             endif


             if (zalt(ik) .gt. 16.  .and.  zalt(ik) .le. 17.) then
                kyyav = .5*kyyall(ij,ik,im10+1) + .5*kyyot(ij,ik,im10)

                if (ABS(lateg4(ij)) .lt. 28.) kyyav = kyyot(ij,ik,im10)

                if (ABS(lateg4(ij)) .lt. 15.)
     >            kyyav= .5*kyyall(ij,ik,im10+1) + .5*kyyot(ij,ik,im10)
             endif


             if (zalt(ik) .gt. 15.  .and.  zalt(ik) .le. 16.) then
                kyyav = .5*kyyall(ij,ik,im10+1) + .5*kyyot(ij,ik,im10)

                if (ABS(lateg4(ij)) .lt. 28.) kyyav = kyyot(ij,ik,im10)

                if (ABS(lateg4(ij)) .lt. 15.) 
     >            kyyav= .2*kyyall(ij,ik,im10+1) + .8*kyyot(ij,ik,im10)
             endif


             if (zalt(ik) .le. 15.) then
                kyyav = .75*kyyall(ij,ik,im10+1) + .25*kyyot(ij,ik,im10)

                if (ABS(lateg4(ij)) .le. 28.) 
     >            kyyav = .5*kyyall(ij,ik,im10+1) + .5*kyyot(ij,ik,im10)

                if (ABS(lateg4(ij)) .le. 24.) 
     >          kyyav = .25*kyyall(ij,ik,im10+1) + .75*kyyot(ij,ik,im10)

                if (ABS(lateg4(ij)) .le. 20.) kyyav = kyyot(ij,ik,im10)
             endif


         IF (kyyall(ij,ik,im10+1) .le. kyyav) kyyall(ij,ik,im10+1)=kyyav


       ENDIF
      ENDIF

 7756  CONTINUE


C
C  Specify minimum Kyy below 50 km:  
C      KYYALL(L$+1,Z$X,74), KYYCL(L$+1, Z$X, 72)  - all are in cm2/sec 
C      make sure Kyy is greater than the specified minimum,   LATEG4(L$+1);  ZALT(Z$X)
C

      do 730 ik=1,Z$

         if (zalt(ik) .le. 50.) then

      do 732 ij=1,L$+1
C
c    first set tropics/midlats - go up to 50N,S:  
C    Note that Kyy mins and maxs have been set appropriately at each level, 
C           but not necessarily consistent from level to level since ACTUAL Kyy changes from level to level
C           
C   for 9FQ run, set min Kyy = 7.e8 for 29-40 km (instead of 10.e8)
C

       IF (ABS(lateg4(ij)) .le. 50.) then

        if (zalt(ik) .ge. 39.7) then
           kyyf = 10.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
        endif


        if (zalt(ik) .ge. 29.7  .and.  zalt(ik) .lt. 39.7) then
           kyyf = 7.e8
cccccc           kyyf = 10.e8
cccccc           if (ABS(lateg4(ij)) .le. 10.) kyyf = 7.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
        endif


        if (zalt(ik) .ge. 28.7  .and.  zalt(ik) .lt. 29.7) then
           kyyf = 7.e8
ccccccc           kyyf = 10.e8
ccccccc           if (ABS(lateg4(ij)) .le. 15.) kyyf = 7.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
        endif


        if (zalt(ik) .ge. 27.7  .and.  zalt(ik) .lt. 28.7) then
           kyyf = 7.e8
ccccccc          if (ABS(lateg4(ij)) .le. 15.) kyyf = 5.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
        endif


        if (zalt(ik) .ge. 26.7  .and.  zalt(ik) .lt. 27.7) then
           kyyf = 5.e8
           if (ABS(lateg4(ij)) .le. 20.) kyyf = 3.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
                                                                            ! also set Max Kyy in tropics
ccx           kyyx = 1000.e8
ccx           if (ABS(lateg4(ij)) .le. 10.) kyyx = 10.e8
ccx           if (kyyall(ij,ik,im10+1) .gt. kyyx) kyyall(ij,ik,im10+1)=kyyx
        endif


        if (zalt(ik) .ge. 21.7  .and.  zalt(ik) .lt. 26.7) then
           kyyf = 5.e8
           if (ABS(lateg4(ij)) .le. 20.) kyyf = 1.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
                                                                            ! also set Max Kyy in tropics
ccx           kyyx = 1000.e8
ccx           if (ABS(lateg4(ij)) .le. 10.) kyyx = 5.e8
ccx           if (kyyall(ij,ik,im10+1) .gt. kyyx) kyyall(ij,ik,im10+1)=kyyx
        endif


        if (zalt(ik) .ge. 20.7  .and.  zalt(ik) .lt. 21.7) then
           kyyf = 5.e8
           if (ABS(lateg4(ij)) .le. 20.) kyyf =  1.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
                                                                            ! also set Max Kyy in tropics
ccx           kyyx = 1000.e8
ccx           if (ABS(lateg4(ij)) .le. 10.) kyyx = 7.e8
ccx           if (kyyall(ij,ik,im10+1) .gt. kyyx) kyyall(ij,ik,im10+1)=kyyx
        endif


        if (zalt(ik) .ge. 19.7  .and.  zalt(ik) .lt. 20.7) then
           kyyf = 5.e8
           if (ABS(lateg4(ij)) .le. 20.) kyyf = 1.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
                                                                            ! also set Max Kyy in tropics
ccx           kyyx = 1000.e8
ccx           if (ABS(lateg4(ij)) .le. 10.) kyyx = 10.e8
ccx           if (kyyall(ij,ik,im10+1) .gt. kyyx) kyyall(ij,ik,im10+1)=kyyx
        endif


        if (zalt(ik) .ge. 18.7 .and. zalt(ik) .lt. 19.7) then
           kyyf = 5.e8
           if (ABS(lateg4(ij)) .le. 20.) kyyf =  1.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
                                                                            ! also set Max Kyy in tropics
ccx           kyyx = 1000.e8
ccx           if (ABS(lateg4(ij)) .le. 10.) kyyx = 20.e8
ccx           if (kyyall(ij,ik,im10+1) .gt. kyyx) kyyall(ij,ik,im10+1)=kyyx
        endif


        if (zalt(ik) .ge. 17.7 .and. zalt(ik) .lt. 18.7) then
           kyyf = 5.e8
           if (ABS(lateg4(ij)) .le. 20.) kyyf =  1.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
                                                                            ! also set Max Kyy in tropics
ccx           kyyx = 1000.e8
ccx           if (ABS(lateg4(ij)) .le. 10.) kyyx = 30.e8
ccx           if (kyyall(ij,ik,im10+1) .gt. kyyx) kyyall(ij,ik,im10+1)=kyyx
        endif

                                                                            
        if (zalt(ik) .ge. 15.7 .and. zalt(ik) .lt. 17.7) then
           kyyf = 5.e8
           if (ABS(lateg4(ij)) .le. 20.) kyyf = 1.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
        endif

                                                                        
        if (zalt(ik) .ge. 12.7 .and. zalt(ik) .lt. 15.7) then
           kyyf = 5.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
        endif


        if (zalt(ik) .ge. 10.7  .and.  zalt(ik) .lt. 12.7) then
           kyyf = 7.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
        endif

        if (zalt(ik) .ge. 7.7  .and.  zalt(ik) .lt. 10.7) then
           kyyf = 10.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
        endif

        if (zalt(ik) .ge. 6.7  .and.  zalt(ik) .lt. 7.7) then
           kyyf = 10.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
        endif

        if (zalt(ik) .ge. 5.7  .and.  zalt(ik) .lt. 6.7) then
           kyyf = 20.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
        endif

        if (zalt(ik) .ge. 4.7  .and.  zalt(ik) .lt. 5.7) then
           kyyf = 50.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
        endif

        if (zalt(ik) .ge. 3.7 .and. zalt(ik) .lt. 4.7) then
           kyyf = 100.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
        endif

         if (zalt(ik) .lt. 3.7) then
           kyyf = 200.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
        endif

      ENDIF
                              ! end tropics/subtropics loop

ccmid 
ccmid   -  MIDLATITUDE LOOP - just use as transition zone, 51-55N,S
 
      IF (ABS(lateg4(ij)) .gt. 50.  .and. ABS(lateg4(ij)) .le. 55.) then
 
         if (zalt(ik) .ge. 12.7) then
           kyyn = 5.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyn) kyyall(ij,ik,im10+1)=kyyn
         endif

         if (zalt(ik) .ge. 10.7  .and.  zalt(ik) .lt. 12.7) then
            kyyn = 7.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyn) kyyall(ij,ik,im10+1)=kyyn
         endif

         if (zalt(ik) .ge. 7.7 .and. zalt(ik) .lt. 10.7) then
           kyyn = 10.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyn) kyyall(ij,ik,im10+1)=kyyn
         endif
 
         if (zalt(ik) .ge. 6.7  .and.  zalt(ik) .lt. 7.7) then
           kyyf = 10.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
         endif

         if (zalt(ik) .ge. 5.7 .and. zalt(ik) .lt. 6.7) then
           kyyn = 20.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyn) kyyall(ij,ik,im10+1)=kyyn
         endif

         if (zalt(ik) .ge. 4.7 .and. zalt(ik) .lt. 5.7) then
           kyyn = 50.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyn) kyyall(ij,ik,im10+1)=kyyn
         endif
 
         if (zalt(ik) .ge. 3.7 .and. zalt(ik) .lt. 4.7) then
           kyyn = 100.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyn) kyyall(ij,ik,im10+1)=kyyn
         endif
 
         if (zalt(ik) .lt. 3.7) then
           kyyn = 200.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyn) kyyall(ij,ik,im10+1)=kyyn
        endif
 
       ENDIF
                                ! end mid-latitude loop


C  for polar region, ramp down Kyy faster w/ height than at midlatitudes to avoid blowing up the vortex. 
C     Note that this only matters for the SH, since the NH Kyy is soooo large most of the time, except
C     in summer when it doesn't matter anyway (things are already well mixed - not much going on.....)
C     and set smaller minimum in middle world and above
C     specify seasonally dependent min Kyy for 10-20 km (defined at top), ramp down at 18-19 km
C
       IF (ABS(lateg4(ij)) .gt. 55.) then

        if (zalt(ik) .ge. 15.7) then
           kyyn = 1.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyn) kyyall(ij,ik,im10+1)=kyyn
        endif

        if (zalt(ik) .ge. 12.7 .and. zalt(ik) .lt. 15.7) then
           kyyn = 3.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyn) kyyall(ij,ik,im10+1)=kyyn
        endif

        if (zalt(ik) .ge. 9.7 .and. zalt(ik) .lt. 12.7) then
           kyyn = 5.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyn) kyyall(ij,ik,im10+1)=kyyn
        endif

        if (zalt(ik) .ge. 7.7 .and. zalt(ik) .lt. 9.7) then
           kyyn = 10.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyn) kyyall(ij,ik,im10+1)=kyyn
        endif

         if (zalt(ik) .ge. 6.7  .and.  zalt(ik) .lt. 7.7) then
           kyyf = 10.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyf) kyyall(ij,ik,im10+1)=kyyf
         endif

         if (zalt(ik) .ge. 5.7 .and. zalt(ik) .lt. 6.7) then
           kyyn = 20.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyn) kyyall(ij,ik,im10+1)=kyyn
         endif

         if (zalt(ik) .ge. 4.7 .and. zalt(ik) .lt. 5.7) then
           kyyn = 50.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyn) kyyall(ij,ik,im10+1)=kyyn
         endif

        if (zalt(ik) .ge. 3.7 .and. zalt(ik) .lt. 4.7) then
           kyyn = 100.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyn) kyyall(ij,ik,im10+1)=kyyn
        endif

        if (zalt(ik) .lt. 3.7) then
           kyyn = 200.e8
           if (kyyall(ij,ik,im10+1) .lt. kyyn) kyyall(ij,ik,im10+1)=kyyn
        endif

      ENDIF
                             ! end high latitude loop

 732       CONTINUE

        endif      

c                               ! 50-55km, set MIN Kyy=15.e8  AT ALL LATITUDES, ramp up to 20.e8 above 55 km

       if (zalt(ik) .gt. 50.  .and.  zalt(ik) .le. 55.) then
         do 742 ij=1,L$+1
         if (kyyall(ij,ik,im10+1) .lt. 15.e8) kyyall(ij,ik,im10+1)=15.e8
 742     continue
       endif


       if (zalt(ik) .gt. 55.) then
         do 744 ij=1,L$+1
         if (kyyall(ij,ik,im10+1) .lt. 20.e8) kyyall(ij,ik,im10+1)=20.e8
 744     continue
       endif


 730   CONTINUE


C  ***********  END  KYY   ADJUSTMENTS  ***********************************************
C
C
C
C
C
C  ****************   SET UP  PRESCRIBED  TROPOSPHERIC/LOWER  STRATOSPHERIC  KZZ's  *********************
C         NOTE: KZZs are all defined at BOX EDGES, use ZALTE(Z$X1)
C
C   Initialize  KZZTROP(L$,Z$X1) in cm2/sec  IN COMMON   (defined at top/bottom BOX EDGES)
C
	DO 200 IK=1,Z$X1
	DO 200 IJ=1,L$
 200       KZZTROP(IJ,IK) = 0.0


c  determine 30 km and 55 km indicies (box edges),  ZALTE(Z$X1)

       do 644 ik=1,Z$X
          if (zalte(ik) .le. 30. .and. zalte(ik+1) .gt. 30.) i30=ik
          if (zalte(ik) .le. 40. .and. zalte(ik+1) .gt. 40.) i40=ik
          if (zalte(ik) .le. 55. .and. zalte(ik+1) .gt. 55.) i55=ik
 644   CONTINUE



C
C   Sept 2006, include Kzz proportional to N^2, du/dz 
C      -  use TEMPALLS(L$S,Z$S) - REAL*4,   and   UBAR1(L$S,Z$S) - REAL*8,  both in COMMON
C  
C     need to interpolate to (L$, Z$X1) -->  ubarkz8(L$,Z$X1) in REAL*8, convert to REAL*4 
C         ubarkz(L$,Z$X1)  is REAL*4 defined above
C
C     LATST(L$S), ZSTR(Z$S), LAT(L$) , ZALTE8(Z$X1)  are all REAL*8 in COMMON
C

            CALL BINTERP8(LATST, L$S, ZSTR, Z$S, UBAR1,
     >                    LAT, L$, ZALTE8, Z$X1, 0, 0, UBARKZ8)

            do 556 ik=1,Z$X1
            do 556 ij=1,L$
 556           ubarkz(ij,ik) = ubarkz8(ij,ik)



C  LATST4(L$S), ZSTR4(Z$S), LAT4(L$), ZALTE(Z$X1), TEMPALLS(L$S,Z$S)  are REAL*4 in COMMON
C        TEMPKZ(L$,Z$X1)  is REAL*4  defined above


            CALL BINTERP(LATST4, L$S, ZSTR4, Z$S, TEMPALLS,
     >                   LAT4, L$, ZALTE, Z$X1, 0, 0, TEMPKZ)

C
C
c new: for Kzz, just base it on dtheta/dz/theta everywhere up to 30 km (similar to dT/dz, but a bit better)
C
C   for Sept 2006, now base it on N^2 - very similar spatial distribution to dtheta/dz/theta 
C        see also Andrews et al., 1987 for static stability proportional to dln(theta)/dz
C
C  KZZTROP(L$,Z$X1) in cm2/sec  IN COMMON,  ZALTE(Z$X1);  use theta, dtheta/dz; PRESS(Z$X)   
C     and set min Kzztrop here to .006e4 (not .01e4) for proper tape recorder propagation
C
C   NOW:  just use B-V frequency - use leap-frog differencing;    TEMPKZ(L$,Z$X1); ZALTE(Z$X1)
C
      do 1730 ik=1,i30
      do 1730 ij=1,L$

            if (ik .eq. 1) then
              th0 = TEMPKZ(ij,1)                         ! *(1013./press(1))**.286
              th1 = TEMPKZ(ij,2)                         ! *(1013./press(2))**.286
              dzalt = (zalte(2) - zalte(1))*1000.

              pr0 = presse8(1) 
              pr1 = presse8(2) 
            endif

            if (ik .ge. 2) then
              th0 = TEMPKZ(ij,ik-1)                      ! *(1013./press(ik-1))**.286
              th1 = TEMPKZ(ij,ik+1)                      ! *(1013./press(ik))**.286
              dzalt = (zalte(ik+1) - zalte(ik-1))*1000.

              pr0 = presse8(ik-1)
              pr1 = presse8(ik+1)
            endif

C
C   also include theta-E contribution:   PRESSE8(Z$X1) is REAL*8
C
        vps0 = 1013.d0/760.d0*DEXP(-88.44375824d0 + 0.7641223392d0*th0 - 
     >           0.002285011125d0*(th0**2) + 2.539292763D-6*(th0**3))

        vps1 = 1013.d0/760.d0*DEXP(-88.44375824d0 + 0.7641223392d0*th1 - 
     >           0.002285011125d0*(th1**2) + 2.539292763D-6*(th1**3))

            wsat0 = 0.622d0*vps0/(pr0 - vps0)
            wsat1 = 0.622d0*vps1/(pr1 - vps1)

            ww0 = 2.5d6*wsat0/(1004.d0*th0)
            ww1 = 2.5d6*wsat1/(1004.d0*th1)

            wz = (ww1 - ww0)/dzalt


            thavg = TEMPKZ(ij,ik)
       xnt = rr1/hh*((th1 - th0)/dzalt + kap*thavg/hh + thavg*wz)*86400.
                                                                         ! w/ theta-E, set upper limit on N^2
            if (xnt .ge. 50.) xnt = 50.
                                                       
C                                                     !! new method,  upper/lower limits
C                                                     !! also include latitude, altitude dependence
            ykf = SQRT(COS(lat4(ij)*DTR))
            if (ABS(lat4(ij)) .ge. 55.) ykf = COS(lat4(ij)*DTR)
            if (ABS(lat4(ij)) .ge. 60.) ykf = (COS(lat4(ij)*DTR))**2
            if (ykf .le. .1) ykf = .1
                                                    ! include adjustment for SH
            yka = (SIN(lat4(ij)*2.*DTR))**2
            if (lat4(ij) .ge. 0.) yka = 0.

            zkf = (22. - zalte(ik))/5.   ! - .1
            if (zkf .ge. 1.) zkf = 1.
            if (zkf .le. .1) zkf = .1


            kz0 = .25*1.e4
            kz1 = 30.*1.e4*ykf

            tz0 = 30.
            tz1 = 0. + 10.*yka


       kzzn2= exp(alog(kz1/kz0)/(tz1-tz0)*(xnt-tz0) + alog(kz0))*zkf



C  need to ramp down Kzz faster for equatorial region where xnt > tz0, blend w/ latitude at 20-30N,S

            kzs = .01*1.e4
            tzs = 38.
            kzztt = exp(alog(kz0/kzs)/(tz0-tzs)*(xnt-tzs) + alog(kzs))

            if (xnt .ge. tz0) then
              kzz0 = kzzn2

              if (ABS(lat4(ij)) .le. 20.) kzzn2 = kzztt

              if (ABS(lat4(ij)) .gt. 20.  .and.  ABS(lat4(ij)) .le. 25.) 
     >            kzzn2 = (2.*kzztt + 1.*kzz0)/3.

              if (ABS(lat4(ij)) .gt. 25.  .and.  ABS(lat4(ij)) .le. 30.)
     >            kzzn2 = (1.*kzztt + 2.*kzz0)/3.
            endif



C  now include Kzz contribution from du/dz, use dimensionless Richardson Number-like quantity
C     ubarkz(L$,Z$X1), use dzalt (in meters) defined above from ZALTE
         

         if (ik .eq. 1) then
            uu0 = UBARKZ(ij,1)
            uu1 = UBARKZ(ij,2)
         endif

         if (ik .ge. 2) then
            uu0 = UBARKZ(ij,ik-1)
            uu1 = UBARKZ(ij,ik+1)
         endif
C
C
C   take dudz/(N^2)^2,  so rii is also in 1/s^2 as is B-V freq defined above, (convert BV to 1/s^2)
C                                                                 divide by 100. for easier use

         rii = (((uu1 - uu0)/dzalt)**2)/((xnt/86400.)**2)/100.

C                                                                    ramp down contribution at 23-25 km
         if (zalte(ik) .ge. 23.) rii = rii/2.
         if (zalte(ik) .ge. 24.) rii = rii/2.
         if (zalte(ik) .ge. 25.) rii = 0.
                                                            !  also ramp down contribution in tropics
         yrf = 1.
         if (ABS(lat4(ij)) .le. 30.)yrf= SIN((lat4(ij)-15.)*DTR*90./15.)
         if (yrf .le. 0.) yrf = 0.
         if (yrf .ge. 1.) yrf = 1.
        

         kzzri = rii*.275*1.e4*yrf

C                                                             ! set larger limit on kzzri everywhere
         kzzxe = .75*1.e4
         if (kzzri .ge. kzzxe) kzzri = kzzxe


                                                             ! add total Kzz contribution 
         kzztrop(ij,ik) = kzzn2 + kzzri



C   Now limit total Kzztrop in polar region at 12-30 km, ramp down to .5 m2/sec at 12 km
 
       if (zalte(ik) .ge. 11.) then
          if (ABS(lat4(ij)) .ge. 50.) then 
             kzzx = ((.5 - .1)/(11.-19.)*(zalte(ik)-19.) + .1)*1.e4
             if (kzzx .le. .1e4) kzzx = .1*1.e4
             if (kzztrop(ij,ik) .ge. kzzx) kzztrop(ij,ik) = kzzx
          endif
       endif



C   set absolute Min and Max on Kzztrop
C   for Max, set to 20.e4 for 30S-30N, 10.e4 poleward of 40N,S;   with simple interpolation in between
C
         kzmin = .005*1.e4
         if (kzztrop(ij,ik) .lt. kzmin) kzztrop(ij,ik) = kzmin

         kzmax = 30.*1.e4
         if (kzztrop(ij,ik) .gt. kzmax) kzztrop(ij,ik) = kzmax



c     set min Kzz=10 m2/sec in Boundary layer, below 2 km  -  KZZTROP(L$,Z$X1) in cm2/sec  IN COMMON

            if (zalte(ik) .le. 1.) then 
              if (kzztrop(ij,ik) .lt. 10.e4) kzztrop(ij,ik) = 10.e4
            endif

            if (zalte(ik) .gt. 1.  .and.  zalte(ik) .le. 2.) then 
              if (kzztrop(ij,ik) .lt. 5.e4) kzztrop(ij,ik) = 5.e4
            endif

 1730   CONTINUE
                           !! end 0-30 km loop



                           !! loop for 30 km - top levels
      do 2730 ij=1,L$

c                                   30 - 55 km  -   ! note, the 8. here is ALOG(3000.)
        do 623 ik=i30+1, i55
           kzztrop(ij,ik) = exp(alog(3000./kzztrop(ij,i30))/
     >                             (i55+1-i30)*(ik - (i55+1)) + 8.)
 623      continue

c                                   55 km and above, KZZTROP(L$,Z$X1)
        do 626 ik=i55+1, Z$X1
 626       kzztrop(ij,ik) = 3000.

 2730   CONTINUE


C   
C  for ERA-40 data, need to set Kzz limit of 1 m2/sec at 80S-88S, 5-11 km 
C     where it gets up to 100 m2/sec in Nov. 1995, 1996, 1997
C     - this prevents numerical instability at high res when Kyy is large;
C     LAT(L$); ZALTE(Z$X1);  KZZTROP(L$,Z$X1) is in cm2/sec
C
       do 1700 ik=1,Z$X1
       do 1700 ij=1,L$
          if (zalte(ik) .ge. 5.  .and.  zalte(ik) .le. 11.  
     >                           .and.  lat(ij) .le. -80.) then 
              if (kzztrop(ij,ik) .gt. 1.e4) kzztrop(ij,ik) = 1.e4
          endif
1700    CONTINUE



c  now add gravity wave and tropospheric Kzzs together into KZZTOT(L$,Z$X1,74),  
C                                      KZZGW(L$,Z$X1), KZZTROP(L$,Z$X1) are in cm^2/sec

       do 717 ik=1,Z$X1
       do 717 ij=1,L$
 717    KZZTOT(ij,ik,im10+1) = KZZGW(ij,ik) + KZZTROP(ij,ik)


C  ****************  END  TROPOSPHERIC/LOWER  STRATOSPHERIC  KZZ's  *********************
C
C
C
C  *******************************     KYZ   CALCULATION   *********************************************
C
C   call separate routine GETKYZ, first load in proper TEMPs, Kyys, Kzzs
C
C       TEMP(L$,Z$X),          EKYY(L$+1,Z$X),         EKZZ(L$,Z$X1) 
C
C       TEMPALL(L$,Z$X,74), KYYALL(L$+1,Z$X,74), KZZTOT(L$,Z$X1,74)
C

        DO 8015 ik=1,Z$X
        DO 8015 ij=1,L$
         TEMP(ij,ik) = TEMPALL(ij,ik,im10+1)
 8015   CONTINUE

 
c    also set upper Kyy limit of 1.e11 everywhere before computing Kyz, then do again below

        DO 8016 ik=1,Z$X
        DO 8016 ij=1,L$+1
       if (kyyall(ij,ik,im10+1) .gt. 1.d11) kyyall(ij,ik,im10+1) = 1.d11
         EKYY(ij,ik) = KYYALL(ij,ik,im10+1)
 8016   CONTINUE


       DO 8017 ik=1,Z$X1
       DO 8017 ij=1,L$
          EKZZ(ij,ik) = KZZTOT(ij,ik,im10+1)
 8017  CONTINUE

c   
c  also adjust the slope, EKYY and/or EKZZ as necessary to ensure the diffusion matrix is positive
C

         CALL GETKYZ


C
C   load into full array for output, EKYZ(L$+1,Z$X1), KYZALL(L$+1,Z$X1)
C                     also do slope, THSLP(L$+1,Z$X1), THSLPALL(L$+1,Z$X1)
C
       DO 4242 IK=1,Z$X1
       DO 4242 IJ=1,L$+1
         KYZALL(ij,ik) = EKYZ(ij,ik)
         THSLPALL(ij,ik) = THSLP(ij,ik)
 4242   CONTINUE


C                               also,  reload Kyy arrays in case adjustments were made
C                               and re-check for lower (use 1.e8) and upper limits
      do 567 ik=1,Z$X
      do 567 ij=1,L$+1
         if (EKYY(ij,ik) .le. 1.e8) EKYY(ij,ik) = 1.e8
         if (EKYY(ij,ik) .ge. 1000.e8) EKYY(ij,ik) = 1000.e8

         KYYALL(ij,ik,im10+1) = EKYY(ij,ik) 
 567  CONTINUE


c                                                 only Kyy is adjusted, NOT KZZ
cc       do 211 ik=1,Z$X1
cc       do 211 ij=1,L$
cc         KZZTOT(ij,ik,im10+1) = EKZZ(ij,ik)
cc 211  CONTINUE



C  now write out extra quantities on model grid which are NOT stored for every time period
C  in COMMON:   VALL(L$+1,Z$X), KZZTROP(L$,Z$X1), KZZGW(L$,Z$X1), KYZALL(L$+1,Z$X1), THSLPALL(L$+1,Z$X1)

          WRITE (57) vall
          WRITE (57) kzztrop
          WRITE (57) kzzgw
          WRITE (57) kyzall
          WRITE (57) thslpall


Ckyy   
C      for KYY testing, write out everything here on fort.37: 
C            kyy0(L$S,61), delf(L$S,61), qy(L$S,61), wout(L$S,61),  and  EKYY(L$+1,Z$X)

                                                 ! just write out dimensions, etc. 1st time through
        if (im324 .eq. 1) then 
           WRITE (37) L$, L$S, Z$S, Z$X, itdyn, NMON$
           WRITE (37) LAT
           WRITE (37) ZALT
           WRITE (37) PRESS
           WRITE (37) ZSTR
           WRITE (37) PRST
           WRITE (37) LATST
           WRITE (37) ZALTE
           WRITE (37) PRESSE
           WRITE (37) LATEG
         endif

           WRITE (37) delf
           WRITE (37) qy
           WRITE (37) kyy0
           WRITE (37) ekyy
           WRITE (37) wout


	RETURN
	END

