C
        SUBROUTINE XCOUPLED(DTMARK, IYR0, IIDAY360, IYRCT0, Z$XC,
     >           LATC, LATCE, ZALT90C, ZALTC, ZALTEC, EPFX, CNC8, MMC8,
     >       tbaro, ubaro, kyyo, kzzo, wbaro, HEATOUT, COOLOUT, DRAGOUT)


C
C    THIS IS THE COUPLED MODEL DRIVER ROUTINE
C
C
!        USE degree_trig

        INTEGER IYR0, IIDAY360, L$C, Z$C, Z$XC, IZ96, IYRCT0
        INTEGER N$, M$, LC1, ZC1, N11, M11

        include "comcfg.h"
        include "timer.inc"
        include "PARAM.INC"


C                           L$C, Z$C now just set from NS$, MS$ in PARAM.INC; L$C=L$; Z$C=Z$
        PARAMETER(L$C=NS$)
        PARAMETER(Z$C=MS$)


        COMMON /SWITC0/ ISWD(NSWTCHS), LH_SW(12), KZZ_SW(15)

C  
C  COMMONs for COUPLED MODEL X* BBC, fixed model Kzz, Kyy;   latent heating from GEOS4 and WACCM
C      read in TEMPIN,  BBCOUP(91,21,360) is in m2/sec;  KYY9FR, KZZ9FR(45,89,360) are in cm2/sec

        COMMON/CCHI/BBCOUP(91,21,360), XLATBBC(91), ZBBC(21), 
     >              XKZZ9FR(45,89,360), xlatkzz(45), zzkzz(89),
     >              XKYY9FR(46,88,360), xlatkyy(46), zzkyy(88)

C                                                                        ! WLH, WHACK are in K/sec

        COMMON/CCLH/WLH(91,117,360)
ccwconv175 - OLD         COMMON/CCLH/XLHC(91,117,360), WLH(91,117,360), WHACK(91,117,360)


        COMMON/CCHEM/CNC(S$C,NS$,MS$), HEATCHM(NS$,MS$),COOLCHM(NS$,MS$)
     >            , HEATCHML(NS$,MS$), DRAGCHM(7,NS$,MS$)


C  common for UV+VIS heating from PHOTHEAT in K/sec

        COMMON/CPHEAT/PHEAT(5,L$C,Z$C), SORADHEAT(15,L$C,Z$C)


C                                                         common for fixed model clim TEMP, EHF (K/sec)
        COMMON/CEHFX/tempfx(91,117,74), ehfx(91,117,74),  
     >               latfx(91), zzfx(117)



c  common for TD GEOS 4 surface temperatures used in COUPLED model, 1950-2100, g4sens - sensitivity to CO2, w/ seas cycle

        COMMON/CTEMPG4/ TEMPG4(1800,91,3), latg4(91), timeg4(1800), 
     >                  g4sens(15,91,3), timesf(14)


c  common for TD GEOS 4 latent heating xlhsens is sensitivity to CO2, w/ seas cycle,
C     G4LHDAY(91,117) is for current day (K/day);  xlhsens5 is GEOS5 latent heating

        COMMON/CG4LH/ xlhsens(15,91,117), G4LHDAY(91,117)
        COMMON/CG5LH/ xlhsens5(15,91,117)


c  common for coupled model sfc albedo from TOMS (season-lat), 1-14 are months, 15=latitude

        COMMON/CTALBEDO/ TALBEDO(15,180), ALBDAY(180), ALBLAT(180)


c  common for UBAR climatology, used in both fixed and COUPLED models, read in TEMPIN
C                                         and UBAR-CO2 Sensitivity from GEOS5 (in m/sec/ppmv)
        COMMON/CUBAR/ UBARCL(91,117,74)
        COMMON/CUBARCO2/ UBARCO2(91,117)


c  common for coupled model TBAR/UBAR - NCEP differences (lat-hgt-season) from 1979-2006: from TEMPIN
C
        COMMON/CUTDIFF/ TDIFFR(14,45,76), UDIFFR(14,45,76), xlatf(45), 
     >                  zzff(76)


c  common for NCEP DELF for baroclinic waves ONLY - for coupled model (in m/sec^2) from TEMPIN
C                                                                    not currently used
ccbaro        COMMON/CDFBARO/ DFBARO(73,33,360), BALAT(73), BZZ(33)


C   common for PW DRAG computed from NCEP UBAR for coupled model
C
ccccccc        COMMON/CDELF/ IYRCT0


C
C  common for the vertical gradients of the DU/km ozone climatology, for Kzz adjustment
C
        COMMON/CDUKZ/dukz(14,45,89)


C     txx is from TEMPLON called from JPWAVE; tcc is for SETDAILY to generate PDFs
C
        COMMON/CTXX/ txx(72,N$,M$)
        COMMON/CTPDF/tcc(72,L$C,Z$C)



        REAL*8 CNC8(S$C,L$C,Z$C), MMC8(L$C,Z$XC)

        REAL*4 DTMARK, LATC(L$C), LATCE(L$C+1), latfx, latg4
        REAL*4 ZALT90C(Z$C), ZALTC(Z$XC), ZALTEC(Z$XC+1)
        REAL*4 h2oin(L$C,Z$C), co2in(L$C,Z$C), o3in(L$C, Z$C)
        REAL*4 BBCDAY(91,21), KYYDAY(46,88)
        REAL*4 LHDAY(91,117), WLHDAY(91,117), HACKDAY(91,117)
        REAL*4 TEMPDAY(91,117), TEMPDAYC(L$C,Z$XC), TEMPG4DAY(91,3)
        REAL*4 IDAYU(74), u74(74), UBARDAY(91,117), EHFDAY(91,117)
        REAL*4 WLHTOT(91,117), EPFX(L$C,Z$XC,74), EPDAY(L$C,Z$XC)
        REAL*4 TDIFFD(45,76), UDIFFD(45,76), BARODAY(73,33)
        REAL*4 txx0(N$,M$), ccx(L$C,Z$C), ttz(N$+1), ttch(L$C)
        REAL*4 BALAT(73), BZZ(33), ttdiff(N$,M$), ttdiffo(L$C,Z$C)
        REAL*4 KZZ0(45,89), KZZTH(45,89), KZZDAY(45,89)


C  these are the arrays from the coupled model dynamics routines:      
C
        REAL*4 ypf(N$), zpf(M$), yppf(N$+1), zppf(M$+1), kzzc(N$+1,M$+1)
        REAL*4 tempc(N$+1,M$+1), ubarc(N$+1,M$+1),kyyc(N$,M$),wsc(N$,M$)


c   these are the dynamics arrays interpolated to the chemistry grid (as in SETDAILY) for output back to MAIN
c
        REAL*4 tbaro(L$C,Z$XC), ubaro(L$C,Z$XC), kyyo(L$C+1,Z$XC)
        REAL*4 kzzo(L$C,Z$XC+1), wbaro(L$C,Z$XC+1)
        REAL*4 kyyminc, kyymaxc, kzzminc, kzzmaxc
        REAL*4 ttime(1), ttout(1), g4in(1800), g4ss(14), zzg4(3)

        REAL*4 HEATOUT(2,L$C,Z$C), COOLOUT(L$C,Z$C), DRAGOUT(7,L$C,Z$C)
  

c jer common for dynamics writeout
        common /dump_2/ iwrite,yrjer

cjer  set flag for dynamics writeout to fort.10

        iwrite=0
c       if(iyr .ge. 44) then 
           if( mod(iiday360-15,30) .eq. 0) iwrite=1
c       endif



C  Update new timing for coupled model only (timer.f), do this BEFORE call to coupled model dynamics
C     to be consistent w/ fixed model - ie, since IDAY360 is updated at the BEGINNING of the DO 100 loop

        CALL TIME_TICK



C  load in REAL*4 CNC array (mixing ratios) in COMMON for use in COUPLED model
C      REAL*8 CNC8(S$C,L$C=NS$,Z$C=MS$), MMC8(L$C=NS$,Z$XC) -> REAL*4 CNC(S$C,NS$,MS$)
C
        do 133 ik=1,MS$
        do 133 ij=1,NS$
        do 133 isss=1,S$C
 133       CNC(isss,ij,ik) = CNC8(isss,ij,ik)/MMC8(ij,ik)



C
C  interpolate time dependent GEOS4 surface temps to current day:  use IIDAY360, IYR0 = IYR (ie, 1950 = 15)
C   TEMPG4(1800,91,3), latg4(91), timeg4(1800) in COMMON above for 1950-2100  -> TEMPG4DAY(91,3)
C   use LINTERP - works fine - OK to interpolate before first time index/after last time index, keeps constant value
c                                                       ttime(1), ttout(1), g4in(1800)
ccsurf
ccsurf         ttime(1) = (IYR0 + 1935.) + IIDAY360*1./360.
ccsurf
ccsurf         DO 488 ik=1,3
ccsurf         DO 488 ij=1,91
ccsurf
ccsurf            do 489 iii=1,1800
ccsurf 489            g4in(iii) = tempg4(iii,ij,ik)
ccsurf 
ccsurf            CALL LINTERP(timeg4, 1800, g4in, ttime, 1, 0, ttout)
ccsurf
ccsurf            tempg4day(ij,ik) = ttout(1)
ccsurf 488   CONTINUE
ccsurf


C  load in surface temps here BASED ON SENSITIVITY FACTORS AND CO2 MIXING RATIO --> load into TEMPG4DAY(91,3)
C   g4sens(15,91,3), 1st index, i=1-14 are 1-14 months, 15=sensitivity factor (K/ppmv)
C   interpolate seasonal cycle to current day (IIDAY360), timesf(14) ->  ttime(1), ttout(1), g4ss(14)
C       ie, ttout(1) is the 1950 temperature interpolated to current IIDAY360
C
C   CNC(S$C,NS$,MS$) is mixing ratio, just use CO2 BC (which is for current day)
C   ie, use level 1 for any latitude (use 9 so it won't blow up w low res), and get change from 1950 (307.7ppm)
C
         ddco2 = CNC(20,9,1)*1.e6 - 307.7

         ttime(1) = IIDAY360*1.

         DO 588 ik=1,3
         DO 588 ij=1,91

            do 589 iii=1,14
 589            g4ss(iii) = g4sens(iii,ij,ik)
 
            CALL LINTERP(timesf, 14, g4ss, ttime, 1, 0, ttout)

            tempg4day(ij,ik) = ttout(1) + ddco2*g4sens(15,ij,ik)
 588   CONTINUE


C  pw51 - set sfc temps for 90S-86S = 84S (so they're not so cold) - TEMPG4DAY(91,3) 
C
         DO 591 ik=1,3
         DO 591 ij=1,3
 591        tempg4day(ij,ik) = tempg4day(4,ik)



C  
C  also interpolate UBARCL(91,117,74) (m/sec) to current day: UBARDAY(91,117), IDAYU(74), u74(74) - use LINTERP
C      also interpolate clim temps for current day TEMPDAY(91,117) from tempfx(91,117,74) read in TEMPIN
C                                              and interpolate EHF (K/sec)   ehfx(91,117,74) -> EHFDAY(91,117)
C
C   also include CO2 sensitivity factor for UBAR in lower troposphere only
C         UBARCO2(91,117) is in m/sec/ppmv; use avg CO2 for 1979-2006 as basis (356.68 ppmv)
C

         do 9772 ijk=1,74
9772     IDAYU(ijk) = (ijk-1)*5.- 2.


         DO 7700 IK=1,117
         DO 7700 IJ=1,91

             ubco2 = (CNC(20,9,1)*1.e6 - 356.68)*ubarco2(ij,ik)
             if (zzfx(ik) .ge. 5.) ubco2 = 0.


             do 7702 iim=1,74
 7702           u74(iim) = UBARCL(ij,ik,iim)

             CALL LINTERP(idayu, 74, u74, ttime, 1, 0, ttout)

             UBARDAY(IJ,IK) = ttout(1) + ubco2


             do 7703 iim=1,74
 7703           u74(iim) = TEMPFX(ij,ik,iim)

             CALL LINTERP(idayu, 74, u74, ttime, 1, 0, ttout)

             TEMPDAY(IJ,IK) = ttout(1)


             do 7704 iim=1,74
 7704           u74(iim) = EHFX(ij,ik,iim)

             CALL LINTERP(idayu, 74, u74, ttime, 1, 0, ttout)

             EHFDAY(IJ,IK) = ttout(1)
 7700    CONTINUE


C  
C   Interpolate Fixed model DELF to current day: EPFX(L$C,Z$XC,74)-> EPDAY(L$C,Z$XC) in m/sec/day
C                                      IDAYU(74), u74(74), ttime(1), ttout(1) - use LINTERP

         DO 7750 IK=1,Z$XC
         DO 7750 IJ=1,L$C

             do 7751 iim=1,74
 7751           u74(iim) = EPFX(ij,ik,iim)

             CALL LINTERP(idayu, 74, u74, ttime, 1, 0, ttout)

             EPDAY(IJ,IK) = ttout(1)
 7750    CONTINUE


C
C  load in daily GCM latent heating and X* BBC on STREAMF grid to current day for use in TWODR
C    BBCOUP(91,21,360) -> BBCDAY(91,21), XLATBBC(91), ZBBC(21) -  BBCOUP/BBCDAY are in m2/sec
C    XLHC(91,117,360) -> LHDAY(91,117) - NO LONGER USED (WCONV175)
C
C    also load in WACCM latent heating, WLH, WHACK are in K/sec, convert to K/day
C     -> WLHDAY(91,117), HACKDAY(91,117) are in K/day,  also get total WLHTOT(91,117) in K/day
C
C  add in extra latent heating if needed, weight by sin(lat)^4, set to 0 above 25 km
C                                   use latfx(91), zzfx(117) on STREAMF grid
C
C   WCONV147 - WLH and HACK combined into one array (WLH, offline) and cleaned up
C
C    WCONV215 - put seasonal mask on NH LHeat enhancement (much weaker in NH summer)
C
            seasf = 1.
            if (IIDAY360 .ge. 131  .and.  IIDAY360 .le. 251)
     >          seasf = SIND(IIDAY360*1.45 + 84.)**4
            if (seasf .le. 0.) seasf = 0.
            if (seasf .ge. 1.) seasf = 1.

            DO 8505 IK=1,117
            DO 8505 IJ=1,91
               LHDAY(IJ,IK) = 0.   ! wconv175 - XLHC(IJ,IK,IIDAY360)

               WLHDAY(IJ,IK) = WLH(IJ,IK,IIDAY360)*86400.
               HACKDAY(IJ,IK) = 0.     ! WHACK(IJ,IK,IIDAY360)*86400.

               xlt0 = WLHDAY(IJ,IK)    !  + HACKDAY(IJ,IK)   ! don't reduce HACK

               xltfac = 1.05 + 1.25*((SIND(latfx(ij)))**4)
               if (latfx(ij) .ge. 0.)
     >             xltfac = 1.05 + seasf*1.25*((SIND(latfx(ij)))**4)

               if (zzfx(ik) .le. 2.) xltfac = (1. + xltfac)/2.
               if (zzfx(ik) .le. 1.) xltfac = 1.
C                                                   reduce xltfac at 2km, set to 1 at 0-1 km
cpw92/93                             xltfac = 1.
               xlt = xlt0*xltfac
C                                                                 ! set limit above 5 km, and > 50N,S
               if (zzfx(ik) .ge. 5.  .and.  xlt .gt. 2.8) xlt = 2.8
  
               xltx = .5 + 2.4*COSD(latfx(ij))
               if (ABS(latfx(ij)) .ge. 50. .and. xlt .gt. xltx) xlt=xltx

               if (zzfx(ik) .GE. 20.) xlt = 0.0
               if (xlt .le. 0.) xlt = 0.

               WLHTOT(IJ,IK) = xlt
 8505       CONTINUE




C
C  load in GEOS 4 TD latent heating here BASED ON SENSITIVITY FACTORS AND CO2 MIXING RATIO --> load into G4LHDAY(91,117)
C   xlhsens(15,91,117), 1st index, i=1-14 are 1-14 months, 15=sensitivity factor (K/day/ppmv)
C       just use GEOS 4 sensitivity factor w/ WACCM seasonal cycle, not worth scaling to WACCM
C   interpolate seasonal cycle to current day (IIDAY360), timesf(14) ->  ttime(1), ttout(1), g4ss(14)
C
C  CNC(S$C,NS$,MS$) is mixing ratio, just use CO2 BC (which is for current day)
C  use level 1 for any latitude (use 9 so it won't blow up w low res), and get change from 1950-2100 avg (468.35 ppm)
C
C   NOTE: GEOS 4 Latent heating seasonal cycle  DOESN'T WORK - too COLD! 
C   so use WACCM seasonal cycle (1999-2003 avg) + CO2 sensitivity from GEOS 4 (use deviation from 2001 CO2, ie, 370 ppmv)
C         WLHTOT(91,117) is in K/day;   G4LHDAY(91,117) is in K/day
C
C  WCONV190 - use GEOS 5 latent heat sensitivity, xlhsens5(15,91,117)
C     as w/ GEOS4, use deviation from 2001 CO2, ie, 370 ppmv)
C
C   both GEOS4 and GEOS5 latent heating sensitivity have been cleaned up (WCONV147)
C                                 also ensure that total latent heating is POSITIVE

         ddco2 = CNC(20,9,1)*1.e6 - 370.

         DO 788 ik=1,117
         DO 788 ij=1,91
              g4lhday(ij,ik) = WLHTOT(ij,ik) + ddco2*xlhsens5(15,ij,ik)
              if (g4lhday(ij,ik) .lt. 0.) g4lhday(ij,ik) = 0.
 788     CONTINUE


c
Cg4lh  to use GEOS 4 latent heating - temps get too cold
Cg4lh
Cg4lh         ddco2 = CNC(20,9,1)*1.e6 - 468.35
Cg4lh
Cg4lh         DO 792 ik=1,117
Cg4lh         DO 792 ij=1,91
Cg4lh
Cg4lh            do 1589 iii=1,14
Cg4lh 1589            g4ss(iii) = xlhsens(iii,ij,ik)
Cg4lh 
Cg4lh            CALL LINTERP(timesf, 14, g4ss, ttime, 1, 0, ttout)
Cg4lh
Cg4lh            g4lhday(ij,ik) = ttout(1) + ddco2*xlhsens(15,ij,ik)
Cg4lh 792   CONTINUE





c   INTERPOLATE TOMS sfc albedo to current day TALBEDO(15,180) (1-14=month, 15=latitude) -> ALBDAY, ALBLAT(180)
c                                             g4ss(14), timesf(14) -> ttout(1) , also load latitude 
         DO 7709 IJ=1,180

             do 7712 iim=1,14
 7712           g4ss(iim) = TALBEDO(iim,ij)

             CALL LINTERP(timesf, 14, g4ss, ttime, 1, 0, ttout)

             ALBDAY(ij) = ttout(1)
             ALBLAT(ij) = TALBEDO(15,ij)
 7709    CONTINUE



c  INTERPOLATE coupled model TBAR/UBAR - NCEP differences to current day; use timesf(14) (-15 -> 375 days)
C    TDIFFR(14,45,76), UDIFFR(14,45,76) ->  TDIFFD(45,76), UDIFFD(45,76)
C        xlatf(45), zzff(76);   g4ss(14),  ttime(1), ttout(1)
C

         DO 790 ik=1,76
         DO 790 ij=1,45

            do 791 iii=1,14
 791            g4ss(iii) = tdiffr(iii,ij,ik)
 
            CALL LINTERP(timesf, 14, g4ss, ttime, 1, 0, ttout)

            tdiffd(ij,ik) = ttout(1)


            do 793 iii=1,14
 793            g4ss(iii) = udiffr(iii,ij,ik)

            CALL LINTERP(timesf, 14, g4ss, ttime, 1, 0, ttout)

            udiffd(ij,ik) = ttout(1)
 790   CONTINUE


C  load in NCEP DELF for baroclinic waves ONLY - for current day -> BARODAY(73,33) - not currently used
C       DFBARO(73,33,360) is in m/sec^2;   BALAT(73), BZZ(33) - need to define for call to TWODS

          DO 797 ik=1,33
              bzz(ik) = ik - 1.
          DO 797 ij=1,73
              balat(ij) = (ij-1)*2.5 - 90.
 797          baroday(ij,ik) = 0.0   !  dfbaro(ij,ik,iiday360)




            DO 8507 IK=1,21
            DO 8507 IJ=1,91
 8507          BBCDAY(IJ,IK) = BBCOUP(IJ,IK,IIDAY360)


C
C   Kzz adjustment for current day done in separate routine:
C

            DO 7517 IK=1,89
            DO 7517 IJ=1,45
 7517           kzz0(IJ,IK) = XKZZ9FR(IJ,IK,IIDAY360)


       CALL KZZADJ(IIDAY360,TIMESF, xlatkzz, zzkzz, KZZ0, KZZTH, KZZDAY)



C
C  load in fixed model Kyy:  XKYY9FR(46,88,360) -> KYYDAY(46,88), xlatkyy(46), zzkyy(88)
C  for use in coupled model dynamics,
C  then after calling coupled model dynamics, load KYYDAY(46,88) array into KYYO(L$C+1,Z$XC)
C  for return to chemistry (this keeps the original high resolution of Kyy)
C   -  DON'T do this when BLENDING fixed/coupled model Kyys
C
      DO 7617 IK=1,88
      DO 7617 IJ=1,46
 7617    KYYDAY(IJ,IK) = XKYY9FR(IJ,IK,IIDAY360)



C  call main coupled model driver:
C
c  NOTE: lower case are returned from TWODR - ie, output from coupled model dynamics:  
C  REAL*4 ypf(N$), zpf(M$), yppf(N$+1), zppf(M$+1)
c  kyyc(N$,M$), wsc(N$,M$), kzzc(N$+1,M$+1), tempc(N$+1,M$+1), ubarc(N$+1,M$+1), output to MAIN for diagnostics
C

       IF (iswfixed .eq. 0) CALL TWODR(DTMARK, L$C, Z$C, LATC, ZALT90C, 
     >     LHDAY, WLHDAY, HACKDAY, TEMPDAY, UBARDAY, TEMPG4DAY, EHFDAY,
     >     BBCDAY, XLATBBC, ZBBC, KZZTH, XLATKZZ, ZZKZZ,
     >                            KYYDAY, XLATKYY, ZZKYY, EPDAY, ZALTC,
     >     TDIFFD, UDIFFD, XLATF, ZZFF, BARODAY, BALAT, BZZ,
     >     ypf, zpf, yppf, zppf, tempc, ubarc, kyyc, kzzc, wsc, ttdiff)


C                             ! convert ubarc(N$+1,M$+1) from cm/sec to m/sec
          do 507 ik=1,M$+1
          do 507 ij=1,N$+1
 507         ubarc(ij,ik) = ubarc(ij,ik)/100.

C
C
C  interpolate coupled model dynamics arrays to chemistry grid as in SETDAILY
C     better to use BINTERP (linear interpolation) here instead of REGRID (polynomial interpolation)
C
C    also, BINTERP automatically takes care of interpolating grid values BEYOND the NATIVE grid 
C        - it just sets to a constant value using the last grid point available
C
c    => tbaro(L$C,Z$XC), ubaro(L$C,Z$XC), kyyo(L$C+1,Z$XC), kzzo(L$C,Z$XC+1), wbaro(L$C,Z$XC+1)
c         LATC(L$C);  ZALTC(Z$XC);    LATCE(L$C+1);  ZALTEC(Z$XC+1)
C                         tempc(N$+1,M$+1), ubarc(N$+1,M$+1), yppf(N$+1), zppf(M$+1)

           n11 = N$+1
           m11 = M$+1

           CALL BINTERP(yppf, N11, zppf, M11, TEMPC,
     >                  LATC, L$C, ZALTC, Z$XC, 0, 0, TBARO)

           CALL BINTERP(yppf, N11, zppf, M11, UBARC,
     >                  LATC, L$C, ZALTC, Z$XC, 0, 0, UBARO)



C  also interpolate ttdiff(N$,M$) tropospheric heating correction to proper grid for output
C      use LATC(L$C), ZALT90C(Z$C) --> ttdiffo(L$C,Z$C) is in K/sec (temperature heating)

           CALL BINTERP(ypf, N$, zpf, M$, TTDIFF,
     >                  LATC, L$C, ZALT90C, Z$C, 0, 0, TTDIFFO)


C
C  for output purposes (visually), adjust high latitude UBAR (>85N,S) as done in TWODLIB
C    do this if interpolation is beyond available points;  ubarc(N$+1,M$+1) ;  ubaro(L$C,Z$XC)
C       ypf(N$);   LATC(L$C)
C
Ccpw17  - DONT NEED TO DO THIS w/ UBARC now from original extended dynamics grid
ccpw17
ccpw17      do 607 ij=L$C,1,-1
ccpw17        if (latc(ij) .lt. ypf(1)) then
ccpw17          cfac = COSD(latc(ij))/COSD(latc(ij+1))
ccpw17
ccpw17          do 617 ik=1,Z$XC
ccpw17 617       ubaro(ij,ik) = ubaro(ij+1,ik)*cfac
ccpw17        endif
ccpw17 607  CONTINUE
ccpw17
ccpw17      do 707 ij=1,L$C
ccpw17        if (latc(ij) .gt. ypf(N$)) then
ccpw17          cfac = COSD(latc(ij))/COSD(latc(ij-1))
ccpw17
ccpw17          do 717 ik=1,Z$XC
ccpw17 717       ubaro(ij,ik) = ubaro(ij-1,ik)*cfac
ccpw17        endif
ccpw17 707  CONTINUE



      LC1 = L$C+1
           CALL BINTERP(ypf, N$, zpf, M$, KYYC,
     >                  LATCE, LC1, ZALTC, Z$XC, 0, 0, KYYO)

      ZC1 = Z$XC+1
           CALL BINTERP(ypf, N$, zpf, M$, WSC,
     >                  LATC, L$C, ZALTEC, ZC1, 0, 0, WBARO)


C
c  for KZZ use the original resolution KZZDAY(45,89) array loaded above here -  xlatkzz(45), zzkzz(89)
C   - DON'T use KZZC which has  been interpolated to a lower resolution
C
C
           CALL BINTERP(xlatkzz, 45, zzkzz, 89, KZZDAY,
     >                  LATC, L$C, ZALTEC, ZC1, 1, 1, KZZO)

c
cccc      n11 = N$+1
cccc      m11 = M$+1
cccc           CALL BINTERP(yppf, n11, zppf, m11, KZZC,
cccc     >                  LATC, L$C, ZALTEC, ZC1, 1, 1, KZZO)




c  reset limits on Kyy, Kzz:  kyyo(L$C+1,Z$XC) ; kzzo(L$C,Z$XC+1),  (in cm2/sec), also done in GETMOM

          kyyminc = 1.e8
          kyymaxc = 1.e11

          do 501 ik=1,Z$XC
          do 501 ij=1,L$C+1
             if (kyyo(ij,ik) .le. kyyminc) kyyo(ij,ik) = kyyminc
             if (kyyo(ij,ik) .ge. kyymaxc) kyyo(ij,ik) = kyymaxc
 501	  CONTINUE


          kzzminc = .001*1.e4
          kzzmaxc = 100.*1.e4

        do 504 ik=1,Z$XC+1
        do 504 ij=1,L$C
           if (kzzo(ij,ik) .le. kzzminc) kzzo(ij,ik) = kzzminc
           if (kzzo(ij,ik) .ge. kzzmaxc) kzzo(ij,ik) = kzzmaxc
 504	CONTINUE

 

c  ramp down w* : coupled model goes 1-2 km -> 93.38 km, so ramp down over 96-114 km, set 116 km = 0.
C  set 0 km = 0.,   if lowest coupled model level is > the 2nd level of WBARO by .2 km,
C  then ramp down the 2nd level, otherwise use as is;  wbaro(L$C,Z$XC+1), ZALTEC(Z$XC+1), iz96=79 (currently)

          iz96 = Z$XC+1-10

          do 300 ij=1,L$C

             do 302 ik=IZ96,Z$XC
      wbaro(ij,ik)=(0.-wbaro(ij,IZ96-1))/(zaltec(Z$XC+1)-zaltec(IZ96-1))
     >                * (zaltec(ik) - zaltec(IZ96-1)) + wbaro(ij,IZ96-1)
 302	     CONTINUE

             wbaro(ij,Z$XC+1) = 0.0
             wbaro(ij,1) = 0.0
             if (zpf(1) .gt. zaltec(2)+.2) wbaro(ij,2) = wbaro(ij,3)/2.
 300      CONTINUE

     
C  load in HEATOUT and COOLOUT arrays for current day for output back to MAIN, convert to K/day 
C  HEATCHM, HEATCHML(NS$,MS$), COOLCHM(NS$,MS$) are in COMMON -> HEATOUT, COOLOUT(L$C=NS$,Z$C=MS$)
C
C  NOTE: HEATCHM here ise ONLY for the IR bands from the GEOS5 SORAD.f
C  ADD in diurnal avg UV+VIS bands from the photolysis calc - PHEAT(5,L$C,Z$C) in K/sec, 5=TOTAL heating
C
C  also load in DRAG terms DRAGCHM(7,NS$,MS$) in COMMON -> DRAGOUT(7,L$C,Z$C), convert from cm/sec^2 to m/sec/day
C
C  DRAGCHM(1) is QBARY in 1/(cm-sec), convert to 1/(m-sec), multiply by 1.e11 for output
C                                            HEATCHML already in K/day in TWODS
C
        do 177 ik=1,Z$C
        do 177 ij=1,L$C
        HEATOUT(1,ij,ik) = HEATCHM(ij,ik)*86400. + PHEAT(5,ij,ik)*86400.
           HEATOUT(2,ij,ik) = HEATCHML(ij,ik)   ! *86400.
           COOLOUT(ij,ik) = COOLCHM(ij,ik)*86400.

           do 277 iig=2,6
 277          DRAGOUT(iig,ij,ik) = DRAGCHM(iig,ij,ik)/100.*86400.

           DRAGOUT(1,ij,ik) = DRAGCHM(1,ij,ik)*100.*1.E11

                                        ! load in tropospheric temp correction for output;  ttdiffo is in K/sec
           DRAGOUT(6,ij,ik) = ttdiffo(ij,ik)*86400.
 177    CONTINUE



C
C  interpolate COUPLED MODEL LONGITUDINAL TEMPERATURE field 
C     to chemistry grid for current day, for PDF use in SETDAILY; 
C     txx(72,N$,M$) is from JPWAVE;  txx0(N$,M$), ccx(L$C,Z$C) -> tcc(72,L$C,Z$C)
C
        do 1200 ix=1,72        
  
           do 1205 ik=1,M$
           do 1205 ij=1,N$
 1205	      txx0(ij,ik) = txx(ix,ij,ik)

           CALL BINTERP(ypf, N$, zpf, M$, TXX0,
     >                  LATC, L$C, ZALT90C, Z$C, 0, 0, ccx)

           do 1225 ik=1,Z$C
           do 1225 ij=1,L$C
 1225	      tcc(ix,ij,ik) = ccx(ij,ik)

 1200    CONTINUE



c
c   for PW35, the lowest level of the coupled model temperatures, tempc(N$+1,M$+1) 
C     ie, ik=1 = 1.015 km, are just extrapolated to the 0.5 km level of the chemistry grid
C     in the BINTERP routine above, SO THIS IS WHAT'S USED IN THE OUTPUT
C     no additional calculations are necessary...... (ie, below)
C


C
C  for coupled model, interpolate TEMPG4DAY(91,3), (ie, surface temps based on GEOS 4/sensitivity to CO2 from above)
C  to chemistry grid here where levels 1-3 are GEOS 4 sfc, 1km, 2km temps for current day.  - latg4(91), zzg4(3)
C
C  then load TEMPDAYC(L$C,Z$XC) -> tbaro(L$C,Z$XC) wherever chemistry grid is 
C           BELOW bottom level of coupled model,  zppf(M$+1),   ZALTC(Z$XC), eg, 0.5 km
C           NEED to do this since there's a big (~5-10K) difference between 0.5 and 1km

      do 888 ik=1,3
 888       zzg4(ik) = ik-1.
    
           CALL BINTERP(latg4, 91, zzg4, 3, TEMPG4DAY,
     >                  LATC, L$C, ZALTC, Z$XC, 0, 0, TEMPDAYC)

        do 577 ik=1,Z$XC
           if (zaltc(ik) .lt. zppf(1) ) then
              do 477 ij=1,L$C
 477             tbaro(ij,ik) = TEMPDAYC(ij,ik)
           endif
 577    CONTINUE



C
C  interpolate lowest coupled model level to chemistry latitude grid for output:
C    tempc(N$+1,M$+1), yppf(N$+1), ttz(N$+1) => tbaro(L$C,Z$XC), LATC(L$C); ttch(L$C)
C    then load into lowest chemistry grid point (eg, 1.015 km -> 0.5 km)

cc35         do 675 ij=1,N$+1
cc35 675        ttz(ij) = tempc(ij,1)
cc35
cc35           CALL LINTERP(yppf, N11, ttz, LATC, L$C, 0, TTCH)
cc35
cc35         do 677 ij=1,L$C
cc35 677        tbaro(ij,1) = ttch(ij)
cc35

        RETURN
	END

!        function cosd(dgr_argument)
!        real(4) cosd
!        real(4) dgr_argument
!        real(16), parameter ::
!     &     quadpi = 3.141592653589793238462643383279502884197Q0
!        real(16), parameter :: dgr_to_rad = (quadpi/180Q0)
!        cosd = cos(dgr_to_rad * dgr_argument)
!        end function
!
!        function sind(dgr_argument)
!        real(4) sind
!        real(4) dgr_argument
!        real(16), parameter ::
!     &     quadpi = 3.141592653589793238462643383279502884197Q0
!        real(16), parameter :: dgr_to_rad = (quadpi/180Q0)
!        sind = sin(dgr_to_rad * dgr_argument)
!        end function

