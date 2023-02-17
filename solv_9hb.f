C
C  /misc/kah02/9g> em solv_9gl.f &
C  [1] 7177
C

        SUBROUTINE SOLVER(DTIMEB, ISNTIME)

C
C   ROUTINE TO SOLVE FOR ALL SLOW SPECIES (TRANSPORTED) CONCENTRATIONS AT POINT(IJ,IK)
C   AFTER A TIME STEP. INCLUDES LOGIC TO SWITCH ON BOUNDARIES.
C
C   Now uses values of diurnally averaged radicals and Js (from AER code)
C      and transported slow species, all from current time step                - EF  (2/24/03)
C      
C      NO MORE NCB!!!
C
C  NOTE: THe C array for the slow species here has been transported, so update 
C        the C array with chemistry to become the CN array
C
C   
C   FOR J and K reactions which are products of 2 or more diurnally varying things, we need to 
C      sum up over the 24-hour cycle first, and then average - use JDC(PH$,18,L$,Z$) and CNDC(S$,18+3,L$,Z$)
C
C      and also ALL 7 HET reactions are now DIURNALLY VARYING - KH(RH$,18,L$,Z$) 
C
C
C   TIME DEPENDENT BCs + outputs for Lifetime calcs (9HA) - Jan 2012
C
C
C

        include "com2d.h"



C   common for GMI sfc deposition, emissions, interpolated to L$ latitude grid (from TEMPIN):
C       for CURRENT DAY -- DEPDAY in 1/sec;  EMISDAY in #/cm3/sec
C
C   DEPOSITION: 1=CH2O  2=HNO3  3=H2O2  4=CH3OOH (MP)  5=N2O5  6=NO2   7=O3
C
C   EMISSIONS : 1=CH2O  2=CO    3=NO (NOx)
C

        COMMON/CDEP/tdep(14), DEP0(14,L$,7), EMIS0(14,L$,3), 
     >              DEPDAY(L$,7), EMISDAY(L$,3)



C  common for CGCM-COMBO T1R8 sfc deposition + wet scavenging
C
C  SCAVDAY(20,L$,Z$) is in 1/sec : the combined surface dep + scavenging for current day
C          SCAVDAY(20,L$,Z$) is in 1/sec:
C
C  SCAVDAY:  1=ch2o  2=hno2     3=hno3   4=ho2no2  5=ho2   6=h2o2  7=mp    8=no2  9=no3  10=n2o5
C           11=o3   12=brono2  13=hbr   14=br     15=brcl 16=hobr 17=clo  18=clono2  19=hcl  20=hocl
C
        COMMON/CDEPR8/DEPR8(10,14,L$), SCAVR8(20,14,L$,Z$),
     >                SCAVDAY(20,L$,Z$)


C  wavelength dependent J-arrays (for diagnostics ONLY), from SPDR and PDFINSOL, TIMEAVJ

        COMMON/CJX39/JX39(PH$,5,L$,Z$), JZ39(18,PH$,5,L$,Z$),
     >                J39(PH$,5,L$,Z$)


C   COMMON for SPARC OH field in (#/cm3)
C
        COMMON/CCOH/ohm(14,L$,Z$), ohday(L$,Z$)



        INTEGER ISNTIME

        DOUBLE PRECISION DTIMEB(L$,ISNTIME)
        REAL*8 jsum1, jsum2, totna
        REAL*8 ksum1, ksum2, ksum3, ksum4, ksum5, ksum6, ksum7, ksum8
        REAL*8 ksum9, ksum0, ksum59, ksum69, ksum10, ksum11
        REAL*8 ksumh(14)

        REAL*8 fnrat, fcrat, fbrat, J39, JX39, JZ39
        REAL*8 noyprod, noyloss, clyprd, clyl, bryprd, bryl, hfp, hfl
        REAL*8 cf2op, cf2ol, h2op, h2ol, h2prod, h2loss, coprod, coloss
        REAL*8 co2prod, co2loss, ch4prod, ch4loss, pn2o, n2oloss
        REAL*8 f11p, f11l, f12p, f12l, ccl4p, ccl4l, ch3clp, ch3cl1
        REAL*8 CH3CCL3P, CH3CCL3L, ch3brp, ch3brl, cbrf3l, cbrclf2l
        REAL*8 chclf2l, c2cl3F3l, c2cl2f4l, c2clf5l, c141l, c142l, c123l
        REAL*8 c2402l, c1202l, cmebrl, cbrmfl, BCSH, BCNH

        REAL*8 xspecn, xspecc, xspecb, NOX0, clx0, brx0
        REAL*8 noyo, xspecno, noxo, clyo, xspecco, clxo, n2o5dep
        REAL*8 bryo, xspecbo, brxo, coemis, noxemis, no2dep, hno3dep
        REAL*8 tcly, tbry, clyrat, bryrat

        REAL*4 xlifec(5,L$,Z$)
C
C  xlife(30,6,L$,Z$), xlifej(30,5,L$,Z$) are in COMMON
C
C  xlife:  1st index:  
C
C   1=CFC-11;  2=CFC-12;  3=Carb Tet;  4=Meth Chlor;   5=HCFC-22;  6=N2O;  7=CH4
C   8=Hal 1211  9=Hal 1301  10=CFC113   11=CFC115   12=HFC-134a  13=HFC-143a  14=HFC-23
C  15=CFC114  16=HCFC-141b   17=HCFC-142b   18=CH3Cl  19=CH3Br   20=Hal1202   21=Hal2402
C  22=HFC-32   23=HFC-125    24=HFC-152a   25=HFC-227ea   26=HFC-245fa
C
C   2nd index (loss):   1=total loss;  2=J (tot);  3=O1D;   4=OH;   5=Cl;   6=number density
C
C
C  xlifej - J's for 5 wavelength bins: J39(PH$,5,L$,Z$)
C     1=Lyman alpha;   2=1695-1905A;   3=1905-2299A;  4=2299-2857A;   5=2857-7350A
C
C
C  xlifec - extra constituents:
C     1=O1D     2=OH     3=Cl     4=M     5=Temperature
C

C  arrays for setting minimum surface BC:  IBCMIN = 0/1 - for well behaved lifetimes

        INTEGER ISBC(25)
        DATA ISBC/34, 35, 36, 37, 40, 49, 50, 51, 52, 53, 54, 55, 69,
     >            70, 71, 72, 75, 81, 82, 83, 84, 85, 86, 87, 88/

        DATA IBCMIN/1/


        SAVE



C order:   for 9GL, HNO3 solved as a FAST SPECIES in AER CHEMISTRY, also update FROM TRANSPORTED SPECIES 
C
C families: 
C  1 Ox (no chemistry)
C  2 CHx (no chemistry)
C  3 HOx (no chemistry)
C  4 NOy
C  5 NOz = NOy - SHNO3
c  6 Cly
c  7 Bry
c  8 HF (Fx)

C  reservoirs 
C  9 CF2O
C 10 H2O
C 11 H2
C 12 CO
c 13 CO2

C  source gases (non - CFCs)
C 14 CH4
C 15 N2O
          
C CFCs (19 total)
C 16-34

C other tracers (age, C-14, etc.)

c

C  ***********************************************************************************************


       DO 1000 ik=1,Z$
       DO 1000 ij=1,L$

C
C
C           
C  DEPOSTION COMPUTED BY DIVIDING .1 CM/SEC BY 2X10+5 CM (SIKE OF GRID)  1/5/88
C  OXDEP(L$,360) from HARVARD model are daily values in m/sec, convert to 1/sec using depth of grid box
C
C
C  Sept. 2006, set to ozone climatology for 0.5-1.5 km:
C       IF (zalt(ik) .le. 2.) c(4,IJ,ik) = OZBC(iday360,ij,ik)*1.E-6*M(ij,ik)
C
C
C  DELTAZ(L$,Z$X) in KM (in COMMON)
c      IF (IK .eq. 1) then
c         deploss = OXDEP(ij,iday360)/(DELTAZ(ij,1)*1000.)
c         cn(39,IJ,1) = cn(39,IJ,1)/(1. + deploss*dt)
c      ENDIF
C
C
CC WCONV200, 277
C  surface deposition and lower tropospheric production of ozone now included in JACOB - March/July 2011
C

C  Ox  -  C(39)  (O + O(1D) + O3)

         C(39,ij,ik) = c(1,ij,ik) + c(2,ij,ik) + c(4,ij,ik) 
         CN(39,ij,ik) = C(39,ij,ik) 

         CN(1,ij,ik) = C(1,ij,ik)
         CN(2,ij,ik) = C(2,ij,ik)
         CN(4,ij,ik) = C(4,ij,ik)
         CN(41,ij,ik) = C(41,ij,ik)



C  **************************************************************************************
C
CWCONV200
C  CHx   -  CN(73)  (CH3O2 + CH2O + CH3OOH) - update from transported value, NO CHEMISTRY
C     CH2O sfc emissions NOW done in JACOB (9HA - Sept 2011)

          CN(22,ij,ik) = C(22,ij,ik)
          CN(23,ij,ik) = C(23,ij,ik)
          CN(24,ij,ik) = C(24,ij,ik)

          CN(73,ij,ik) = cn(22,ij,ik) + cn(23,ij,ik) + cn(24,ij,ik) 
          C(73,ij,ik) = CN(73,ij,ik)



C  ****************************************************************************************
c
c
C  HOx - C(74)  (H + OH + HO2 + 2*H2O2) - update from transported value, NO CHEMISTRY
C
C
C  TO RE-SET tropospheric OH (which has been transported) with TRANSCOM climatology:
C    (IFIXOH = 1, set in CONTROL.dat)
C
C    OHDAY(L$,Z$) is in #/cm3 for current day;     ITROP360(L$,360) (in COMMON) are the
C    tropopause FORTRAN indicies for the current Z$ grid, defined in TROPKZZ
C
C   get tropopause index for current latitude,  do smooth blend around tropopause

          iktrop = ITROP360(ij,iday360)
          if (ABS(LAT4(ij)) .ge. 27.  .and.  iktrop .gt. 14) iktrop = 14
          if (ABS(LAT4(ij)) .ge. 30.  .and.  iktrop .gt. 13) iktrop = 13
          if (ABS(LAT4(ij)) .ge. 35.  .and.  iktrop .gt. 12) iktrop = 12
          if (ABS(LAT4(ij)) .ge. 65.) iktrop = 11
ccccccccc       if (ABS(LAT4(ij)) .ge. 75.  .and.  iktrop .gt. 10) iktrop = 10
ccccccccc       if (ABS(LAT4(ij)) .ge. 80.  .and.  iktrop .gt.  9) iktrop = 9

      IF (IFIXOH .eq. 1) then
       if (ik .eq. iktrop) c(13,ij,ik) = .5*(OHDAY(ij,ik) + c(13,ij,ik))
       if (ik .lt. iktrop) c(13,ij,ik) = OHDAY(ij,ik)
      ENDIF


          C(74,ij,ik) = c(12,ij,ik) + c(13,ij,ik) + c(14,ij,ik) 
     >                + 2.*c(16,ij,ik)
          CN(74,ij,ik) = C(74,ij,ik) 

          CN(12,ij,ik) = C(12,ij,ik) 
          CN(13,ij,ik) = C(13,ij,ik) 
          CN(14,ij,ik) = C(14,ij,ik) 
          CN(16,ij,ik) = C(16,ij,ik) 



C  **************************************************************************************
C
c
c   also update solid H2O from transported value, NO CHEMISTRY
C        
          CN(67,ij,ik) = C(67,ij,ik) 



C********************************************************************************************
C
C
C  also update the CN array for the x-species from transported values only 
C  (these are NOT changed by the family updates):   ClONO2, ClNO2, BrONO2, BrNO2, BrCl
C

          CN(30,ij,ik) = C(30,ij,ik) 
          CN(63,ij,ik) = C(63,ij,ik) 
          CN(47,ij,ik) = C(47,ij,ik) 
          CN(79,ij,ik) = C(79,ij,ik) 
          CN(60,ij,ik) = C(60,ij,ik) 



C
C  ****************************************************************************************
C
C    do all het reactions here:  some are NOT used,  KH(RH$,18,L$,Z$) - now DIURNALLY VARYING (REAL*8)
C    
C    
      do 214 ikh=1,14
 214     ksumh(ikh) = 0.D0

      do 215 it=1,isntime-1
        ksumh(1) = ksumh(1)
     >    + 0.5D0*(KH(1,it,ij,ik)*CNDC(30,it,ij,ik) 
     >    +    KH(1,it+1,ij,ik)*CNDC(30,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh(2) = ksumh(2) + 0.5D0*(KH(2,it,ij,ik)*CNDC(30,it,ij,ik) + 
     >                         KH(2,it+1,ij,ik)*CNDC(30,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh(3) = ksumh(3) + 0.5D0*(KH(3,it,ij,ik)*CNDC(8,it,ij,ik) + 
     >                         KH(3,it+1,ij,ik)*CNDC(8,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh(4) = ksumh(4)
     >    + 0.5D0*(KH(4,it,ij,ik)*CNDC(8,it,ij,ik) 
     >    +    KH(4,it+1,ij,ik)*CNDC(8,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh(5) = ksumh(5)
     >    + 0.5D0*(KH(5,it,ij,ik)*CNDC(25,it,ij,ik)                   
     >    +    KH(5,it+1,ij,ik)*CNDC(25,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh(6) = ksumh(6) + 0.5D0*(KH(6,it,ij,ik)*CNDC(47,it,ij,ik) + 
     >                         KH(6,it+1,ij,ik)*CNDC(47,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh(7) = ksumh(7)
     >    + 0.5D0*(KH(7,it,ij,ik)*CNDC(68,it,ij,ik)
     >    +    KH(7,it+1,ij,ik)*CNDC(68,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh(8) = ksumh(8)
     >    + 0.5D0*(KH(8,it,ij,ik)*CNDC(7,it,ij,ik)
     >    +        KH(8,it+1,ij,ik)*CNDC(7,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh(9) = ksumh(9)
     >    + 0.5D0*(KH(9,it,ij,ik)*CNDC(64,it,ij,ik)
     >    +        KH(9,it+1,ij,ik)*CNDC(64,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh(10) = ksumh(10)
     >    + 0.5D0*(KH(10,it,ij,ik)*CNDC(8,it,ij,ik)
     >    +        KH(10,it+1,ij,ik)*CNDC(8,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh(11) = ksumh(11)
     >    + 0.5D0*(KH(11,it,ij,ik)*CNDC(25,it,ij,ik)
     >    +        KH(11,it+1,ij,ik)*CNDC(25,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh(12) = ksumh(12)
     >    + 0.5D0*(KH(12,it,ij,ik)*CNDC(30,it,ij,ik)
     >    +        KH(12,it+1,ij,ik)*CNDC(30,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh(13) = ksumh(13)
     >    + 0.5D0*(KH(13,it,ij,ik)*CNDC(68,it,ij,ik)
     >    +        KH(13,it+1,ij,ik)*CNDC(68,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksumh(14) = ksumh(14)
     >    + 0.5D0*(KH(14,it,ij,ik)*CNDC(47,it,ij,ik)
     >    +        KH(14,it+1,ij,ik)*CNDC(47,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0
 215  CONTINUE



C  *************************************************************************************
C
C
C   first sum up NOx species from transported values, use C ARRAY, 
C      these will be updated with updates to NOy, include HOONO
C      include SOLID HNO3 ;  DO NOT include x-species (ClONO2, BrONO2, etc) here
C
C      sum up x-species separately (ClONO2, BrONO2, ClNO2, BrNO2)
C
C
       NOX0 =  c(5,ij,ik) + c(6,ij,ik) + c(7,ij,ik)
     >    + 2.*c(8,ij,ik) + c(9,ij,ik) + c(10,ij,ik) + c(66,ij,ik)
     >    +   c(38,ij,ik) + c(42,ij,ik) + c(89,ij,ik)


       xspecn= c(30,ij,ik) + c(63,ij,ik) + c(47,ij,ik) + c(79,ij,ik)


C
C   NOy - CN(31) - first sum up family from transported values, use C ARRAY 
C     - this is the total ODD NITROGEN in model, 
C       so INCLUDE ALL ODD NITROGEN SPECIES here (including SHNO3)
C

       C(31,ij,ik) = c(5,ij,ik) + c(6,ij,ik) + c(7,ij,ik) 
     >          + 2.*c(8,ij,ik) + c(9,ij,ik) + c(10,ij,ik) + c(66,ij,ik)
     >            + c(38,ij,ik) + c(42,ij,ik) + c(89,ij,ik)
     >           + c(30,ij,ik) + c(63,ij,ik) + c(47,ij,ik) + c(79,ij,ik) 


C                     account for diurnally varying radicals
      ksum59 = 0.D0  
      ksum69 = 0.D0   
      do 218 it=1,isntime-1
        ksum59 = ksum59 + 0.5D0*(CNDC(5,it,ij,ik)*CNDC(9,it,ij,ik) +
     >                           CNDC(5,it+1,ij,ik)*CNDC(9,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

        ksum69 = ksum69 + 0.5D0*(CNDC(6,it,ij,ik)*CNDC(9,it,ij,ik) +
     >                           CNDC(6,it+1,ij,ik)*CNDC(9,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0
 218  CONTINUE


C ----------------------------------------------------------------------


C ADD LIGHTNING PRODUCTION OF NOX  12/15/88  H2 SOURCE OF KO ET AL. 1986,  
C    updated for new high res model -  ZALT(Z$X) add LIGHTNING source at and below 16 km, 25S-25N, LAT(L$)
C    also, ADD GCR PRODUCTION OF NOX  12/13/88,   NOy LOSS also include HNO3 rainout
C
CGMI       IF (IK .GE. 2  .AND.  zalt(IK) .LE. 16.) THEN
CGMI          IF (lat(ij) .GE. -25.  .AND.  lat(ij) .LE. 25.) PLIGHTN = 1.E3
CGMI       ENDIF

CGMI - use GMI lightning here - DLIGHT(L$,Z$) (in COMMON) interpolated to current day in BOUNDC
C  this is production of N  in #/cm^3/sec;  need to reduce to agree w/ GMI NOy -  (EF, July, 2009)
C
       PLIGHTN = 0.0
       if (zalt(ik) .le. 20.) then
          PLIGHTN = DLIGHT(ij,ik)/4.
       endif
  
! AYelland 10Oct2020 - introduced "ionsflag_solv" & "lightnflag_solv" to alter calculation in run

C       print *, "SOLVER: ionsflag_solv = ", ionsflag_solv
C       print *, "SOLVER: lightnflag_solv = ", lightnflag_solv

C         Checking ionSource Read-in:  [UNITS: ((ions)/(cm^3*s))]
C           print *, "SOLVER: ionSource = "
C           print *, ionSource

C         Checking GCRions Read-in:  [UNITS: ((ions)/(cm^3*s))]
C           print *, "SOLVER: GCRions = "
C           print *, GCRions

C         Checking ionsRatio Read-in
C           print *, "SOLVER: ionsRatio = "
C           print *, ionsRatio

      ! We use "ionsRatio" ["ionSource" (max SNCR ions) / "GCR" (background ions)] as the multiplier for the lightning effect
      if (ionsflag_solv .eq. 1) then

        if (lightnflag_solv .eq. 1) then
            PLIGHTN = ionsRatio(IJ,IK)*PLIGHTN
         end if
!!%!BThomas 8Nov2022 - restricting ion input to altitude above 70 km
!!%         if (ik.lt.65) then
!!%            ionSource(ij,ik) = 0.0
!!%            print *, "0 ions below ik=65"
!!%            print *, "ik, press(ik), zalt90arr(ik): "
!!%            print *, "   " ik, press(ik), zalt90arr(ik)
!!%         end if
!!%!     done as a test of O3 effect of ionization only input at high altitude
!!%!end BThomas 8Nov2022
        noyprod = k(45,ij,ik)*c(2,ij,ik)*c(11,ij,ik)*2.
     >            + 1.25*GCR(ij,ik) + 1.25*ionSource(ij,ik) + PLIGHTN
        !BThomas 24Jan2019 - adding ionSource here, analogous to GCR - see ionSourceIn.f

      else

        noyprod = k(45,ij,ik)*c(2,ij,ik)*c(11,ij,ik)*2.
     >            + 1.25*GCR(ij,ik) + PLIGHTN
        !AYelland 10Oct2020 - If ionization isn't being factored into equation (ionSource), lightning is not included

      endif

C ----------------------------------------------------------------------


C  HNO3 enhanced washout - HNO3W(L$,Z$) is in 1/sec, in COMMON
C
C  now use GMI COMBO sfc deposition + scavenging: SCAVDAY(20,L$,Z$) is in 1/sec (June 2011)
C                                    include HNO3 and HOONO isomer - C(89)

       noyloss = (K(47,ij,ik)*ksum59*2.D0
     >         +  K(81,ij,ik)*ksum69*2.D0 

     >         +  SCAVDAY(2,ij,ik)*c(42,ij,ik)
     >         +  SCAVDAY(3,ij,ik)*(c(10,ij,ik) + c(89,ij,ik))
     >         +  SCAVDAY(4,ij,ik)*c(38,ij,ik)
     >         +  SCAVDAY(8,ij,ik)*c(6,ij,ik)
     >         +  SCAVDAY(9,ij,ik)*c(7,ij,ik)
     >         +  SCAVDAY(10,ij,ik)*c(8,ij,ik)
     >         +  SCAVDAY(12,ij,ik)*c(47,ij,ik)
     >         +  SCAVDAY(18,ij,ik)*c(30,ij,ik) )/c(31,ij,ik)


       CN(31,ij,ik) = (C(31,ij,ik) + noyprod*DT)/(1. + noyloss*DT) 


C   use GMI BC for NOy here;  BVALG(S$,L$,10) is in COMMON (ppp); ik: 0.5 - 9.5 km
CWCONV200  - do individual species BCs below
CGMI:
C       IF (zalt(ik) .le. 5.) CN(31,ij,ik) = BVALG(31,ij,ik)*M(ij,ik)
C


       CTPROD(2,IJ,IK) = noyprod/m(ij,ik)       
       CTLOSS(2,IJ,IK) = noyloss*c(31,ij,ik)/m(ij,ik)



C  ********************************************
C
C
C  NOW update NOx radicals and SHNO3 using new NOy family - update BOTH the C and CN arrays
C       this puts all the NOy increase into the NOx species,  cross-species are NOT changed
C
C   NOTE: FOR THE CURRENT GRID POINT: if all the NOy is in the XSPECS, ie, fnrat LE 0, then set fnrat=1, 
C         and the NOy family and the NOx radicals are reset to the transported values (no chemical updates)
C
        fnrat = (CN(31,ij,ik) - xspecn)/NOX0

        if (fnrat .le. 0.) then
          noyo    = CN(31,ij,ik)/M(IJ,IK)*1.d9
          xspecno = xspecn/M(IJ,IK)*1.d9
          noxo    = nox0/M(IJ,IK)*1.d9

          write(1581,221) iyr, iday360, ij, ik,
     >            fnrat, noyo, xspecno, noxo
 221      format('IN SOLVER, NOy-fnrat le 0 : ', 4I4, 2x, 1P4D13.4)

          fnrat = 1.d0
          CN(31,ij,ik) = C(31,ij,ik)
        endif


        CN(5,ij,ik) = C(5,ij,ik)*fnrat
        CN(6,ij,ik) = C(6,ij,ik)*fnrat
        CN(7,ij,ik) = C(7,ij,ik)*fnrat 
        CN(8,ij,ik) = C(8,ij,ik)*fnrat
        CN(9,ij,ik) = C(9,ij,ik)*fnrat 
        CN(10,ij,ik) = C(10,ij,ik)*fnrat 
        CN(89,ij,ik) = C(89,ij,ik)*fnrat 
        CN(38,ij,ik) = C(38,ij,ik)*fnrat 
        CN(42,ij,ik) = C(42,ij,ik)*fnrat 
        CN(66,ij,ik) = C(66,ij,ik)*fnrat 



        C(5,ij,ik) = CN(5,ij,ik)
        C(6,ij,ik) = CN(6,ij,ik)
        C(7,ij,ik) = CN(7,ij,ik)
        C(8,ij,ik) = CN(8,ij,ik)
        C(9,ij,ik) = CN(9,ij,ik)
        C(10,ij,ik) = CN(10,ij,ik)
        C(89,ij,ik) = CN(89,ij,ik)
        C(38,ij,ik) = CN(38,ij,ik)
        C(42,ij,ik) = CN(42,ij,ik)
        C(66,ij,ik) = CN(66,ij,ik)


CWCONV200 - NOT DONE
C  re do NOy surface value, do HNO3, NO, NO2, N2O5 separately, then re-sum total NOy family
C
C  NOx surface emissions for current day - EMISDAY(L$,3) is in #/cm3/sec;  1=CH2O  2=CO  3=NO (NOx)
C
C  surface deposition loss - DEPDAY(L$,7) is in 1/sec
C  1=CH2O  2=HNO3  3=H2O2  4=CH3OOH (MP)  5=N2O5  6=NO2   7=O3
C
c       no2dep = 0.d0
c       hno3dep = 0.d0
c       n2o5dep = 0.d0
c       noxemis = 0.d0
c
c       IF (ik .eq. 1) THEN
c           hno3dep = DEPDAY(ij,2)
c           n2o5dep = DEPDAY(ij,5)
c           no2dep  = DEPDAY(ij,6)
c
c           noxemis = EMISDAY(ij,3)
c           ynrat = 1.  ! amount of NOx emis in NO2
c
C NO
c           CN(5,ij,1) = CN(5,ij,1) + (1.-ynrat)*noxemis*DT
C NO2
c           CN(6,ij,1) = (CN(6,ij,1) + ynrat*noxemis*DT)/(1. + no2dep*DT)
C N2O5
c           CN(8,ij,1) = CN(8,ij,1)/(1. + n2o5dep*DT)
C HNO3
c           CN(10,ij,1) = CN(10,ij,1)/(1. + hno3dep*DT)
C NOy
c           CN(31,ij,1) = cn(5,ij,1) + cn(6,ij,1) + cn(7,ij,1)
c     >         + 2.*cn(8,ij,1) + cn(9,ij,1) + cn(10,ij,1) + cn(66,ij,1)
c     >           + cn(38,ij,1) + cn(42,ij,1) + cn(89,ij,1)
c     >           + cn(30,ij,1) + cn(63,ij,1) + cn(47,ij,1) + cn(79,ij,1)
c
c           C(5,ij,1) = CN(5,ij,1)
c           C(6,ij,1) = CN(6,ij,1)
c           C(8,ij,1) = CN(8,ij,1)
c           C(10,ij,1) = CN(10,ij,1)
c           C(89,ij,1) = CN(89,ij,1)
c           C(31,ij,1) = CN(31,ij,1)
c       ENDIF


CWCONV200 - update NOx BCs w/ GMI below 5 km (HNO2, HNO3, HNO4, NO, NO2, NO3, N2O5), re-sum total NOy
C
       IF (zalt(ik) .le. 5.) then
           CN(5,ij,ik)  = BVALG(5,ij,ik)*M(ij,ik)
           CN(6,ij,ik)  = BVALG(6,ij,ik)*M(ij,ik)
           CN(7,ij,ik)  = BVALG(7,ij,ik)*M(ij,ik)
           CN(8,ij,ik)  = BVALG(8,ij,ik)*M(ij,ik)
           CN(10,ij,ik) = BVALG(10,ij,ik)*M(ij,ik)
           CN(38,ij,ik) = BVALG(38,ij,ik)*M(ij,ik)
           CN(42,ij,ik) = BVALG(42,ij,ik)*M(ij,ik)

           CN(31,ij,ik) = cn(5,ij,ik) + cn(6,ij,ik) + cn(7,ij,ik) 
     >      + 2.*cn(8,ij,ik) + cn(9,ij,ik) + cn(10,ij,ik) + cn(66,ij,ik)
     >       + cn(38,ij,ik) + cn(42,ij,ik) + cn(89,ij,ik)
     >       + cn(30,ij,ik) + cn(63,ij,ik) + cn(47,ij,ik) + cn(79,ij,ik) 

           C(5,ij,ik) = CN(5,ij,ik)
           C(6,ij,ik) = CN(6,ij,ik)
           C(7,ij,ik) = CN(7,ij,ik)
           C(8,ij,ik) = CN(8,ij,ik)
           C(10,ij,ik) = CN(10,ij,ik)
           C(89,ij,ik) = CN(89,ij,ik)
           C(38,ij,ik) = CN(38,ij,ik)
           C(42,ij,ik) = CN(42,ij,ik)
           C(31,ij,ik) = CN(31,ij,ik)
       ENDIF


C  ********************************************************************************
C
C
C   first sum up Clx species from transported values, use C ARRAY, these will be updated with updates to Cly
C     DO NOT include x-species (ClONO2, ClNO2, BrCl) here
C
C     sum up x-species separately
C
C
CC  tfam - first sum up individual transported members, then apply the ratio of the 
CC         transported family, C(32), to each member, to adjust for discontinuities
C
C       include BrCl in Bry below; OK to include ClONO2, ClNO2 in Cly
C

       tcly = c(25,ij,ik) + c(26,ij,ik) + c(27,ij,ik)  
     >      + c(28,ij,ik) + c(29,ij,ik) + c(30,ij,ik) 
     >      + 2.*c(61,ij,ik) + c(62,ij,ik) + c(63,ij,ik)
     >      + 2.*c(64,ij,ik) + c(65,ij,ik) 

       clyrat = c(32,ij,ik)/tcly

       c(25,ij,ik) = c(25,ij,ik)*clyrat
       c(26,ij,ik) = c(26,ij,ik)*clyrat
       c(27,ij,ik) = c(27,ij,ik)*clyrat
       c(28,ij,ik) = c(28,ij,ik)*clyrat
       c(29,ij,ik) = c(29,ij,ik)*clyrat
       c(30,ij,ik) = c(30,ij,ik)*clyrat

       c(61,ij,ik) = c(61,ij,ik)*clyrat
       c(62,ij,ik) = c(62,ij,ik)*clyrat
       c(63,ij,ik) = c(63,ij,ik)*clyrat
       c(64,ij,ik) = c(64,ij,ik)*clyrat
       c(65,ij,ik) = c(65,ij,ik)*clyrat

C                                 also update ClONO2, ClONO here, and CN(32):
       CN(30,ij,ik) = C(30,ij,ik) 
       CN(63,ij,ik) = C(63,ij,ik) 
       CN(32,ij,ik) = C(32,ij,ik) 

CC  end TFAM



       CLX0 = c(25,ij,ik) + c(26,ij,ik) + c(27,ij,ik) + c(28,ij,ik) 
     >      + c(29,ij,ik) + 2.*c(61,ij,ik) + c(62,ij,ik) 
     >                    + 2.*c(64,ij,ik) + c(65,ij,ik) 


       xspecc = c(30,ij,ik) + c(63,ij,ik) + c(60,ij,ik)



C  Cly - CN(33) - sum up TOTAL Cly from transported/adjusted radicals, use C ARRAY, include X-species


       C(33,ij,ik) = c(25,ij,ik) + c(26,ij,ik) + c(27,ij,ik)  
     >             + c(28,ij,ik) + c(29,ij,ik) + c(30,ij,ik) 
     >             + c(60,ij,ik) + 2.*c(61,ij,ik) + c(62,ij,ik) 
     >             + c(63,ij,ik) + 2.*c(64,ij,ik) + c(65,ij,ik) 


       clyprd = J(19,ij,ik)*c(36,ij,ik)*4. + J(20,ij,ik)*c(37,ij,ik) 
     >        + J(21,ij,ik)*c(34,ij,ik)*3. + J(22,ij,ik)*c(35,ij,ik)*2.
     >        + J(26,ij,ik)*c(40,ij,ik)*3. + J(31,ij,ik)*c(51,ij,ik)
     >        + J(32,ij,ik)*c(52,ij,ik)    + J(33,ij,ik)*c(53,ij,ik)*3.
     >        + J(34,ij,ik)*c(54,ij,ik)*2. + J(35,ij,ik)*c(55,ij,ik)  
     >        + J(65,IJ,IK)*C(69,IJ,IK)*2. + J(66,ij,ik)*c(70,ij,ik) 
     >        + J(67,IJ,IK)*C(71,IJ,IK)*2.  
     >        + K(16,ij,ik)*c(13,ij,ik)*c(37,ij,ik) 
     >        + K(54,ij,ik)*c(36,ij,ik)*c(2,ij,ik)*4.
     >        + K(75,ij,ik)*c(13,ij,ik)*c(40,ij,ik)*3. 
     >        + K(80,ij,ik)*c(2,ij,ik)*c(35,ij,ik)*2. 
     >        + K(83,ij,ik)*c(2,ij,ik)*c(34,ij,ik)*3.
     >        + K(98,ij,ik)*c(13,ij,ik)*c(52,ij,ik)
     >        + K(99,ij,ik)*c(2,ij,ik)*c(53,ij,ik)*3. 
     >        + K(100,ij,ik)*c(2,ij,ik)*c(54,ij,ik)*2.
     >        + K(101,ij,ik)*c(2,ij,ik)*c(55,ij,ik)
     >        + K(112,ij,ik)*c(2,ij,ik)*c(51,ij,ik)
     >        + K(114,ij,ik)*c(2,ij,ik)*c(69,ij,ik)*2.
     >        + K(115,ij,ik)*c(13,ij,ik)*c(69,ij,ik)*2.
     >        + K(116,ij,ik)*c(27,ij,ik)*c(69,ij,ik)*2.
     >        + K(117,ij,ik)*c(70,ij,ik)*c(2,ij,ik)
     >        + K(118,ij,ik)*c(70,ij,ik)*c(13,ij,ik)
     >        + K(119,ij,ik)*c(70,ij,ik)*c(27,ij,ik)
     >        + K(120,ij,ik)*c(71,ij,ik)*c(2,ij,ik)*2.
     >        + K(121,ij,ik)*c(71,ij,ik)*c(13,ij,ik)*2.
     >        + K(122,ij,ik)*c(71,ij,ik)*c(27,ij,ik)*2.
     >        + K(125,ij,ik)*c(52,ij,ik)*c(2,ij,ik)
     >        + K(132,IJ,IK)*c(27,IJ,IK)*c(37,IJ,IK)
     >        + K(136,IJ,IK)*c(2,IJ,IK)*c(37,IJ,IK)
     >        + K(181,IJ,IK)*c(52,ij,ik)*c(27,ij,ik)
     >        + K(182,IJ,IK)*c(40,ij,ik)*c(27,ij,ik)*3.
     >        + K(129,IJ,IK)*c(13,IJ,IK)*c(36,IJ,IK)*4.
     >        + K(130,IJ,IK)*c(13,IJ,IK)*c(34,IJ,IK)*3.
     >        + K(131,IJ,IK)*c(13,IJ,IK)*c(35,IJ,IK)*2.

     >        + K(229,IJ,IK)*c(40,ij,ik)*c(2,ij,ik)*3.
     >        + K(238,IJ,IK)*c(13,IJ,IK)*c(53,IJ,IK)*3.
     >        + K(239,IJ,IK)*c(13,IJ,IK)*c(54,IJ,IK)*2.
     >        + K(240,IJ,IK)*c(13,IJ,IK)*c(55,IJ,IK)
     >        + K(243,IJ,IK)*c(13,IJ,IK)*c(51,IJ,IK)
     >        + K(259,IJ,IK)*c(27,IJ,IK)*c(36,IJ,IK)*4.
     >        + K(260,IJ,IK)*c(27,IJ,IK)*c(34,IJ,IK)*3.
     >        + K(261,IJ,IK)*c(27,IJ,IK)*c(35,IJ,IK)*2.
     >        + K(262,IJ,IK)*c(27,IJ,IK)*c(53,IJ,IK)*3.
     >        + K(263,IJ,IK)*c(27,IJ,IK)*c(54,IJ,IK)*2.
     >        + K(264,ij,ik)*c(27,ij,ik)*c(55,ij,ik)
     >        + K(265,ij,ik)*c(27,ij,ik)*c(51,ij,ik)



C  Cly loss due to rainout of HCl (and others),  normalize total loss to Clx concentration
C  now use GMI COMBO sfc deposition + scavenging: SCAVDAY(20,L$,Z$) is in 1/sec (June 2011)
C

       clyl = (SCAVDAY(15,ij,ik)*c(60,ij,ik)
     >      +  SCAVDAY(17,ij,ik)*c(28,ij,ik)
     >      +  SCAVDAY(18,ij,ik)*c(30,ij,ik)
     >      +  SCAVDAY(19,ij,ik)*c(29,ij,ik)
     >      +  SCAVDAY(20,ij,ik)*c(25,ij,ik) )/c(33,ij,ik)


       CN(33,ij,ik) = (C(33,ij,ik) + clyprd*DT)/(1. + clyl*DT) 


       CTPROD(3,IJ,IK) = clyprd/m(ij,ik)
       CTLOSS(3,IJ,IK) = clyl*c(33,ij,ik)/m(ij,ik)


C
C
C  NOW update Clx radicals using new Cly family - update BOTH the C and CN arrays
C    this puts all the Cly increase into the Clx species,  cross-species are NOT changed
C
C   NOTE: as w/ NOy FOR THE CURRENT GRID POINT: 
C     if all the Cly is in the XSPECS, ie, fcrat LE 0, then set fcrat=1, 
C     and the Cly family and the Clx radicals are reset to the transported values (no chemical updates)


        fcrat = (CN(33,ij,ik) - xspecc)/CLX0

        if (fcrat .le. 0.) then 
          clyo    = CN(33,ij,ik)/M(IJ,IK)*1.d9
          xspecco = xspecc/M(IJ,IK)*1.d9
          clxo    = CLX0/M(IJ,IK)*1.d9

          write(1582,222) iyr, iday360, ij, ik,
     >            fcrat, clyo, xspecco, clxo
 222      format('IN SOLVER, Cly-fcrat le 0 : ', 4I4, 2x, 1P4D13.4)

          fcrat = 1.d0
          CN(33,ij,ik) = C(33,ij,ik)
        endif


        CN(25,ij,ik) = C(25,ij,ik)*fcrat
        CN(26,ij,ik) = C(26,ij,ik)*fcrat
        CN(27,ij,ik) = C(27,ij,ik)*fcrat
        CN(28,ij,ik) = C(28,ij,ik)*fcrat
        CN(29,ij,ik) = C(29,ij,ik)*fcrat
        CN(61,ij,ik) = C(61,ij,ik)*fcrat
        CN(62,ij,ik) = C(62,ij,ik)*fcrat
        CN(64,ij,ik) = C(64,ij,ik)*fcrat
        CN(65,ij,ik) = C(65,ij,ik)*fcrat

        C(25,ij,ik) = CN(25,ij,ik)
        C(26,ij,ik) = CN(26,ij,ik) 
        C(27,ij,ik) = CN(27,ij,ik) 
        C(28,ij,ik) = CN(28,ij,ik) 
        C(29,ij,ik) = CN(29,ij,ik) 
        C(61,ij,ik) = CN(61,ij,ik) 
        C(62,ij,ik) = CN(62,ij,ik) 
        C(64,ij,ik) = CN(64,ij,ik) 
        C(65,ij,ik) = CN(65,ij,ik) 


C  **********************************************************************************
C
C
C   first sum up Brx species from transported values, use C ARRAY , these will be updated with updates to Bry
C     DO NOT include x-species (BrONO2, BrNO2, BrCl) here
C
C     sum up x-species separately
C
CC  tfam - first sum up individual transported members, then apply the ratio of the 
CC         transported family, C(58), to each member, to adjust for discontinuities
C          INCLUDE cross-species BrCl, BrNO2, BrONO2, since these are small in Cl, NOx
C

        tbry = 2.*c(43,ij,ik) + c(44,ij,ik) + c(45,ij,ik)  
     >       + c(46,ij,ik) + c(47,ij,ik) + c(60,ij,ik)
     >       + c(68,ij,ik) + c(79,ij,ik)

        bryrat = c(58,ij,ik)/tbry

        c(43,ij,ik) = c(43,ij,ik)*bryrat
        c(44,ij,ik) = c(44,ij,ik)*bryrat
        c(45,ij,ik) = c(45,ij,ik)*bryrat
        c(46,ij,ik) = c(46,ij,ik)*bryrat
        c(47,ij,ik) = c(47,ij,ik)*bryrat
        c(60,ij,ik) = c(60,ij,ik)*bryrat
        c(68,ij,ik) = c(68,ij,ik)*bryrat
        c(79,ij,ik) = c(79,ij,ik)*bryrat

C                                    also update CN(47), CN(60), CN(79), CN(58) here
        CN(47,ij,ik) = C(47,ij,ik) 
        CN(60,ij,ik) = C(60,ij,ik) 
        CN(79,ij,ik) = C(79,ij,ik) 

        CN(58,ij,ik) = C(58,ij,ik) 

CC  end TFAM



        BRX0 = 2.*c(43,ij,ik) + c(44,ij,ik) + c(45,ij,ik)  
     >          + c(46,ij,ik) + c(68,ij,ik)


        xspecb = c(47,ij,ik) + c(79,ij,ik) + c(60,ij,ik)




c  Bry - CN(48) -  sum up total Bry family from transported/adjusted radicals, use C ARRAY 

 
        C(48,ij,ik) = 2.*c(43,ij,ik) + c(44,ij,ik) + c(45,ij,ik)  
     >                 + c(46,ij,ik) + c(47,ij,ik) + c(60,ij,ik)
     >                 + c(68,ij,ik) + c(79,ij,ik)



        bryprd = J(29,ij,ik)*c(49,ij,ik) + J(30,ij,ik)*c(50,ij,ik)
     >         + J(31,ij,ik)*c(51,ij,ik) + J(68,IJ,IK)*c(72,IJ,IK)*2.
     >         + J(54,IJ,IK)*c(75,IJ,IK)*2.
     >         + J(69,IJ,IK)*c(76,IJ,IK)*2.
     >         + J(70,IJ,IK)*c(77,IJ,IK)*3.
     >         + K(97,ij,ik)*c(49,ij,ik)*c(13,ij,ik)
     >         + K(112,IJ,IK)*c(51,ij,ik)*c(2,ij,ik)
     >         + K(113,IJ,IK)*c(50,ij,ik)*c(2,ij,ik)
     >         + K(123,IJ,IK)*c(72,ij,ik)*c(2,ij,ik)*2.
     >         + K(124,IJ,IK)*c(49,ij,ik)*c(2,ij,ik)
     >         + K(174,IJ,IK)*c(75,ij,ik)*c(2,ij,ik)*2.
     >         + K(176,IJ,IK)*c(49,ij,ik)*c(27,ij,ik)
     >         + K(177,IJ,IK)*c(76,ij,ik)*c(13,ij,ik)*2.
     >         + K(178,IJ,IK)*c(77,ij,ik)*c(13,ij,ik)*3.
     >         + K(179,IJ,IK)*c(76,ij,ik)*c(27,ij,ik)*2.
     >         + K(180,IJ,IK)*c(77,ij,ik)*c(27,ij,ik)*3.
     >         + K(183,IJ,IK)*c(76,ij,ik)*c(2,ij,ik)*2.
     >         + K(184,IJ,IK)*c(77,ij,ik)*c(2,ij,ik)*3.

     >         + K(241,IJ,IK)*c(13,IJ,IK)*c(75,IJ,IK)*2.
     >         + K(242,IJ,IK)*c(13,IJ,IK)*c(50,IJ,IK)
     >         + K(243,IJ,IK)*c(13,IJ,IK)*c(51,IJ,IK)
     >         + K(244,IJ,IK)*c(13,IJ,IK)*c(72,IJ,IK)*2.
     >         + K(265,ij,ik)*c(27,ij,ik)*c(51,ij,ik)
     >         + K(266,ij,ik)*c(27,ij,ik)*c(50,IJ,IK)
     >         + K(267,ij,ik)*c(27,ij,ik)*c(75,IJ,IK)*2.
     >         + K(268,ij,ik)*c(27,ij,ik)*c(72,IJ,IK)*2.



C  Bry loss due to rainout of HBr and others, normalize total loss to Bry concentration
C  now use GMI COMBO sfc deposition + scavenging: SCAVDAY(20,L$,Z$) is in 1/sec (June 2011)

        bryl = (SCAVDAY(12,ij,ik)*c(47,ij,ik)
     >       +  SCAVDAY(13,ij,ik)*c(46,ij,ik)
     >       +  SCAVDAY(14,ij,ik)*c(45,ij,ik)
     >       +  SCAVDAY(15,ij,ik)*c(60,ij,ik)
     >       +  SCAVDAY(16,ij,ik)*c(68,ij,ik) )/c(48,ij,ik)


        CN(48,ij,ik) = (C(48,ij,ik) + bryprd*DT)/(1. + bryl*DT) 


        CTPROD(15,IJ,IK) = bryprd/m(ij,ik)
        CTLOSS(15,IJ,IK) = bryl*c(48,ij,ik)/m(ij,ik)


C
C
C  NOW update Brx radicals using new Bry family - update BOTH the C and CN arrays
C       this puts all the Bry increase into the Brx species,  cross-species are NOT changed
C
C   NOTE: as w/ NOy/Cly FOR THE CURRENT GRID POINT: 
C         if all the Bry is in the XSPECS, ie, fbrat LE 0, then set fbrat=1, 
C         and the Bry family and the Brx radicals are reset to the transported values (no chemical updates)
C

        fbrat = (CN(48,ij,ik) - xspecb)/BRX0

        if (fbrat .le. 0.) then 
          bryo    = CN(48,ij,ik)/M(IJ,IK)*1.d12
          xspecbo = xspecb/M(IJ,IK)*1.d12
          brxo    = BRX0/M(IJ,IK)*1.d12

          write(1583,223) iyr, iday360, ij, ik,
     >            fbrat, bryo, xspecbo, brxo
 223      format('IN SOLVER, Bry-fbrat le 0 : ', 4I4, 2x, 1P4D13.4)

          fbrat = 1.d0
          CN(48,ij,ik) = C(48,ij,ik)
        endif



        CN(43,ij,ik) = C(43,ij,ik)*fbrat
        CN(44,ij,ik) = C(44,ij,ik)*fbrat 
        CN(45,ij,ik) = C(45,ij,ik)*fbrat 
        CN(46,ij,ik) = C(46,ij,ik)*fbrat 
        CN(68,ij,ik) = C(68,ij,ik)*fbrat 

        C(43,ij,ik) = CN(43,ij,ik)
        C(44,ij,ik) = CN(44,ij,ik)
        C(45,ij,ik) = CN(45,ij,ik)
        C(46,ij,ik) = CN(46,ij,ik)
        C(68,ij,ik) = CN(68,ij,ik)



C  ************************************************************************************
C
C
c   HF (Fx) -  CN(56)
C   NOTE: for washout of HF and CF2O, use HCl values from GMI COMBO, June 2011
c                           SCAVDAY(20,L$,Z$) is in 1/sec 

        HFP = J(21,ij,ik)*c(34,ij,ik) + J(30,ij,ik)*c(50,ij,ik)
     >      + J(33,IJ,IK)*c(53,IJ,IK) + J(35,ij,ik)*c(55,IJ,IK) 
     >	    + J(44,IJ,IK)*c(57,IJ,IK)*2. + J(65,IJ,IK)*C(69,IJ,IK)      
     >      + J(67,IJ,IK)*C(71,IJ,IK)
     >      + J(71,IJ,IK)*c(81,IJ,IK)*2.
     >      + J(72,IJ,IK)*c(82,IJ,IK)*3.
     >      + J(73,IJ,IK)*c(83,IJ,IK)
     >      + J(74,IJ,IK)*c(84,IJ,IK)*2.
     >      + J(75,IJ,IK)*c(85,IJ,IK)
     >      + J(76,IJ,IK)*c(86,IJ,IK)*2.
     >      + J(77,IJ,IK)*c(87,IJ,IK)
     >      + J(78,IJ,IK)*c(88,IJ,IK)
     >      + K(66,IJ,IK)*c(2,IJ,IK)*c(57,IJ,IK)*2. 
     >      + K(83,IJ,IK)*c(2,IJ,IK)*c(34,IJ,IK)
     >	    + K(99,IJ,IK)*c(2,IJ,IK)*c(53,IJ,IK) 
     >      + K(101,IJ,IK)*c(2,IJ,IK)*c(55,IJ,IK) 
     >      + K(113,IJ,IK)*c(50,ij,ik)*c(2,ij,ik) 
     >      + K(114,IJ,IK)*c(69,ij,ik)*c(2,ij,ik)
     >      + K(115,IJ,IK)*c(69,ij,ik)*c(13,ij,ik)
     >      + K(116,ij,ik)*c(27,ij,ik)*c(69,ij,ik)
     >      + K(120,IJ,IK)*c(2,IJ,IK)*c(71,IJ,IK)
     >      + K(121,IJ,IK)*c(13,IJ,IK)*c(71,IJ,IK)
     >      + K(122,ij,ik)*c(71,ij,ik)*c(27,ij,ik)

     >      + K(215,IJ,IK)*c(81,ij,ik)*c(2,ij,ik)*2.
     >      + K(224,IJ,IK)*c(82,ij,ik)*c(2,ij,ik)*3.
     >      + K(248,IJ,IK)*c(82,ij,ik)*c(13,ij,ik)*3.
     >      + K(256,IJ,IK)*c(82,ij,ik)*c(27,ij,ik)*3.
     >      + K(203,IJ,IK)*c(83,ij,ik)*c(2,ij,ik)
     >      + K(246,IJ,IK)*c(83,ij,ik)*c(13,ij,ik)
     >      + K(254,IJ,IK)*c(83,ij,ik)*c(27,ij,ik)
     >      + K(201,IJ,IK)*c(84,ij,ik)*c(2,ij,ik)*2.
     >      + K(253,IJ,IK)*c(84,ij,ik)*c(27,ij,ik)*2.
     >      + K(218,IJ,IK)*c(85,ij,ik)*c(2,ij,ik)
     >      + K(250,IJ,IK)*c(85,ij,ik)*c(13,ij,ik)
     >      + K(258,IJ,IK)*c(85,ij,ik)*c(27,ij,ik)
     >      + K(211,IJ,IK)*c(86,ij,ik)*c(2,ij,ik)*2.
     >      + K(247,IJ,IK)*c(86,ij,ik)*c(13,ij,ik)*2.
     >      + K(255,IJ,IK)*c(86,ij,ik)*c(27,ij,ik)*2.
     >      + K(233,IJ,IK)*c(87,ij,ik)*c(2,ij,ik)
     >      + K(252,IJ,IK)*c(87,ij,ik)*c(13,ij,ik)
     >      + K(269,IJ,IK)*c(87,ij,ik)*c(27,ij,ik)
     >      + K(235,IJ,IK)*c(88,ij,ik)*c(2,ij,ik)
     >      + K(251,IJ,IK)*c(88,ij,ik)*c(13,ij,ik)
     >      + K(270,IJ,IK)*c(88,ij,ik)*c(27,ij,ik)*3.

     >      + K(130,IJ,IK)*c(13,IJ,IK)*c(34,IJ,IK)  
     >      + K(240,IJ,IK)*c(13,IJ,IK)*c(55,IJ,IK)
     >      + K(242,IJ,IK)*c(13,IJ,IK)*c(50,IJ,IK)
     >      + K(260,IJ,IK)*c(27,IJ,IK)*c(34,IJ,IK)
     >      + K(262,IJ,IK)*c(27,IJ,IK)*c(53,IJ,IK)
     >      + K(264,ij,ik)*c(27,ij,ik)*c(55,ij,ik)
     >      + K(266,ij,ik)*c(27,ij,ik)*c(50,IJ,IK)



	hfl = SCAVDAY(19,ij,ik)    ! c(59,IJ,IK)


	cn(56,IJ,IK) = (c(56,IJ,IK) + hfp*dt)/(1. + hfl*dt)


        CTPROD(23,IJ,IK) = hfp/m(ij,ik)
        CTLOSS(23,IJ,IK) = hfl*c(56,ij,ik)/m(ij,ik)

C  ***********************************************************************************


c  CF2O   -   CN(57)   


        CF2Op = J(22,ij,ik)*c(35,ij,ik) + J(31,ij,ik)*c(51,ij,ik)
     >        + J(32,ij,ik)*c(52,ij,ik) + J(33,IJ,IK)*c(53,IJ,IK) 
     >	      +	J(34,IJ,IK)*c(54,IJ,IK)*2. + J(35,IJ,IK)*c(55,IJ,IK)*2. 
     >	      + J(30,IJ,IK)*c(50,IJ,IK) + J(66,ij,ik)*c(70,ij,ik)  
     >        + J(67,ij,ik)*c(71,ij,ik) + J(68,ij,ik)*c(72,ij,ik)*2.
     >        + J(54,IJ,IK)*c(75,IJ,IK)
     >        + J(71,IJ,IK)*c(81,IJ,IK)
     >        + J(73,IJ,IK)*c(83,IJ,IK)
     >        + J(75,IJ,IK)*c(85,IJ,IK)*2.
     >        + J(77,IJ,IK)*c(87,IJ,IK)*3.
     >        + J(78,IJ,IK)*c(88,IJ,IK)*2.
     >        + K(80,IJ,IK)*c(35,ij,ik)*c(2,ij,ik)
     >	      +	K(98,IJ,IK)*c(13,IJ,IK)*c(52,IJ,IK) 
     >	      +	K(99,IJ,IK)*c(2,IJ,IK)*c(53,IJ,IK) 
     >	      +	K(100,IJ,IK)*c(2,IJ,IK)*c(54,IJ,IK)*2. 
     >	      + K(101,IJ,IK)*c(2,IJ,IK)*c(55,IJ,IK)*2.  
     >        + K(112,IJ,IK)*c(51,ij,ik)*c(2,ij,ik)
     >        + K(113,IJ,IK)*c(50,ij,ik)*c(2,ij,ik) 
     >        + K(117,IJ,IK)*c(2,IJ,IK)*c(70,IJ,IK)
     >        + K(118,IJ,IK)*c(13,IJ,IK)*c(70,IJ,IK)
     >        + K(119,ij,ik)*c(70,ij,ik)*c(27,ij,ik)
     >        + K(120,IJ,IK)*c(2,IJ,IK)*c(71,IJ,IK)
     >        + K(121,IJ,IK)*c(13,IJ,IK)*c(71,IJ,IK)
     >        + K(122,ij,ik)*c(71,ij,ik)*c(27,ij,ik)
     >        + K(123,IJ,IK)*c(72,ij,ik)*c(2,ij,ik)*2.
     >        + K(125,IJ,IK)*c(52,IJ,IK)*c(2,IJ,IK)
     >        + K(174,IJ,IK)*c(75,ij,ik)*c(2,ij,ik)
     >        + K(181,IJ,IK)*c(52,ij,ik)*c(27,ij,ik)

     >        + K(215,IJ,IK)*c(81,ij,ik)*c(2,ij,ik)
     >        + K(249,IJ,IK)*c(81,ij,ik)*c(13,ij,ik)*2.
     >        + K(257,IJ,IK)*c(81,ij,ik)*c(27,ij,ik)*2.
     >        + K(203,IJ,IK)*c(83,ij,ik)*c(2,ij,ik)
     >        + K(246,IJ,IK)*c(83,ij,ik)*c(13,ij,ik)
     >        + K(254,IJ,IK)*c(83,ij,ik)*c(27,ij,ik)
     >        + K(245,IJ,IK)*c(84,ij,ik)*c(13,ij,ik)
     >        + K(218,IJ,IK)*c(85,ij,ik)*c(2,ij,ik)*2.
     >        + K(250,IJ,IK)*c(85,ij,ik)*c(13,ij,ik)*2.
     >        + K(258,IJ,IK)*c(85,ij,ik)*c(27,ij,ik)*2.
     >        + K(233,IJ,IK)*c(87,ij,ik)*c(2,ij,ik)*3.
     >        + K(252,IJ,IK)*c(87,ij,ik)*c(13,ij,ik)*3.
     >        + K(269,IJ,IK)*c(87,ij,ik)*c(27,ij,ik)*3.
     >        + K(235,IJ,IK)*c(88,ij,ik)*c(2,ij,ik)*2.
     >        + K(251,IJ,IK)*c(88,ij,ik)*c(13,ij,ik)*2.
     >        + K(270,IJ,IK)*c(88,ij,ik)*c(27,ij,ik)
     >        + K(131,IJ,IK)*c(13,IJ,IK)*c(35,IJ,IK)

     >        + K(238,IJ,IK)*c(13,IJ,IK)*c(53,IJ,IK)
     >        + K(239,IJ,IK)*c(13,IJ,IK)*c(54,IJ,IK)*2.
     >        + K(240,IJ,IK)*c(13,IJ,IK)*c(55,IJ,IK)*2.
     >        + K(241,IJ,IK)*c(13,IJ,IK)*c(75,IJ,IK)
     >        + K(242,IJ,IK)*c(13,IJ,IK)*c(50,IJ,IK)
     >        + K(243,IJ,IK)*c(13,IJ,IK)*c(51,IJ,IK)
     >        + K(244,IJ,IK)*c(13,IJ,IK)*c(72,IJ,IK)*2.
     >        + K(261,IJ,IK)*c(27,IJ,IK)*c(35,IJ,IK)
     >        + K(262,IJ,IK)*c(27,IJ,IK)*c(53,IJ,IK)
     >        + K(263,IJ,IK)*c(27,IJ,IK)*c(54,IJ,IK)*2.
     >        + K(264,ij,ik)*c(27,ij,ik)*c(55,ij,ik)*2.
     >        + K(265,ij,ik)*c(27,ij,ik)*c(51,ij,ik)
     >        + K(266,ij,ik)*c(27,ij,ik)*c(50,IJ,IK)
     >        + K(267,ij,ik)*c(27,ij,ik)*c(75,IJ,IK)
     >        + K(268,ij,ik)*c(27,ij,ik)*c(72,IJ,IK)*2.



	CF2Ol = J(44,IJ,IK) + K(66,IJ,IK)*c(2,IJ,IK) + SCAVDAY(19,ij,ik)    ! c(59,IJ,IK)   


        cn(57,IJ,IK) = (c(57,IJ,IK) + cf2op*dt)/(1. + CF2Ol*dt)


        CTPROD(24,IJ,IK) = CF2Op/m(ij,ik)
        CTLOSS(24,IJ,IK) = CF2Ol*c(57,ij,ik)/m(ij,ik)

C  **********************************************************************
C
C
C   H2O   -   CN(15)   
c
c      first loop through diurnally varying production, CNDC(13 (OH) is already adjusted
c
      ksum1 = 0.D0   
      ksum2 = 0.D0   
      ksum3 = 0.D0   
      ksum4 = 0.D0   
      ksum5 = 0.D0   
      ksum6 = 0.D0   
      ksum7 = 0.D0   
      ksum8 = 0.D0   
      ksum9 = 0.D0   
      ksum0 = 0.D0   
      ksum10 = 0.D0   
      ksum11 = 0.D0
      do 217 it=1,isntime-1
        ksum1 = ksum1 + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(29,it,ij,ik) +
     >                         CNDC(13,it+1,ij,ik)*CNDC(29,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksum2 = ksum2 + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(16,it,ij,ik) +
     >                         CNDC(13,it+1,ij,ik)*CNDC(16,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksum3 = ksum3 + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(14,it,ij,ik) +
     >                         CNDC(13,it+1,ij,ik)*CNDC(14,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksum4 = ksum4 + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(23,it,ij,ik) +
     >                         CNDC(13,it+1,ij,ik)*CNDC(23,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksum5 = ksum5 + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(38,it,ij,ik) +
     >                         CNDC(13,it+1,ij,ik)*CNDC(38,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksum6 = ksum6 + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(24,it,ij,ik) +
     >                         CNDC(13,it+1,ij,ik)*CNDC(24,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksum7 = ksum7 + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(13,it,ij,ik) +
     >                         CNDC(13,it+1,ij,ik)*CNDC(13,it+1,ij,ik))
     >                               *DTIMEB(ij,it)/86400.D0

        ksum8 = ksum8 + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(25,it,ij,ik) +
     >                         CNDC(13,it+1,ij,ik)*CNDC(25,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksum9 = ksum9 + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(46,it,ij,ik) +
     >                         CNDC(13,it+1,ij,ik)*CNDC(46,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksum0 = ksum0 + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(42,it,ij,ik) +
     >                         CNDC(13,it+1,ij,ik)*CNDC(42,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksum10 = ksum10 + 0.5D0*(CNDC(12,it,ij,ik)*CNDC(14,it,ij,ik) +
     >                         CNDC(12,it+1,ij,ik)*CNDC(14,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksum11 = ksum11 + 0.5D0*(CNDC(13,it,ij,ik)*CNDC(89,it,ij,ik) +
     >                         CNDC(13,it+1,ij,ik)*CNDC(89,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0
 217  CONTINUE


       h2op = J(20,ij,ik)*c(37,ij,ik) + J(26,ij,ik)*c(40,ij,ik) 
     >      + J(65,IJ,IK)*C(69,IJ,IK) + J(66,ij,ik)*c(70,ij,ik) 
     >      + J(76,IJ,IK)*c(86,IJ,IK)
     >      + J(78,IJ,IK)*c(88,IJ,IK)
     >      + K(14,ij,ik)*c(13,ij,ik)*c(18,ij,ik)
     >      + K(27,ij,ik)*ksum1                                                ! c(13,ij,ik)*c(29,ij,ik)
     >      + K(29,ij,ik)*ksum2                                                ! c(13,ij,ik)*c(16,ij,ik)
     >      + K(30,ij,ik)*c(13,ij,ik)*c(17,ij,ik)
     >      + K(37,ij,ik)*c(13,ij,ik)*c(10,ij,ik)
     >      + K(40,ij,ik)*ksum3                                                ! c(13,ij,ik)*c(14,ij,ik)
     >      + K(51,ij,ik)*ksum4                                                ! c(13,ij,ik)*c(23,ij,ik)
     >      + K(55,ij,ik)*ksum5                                                ! c(13,ij,ik)*c(38,ij,ik)
     >      + K(58,ij,ik)*ksum6                                                ! c(13,ij,ik)*c(24,ij,ik)
     >      + K(59,ij,ik)*ksum7                                                ! c(13,ij,ik)*c(13,ij,ik)
     >      + K(62,ij,ik)*ksum8                                                ! c(13,ij,ik)*c(25,ij,ik)
     >      + K(73,IJ,IK)*ksum10                                               ! c(12,IJ,IK)*c(14,IJ,IK)
     >      + K(75,ij,ik)*c(13,ij,ik)*c(40,ij,ik)*2.
     >      + K(95,ij,ik)*ksum9                                                ! c(13,ij,ik)*c(46,ij,ik)
     >      + K(97,IJ,IK)*c(13,IJ,IK)*c(49,ij,ik)
     >	    + K(98,IJ,IK)*c(13,IJ,IK)*c(52,IJ,IK)
     >      + K(114,IJ,IK)*c(69,ij,ik)*c(2,ij,ik)
     >      + K(115,IJ,IK)*c(69,ij,ik)*c(13,ij,ik)*2.
     >      + K(117,IJ,IK)*c(70,ij,ik)*c(2,ij,ik)
     >      + K(118,IJ,IK)*c(70,ij,ik)*c(13,ij,ik)*2.
     >      + K(121,IJ,IK)*c(13,IJ,IK)*c(71,IJ,IK)
     >      + K(177,IJ,IK)*c(76,ij,ik)*c(13,ij,ik)
     >      + K(178,IJ,IK)*c(77,ij,ik)*c(13,ij,ik)

     >      + K(248,IJ,IK)*c(82,ij,ik)*c(13,ij,ik)
     >      + K(249,IJ,IK)*c(81,ij,ik)*c(13,ij,ik)
     >      + K(246,IJ,IK)*c(83,ij,ik)*c(13,ij,ik)
     >      + K(245,IJ,IK)*c(84,ij,ik)*c(13,ij,ik)
     >      + K(250,IJ,IK)*c(85,ij,ik)*c(13,ij,ik)
     >      + K(211,IJ,IK)*c(86,ij,ik)*c(2,ij,ik)
     >      + K(247,IJ,IK)*c(86,ij,ik)*c(13,ij,ik)
     >      + K(252,IJ,IK)*c(87,ij,ik)*c(13,ij,ik)
     >      + K(235,IJ,IK)*c(88,ij,ik)*c(2,ij,ik)
     >      + K(251,IJ,IK)*c(88,ij,ik)*c(13,ij,ik)
     >      + K(229,IJ,IK)*c(40,ij,ik)*c(2,ij,ik)

     >      + K(148,IJ,IK)*ksum6                                               ! c(13,ij,ik)*c(24,ij,ik)
     >      + K(150,IJ,IK)*ksum0                                               ! c(13,ij,ik)*c(42,ij,ik)
     >      + K(273,IJ,IK)*ksum11                                              ! c(13,IJ,IK)*c(89,IJ,IK)

     >      + ksumh(5)                                                         ! HET, HOCl+HCl
     >      + ksumh(7)                                                         ! HET, HOBr+HCl
     >      + ksumh(11)                                                        ! HET, HOCl+HBr
     >      + ksumh(13)                                                        ! HET, HOBr+HBr


ccccccccc     >      + K(173,IJ,IK)*c(25,ij,ik)*c(46,ij,ik)   - currently = 0.0 

cccccc     >      + c(25,ij,ik)*c(29,ij,ik)*kh(5,ij,ik)
cccccc     >      + c(68,ij,ik)*c(29,ij,ik)*kh(7,ij,ik)


       h2ol = J(4,ij,ik) + J(25,ij,ik) 
     >      + K(39,ij,ik)*c(2,ij,ik)
     >      + K(276,IJ,IK)*c(6,IJ,IK) 
     >      + ksumh(2)/c(15,ij,ik)                                                  ! HET ClONO2 + H2O
     >      + ksumh(3)/c(15,ij,ik)                                                  ! HET N2O5 + H2O
     >      + ksumh(6)/c(15,ij,ik)                                                  ! HET BrONO2 + H2O
     >      + ksumh(8)/c(15,ij,ik)                                                  ! HET NO3 + H2O


ccccccc     >      + K(79,ij,ik)*c(8,ij,ik)                               ! gas phase N2O5+H2O = 0.0
ccccccc     >      + c(30,ij,ik)*kh(2,ij,ik) + c(8,ij,ik)*kh(3,ij,ik)
ccccccc     >      + c(47,ij,ik)*kh(6,ij,ik)


       cn(15,ij,ik) = (c(15,ij,ik) + h2op*dt)/(1. + h2ol*dt)


       CTPROD(25,IJ,IK) = h2op/m(ij,ik)
       CTLOSS(25,IJ,IK) = h2ol*c(15,ij,ik)/m(ij,ik)


C  ****************************************************************************************


C  H2   -   CN(17) 


      jsum1 = 0.D0   
      ksum1 = 0.D0   
      do 212 it=1,isntime-1
        jsum1 = jsum1 + 0.5D0*(JDC(11,it,ij,ik)*CNDC(23,it,ij,ik) +
     >                         JDC(11,it+1,ij,ik)*CNDC(23,it+1,ij,ik))
     >                                    *DTIMEB(ij,it)/86400.D0
 212  CONTINUE


        h2prod = jsum1 + J(25,IJ,IK)*c(15,IJ,IK)                     ! J(11,IJ,IK)*C(23,IJ,IK)
     >     + (J(59,ij,ik) + J(60,ij,ik) + J(61,ij,ik))*c(18,ij,ik)*2.
     >         + K(50,IJ,IK)*c(2,IJ,IK)*c(18,IJ,IK)
     >         + K(71,IJ,IK)*ksum10                                  ! c(12,IJ,IK)*c(14,IJ,IK)

        h2loss = K(23,IJ,IK)*c(27,IJ,IK) + K(30,IJ,IK)*c(13,IJ,IK)
     >         + K(48,IJ,IK)*c(2,IJ,IK)


       	cn(17,IJ,IK) = (c(17,IJ,IK) + h2prod*dt)/(1. + h2loss*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF(IK.EQ.1)THEN
           IF(LBCMRSS(17)) CN(17,IJ,IK) = BVAL(17,IJ)*M(IJ,IK)
	ENDIF


        CTPROD(10,IJ,IK) = H2PROD/m(ij,ik)
        CTLOSS(10,IJ,IK) = H2LOSS*c(17,ij,ik)/m(ij,ik)
        RLOSS(28,IJ,IK)=H2LOSS

C  ****************************************************************************************
C
c
C   CO  -   CN(19)   
C
C   first sum up diurnal variation of JH2CO; JDC(PH$,18,L$,Z$), CNDC(S$,18+3,L$,Z$), DTIMEB(L$,ISNTIME)
C
      jsum1 = 0.D0   
      jsum2 = 0.D0   
      ksum1 = 0.D0   
      ksum2 = 0.D0   
      ksum3 = 0.D0   
      ksum4 = 0.D0   
      do 211 it=1,isntime-1
        jsum1 = jsum1 + 0.5D0*(JDC(10,it,ij,ik)*CNDC(23,it,ij,ik) +
     >                         JDC(10,it+1,ij,ik)*CNDC(23,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

        jsum2 = jsum2 + 0.5D0*(JDC(11,it,ij,ik)*CNDC(23,it,ij,ik) +
     >                         JDC(11,it+1,ij,ik)*CNDC(23,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

        ksum1 = ksum1 + 0.5D0*(CNDC(23,it,ij,ik)*CNDC(1,it,ij,ik) +
     >                         CNDC(23,it+1,ij,ik)*CNDC(1,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

        ksum2 = ksum2 + 0.5D0*(CNDC(23,it,ij,ik)*CNDC(13,it,ij,ik) +
     >                         CNDC(23,it+1,ij,ik)*CNDC(13,it+1,ij,ik))
     >                                     *DTIMEB(ij,it)/86400.D0

        ksum3 = ksum3 + 0.5D0*(CNDC(23,it,ij,ik)*CNDC(27,it,ij,ik) +
     >                         CNDC(23,it+1,ij,ik)*CNDC(27,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0

        ksum4 = ksum4 + 0.5D0*(CNDC(23,it,ij,ik)*CNDC(45,it,ij,ik) +
     >                         CNDC(23,it+1,ij,ik)*CNDC(45,it+1,ij,ik))
     >                                      *DTIMEB(ij,it)/86400.D0
 211  CONTINUE


CWCONV200 - CO surface emissions didn't match GMI very well (needed to increase 2x) - March 2011
CWCONV200   EMISDAY(L$,3) is in #/cm3/sec;  1=CH2O  2=CO    3=NO (NOx)
 
ccc        coemis = 0.d0
ccc        if (ik .eq. 1) coemis = 2.*EMISDAY(ij,2)


	COPROD = jsum1 + jsum2                          ! J(10,IJ,IK)*C(23,IJ,IK) + J(11,IJ,IK)*C(23,IJ,IK)
     >         + (J(12,IJ,IK) + J(41,IJ,IK))*C(20,IJ,IK)
     >         + J(19,ij,ik)*c(36,ij,ik) + J(20,IJ,IK)*c(37,ij,ik) 
     >         + J(21,IJ,IK)*c(34,ij,ik) + J(26,IJ,IK)*c(40,ij,ik)*2.
     >         + J(33,IJ,IK)*c(53,ij,ik) + J(44,IJ,IK)*c(57,IJ,IK)
     >         +(J(59,IJ,IK) + J(60,IJ,IK) + J(61,IJ,IK))*C(18,IJ,IK)
     >         + J(65,IJ,IK)*c(69,IJ,IK)*2.
     >         + J(66,IJ,IK)*c(70,IJ,IK) + J(67,IJ,IK)*c(71,IJ,IK)
     >         + J(69,IJ,IK)*c(76,IJ,IK) + J(70,IJ,IK)*c(77,IJ,IK)
     >         + J(71,IJ,IK)*c(81,IJ,IK)
     >         + J(72,IJ,IK)*c(82,IJ,IK)*2.
     >         + J(74,IJ,IK)*c(84,IJ,IK)
     >         + J(76,IJ,IK)*c(86,IJ,IK)*2.
     >         + J(78,IJ,IK)*c(88,IJ,IK)
     >         + K(21,IJ,IK)*ksum1                                      ! c(23,IJ,IK)*c(1,IJ,IK) 
     >         + K(51,IJ,IK)*ksum2                                      ! c(23,IJ,IK)*c(13,IJ,IK)
     >         + K(54,IJ,IK)*c(36,ij,ik)*c(2,ij,ik)
     >         + K(63,IJ,IK)*ksum3                                      ! c(23,IJ,IK)*c(27,IJ,IK)
     >         + K(75,IJ,IK)*c(13,IJ,IK)*c(40,IJ,IK)*2.
     >         + K(83,IJ,IK)*c(2,IJ,IK)*c(34,IJ,IK)
     >         + K(97,IJ,IK)*c(13,IJ,IK)*c(49,IJ,IK)
     >         + K(99,IJ,IK)*c(2,IJ,IK)*c(53,IJ,IK)
     >         + K(109,IJ,IK)*ksum4                                     ! c(23,IJ,IK)*c(45,IJ,IK)
     >         + K(114,IJ,IK)*c(2,IJ,IK)*c(69,IJ,IK)*2.
     >         + K(115,IJ,IK)*c(13,IJ,IK)*c(69,IJ,IK)*2.
     >         + K(116,ij,ik)*c(27,ij,ik)*c(69,ij,ik)*2.
     >         + K(117,IJ,IK)*c(2,IJ,IK)*c(70,IJ,IK)
     >         + K(118,IJ,IK)*c(13,IJ,IK)*c(70,IJ,IK)
     >         + K(119,ij,ik)*c(70,ij,ik)*c(27,ij,ik)
     >         + K(120,IJ,IK)*c(2,IJ,IK)*c(71,IJ,IK)
     >         + K(121,IJ,IK)*c(13,IJ,IK)*c(71,IJ,IK)
     >         + K(122,ij,ik)*c(27,ij,ik)*c(71,ij,ik)
     >         + K(132,IJ,IK)*c(27,IJ,IK)*c(37,IJ,IK)
     >         + K(176,IJ,IK)*c(49,ij,ik)*c(27,ij,ik)
     >         + K(177,IJ,IK)*c(76,ij,ik)*c(13,ij,ik)
     >         + K(178,IJ,IK)*c(77,ij,ik)*c(13,ij,ik)
     >         + K(179,IJ,IK)*c(76,ij,ik)*c(27,ij,ik)
     >         + K(180,IJ,IK)*c(77,ij,ik)*c(27,ij,ik)
     >         + K(182,IJ,IK)*c(40,ij,ik)*c(27,ij,ik)*2.
     >         + K(184,IJ,IK)*c(77,ij,ik)*c(2,ij,ik)

     >         + K(215,IJ,IK)*c(81,ij,ik)*c(2,ij,ik)
     >         + K(224,IJ,IK)*c(82,ij,ik)*c(2,ij,ik)*2.
     >         + K(248,IJ,IK)*c(82,ij,ik)*c(13,ij,ik)*2.
     >         + K(256,IJ,IK)*c(82,ij,ik)*c(27,ij,ik)*2.
     >         + K(201,IJ,IK)*c(84,ij,ik)*c(2,ij,ik)
     >         + K(253,IJ,IK)*c(84,ij,ik)*c(27,ij,ik)
     >         + K(211,IJ,IK)*c(86,ij,ik)*c(2,ij,ik)*2.
     >         + K(247,IJ,IK)*c(86,ij,ik)*c(13,ij,ik)*2.
     >         + K(255,IJ,IK)*c(86,ij,ik)*c(27,ij,ik)*2.
     >         + K(235,IJ,IK)*c(88,ij,ik)*c(2,ij,ik)
     >         + K(251,IJ,IK)*c(88,ij,ik)*c(13,ij,ik)
     >         + K(270,IJ,IK)*c(88,ij,ik)*c(27,ij,ik)*2.
     >         + K(129,IJ,IK)*c(13,IJ,IK)*c(36,IJ,IK) 
     >         + K(130,IJ,IK)*c(13,IJ,IK)*c(34,IJ,IK)  
     >         + K(229,IJ,IK)*c(40,ij,ik)*c(2,ij,ik)*2.
     >         + K(238,IJ,IK)*c(13,IJ,IK)*c(53,IJ,IK)  
     >         + K(259,IJ,IK)*c(27,IJ,IK)*c(36,IJ,IK)
     >         + K(260,IJ,IK)*c(27,IJ,IK)*c(34,IJ,IK)
     >         + K(262,IJ,IK)*c(27,IJ,IK)*c(53,IJ,IK)

ccccccc     >         + coemis


        COLOSS = K(36,IJ,IK)*C(13,IJ,IK) + K(57,IJ,IK)*C(13,IJ,IK)
     >         + K(111,IJ,IK)*C(1,IJ,IK)*M(IJ,IK)

   
	CN(19,IJ,IK) = (C(19,IJ,IK) + COPROD*DT)/(1. + COLOSS*DT)

c
C PUT IN MIXING RATIO BOUNDARY CONDITIONS
ccccccc	        IF(IK.EQ.1)THEN
ccccccc           IF(LBCMRSS(19)) CN(19,IJ,IK) = BVAL(19,IJ)*M(IJ,IK)
ccccccc    	ENDIF
C
C
CWCONV200   use GMI values at surface - BVALG(S$,L$,10) is in COMMON (ppp)
C           IF (zalt(ik) .le. 5.) CN(19,IJ,IK) = BVALG(19,ij,ik)*M(ij,ik)

        IF (ik .eq. 1) CN(19,IJ,1) = BVALG(19,ij,1)*M(ij,ik)


        CTPROD(11,IJ,IK) = COPROD/m(ij,ik)
        CTLOSS(11,IJ,IK) = COLOSS*c(19,ij,ik)/m(ij,ik)
        RLOSS(27,IJ,IK) = COLOSS


C
C   ***************   CO2  computation  (Chemically active ONLY)  *******************************
C
C   *****  Time Dependent CO2 mixing ratio Boundary conditions in PPMV for 1950-2050, daily values
C         outside of 1979-1995 were extrapolated using linear increase of 1.4 ppmv/year
C         CO2BC(L$,101,360) is in COMMON, ij=1-L$ (interpolated from 18 latitudes), 101 years, 360 days,
C
C
C   CO2 = CN(20)  
C

	CO2PROD = K(36,IJ,IK)*C(13,IJ,IK)*C(19,IJ,IK)
     >          + K(57,IJ,IK)*C(13,IJ,IK)*C(19,IJ,IK)
     >          + K(111,IJ,IK)*C(19,IJ,IK)*C(1,IJ,IK)*M(IJ,IK)
     >          + K(66,IJ,IK)*c(2,IJ,IK)*c(57,IJ,IK)                   ! CF2O + O(1D) -> 2HF + CO2 (now)


	CO2LOSS = J(12,IJ,IK) + J(41,IJ,IK)

     
	CN(20,IJ,IK) = (C(20,IJ,IK) + CO2PROD*DT)/(1. + CO2LOSS*DT)

C
C  for CO2 BC, choose either daily BC (for 1950-2050, above), or WMO yearly values from BC.DAT file 
C                                                            for 1950-2099 defined in BOUNDC
        IF (IK .EQ. 1) then 
ccccccc           CN(20,IJ,IK) = co2bc(ij,iyr-14,iday360)*1.E-6*M(IJ,IK)
           CN(20,IJ,IK) = BVAL(20,IJ)*M(IJ,IK)
ccccccc           CN(20,IJ,IK) = 389.*1.e-6*M(IJ,IK)
        ENDIF



        CTPROD(12,IJ,IK) = CO2PROD/m(ij,ik)
        CTLOSS(12,IJ,IK) = CO2LOSS*c(20,ij,ik)/m(ij,ik)
        RLOSS(32,IJ,IK) = CO2LOSS

C   ***************   END   CO2  computation   **************************************************



C  CH4   -    CN(18)  

	ch4prod = 0.0

        ch4loss = J(59,IJ,IK) + J(60,ij,ik) + J(61,ij,ik)
     >          + K(14,IJ,IK)*c(13,IJ,IK) + K(26,IJ,IK)*c(27,IJ,IK)
     >          + (K(49,IJ,IK) + K(50,IJ,IK))*c(2,IJ,IK)


	cn(18,IJ,IK) = (c(18,IJ,IK) + ch4prod*dt)/(1. + ch4loss*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS

      IF (IK .EQ. 1) THEN
        IF(LBCMRSS(18).OR.LBCMRTD(18)) CN(18,IJ,IK)=BVAL(18,IJ)*M(IJ,IK)
ccccccc       CN(18,IJ,IK) = 1802.*1.e-9*M(IJ,IK)
      ENDIF


        CTPROD(9,IJ,IK) = CH4PROD/m(ij,ik)
        CTLOSS(9,IJ,IK) = CH4LOSS*c(18,ij,ik)/m(ij,ik)
        RLOSS(6,IJ,IK) = CH4LOSS


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 7

        xlife(ipf,1,ij,ik) = ch4loss
        xlife(ipf,2,ij,ik) = J(59,IJ,IK) + J(60,ij,ik) + J(61,ij,ik)
        xlife(ipf,3,ij,ik) = (K(49,IJ,IK) + K(50,IJ,IK))*c(2,IJ,IK)
        xlife(ipf,4,ij,ik) = K(14,IJ,IK)*c(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(26,IJ,IK)*c(27,IJ,IK)
        xlife(ipf,6,ij,ik) = cn(18,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = 
     >      J39(59,ijj0,ij,ik) + J39(60,ijj0,ij,ik) + J39(61,ijj0,ij,ik)
        end do


C   *********************************************************************************************


C   N2O    -   CN(11)   

                                                      
        pn2o = K(81,IJ,IK)*ksum69                       ! c(6,IJ,IK)*c(9,IJ,IK) 
     >       + K(82,IJ,IK)*C(2,IJ,IK)*N2(IJ,IK)         ! K(82) NOW DOES include M, also N2 array used for N2

        N2OLOSS = J(14,ij,ik) + (K(45,ij,ik) + K(78,ij,ik))*c(2,ij,ik)
     >                        +  K(237,ij,ik)*c(13,ij,ik)


	cn(11,IJ,IK) = (c(11,IJ,IK) + pn2o*dt)/(1. + N2OLOSS*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
      IF (IK .EQ. 1) THEN
        IF(LBCMRSS(11).OR.LBCMRTD(11)) CN(11,IJ,IK)=BVAL(11,IJ)*M(IJ,IK)
ccccccc      CN(11,IJ,IK) = 323.5*1.e-9*M(IJ,IK)
      ENDIF


        CTPROD(4,IJ,IK) = pn2o/m(ij,ik)
        CTLOSS(4,IJ,IK) = n2oloss*c(11,ij,ik)/m(ij,ik)
        RLOSS(1,IJ,IK)=N2OLOSS


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 6

        xlife(ipf,1,ij,ik) = n2oloss
        xlife(ipf,2,ij,ik) = J(14,ij,ik)
        xlife(ipf,3,ij,ik) = (K(45,ij,ik) + K(78,ij,ik))*c(2,ij,ik)
        xlife(ipf,4,ij,ik) = K(237,ij,ik)*c(13,ij,ik)
        xlife(ipf,5,ij,ik) = 0.
        xlife(ipf,6,ij,ik) = cn(11,IJ,IK)

C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(14,ijj0,ij,ik)
        end do


C   *********************************************************************************************
C
C
C  CN(34)  IS CFCL3   


	f11p=0.

C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
	IF(.NOT.LBCMRSS(34) .OR. .NOT.LBCMRTD(34))
     *    F11P=BVAL(34,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	F11L = J(21,IJ,IK) + C(2,IJ,IK)*K(83,IJ,IK)
     >                     + K(130,IJ,IK)*c(13,IJ,IK)
     >                     + K(260,IJ,IK)*c(27,IJ,IK)


	cn(34,IJ,IK) = (c(34,IJ,IK) + f11p*dt)/(1. + f11l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
      IF (IK .EQ. 1) THEN
        IF(LBCMRSS(34).OR.LBCMRTD(34)) CN(34,IJ,IK)=BVAL(34,IJ)*M(IJ,IK)
ccccccc          CN(34,IJ,IK) = .2409*1.e-9*M(IJ,IK)
      ENDIF


        CTPROD(5,IJ,IK) = f11p/m(ij,ik)
        CTLOSS(5,IJ,IK) = F11L*c(34,ij,ik)/m(ij,ik)
        RLOSS(2,IJ,IK)=F11L


C
C   xlife(30,6,L$,Z$)
C
C   xlife:  1st index:  1=CFC-11;  2=CFC-12;  3=Carb Tet;  4=Meth Chlor;   5=HCFC-22;  6=N2O;  7=CH4
C           2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 1

        xlife(ipf,1,ij,ik) = F11L
        xlife(ipf,2,ij,ik) = J(21,IJ,IK)
        xlife(ipf,3,ij,ik) = C(2,IJ,IK)*K(83,IJ,IK)
        xlife(ipf,4,ij,ik) = K(130,IJ,IK)*c(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(260,IJ,IK)*c(27,IJ,IK)
        xlife(ipf,6,ij,ik) = cn(34,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(21,ijj0,ij,ik)
        end do


C   *********************************************************************************************
C
C
C  CN(35)  IS CF2CL2   


	f12p=0.

C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
	IF(.NOT.LBCMRSS(35) .OR. .NOT.LBCMRTD(35))
     *  F12P=BVAL(35,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	F12L = J(22,IJ,IK) + C(2,IJ,IK)*K(80,IJ,IK)
     >                     + K(131,IJ,IK)*c(13,IJ,IK)
     >                     + K(261,IJ,IK)*c(27,IJ,IK)


	cn(35,IJ,IK) = (c(35,IJ,IK) + f12p*dt)/(1. + f12l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF (IK.EQ.1) THEN
	IF(LBCMRSS(35).OR.LBCMRTD(35))
     *    CN(35,IJ,IK)=BVAL(35,IJ)*M(IJ,IK)
ccccccc          CN(35,IJ,IK) = .5325*1.e-9*M(IJ,IK)
	ENDIF


        CTPROD(6,IJ,IK) = F12p/m(ij,ik)
        CTLOSS(6,IJ,IK) = F12L*c(35,ij,ik)/m(ij,ik)
        RLOSS(3,IJ,IK)=F12L


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 2

        xlife(ipf,1,ij,ik) = F12L
        xlife(ipf,2,ij,ik) = J(22,IJ,IK)
        xlife(ipf,3,ij,ik) = C(2,IJ,IK)*K(80,IJ,IK)
        xlife(ipf,4,ij,ik) = K(131,IJ,IK)*c(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(261,IJ,IK)*c(27,IJ,IK)
        xlife(ipf,6,ij,ik) = cn(35,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(22,ijj0,ij,ik)
        end do


C   *********************************************************************************************


C  CN(36)  IS CCL4    


	ccl4p=0.

C PUT IN FLUX BOUNDARY CONDITIONS  
	IF (IK .EQ. 1) THEN
          IF(.NOT.LBCMRSS(36) .OR. .NOT.LBCMRTD(36))
     >        CCL4P = BVAL(36,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	CCL4L = J(19,IJ,IK) + K(54,IJ,IK)*C(2,IJ,IK)
     >                      + K(129,IJ,IK)*c(13,IJ,IK)
     >                      + K(259,IJ,IK)*c(27,IJ,IK)


	cn(36,IJ,IK) = (c(36,IJ,IK) + ccl4p*dt)/(1. + ccl4l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF (IK.EQ.1) THEN
	IF(LBCMRSS(36).OR.LBCMRTD(36))
     >    CN(36,IJ,IK) = BVAL(36,IJ)*M(IJ,IK)
ccccccc            CN(36,IJ,IK) = .0876*1.e-9*M(IJ,IK)
	ENDIF


        CTPROD(7,IJ,IK) = ccl4p/m(ij,ik)
        CTLOSS(7,IJ,IK) = CCL4L*c(36,ij,ik)/m(ij,ik)
        RLOSS(4,IJ,IK)=CCL4L


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 3

        xlife(ipf,1,ij,ik) = CCL4L
        xlife(ipf,2,ij,ik) = J(19,IJ,IK)
        xlife(ipf,3,ij,ik) = K(54,IJ,IK)*C(2,IJ,IK)
        xlife(ipf,4,ij,ik) = K(129,IJ,IK)*c(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(259,IJ,IK)*c(27,IJ,IK)   
        xlife(ipf,6,ij,ik) = cn(36,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(19,ijj0,ij,ik)
        end do


C   *********************************************************************************************


C  CN(37)  is CH3Cl   


	ch3clp=0.

	CH3CLl = J(20,IJ,IK) + K(16,IJ,IK)*C(13,IJ,IK)
     >                       + K(132,IJ,IK)*c(27,IJ,IK)
     >                       + K(136,IJ,IK)*c(2,IJ,IK)


	CN(37,ij,ik) = (C(37,ij,ik) + ch3clP*dt)/(1. + ch3cll*dt) 


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF (IK .EQ. 1) THEN
          IF (LBCMRSS(37) .OR. LBCMRTD(37))
     >      CN(37,IJ,IK) = BVAL(37,IJ)*M(IJ,IK)
ccccccc            CN(37,IJ,IK) = .5500*1.e-9*M(IJ,IK)
        ENDIF


        CTPROD(8,IJ,IK) = ch3clp/m(ij,ik)
        CTLOSS(8,IJ,IK) = CH3CLL*c(37,ij,ik)/m(ij,ik)
        RLOSS(5,IJ,IK)=CH3CLL


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 18

        xlife(ipf,1,ij,ik) = CH3CLL
        xlife(ipf,2,ij,ik) = J(20,IJ,IK)
        xlife(ipf,3,ij,ik) = K(136,IJ,IK)*c(2,IJ,IK)
        xlife(ipf,4,ij,ik) = K(16,IJ,IK)*c(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(132,IJ,IK)*c(27,IJ,IK)
        xlife(ipf,6,ij,ik) = cn(37,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(20,ijj0,ij,ik)
        end do


C   *********************************************************************************************


C  CN(40)  IS CH3CCL3  


	CH3CCL3P=0.0E0

C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
           IF(.NOT.LBCMRSS(40) .OR. .NOT.LBCMRTD(40))
     >         CH3CCL3P = BVAL(40,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	CH3CCL3L = J(26,IJ,IK) + K(75,IJ,IK)*C(13,IJ,IK)
     >                         + K(182,IJ,IK)*c(27,ij,ik)
     >                         + K(229,IJ,IK)*c(2,ij,ik)


	CN(40,IJ,IK) = (C(40,IJ,IK) + CH3CCL3P*DT)/(1. + CH3CCL3L*DT)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
      IF (IK .EQ. 1) THEN
        IF(LBCMRSS(40).OR.LBCMRTD(40)) CN(40,IJ,IK)=BVAL(40,IJ)*M(IJ,IK)
ccccccc             CN(40,IJ,IK) = .0083*1.e-9*M(IJ,IK)
      ENDIF


        CTPROD(13,IJ,IK) = CH3CCL3P/m(ij,ik)
        CTLOSS(13,IJ,IK) = CH3CCL3L*c(40,ij,ik)/m(ij,ik)
        RLOSS(7,IJ,IK)=CH3CCL3L


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 4

        xlife(ipf,1,ij,ik) = CH3CCL3L
        xlife(ipf,2,ij,ik) = J(26,IJ,IK)
        xlife(ipf,3,ij,ik) = K(229,IJ,IK)*c(2,ij,ik)
        xlife(ipf,4,ij,ik) = K(75,IJ,IK)*C(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(182,IJ,IK)*c(27,ij,ik)
        xlife(ipf,6,ij,ik) = CN(40,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(26,ijj0,ij,ik)
        end do


C   *********************************************************************************************
C
c
c  CH3Br   -   CN(49) 
    

	ch3brp=0.0


	ch3brl = J(29,ij,ik) + K(97,IJ,IK)*c(13,IJ,IK)       
     >                       + K(124,IJ,IK)*c(2,IJ,IK)
     >                       + K(176,IJ,IK)*c(27,ij,ik)


	cn(49,IJ,IK) = (c(49,IJ,IK) + ch3brp*dt)/(1. + ch3brl*dt)

C PUT IN MIXING RATIO BOUNDARY CONDITIONS
      IF (IK .EQ. 1) THEN
       IF (LBCMRSS(49).OR.LBCMRTD(49)) CN(49,IJ,IK)=BVAL(49,IJ)*M(IJ,IK)
ccccccc            CN(49,IJ,IK) = .0072*1.e-9*M(IJ,IK)
      ENDIF


        CTPROD(16,IJ,IK) = ch3brp/m(ij,ik)
        CTLOSS(16,IJ,IK) = CH3BRL*c(49,ij,ik)/m(ij,ik)
        RLOSS(10,IJ,IK)=CH3BRL


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 19

        xlife(ipf,1,ij,ik) = CH3BRL
        xlife(ipf,2,ij,ik) = J(29,IJ,IK)
        xlife(ipf,3,ij,ik) = K(124,IJ,IK)*c(2,IJ,IK)
        xlife(ipf,4,ij,ik) = K(97,IJ,IK)*c(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(176,IJ,IK)*c(27,IJ,IK)
        xlife(ipf,6,ij,ik) = cn(49,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(29,ijj0,ij,ik)
        end do



C   *********************************************************************************************
C
C
c  CBrF3  -     CN(50)    -  Halon-1301

	cbrf3p=0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF (IK .EQ. 1 )THEN
          IF (.NOT.LBCMRSS(50) .OR. .NOT.LBCMRTD(50))
     >        CBRF3P = BVAL(50,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	cbrf3l = J(30,ij,ik) + C(2,IJ,IK)*K(113,IJ,IK)
     >                       + K(242,IJ,IK)*c(13,IJ,IK)
     >                       + K(266,ij,ik)*c(27,ij,ik)



	cn(50,IJ,IK) = (c(50,IJ,IK) + cbrf3p*dt)/(1. + cbrf3l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
      IF (IK .EQ. 1) THEN
       IF (LBCMRSS(50).OR.LBCMRTD(50)) CN(50,IJ,IK)=BVAL(50,IJ)*M(IJ,IK)
ccccccc             CN(50,IJ,IK) = .0032*1.e-9*M(IJ,IK)
      ENDIF


        CTPROD(22,IJ,IK) = cbrf3p/m(ij,ik)
        CTLOSS(22,IJ,IK) = CBRF3L*c(50,ij,ik)/m(ij,ik)
        RLOSS(16,IJ,IK)=CBRF3L

C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 9

        xlife(ipf,1,ij,ik) = cbrf3l
        xlife(ipf,2,ij,ik) = J(30,IJ,IK)
        xlife(ipf,3,ij,ik) = K(113,IJ,IK)*c(2,IJ,IK)
        xlife(ipf,4,ij,ik) = K(242,IJ,IK)*c(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(266,ij,ik)*c(27,ij,ik)
        xlife(ipf,6,ij,ik) = cn(50,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(30,ijj0,ij,ik)
        end do



C   *********************************************************************************************
C
C
c  CBrClF2   -    CN(51)   -   Halon-1211
C

	cbrclf2p = 0.0


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF (IK .EQ. 1) THEN 
          IF (.NOT.LBCMRSS(51) .OR. .NOT.LBCMRTD(51))
     >         CBRCLF2P = BVAL(51,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


        cbrclf2l = J(31,ij,ik) + C(2,IJ,IK)*K(112,IJ,IK)
     >                         + K(243,IJ,IK)*c(13,IJ,IK)
     >                         + K(265,ij,ik)*c(27,ij,ik)


	cn(51,IJ,IK) = (c(51,IJ,IK) + cbrclf2p*dt)/(1. + cbrclf2l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF (IK .EQ. 1) THEN
         IF(LBCMRSS(51).OR.LBCMRTD(51))CN(51,IJ,IK)=BVAL(51,IJ)*M(IJ,IK)
ccccccc              CN(51,IJ,IK) = .00407*1.e-9*M(IJ,IK)
	ENDIF


        CTPROD(21,IJ,IK) = cbrclf2p/m(ij,ik)
        CTLOSS(21,IJ,IK) = CBRCLF2L*c(51,ij,ik)/m(ij,ik)
        RLOSS(15,IJ,IK) = CBRCLF2L


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 8

        xlife(ipf,1,ij,ik) = cbrclf2l
        xlife(ipf,2,ij,ik) = J(31,IJ,IK)
        xlife(ipf,3,ij,ik) = K(112,IJ,IK)*c(2,IJ,IK)
        xlife(ipf,4,ij,ik) = K(243,IJ,IK)*c(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(265,ij,ik)*c(27,ij,ik)
        xlife(ipf,6,ij,ik) = cn(51,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(31,ijj0,ij,ik)
        end do



C   *********************************************************************************************
C
C
c  CHClF2    -    CN(52)  -  HCFC-22


	chclf2p=0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF (IK .EQ. 1) THEN
 	IF(.NOT.LBCMRSS(52) .OR. .NOT.LBCMRTD(52))
     *     CHCLF2P = BVAL(52,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	chclf2l = J(32,IJ,IK) + K(98,IJ,IK)*c(13,IJ,IK)
     >                        + K(125,IJ,IK)*c(2,IJ,IK)
     >                        + K(181,IJ,IK)*c(27,ij,ik)


	cn(52,IJ,IK) = (c(52,IJ,IK) + chclf2p*dt)/(1. + chclf2l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
      IF (IK .EQ. 1) THEN
       IF (LBCMRSS(52).OR.LBCMRTD(52)) CN(52,IJ,IK)=BVAL(52,IJ)*M(IJ,IK)
ccccccc             CN(52,IJ,IK) = .20680*1.e-9*M(IJ,IK)
      ENDIF


        CTPROD(17,IJ,IK) = chclf2p/m(ij,ik)
        CTLOSS(17,IJ,IK) = CHCLF2L*c(52,ij,ik)/m(ij,ik)
        RLOSS(11,IJ,IK)=CHCLF2L


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 5

        xlife(ipf,1,ij,ik) = chclf2l
        xlife(ipf,2,ij,ik) = J(32,IJ,IK)
        xlife(ipf,3,ij,ik) = K(125,IJ,IK)*c(2,IJ,IK)
        xlife(ipf,4,ij,ik) = K(98,IJ,IK)*c(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(181,IJ,IK)*c(27,ij,ik)
        xlife(ipf,6,ij,ik) = cn(52,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(32,ijj0,ij,ik)
        end do


C   *********************************************************************************************


c  C2Cl3F3     -   CN(53)


	c2cl3f3p=0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
	  IF(.NOT.LBCMRSS(53) .OR. .NOT.LBCMRTD(53))
     >       C2CL3F3P=BVAL(53,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	c2cl3F3l = J(33,IJ,IK) + C(2,IJ,IK)*K(99,IJ,IK)
     >                         + K(238,IJ,IK)*c(13,IJ,IK)
     >                         + K(262,IJ,IK)*c(27,IJ,IK)


	cn(53,IJ,IK) = (c(53,IJ,IK) + c2cl3f3p*dt)/(1. + c2cl3f3l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
      IF (IK.EQ.1) THEN
        IF(LBCMRSS(53).OR.LBCMRTD(53)) CN(53,IJ,IK)=BVAL(53,IJ)*M(IJ,IK)
ccccccc              CN(53,IJ,IK) = .0756*1.e-9*M(IJ,IK)
      ENDIF


        CTPROD(18,IJ,IK) = C2CL3F3p/m(ij,ik)
        CTLOSS(18,IJ,IK) = C2CL3F3L*c(53,ij,ik)/m(ij,ik)
        RLOSS(12,IJ,IK)=C2CL3F3L


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 10

        xlife(ipf,1,ij,ik) = c2cl3F3l
        xlife(ipf,2,ij,ik) = J(33,IJ,IK)
        xlife(ipf,3,ij,ik) = K(99,IJ,IK)*c(2,IJ,IK)
        xlife(ipf,4,ij,ik) = K(238,IJ,IK)*c(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(262,IJ,IK)*c(27,IJ,IK)
        xlife(ipf,6,ij,ik) = cn(53,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(33,ijj0,ij,ik)
        end do



C   *********************************************************************************************
C
C
c  C2Cl2F4     -   cn(54)


	c2cl2f4p=0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
	  IF(.NOT.LBCMRSS(54) .OR. .NOT.LBCMRTD(54))
     >       C2CL2F4P = BVAL(54,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	c2cl2f4l = J(34,ij,ik) + C(2,IJ,IK)*K(100,IJ,IK)
     >                         + K(239,IJ,IK)*c(13,IJ,IK)
     >                         + K(263,IJ,IK)*c(27,IJ,IK)


	cn(54,IJ,IK) = (c(54,IJ,IK) + c2cl2f4p*dt)/(1. + c2cl2f4l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
       IF (IK.EQ.1) THEN
        IF(LBCMRSS(54).OR.LBCMRTD(54)) CN(54,IJ,IK)=BVAL(54,IJ)*M(IJ,IK)
ccccccc             CN(54,IJ,IK) = .0164*1.e-9*M(IJ,IK)
       ENDIF


        CTPROD(19,IJ,IK) = c2cl2f4p/m(ij,ik)
        CTLOSS(19,IJ,IK) = C2CL2F4L*c(54,ij,ik)/m(ij,ik)
        RLOSS(13,IJ,IK)=C2CL2F4L


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 15

        xlife(ipf,1,ij,ik) = c2cl2f4l
        xlife(ipf,2,ij,ik) = J(34,IJ,IK)
        xlife(ipf,3,ij,ik) = K(100,IJ,IK)*c(2,IJ,IK)
        xlife(ipf,4,ij,ik) = K(239,IJ,IK)*c(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(263,IJ,IK)*c(27,IJ,IK)
        xlife(ipf,6,ij,ik) = cn(54,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(34,ijj0,ij,ik)
        end do


C   *********************************************************************************************
C
C
c    C2ClF5  -   CN(55)


	c2clf5p=0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
	IF(.NOT.LBCMRSS(55) .OR. .NOT.LBCMRTD(55))
     >    C2CLF5P = BVAL(55,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	c2clf5l = J(35,IJ,IK) + C(2,IJ,IK)*K(101,IJ,IK)
     >                        + K(240,IJ,IK)*c(13,IJ,IK)
     >                        + K(264,ij,ik)*c(27,ij,ik)


	cn(55,IJ,IK) = (c(55,IJ,IK) + c2clf5p*dt)/(1. + c2clf5l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
       IF (IK.EQ.1) THEN
        IF(LBCMRSS(55).OR.LBCMRTD(55)) CN(55,IJ,IK)=BVAL(55,IJ)*M(IJ,IK)
ccccccc                CN(55,IJ,IK) = .0084*1.e-9*M(IJ,IK)
       ENDIF


        CTPROD(20,IJ,IK) = c2clf5p/m(ij,ik)
        CTLOSS(20,IJ,IK) = C2CLF5L*c(55,ij,ik)/m(ij,ik)
        RLOSS(14,IJ,IK)=C2CLF5L


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 11

        xlife(ipf,1,ij,ik) = c2clf5l
        xlife(ipf,2,ij,ik) = J(35,IJ,IK)
        xlife(ipf,3,ij,ik) = K(101,IJ,IK)*c(2,IJ,IK)
        xlife(ipf,4,ij,ik) = K(240,IJ,IK)*c(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(264,ij,ik)*c(27,ij,ik)
        xlife(ipf,6,ij,ik) = cn(55,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(35,ijj0,ij,ik)
        end do


C   *********************************************************************************************


c  HCFC-141b  -   CN(69)


	c141p = 0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
          IF(.NOT.LBCMRSS(69) .OR. .NOT.LBCMRTD(69))
     >      C141P = BVAL(69,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


        c141l = J(65,ij,ik) + K(114,IJ,IK)*c(2,IJ,IK) 
     >                      + K(115,IJ,IK)*c(13,IJ,IK)
     >                      + K(116,IJ,IK)*c(27,IJ,IK) 


	cn(69,IJ,IK) = (c(69,IJ,IK) + c141p*dt)/(1. + c141l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF(IK.EQ.1)THEN
          IF(LBCMRSS(69).OR.LBCMRTD(69))
     >      CN(69,IJ,IK) = BVAL(69,IJ)*M(IJ,IK)
ccccccc             CN(69,IJ,IK) = .0203*1.e-9*M(IJ,IK)
	ENDIF


        CTPROD(29,IJ,IK) = c141p/m(ij,ik)
        CTLOSS(29,IJ,IK) = c141l*c(69,ij,ik)/m(ij,ik)
        RLOSS(18,IJ,IK) = c141l


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 16

        xlife(ipf,1,ij,ik) = c141l
        xlife(ipf,2,ij,ik) = J(65,IJ,IK)
        xlife(ipf,3,ij,ik) = K(114,IJ,IK)*c(2,IJ,IK)
        xlife(ipf,4,ij,ik) = K(115,IJ,IK)*c(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(116,IJ,IK)*c(27,IJ,IK) 
        xlife(ipf,6,ij,ik) = cn(69,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(65,ijj0,ij,ik)
        end do


C   *********************************************************************************************


C HCFC-142b   -    CN(70)


	c142p = 0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
          IF(.NOT.LBCMRSS(70) .OR. .NOT.LBCMRTD(70))
     >      C142P = BVAL(70,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


        c142l = J(66,ij,ik) + K(117,IJ,IK)*c(2,IJ,IK) 
     >                      + K(118,IJ,IK)*c(13,IJ,IK) 
     >                      + K(119,IJ,IK)*c(27,IJ,IK) 


	cn(70,IJ,IK) = (c(70,IJ,IK) + c142p*dt)/(1. + c142l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF(IK.EQ.1)THEN
          IF(LBCMRSS(70).OR.LBCMRTD(70))
     >      CN(70,IJ,IK) = BVAL(70,IJ)*M(IJ,IK)
ccccccc               CN(70,IJ,IK) = .0205*1.e-9*M(IJ,IK)
	ENDIF


        CTPROD(30,IJ,IK) = c142p/m(ij,ik)
        CTLOSS(30,IJ,IK) = c142l*c(70,ij,ik)/m(ij,ik)
        RLOSS(19,IJ,IK) = c142l


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 17

        xlife(ipf,1,ij,ik) = c142l
        xlife(ipf,2,ij,ik) = J(66,IJ,IK)
        xlife(ipf,3,ij,ik) = K(117,IJ,IK)*c(2,IJ,IK)
        xlife(ipf,4,ij,ik) = K(118,IJ,IK)*c(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(119,IJ,IK)*c(27,IJ,IK) 
        xlife(ipf,6,ij,ik) = cn(70,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(66,ijj0,ij,ik)
        end do



C   *********************************************************************************************


C HCFC-123  -   CN(71)


	c123p = 0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
          IF(.NOT.LBCMRSS(71) .OR. .NOT.LBCMRTD(71))
     >      C123P = BVAL(71,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	c123l = J(67,ij,ik) + K(120,IJ,IK)*c(2,IJ,IK) 
     >                      + K(121,IJ,IK)*c(13,IJ,IK)
     >                      + K(122,ij,ik)*c(27,ij,ik)

	cn(71,IJ,IK) = (c(71,IJ,IK) + c123p*dt)/(1. + c123l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF(IK.EQ.1)THEN
          IF(LBCMRSS(71).OR.LBCMRTD(71))
     >      CN(71,IJ,IK) = BVAL(71,IJ)*M(IJ,IK)
ccccccc               CN(71,IJ,IK) = .0025*1.e-9*M(IJ,IK)
	ENDIF


        CTPROD(31,IJ,IK) = c123p/m(ij,ik)
        CTLOSS(31,IJ,IK) = c123l*c(71,ij,ik)/m(ij,ik)
        RLOSS(29,IJ,IK) = c123l


C   *********************************************************************************************


C  Halon-2402    -   CN(72)


	c2402p = 0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
          IF(.NOT.LBCMRSS(72) .OR. .NOT.LBCMRTD(72))
     >      C2402P = BVAL(72,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	c2402l = J(68,ij,ik) + K(123,IJ,IK)*c(2,IJ,IK) 
     >                       + K(244,IJ,IK)*c(13,IJ,IK)
     >                       + K(268,ij,ik)*c(27,ij,ik)


	cn(72,IJ,IK) = (c(72,IJ,IK) + c2402p*dt)/(1. + c2402l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF(IK.EQ.1)THEN
          IF(LBCMRSS(72).OR.LBCMRTD(72))
     >      CN(72,IJ,IK) = BVAL(72,IJ)*M(IJ,IK)
ccccccc                CN(72,IJ,IK) = .00046*1.e-9*M(IJ,IK)
	ENDIF


        CTPROD(32,IJ,IK) = c2402p/m(ij,ik)
        CTLOSS(32,IJ,IK) = c2402l*c(72,ij,ik)/m(ij,ik)
        RLOSS(21,IJ,IK) = c2402l


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 21

        xlife(ipf,1,ij,ik) = c2402l
        xlife(ipf,2,ij,ik) = J(68,IJ,IK)
        xlife(ipf,3,ij,ik) = K(123,IJ,IK)*c(2,IJ,IK)
        xlife(ipf,4,ij,ik) = K(244,IJ,IK)*c(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(268,ij,ik)*c(27,ij,ik)
        xlife(ipf,6,ij,ik) = cn(72,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(68,ijj0,ij,ik)
        end do


C   *********************************************************************************************


C  Halon-1202  (CF2Br2)  -   CN(75)  (Set BC to .00001 ppbv)


	c1202p = 0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
          IF(.NOT.LBCMRSS(75) .OR. .NOT.LBCMRTD(75))
     >      C1202P = BVAL(75,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	c1202l = J(54,ij,ik) + K(174,IJ,IK)*c(2,ij,ik)
     >                       + K(241,IJ,IK)*c(13,IJ,IK)
     >                       + K(267,ij,ik)*c(27,ij,ik)


	cn(75,IJ,IK) = (c(75,IJ,IK) + c1202p*dt)/(1. + c1202l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF(IK.EQ.1)THEN
          IF(LBCMRSS(75).OR.LBCMRTD(75))
     >      CN(75,IJ,IK) = BVAL(75,IJ)*M(IJ,IK)
ccccccc                CN(75,IJ,IK) = .00001*1.e-9*M(IJ,IK)
	ENDIF


        RLOSS(22,IJ,IK) = c1202l           


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 20

        xlife(ipf,1,ij,ik) = c1202l
        xlife(ipf,2,ij,ik) = J(54,IJ,IK)
        xlife(ipf,3,ij,ik) = K(174,IJ,IK)*c(2,IJ,IK)
        xlife(ipf,4,ij,ik) = K(241,IJ,IK)*c(13,IJ,IK)
        xlife(ipf,5,ij,ik) = K(267,ij,ik)*c(27,ij,ik)
        xlife(ipf,6,ij,ik) = cn(75,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(54,ijj0,ij,ik)
        end do



C   *********************************************************************************************


C  Dibromomethane (Methylene Bromide)   (CH2Br2)  -   CN(76)


	cmebrp = 0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
          IF(.NOT.LBCMRSS(76) .OR. .NOT.LBCMRTD(76))
     >      cmebrp = BVAL(76,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	cmebrl = J(69,ij,ik) + K(177,IJ,IK)*c(13,ij,ik)
     >                       + K(179,IJ,IK)*c(27,ij,ik)
     >                       + K(183,IJ,IK)*c(2,ij,ik)


	cn(76,IJ,IK) = (c(76,IJ,IK) + cmebrp*dt)/(1. + cmebrl*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF(IK.EQ.1)THEN
ccccccc          IF(LBCMRSS(76).OR.LBCMRTD(76))
ccccccc     >      CN(76,IJ,IK) = BVAL(76,IJ)*M(IJ,IK)
               CN(76,IJ,IK) = .0000*1.e-9*M(IJ,IK)
	ENDIF


        RLOSS(30,IJ,IK) = cmebrl


C   *********************************************************************************************

C  Bromoform   (CHBr3)  -   CN(77)


	cbrmfp = 0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
          IF(.NOT.LBCMRSS(77) .OR. .NOT.LBCMRTD(77))
     >      cbrmfp = BVAL(77,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	cbrmfl = J(70,ij,ik) + K(178,IJ,IK)*c(13,ij,ik)
     >                       + K(180,IJ,IK)*c(27,ij,ik)
     >                       + K(184,IJ,IK)*c(2,ij,ik)


	cn(77,IJ,IK) = (c(77,IJ,IK) + cbrmfp*dt)/(1. + cbrmfl*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF(IK.EQ.1)THEN
ccccccc          IF(LBCMRSS(77).OR.LBCMRTD(77))
ccccccc     >      CN(77,IJ,IK) = BVAL(77,IJ)*M(IJ,IK)
              CN(77,IJ,IK) = .0000*1.e-9*M(IJ,IK)
	ENDIF


        RLOSS(31,IJ,IK) = cbrmfl



C *********************************************************************************************
C
C     HFCs:
C
C *********************************************************************************************
C
C
C  HFC-134a    (CH2FCF3)  -   CN(81)


	ch134p = 0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
          IF(.NOT.LBCMRSS(81) .OR. .NOT.LBCMRTD(81))
     >      ch134p = BVAL(81,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	ch134l = J(71,IJ,IK) + K(215,IJ,IK)*c(2,ij,ik)
     >                       + K(249,IJ,IK)*c(13,ij,ik)
     >                       + K(257,IJ,IK)*c(27,ij,ik)


	cn(81,IJ,IK) = (c(81,IJ,IK) + ch134p*dt)/(1. + ch134l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF(IK.EQ.1)THEN
        IF(LBCMRSS(81).OR.LBCMRTD(81)) CN(81,IJ,IK)=BVAL(81,IJ)*M(IJ,IK)
ccccccc           CN(81,IJ,IK) = 59.45*1.e-12*M(IJ,IK)
	ENDIF


        RLOSS(25,IJ,IK) = ch134l


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 12

        xlife(ipf,1,ij,ik) = ch134l
        xlife(ipf,2,ij,ik) = J(71,IJ,IK)
        xlife(ipf,3,ij,ik) = K(215,IJ,IK)*c(2,ij,ik)
        xlife(ipf,4,ij,ik) = K(249,IJ,IK)*c(13,ij,ik)
        xlife(ipf,5,ij,ik) = K(257,IJ,IK)*c(27,ij,ik)
        xlife(ipf,6,ij,ik) = cn(81,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(71,ijj0,ij,ik)
        end do


C   *********************************************************************************************
C
C  HFC-143a    (CF3CH3)  -   CN(82)


	ch143p = 0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
          IF(.NOT.LBCMRSS(82) .OR. .NOT.LBCMRTD(82))
     >      ch143p = BVAL(82,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	ch143l = J(72,IJ,IK) + K(224,IJ,IK)*c(2,ij,ik)
     >                       + K(248,IJ,IK)*c(13,ij,ik)
     >                       + K(256,IJ,IK)*c(27,ij,ik)


	cn(82,IJ,IK) = (c(82,IJ,IK) + ch143p*dt)/(1. + ch143l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF(IK.EQ.1)THEN
        IF(LBCMRSS(82).OR.LBCMRTD(82)) CN(82,IJ,IK)=BVAL(82,IJ)*M(IJ,IK)
ccccccc           CN(82,IJ,IK) = 10.3*1.e-12*M(IJ,IK)
	ENDIF


        RLOSS(26,IJ,IK) = ch143l


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 13

        xlife(ipf,1,ij,ik) = ch143l
        xlife(ipf,2,ij,ik) = J(72,IJ,IK)
        xlife(ipf,3,ij,ik) = K(224,IJ,IK)*c(2,ij,ik)
        xlife(ipf,4,ij,ik) = K(248,IJ,IK)*c(13,ij,ik)
        xlife(ipf,5,ij,ik) = K(256,IJ,IK)*c(27,ij,ik)
        xlife(ipf,6,ij,ik) = cn(82,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(72,ijj0,ij,ik)
        end do


C   *********************************************************************************************
C
C   HFC-23    (CHF3)  -   CN(83)


	ch23p = 0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
          IF(.NOT.LBCMRSS(83) .OR. .NOT.LBCMRTD(83))
     >      ch23p = BVAL(83,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	ch23l = J(73,IJ,IK) + K(203,IJ,IK)*c(2,ij,ik)
     >                      + K(246,IJ,IK)*c(13,ij,ik)
     >                      + K(254,IJ,IK)*c(27,ij,ik)


	cn(83,IJ,IK) = (c(83,IJ,IK) + ch23p*dt)/(1. + ch23l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF(IK.EQ.1)THEN
        IF(LBCMRSS(83).OR.LBCMRTD(83)) CN(83,IJ,IK)=BVAL(83,IJ)*M(IJ,IK)
ccccccc           CN(83,IJ,IK) = 22.95*1.e-12*M(IJ,IK)
	ENDIF


cccccc        RLOSS(27,IJ,IK) = ch23l


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 14

        xlife(ipf,1,ij,ik) = ch23l
        xlife(ipf,2,ij,ik) = J(73,IJ,IK)
        xlife(ipf,3,ij,ik) = K(203,IJ,IK)*c(2,ij,ik)
        xlife(ipf,4,ij,ik) = K(246,IJ,IK)*c(13,ij,ik)
        xlife(ipf,5,ij,ik) = K(254,IJ,IK)*c(27,ij,ik)
        xlife(ipf,6,ij,ik) = cn(83,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(73,ijj0,ij,ik)
        end do


C   *********************************************************************************************
C
C   HFC-32    (CH2F2)  -   CN(84)


	ch32p = 0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
          IF(.NOT.LBCMRSS(84) .OR. .NOT.LBCMRTD(84))
     >      ch32p = BVAL(84,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	ch32l = J(74,ij,ik) + K(201,IJ,IK)*c(2,ij,ik)
     >                      + K(245,IJ,IK)*c(13,ij,ik)
     >                      + K(253,IJ,IK)*c(27,ij,ik)


	cn(84,IJ,IK) = (c(84,IJ,IK) + ch32p*dt)/(1. + ch32l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF(IK.EQ.1)THEN
        IF(LBCMRSS(84).OR.LBCMRTD(84)) CN(84,IJ,IK)=BVAL(84,IJ)*M(IJ,IK)
ccccccc           CN(84,IJ,IK) = 5.75*1.e-12*M(IJ,IK)
	ENDIF


ccccccc        RLOSS(28,IJ,IK) = ch32l


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 22

        xlife(ipf,1,ij,ik) = ch32l
        xlife(ipf,2,ij,ik) = J(74,IJ,IK)
        xlife(ipf,3,ij,ik) = K(201,IJ,IK)*c(2,ij,ik)
        xlife(ipf,4,ij,ik) = K(245,IJ,IK)*c(13,ij,ik)
        xlife(ipf,5,ij,ik) = K(253,IJ,IK)*c(27,ij,ik)
        xlife(ipf,6,ij,ik) = cn(84,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(74,ijj0,ij,ik)
        end do



C   *********************************************************************************************
C
C   HFC-125   (CHF2CF3)  -   CN(85)


	ch125p = 0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
          IF(.NOT.LBCMRSS(85) .OR. .NOT.LBCMRTD(85))
     >      ch125p = BVAL(85,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	ch125l = J(75,IJ,IK) + K(218,IJ,IK)*c(2,ij,ik)
     >                       + K(250,IJ,IK)*c(13,ij,ik)
     >                       + K(258,IJ,IK)*c(27,ij,ik)


	cn(85,IJ,IK) = (c(85,IJ,IK) + ch125p*dt)/(1. + ch125l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF(IK.EQ.1)THEN
        IF(LBCMRSS(85).OR.LBCMRTD(85)) CN(85,IJ,IK)=BVAL(85,IJ)*M(IJ,IK)
ccccccc           CN(85,IJ,IK) = 8.8*1.e-12*M(IJ,IK)
	ENDIF


ccccccc        RLOSS(29,IJ,IK) = ch125l


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 23

        xlife(ipf,1,ij,ik) = ch125l
        xlife(ipf,2,ij,ik) = J(75,IJ,IK)
        xlife(ipf,3,ij,ik) = K(218,IJ,IK)*c(2,ij,ik)
        xlife(ipf,4,ij,ik) = K(250,IJ,IK)*c(13,ij,ik)
        xlife(ipf,5,ij,ik) = K(258,IJ,IK)*c(27,ij,ik)
        xlife(ipf,6,ij,ik) = cn(85,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(75,ijj0,ij,ik)
        end do



C   *********************************************************************************************
C
C   HFC-152a  (CH3CHF2)  -   CN(86)


	ch152p = 0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
          IF(.NOT.LBCMRSS(86) .OR. .NOT.LBCMRTD(86))
     >      ch152p = BVAL(86,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	ch152l = J(76,IJ,IK) + K(211,IJ,IK)*c(2,ij,ik)
     >                       + K(247,IJ,IK)*c(13,ij,ik)
     >                       + K(255,IJ,IK)*c(27,ij,ik)


	cn(86,IJ,IK) = (c(86,IJ,IK) + ch152p*dt)/(1. + ch152l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF(IK.EQ.1)THEN
        IF(LBCMRSS(86).OR.LBCMRTD(86)) CN(86,IJ,IK)=BVAL(86,IJ)*M(IJ,IK)
ccccccc           CN(86,IJ,IK) = 6.75*1.e-12*M(IJ,IK)
	ENDIF


ccccccc        RLOSS(30,IJ,IK) = ch152l


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 24

        xlife(ipf,1,ij,ik) = ch152l
        xlife(ipf,2,ij,ik) = J(76,IJ,IK)
        xlife(ipf,3,ij,ik) = K(211,IJ,IK)*c(2,ij,ik)
        xlife(ipf,4,ij,ik) = K(247,IJ,IK)*c(13,ij,ik)
        xlife(ipf,5,ij,ik) = K(255,IJ,IK)*c(27,ij,ik)
        xlife(ipf,6,ij,ik) = cn(86,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(76,ijj0,ij,ik)
        end do



C   *********************************************************************************************
C
C   HFC-227ea  (CF3CHFCF3)  -   CN(87)


	ch227p = 0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
          IF(.NOT.LBCMRSS(87) .OR. .NOT.LBCMRTD(87))
     >      ch227p = BVAL(87,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	ch227l = J(77,IJ,IK) + K(233,IJ,IK)*c(2,ij,ik)
     >                       + K(252,IJ,IK)*c(13,ij,ik)
     >                       + K(269,IJ,IK)*c(27,ij,ik)


	cn(87,IJ,IK) = (c(87,IJ,IK) + ch227p*dt)/(1. + ch227l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF(IK.EQ.1)THEN
        IF(LBCMRSS(87).OR.LBCMRTD(87)) CN(87,IJ,IK)=BVAL(87,IJ)*M(IJ,IK)
ccccccc           CN(87,IJ,IK) = 2.*1.e-12*M(IJ,IK)
	ENDIF


ccccccc        RLOSS(31,IJ,IK) = ch227l


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 25

        xlife(ipf,1,ij,ik) = ch227l
        xlife(ipf,2,ij,ik) = J(77,IJ,IK)
        xlife(ipf,3,ij,ik) = K(233,IJ,IK)*c(2,ij,ik)
        xlife(ipf,4,ij,ik) = K(252,IJ,IK)*c(13,ij,ik)
        xlife(ipf,5,ij,ik) = K(269,IJ,IK)*c(27,ij,ik)
        xlife(ipf,6,ij,ik) = cn(87,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(77,ijj0,ij,ik)
        end do



C   *********************************************************************************************
C
C   HFC-245fa  (CHF2CH2CF3)  -   CN(88)


	ch245p = 0.


C PUT IN FLUX BOUNDARY CONDITIONS  
	IF(IK .EQ. 1)THEN
          IF(.NOT.LBCMRSS(88) .OR. .NOT.LBCMRTD(88))
     >      ch245p = BVAL(88,IJ)/DELTAZ(IJ,IK)/1.E5
	ENDIF


	ch245l = J(78,IJ,IK) + K(235,IJ,IK)*c(2,ij,ik)
     >                       + K(251,IJ,IK)*c(13,ij,ik)
     >                       + K(270,IJ,IK)*c(27,ij,ik)


	cn(88,IJ,IK) = (c(88,IJ,IK) + ch245p*dt)/(1. + ch245l*dt)


C PUT IN MIXING RATIO BOUNDARY CONDITIONS
	IF(IK.EQ.1)THEN
        IF(LBCMRSS(88).OR.LBCMRTD(88)) CN(88,IJ,IK)=BVAL(88,IJ)*M(IJ,IK)
ccccccc           CN(88,IJ,IK) = 1.7*1.e-12*M(IJ,IK)
	ENDIF


ccccccc        RLOSS(32,IJ,IK) = ch245l


C
C   xlife(30,6,L$,Z$);  2nd index (loss):   1=J (tot);  2=O1D;   3=OH;   4=Cl;   5=number density
C
        ipf = 26

        xlife(ipf,1,ij,ik) = ch245l
        xlife(ipf,2,ij,ik) = J(78,IJ,IK)
        xlife(ipf,3,ij,ik) = K(235,IJ,IK)*c(2,ij,ik)
        xlife(ipf,4,ij,ik) = K(251,IJ,IK)*c(13,ij,ik)
        xlife(ipf,5,ij,ik) = K(270,IJ,IK)*c(27,ij,ik)
        xlife(ipf,6,ij,ik) = cn(88,IJ,IK)


C   load in wavelength dependent J's into xlifej(30,5,L$,Z$), J39(PH$,5,L$,Z$)

        do ijj0=1,5
           xlifej(ipf,ijj0,ij,ik) = J39(78,ijj0,ij,ik)
        end do



C
C  load in extra constituents into array : xlifec(5,L$,Z$);   TEMP(L$,Z$X) is in COMMON
C
C     1=O1D ;    2=OH ;   3=Cl ;   4=M ;   5=Temperature
C

        xlifec(1,ij,ik) = c(2,ij,ik)
        xlifec(2,ij,ik) = c(13,ij,ik)
        xlifec(3,ij,ik) = c(27,ij,ik)
        xlifec(4,ij,ik) = m(ij,ik)
        xlifec(5,ij,ik) = temp(ij,ik)



C  ***************************  SF6 --  AGE    ***************************************
C
C  for SF6,  no chemistry, just add in production in levels 1-3 (doesn't matter what res), at 35N-55N
c  Use Malcolm Ko's formula for SF6 emission per day (sf6a), key on IYR, where IYR = 0 = 1966 (IYR-31 HERE)
C  NOTE: Malcolm's formula gives a rate of .1 Ggm/yr on day 1, 1966, and .3 Ggm/yr on day 1, 1967,
C  which is 6 months off, so just subtract .1 Ggm/yr from formula, or .1/360. Ggm/day for correct rate
C
C    Then, multiply by 4.771984e25 to convert from giga-grams per day to molecules/sec 
C    then divide by total volume at 35N-55N over the first 3 model levels to get 
C    sf6prod in #molecules/cm3-sec input. To get the total volume, take the area at each
C    latitude, multiply by the total deltaz (im cm) over the first 3 levels at that latitude,
C    and integrate this over 35N-55N (ij35-ij55). I know this is probably more detail than is 
C    needed, but we'll do it to be consistent with the real atmosphere
C    
cccc      if (ik .le. 3  .and.  ij .ge. 12  .and.  ij .le. 15) then               ! put in at 25N-55N - CHJ
cccc        sf6a = .2*(IYR+1)/360. + .2*(iday360 - 361./2)/129600. - .1/360.  ! Malcolm Ko's modified formula
cccc        sf6a = .2*(IYR+1)/360. + .2*(iday360 - 361./2)/129600.            ! Malcolm Ko's original formula
cccc            sf6a = (.2*IYR + .1)/360.                                         ! simple formula
cccc          totvol = 1.336e18*DELTAZ(IJ,IK)*1.E5*3.                   ! Charley's total volume for 25N-60N

C  ij35, ij55 in COMMON, AREA(L$) is the cosine of latitude divided by the total of the cosine of all lats

cage        sf6loss = 0.0E0
cage        sf6prod = 0.0E0

cage      if (iyr .ge. 31) then
cage      if (ik .le. 3  .and.  ij .ge. ij35  .and.  ij .le. ij55) then    ! put in at 35N-55N at lowest 3 levels
                                                                           ! no matter what the resolution
cage      sf6a= .2*(IYR-31+1)/360. + .2*(iday360 - 361./2)/129600. - .1/360.  ! Malcolm Ko's modified formula
C                                                                         ! deltaz is in KM
cage         totvol=0.0
cage         do 447 ijs=ij35,ij55
cage           totvol = totvol + area(ijs)*globe*
cage     >              (deltaz(ijs,1) + deltaz(ijs,2) + deltaz(ijs,3))*1.e5 
cage 447     CONTINUE

cage         sf6prod = sf6a*4.771984e25/totvol
cage      endif
cage      endif
        
cage      cn(77,ij,ik) = (c(77,ij,ik) + sf6prod*dt)/(1. + sf6loss*dt) 
C
C
C   ***************************  END   SF6/AGE  *****************************************
C
C
C  ***************************  Simple time-increasing  AGE      CN(78)  *********************************
C
C  for time increasing mixing ratio at ground, CN(78), just set equal to year, M(L$,Z$58)
C    starting at 1960, day 1 = 1960.0
C
      prod78 = 0.0
      loss78 = 0.0
      cn(78,ij,ik) = (c(78,ij,ik) + prod78*dt)/(1. + loss78*dt) 

      if (ik .eq. 1) CN(78,IJ,IK)=(iyr+1935.+(iday360-1.)/360.)*M(IJ,IK)

C
C   ***************************  END   simple  AGE  *****************************************
C
C
C 
C  ***********  CARBON-14  - CN(21)  adapted from Charley's solvgah.f (same as solvgaj.f) *************
c
C  ******   Carbon-14 initialized in 1992 (INTERANNUAL transport) and 2012 for CLIMATOLOGICAL TRANSPORT

        C14LOSS=0.0E0
        C14PROD=0.0E0

        CN(21,IJ,IK) = (C(21,IJ,IK) + c14prod*dt)/(1. + c14loss*dt) 

c                                                          month counter = 0 for Oct. 1963
CCCC        TMON = IYR*12. + (DAY360 + 14.5)/30. - 10.          

        TMON = TMONC14
c                                    do Time-dependent BC's,   LAT(L$)
        IF(IK .EQ. 1)THEN
C                                 first FOR SOUTHERN HEMISPHERE, set = NH  FOR ANYTIME PAST JUNE 15, 1968
          IF(LAT(IJ) .LT. 0.) THEN 
               IF (TMON .LE. 56.) BCSH = 44.5 + 1.02535*TMON 
     *             - 2.13565E-2*TMON*TMON + 8.61853E-5*TMON*TMON*TMON
C                                                 
               IF (TMON .GT. 56.) BCSH = 73.0 - 0.27823*TMON 
     c            - 3.45648E-3*TMON*TMON + 4.21159E-5*TMON*TMON*TMON

               cn(21,IJ,1) = BCSH*4.82E-18*m(ij,1)
          ENDIF

C                                               FOR NORTHERN HEMISPHERE
         IF(LAT(IJ) .GE. 0.) THEN
            BCNH = 73.0 - 0.27823*TMON - 3.45648E-3*TMON*TMON
     *              + 4.21159E-5*TMON*TMON*TMON

            cn(21,IJ,1) = BCNH*4.82E-18*m(ij,1)
         ENDIF

       ENDIF

C   ***************************   END CARBON-14  *****************************************
C
C
c  Set limit to 1.E-12 number density -- problem with small NO (E-16), large CH3O2 (E+22)
C    at 85S, 212 mb, starting sometime before Nov 20 on 8th year of run
C    (m at 115 km ~1.e12, and ~3.e19 at the ground, so the minimum mixing ratio will be ~3.E-31)
C    BUT DON'T ADJUST RAINOUT PARAMETERS
C
       DO 7382 III=1,S$
          if (iii .ne. 58  .and.  iii .ne. 59) then
             IF (CN(III,IJ,IK) .LT. 1.D-12) CN(III,IJ,IK) = 1.D-12
          endif
7382   CONTINUE


1000   CONTINUE


                               ! increment C-14 month counter 1X/day - OUTSIDE of DO 1000 loop !!
       TMONC14 = TMONC14 + 1./30.



C  if IBCMIN=1, loop through the 25 source gases, defined by ISBC(25)
C     and set MIN BC of 1e-7 ppbv (2000-3000 #/cm^3) for well-behaved lifetimes

        IF (IBCMIN .eq. 1) then

           do 775 ij=1,L$
              bcmin = 1.e-7*1.e-9*M(IJ,1)

              do 776 ijk = 1,25
                 izs = ISBC(ijk)
                 if (CN(izs,ij,1) .lt. bcmin) CN(izs,ij,1) = bcmin
 776          CONTINUE
 775       CONTINUE

        ENDIF



c  write out real*4 xlife(30,6,L$,Z$), xlifej(30,5,L$,Z$) (in COMMON),  xlifec(5,L$,Z$) 
C    every month on final year of run, 
C    based on INYEARS (number of years in run, in COMMON) and IYRCT

       if (IYRCT .eq. INYEARS) then
         if (mod(iday360-15,30) .eq. 0.0) then
            write (2277) iyr, iday360, L$, Z$, Z$X
            write (2277) LAT4, zalt90, pres90, xlife, xlifej, xlifec
         endif
       endif


      RETURN
      END
