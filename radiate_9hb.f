C
c  Radiation Model, Main driver
C
C  this uses the GEOS5 radiation routines adopted from Code obtained from Dude
C
c  Program runs for COUPLED MODEL dynamics lat grid, 
C        uses COUPLED model altitude grid for heating rates (TOP - DOWN)
C
C                                                 SOLDEC is in DEGREES, IYEARC (15=1950)

      SUBROUTINE RADIATE_G5(yp, zp, ypp, zpp, TMOD0, TEMPBX,ychem,zchem,
     > SOLDEC, IDOY360, IYEARC, dayl, TTOTHX, VNC, THGX, HEATG5, COOLG5)


      INCLUDE 'PARAM.INC'


C   number of soundings (MSR), number of aerosols (NAERO), number of sub-grid types (NSURF), for RAD CODE
C                       R$ = number of LAYERS for RADIATION CODE (IRRAD/SORAD)
      INTEGER R$

      PARAMETER (MSR=1)
      PARAMETER (NAERO=0)
      PARAMETER (NSURF=1)
      PARAMETER (R$=56)


C   COMMON for CONSTITS -mixing ratio for current day, and profs, and HEATING, COOLING RATES in CHEM GRID
C       in PARAM:  NS$=45;  MS$=76  (chemistry grid)

C
      COMMON/CCHEM/CNC(S$C,NS$,MS$), HEATCHM(NS$,MS$), COOLCHM(NS$,MS$),
     >             HEATCHML(NS$,MS$), DRAGCHM(7,NS$,MS$)


C  common of current day's NCEP clim temps, TNCEPSFC are temps at 0-1 km (NOT USED)
C     TD GEOS4 temps at SFC, 1, 2, km
C     also temp, ubar at lowest extended grid point (1.015 km) - UT1KM(NP$,2) - 1=TEMP;  2=UBAR

      COMMON/CTNCEP/TNCEP(N$,M$), TNCEPSFC(NP$,2), TEMPG42D(NP$,3), 
     >              UT1KM(NP$,2), U3KM(NP$,2)



c  common of diurnally avged heating rates computed from the photolysis (from previous day), in K/sec
C   SORADHEAT(15,NS$,MS$) are the individual bands from SORAD.f interpolated to the chemistry grid for output in MAIN
c                                                                     NS$=L$  MS$=Z$ (in PARAM.INC)
      COMMON/CPHEAT/PHEAT(5,NS$,MS$), SORADHEAT(15,NS$,MS$)


c  common for coupled model sfc albedo (UV-VIS) from TOMS (season-lat), interpolated to current day in XCOUPLED
C     ALBDAY(180) is albedo (fraction), ALBLAT are the 180 latitudes

      COMMON/CTALBEDO/ TALBEDO(15,180), ALBDAY(180), ALBLAT(180)


C   common for latent heating, added here at the end
C 
      COMMON/CLATH/GLHDAYC(N$,M$), WLHDAYC(N$,M$), HACKDAYC(N$,M$),
     >             G4LHDAYC(N$,M$), LTNHT(N$,M$),  TDIFFDAY(N$,M$),
     >             G4LHDAYCX(NP$,MP$)


c   solar flux data from Judith Lean:
C   index 1 = time;  2-40 are the 39 model wavelength bins (in ph/cm^2/sec); 41=TOA flux (W/m2)
C     SJFLUXDAY(41) are TD values interpolated to the current day
C
C     HJFLUXDAY(12) are for the heating rate intervals (in mW/m2), interpolated to current day 
C
      COMMON/CSJFLUXD/ SJFLUXDAY(41), HJFLUXDAY(12), FLXA(41), ISOLCYC


  
      !INTEGER N$, M$, NP$, MP$, NS$, MS$, MSR, naero, nsurf
      INTEGER IDOY360,IYEARC
      REAL ZP(M$), ZP2(R$), ZPL(R$), ZPE(R$+1), PRMODE(MSR,R$+1)
      REAL YP(N$), ychem(NS$), zchem(MS$), dayl, prlay(R$)
      REAL ypp(NP$), zpp(MP$), ttx(MP$), ttg4(3), zzg4(3)
      REAL ttr(R$), TTOTHR(R$), VNCR(R$), THGR(R$), zzgr(2), ttgr(2)
      REAL TMOD0(N$,M$), TMODR(MSR,R$), TTOTHX(MP$), VNC(M$), THGX(MP$)
      REAL cnc1(NS$,MS$), cnc2(NP$,R$), CNCH(S$C,NP$,R$), cosz(MSR)
      REAL vmro3(R$), oa(MSR,R$), wa(MSR,R$), co2vmr(MSR,R$)
      REAL vmrch4(MSR,R$), vmrn2o(MSR,R$)
      REAL vmrcfc11(MSR,R$), vmrcfc12(MSR,R$), vmrcfc22(MSR,R$)
      REAL cwc(msr,R$,3), fcld(msr,R$), reff(msr,R$,3)
      REAL tg(msr,nsurf), tv(msr,nsurf), eg(msr,nsurf,10), tsurf(MSR)
      REAL ev(msr,nsurf,10), rv(msr,nsurf,10), fs(msr,nsurf), albs(NP$)
      REAL tempbx(NP$,MP$), tempbxs(NP$,MP$)
      REAL heatxx(N$,M$), coolxx(N$,M$)

      REAL*8 ychem8(NS$), zchem8(MS$), ypp8(NP$), zp8(R$)
      REAL*8 cnc18(NS$,MS$), cnc28(NP$,R$)

      character*40 saerosols(naero)
      REAL raero(msr,R$,naero), rhaero(msr,R$)
      REAL rsuvbm(MSR), rsuvdf(MSR), rsirbm(MSR), rsirdf(MSR)

      REAL*8 sjfluxday, hjfluxday, flxa


c-----output parameters from RADRED9 for clouds

      real cwc1(R$,3), cloudi(R$)


c-----output parameters - SORAD

      real flx(msr,R$+1), flc(msr,R$+1), flxb(msr,R$+1,15)
      real fdiruv (msr), fdifuv (msr)
      real fdirpar(msr), fdifpar(msr)
      real fdirir (msr), fdifir (msr)


c-----output parameters - IRRAD

      real flxir(msr,R$+1), flcir(msr,R$+1), flair(msr,R$+1)
      real dfdts(msr,R$+1), sfcem(msr), taudiag(msr,R$,10)

      real flxda(msr,R$+1), HEATT(R$), COOLT(R$), irflux(NP$,R$+1,3)
      real HEATG5E(NP$,R$), COOLG5E(NP$,R$), NETG5E(NP$,R$)
      real HEATG5(N$,M$), COOLG5(N$,M$), PHEATC(N$,M$)

      real flxbda(msr,R$+1,15), heatuv(NP$,R$,15), chc1(M$+1),chc2(MS$)


      logical high, trace, overcast, aerosol

      logical lprint,ascatt,aerfl
C
C
      logical aprint(3), fprint(3), solar, ir
C

      real latcld(72), camt(360,3,72), cpres(3,72), ctau(360,3,72)
      real cwcr(360,72,10,2), clfrac(360,72,10), zzcld(10)
      dimension cloudj(R$), taucld(R$), omgcld(3,R$), gcld(3,R$)
      logical cldflg

      COMMON/CLOUDS1/latcld, camt, cpres, ctau, cwcr, clfrac, zzcld
      COMMON/CLOUDS2/cldflg, taucld, omgcld, gcld, cloudj

C                                   WACCM cloud parameters for current day on radiation grid
      REAL clmr(NP$,R$,3), clfr(NP$,R$), creff(NP$,R$,3)

C                              !  GEOS 5 cloud parameters for current day on radiation grid
      REAL g5clday(5,NP$,R$)



      dtr = 3.1415926/180.
      rtd = 180./3.1415926
                              ! load in SJFLUXDAY(41), index 41 = TOA flux (W/m2) for use below (current day's value)
      sjtoa = SJFLUXDAY(41)

C
C
C  set logical arrays as defined in radiate.cfg:
C                   .f. .f. .f.                   debug print flags: radred, htrdly, irdriv - aprint(3)
C                   .f. .f. .t.                   output flags: ascii,unformatted,DF        - fprint(3)
C                   .t. .t.                       solar heating; longwave cooling           - solar, ir
C                   .t.                           cloud flag: do cloud processing           - cldflg

         aprint(1) = .false.
         aprint(2) = .false.
         aprint(3) = .false.

         fprint(1) = .false.
         fprint(2) = .false.
         fprint(3) = .true.
         lprint =    .false.

         solar  = .true.
            ir  = .true.

         cldflg   = .false.
         overcast = .false.       ! .true. = the layer cloud cover is either 0 or 1. - .true. is faster

                                  ! use table look-up in IR code, and trace gases in IR
         high  = .true.
         trace = .true.



c PW32 - NEW ALTITUDE GRID for RAD CODE, sfc pressure is 1000 mbar (0 km),
C    then use 1.015 km resolution up to 20 km (1/2 of coupled model res)
C    then 2.03 km for 22-93.38 km, as in chemistry grid
C        ZPE(R$+1=57), goes TOP-DOWN

         delzr = 2.03

         do 110 ik=1,21
 110        zpe(R$+1+1-ik) = (ik-1.)*delzr/2.

         do 112 ik=22,R$+1
 112        zpe(R$+1+1-ik) = (ik-22/2.)*delzr


C
C  get layer mid-poiint altitudes - ZPL(R$) goes TOP-DOWN, ZP2(R$) goes BOTTOM-UP
C
         do 105 ik=1,R$
            zpl(ik) = (zpe(ik) + zpe(ik+1))/2.
 105        zp2(R$+1-ik) = zpl(ik)

C                            get lowest 2 levels for relaxation to NCEP below - zzgr(2)
         do 106 ik=1,2
 106        zzgr(ik) = zp2(ik)

         do 107 ik=1,3
 107        zzg4(ik) = ik-1.


C                            ! find 80 km index for CO2 band, use ZP2(R$) (get last index at or below 80km)
         i80 = R$
         do 117 ik=1,R$
 117        if (zp2(ik) .le. 80.) i80 = ik



C  get pressures for RAD CODE - prmode(msr,R$+1=57)  (TOP-DOWN), should be in mbar for both SORAD and IRRAD

        do 111 ik=1,R$+1
 111       prmode(1,ik) = 1000.*EXP(-zpe(ik)/7.)



C  also get LAYER pressures  - prlay(R$)  (TOP-DOWN), so use ZPL(R$)

        do 114 ik=1,R$
 114       prlay(ik) = 1000.*EXP(-zpl(ik)/7.)



C  get pressure level index separating high and middle clouds (ict) -  use 1st index below ~400 mb
C            and level index separating middle and low clouds (icb) - use 1st index below 700 mb 
C                      use PRMODE(R$+1) array (top-down)  -  pressure LEVELS
        ict = 51
        icb = 55

        do 151 ik=R$+1,1,-1
          if (prmode(1,ik) .ge. 400.) ict = ik
          if (prmode(1,ik) .ge. 700.) icb = ik
 151    CONTINUE


C  load in REAL*8 arrays for BINTER8

        do 1701 ij=1,NS$
 1701      ychem8(ij) = ychem(ij)

        do 1702 ik=1,MS$
 1702      zchem8(ik) = zchem(ik)

        do 1703 ij=1,NP$
 1703      ypp8(ij) = ypp(ij)

        do 1705 ik=1,R$
 1705      zp8(ik) = zp2(ik)


C  
C  interpolate TTOTHX(MP$), VNC(M$), THGX(MP$) to radiation layer grid (use logarithmic)
C          ZP(M$); ZP2(R$) all are BOTTOM-UP  -->  TTOTHR(R$), VNCR(R$), THGR(R$)
C
         CALL LINTERP(zpp, MP$, TTOTHX, zp2, R$, 1, TTOTHR)
         CALL LINTERP(zpp, MP$, THGX, zp2, R$, 1, THGR)
         CALL LINTERP(zp, M$, VNC, zp2, R$, 1, VNCR)



C  interpolate surface UV-VIS albedo for current day to coupled model latitude grid: 
C    ALBDAY(180) is albedo (fraction), ALBLAT(180) are the latitudes -> albs(NP$), YPP(NP$)


         CALL LINTERP(ALBLAT, 180, ALBDAY, ypp, NP$, 0, ALBS)



C  initialize aerosol parameters, these are not used, so also set = 0.
C               number of aerosol types (naero) 
c  names of aerosols (string) - saerosols(naero)
c  aerosol mixing ratios   - raero(MSR,R$,naero)
c  relative humidity for aerosol calc. (fraction) - rhaero(MSR,R$)
C

         do 230 ii=1,naero
         do 230 ik=1,R$
         do 230 ins=1,MSR
 230        raero(ins,ik,ii) = 0.0

         do 231 ik=1,R$
         do 231 ins=1,MSR
 231        rhaero(ins,ik) = 0.0

c
CEGF      ncir=17
c
C
C interpolate model constits to heating rate latitude/altitude grid - use BINTERP (do LOG INTERP in ALT)
C    heating rate LAYER altitude grid, ZP2(R$) goes bottom-up
C
C  NOTE: if heating rate grid goes BEYOND constit grid, values are set to last grid point available
C        CNC(S$C,NS$,MS$), ychem8(NS$), zchem8(MS$), cnc18(NS$,MS$) -> cnc28(NP$,R$), YPP8(NP$), ZP8(R$)
C        FOR BINTERP and LINTERP, must have ALTITUDE arrays INCREASING in value!!!!!!!
C        CNCH(S$C,NP$,R$) is then reversed to go top-down;   9HB - do in REAL*8
C
        do 505 is=1,S$C

           do 506 ik=1,MS$
           do 506 ij=1,NS$
 506          cnc18(ij,ik) = cnc(is,ij,ik)

           CALL BINTERP8(ychem8, NS$, zchem8, MS$, cnc18,
     >                   YPP8, NP$, ZP8, R$, 0, 1, CNC28)

           DO 507 ik=1,R$
           DO 507 ij=1,NP$
 507	     cnch(is,ij,R$+1-ik) = cnc28(ij,ik)

 505	 CONTINUE



C  also interpolate TOTAL diurnally avged heating rate from photolysis calc. (K/sec)
c  onto coupled model EXTENDED latitude/altitude grid, goes bottom-up:
C    PHEAT(5,NS$,MS$) -> PHEATC(N$,M$), ZP2(R$)

           do 617 ik=1,MS$
           do 617 ij=1,NS$
 617          cnc1(ij,ik) = PHEAT(5,ij,ik)

           CALL BINTERP(ychem, NS$, zchem, MS$, cnc1, 
     >                  YP, N$, ZP, M$, 0, 0, PHEATC)


C  
C   smooth Temperature field before doing radiation, TEMPBX(NP$,MP$), TEMPBXS(NP$,MP$)
C      SMOOTH54 (below) is for REAL*4 arrays

cccpw82           CALL SMOOTH54(tempbx, tempbxs, NP$, MP$)



C  
C   get WACCM cloud parameters, interpolated to current radiation grid, TOP-DOWN
C                                                      ZP2(R$), ZPL(R$)    
C   CLMR(NP$,R$,3)  -  mixing ratio (kg/kg, or gm/gm)   TOP-DOWN
C   CLFR(NP$,R$)    -  fraction
C   CREFF(NP$,R$,3) -  effective size of cloud particles (micron)
c      index 1 for ice particle
c      index 2 for liquid drops
c      index 3 for rain drops (NOT USED CURRENTLY)


CCWCONV191           CALL WACL(YPP, NP$, ZPL, ZP2, R$, IDOY360, CLMR, CLFR, CREFF)


CCWCONV191 - use GEOS 5 TD cloud parameters;  get CO2 BC for current day, CNC(S$C,NS$,MS$) is mixing ratio (ppp)

           co2day = CNC(20,9,1)*1.e6

           CALL G5CL(YPP, NP$, ZPL, ZP2, R$, IDOY360, CO2DAY, G5CLDAY)



c                                 ...beginning of lat loop
      DO 2000 IJ = 1, NP$


c  define surface reflectivity, function of latitude/season - ALBS(NP$) is for UV-VIS
C    use for both BEAM and DIFFUSE
c
c   as a VERY ROUGH APPROXIMATION, define IR sfc albedo as 75% of UV/VIS albedo
C      eg, as deduced from Roesch et al., 1992;   but compared with using 0 albedo, 
C      this has only a VERY SMALL impact on the resulting SW (a few %) and 
C      temperature (~.01K-.05K) in the lower troposphere
C
c        in the UV+par region:
c           for beam insolation    (rsuvbm)        fraction     m 
c           for diffuse insolation (rsuvdf)        fraction     m 
c        in the near-ir region:
c           for beam insolation    (rsirbm)        fraction     m  
c           for diffuse insolation (rsirdf)        fraction     m
C
         rsuvbm(1) = ALBS(ij)
         rsuvdf(1) = ALBS(ij)
         rsirbm(1) = ALBS(ij)*.75     ! 0.0            ! - changing from 0 has a VERY SMALL effect in SORAD
         rsirdf(1) = ALBS(ij)*.75     ! 0.0


C  define sub-grid parameters:   number of sub-grid surface types - nsurf  
C  land or ocean srface emissivity - eg(msr,nsurf,10) - srface reflectivity in 10 ir bands, set = 1, per JOAN
C  probably best NOT to be consistent with rsirbm/rsirdf, since the SOLIR wavelengths are 
C  mostly much SHORTER than in IRRAD 
C   (at the surface: A(l) + R(l) = 1 ; A=absorptivity, R=reflectivity, l=wavelength,(6.11, p.295,Wallace&H)
C
C  for vegetation emissivity and reflectivity, ev, rv(msr,nsurf,10) just set to 0.0

          do 280 iir = 1,10
          do 280 ik=1,nsurf
          do 280 ins=1,MSR
             eg(ins,ik,iir) = 1.
             ev(ins,ik,iir) = 0.
 280         rv(ins,ik,iir) = 0.



C
C  ****************   FOR  ISCCP  CLOUD  CLIMATOLOGY  *************************************
C
c ...get profiles of clouds, diurnal average quantities, also pass latitude and pressure grids
C      YPP(NP$), PRLAY(R$) (top-down),  cwc1(R$,3) and cloudi(R$) go top down

cccccccc       CALL RADRED9(YPP, ij, PRLAY, msr, lprint, IDOY360, CWC1, CLOUDI)

c
c  load in CLOUD Parameters cwc1(R$,3) and cloudi(R$) defined from CLOUD CLIM in RADRED9
C
C    effective size of cloud particles from ISCCP - reff(msr,R$,3): 
C             1=ice particles (30 um);  2=liquid drops (10 um)
C
C  NOTE: liquid water clouds (cwc(2)) occur ONLY at and below the 4.06 km level (~550 mbar).
C    Since the 2.03 km level (bottom) is specified, the interaction of the radiative effects
C    of the liquid clouds amongst the different levels is NOT properly handled. So it's best
C    to TURN OFF the liquid water clouds. Also, effect of the ice clouds is too large, so reduce by 50%.
C    
C      ALL CLOUDS TURNED OFF  !!!!! - CLDFLG = .false. also turns off clouds
C
C  ****************   END  ISCCP  CLOUD  CLIMATOLOGY  *************************************
C
C
C
C   LOAD IN WACCM CLOUD PARAMETERS:      clmr(NP$,R$,3), creff(NP$,R$,3), clfr(NP$,R$),
C
C  cloud water mixing ratio (gm/gm or kg/kg) - cwc(msr,R$,3),  reff(msr,R$,3) :
C            1=ice particles;  2=liquid drops;  3=rain drops (not used)
C     
C   cloud fraction - fcld(msr,R$),  and   cloud indicies icb, ict are defined above
C
CCWCONV191
c         do 120 ik=1,R$
c         do 120 ins=1,MSR
c            cwc(ins,ik,1) = clmr(ij,ik,1)
c            cwc(ins,ik,2) = clmr(ij,ik,2)
c            cwc(ins,ik,3) = 0.0
c
c           reff(ins,ik,1) = creff(ij,ik,1)        ! 30.
c           reff(ins,ik,2) = creff(ij,ik,2)        ! 10.
c           reff(ins,ik,3) = 0.0
c 120     CONTINUE
c
c         do 121 ik=1,R$
c         do 121 ins=1,MSR
c 121       fcld(ins,ik) = clfr(ij,ik)



CCWCONV191 - GEOS 5 TD cloud parameters, interpolated to current day on radiation grid 
C                    G5CLDAY(5,NP$,R$):
C
C       1 = mass_fraction_of_cloud_liquid_water (kg/kg)
C       2 = mass_fraction_of_cloud_ice_water    (kg/kg)
C       3 = cloud_fraction
C       4 = liquid_cloud_particle_effective_radius (um - microns)
C       5 = ice_cloud_particle_effective_radius (um - microns)
C
C
C  cloud water mixing ratio (gm/gm or kg/kg) - cwc(msr,R$,3),  reff(msr,R$,3) :
C            1=ice particles;  2=liquid drops;  3=rain drops (not used)
C     
C   cloud fraction - fcld(msr,R$),  and   cloud indicies icb, ict are defined above
C
         do 410 ik=1,R$
         do 410 ins=1,MSR
            cwc(ins,ik,1) = G5CLDAY(2,ij,ik)
            cwc(ins,ik,2) = G5CLDAY(1,ij,ik)
            cwc(ins,ik,3) = 0.0

           reff(ins,ik,1) = G5CLDAY(5,ij,ik)
           reff(ins,ik,2) = G5CLDAY(4,ij,ik)
           reff(ins,ik,3) = 0.0

           fcld(ins,ik) = G5CLDAY(3,ij,ik)
 410     CONTINUE



ccwrite - cloudj(R$), taucld(R$), omgcld(3,R$), gcld(3,R$)
ccwrite             write (711, 1000)
ccwrite             write (711, 1001) idoy360, ij, ypp(ij)
ccwrite          
ccwrite             do 555 ik=38,R$
ccwrite 555            write (711, 1002) zpl(ik), prlay(ik), cloudj(ik), 
ccwrite     >            taucld(ik), gcld(1,ik), gcld(2,ik), gcld(3,ik),
ccwrite     >            omgcld(1,ik), omgcld(2,ik), omgcld(3,ik)
ccwrite 1000        format(5x)
ccwrite 1001        format(5x, 2I10, F10.3)
ccwrite 1002        format(F6.2, F9.3, 2F12.4, 6F12.3)
ccwrite


C
C  load in temperatures here, use TEMPBX(NP$,MP$) - coupled model extended grid, interpolate to LAYER TEMPS
C    then reverse to go TOP-DOWN for heating rates   -->    TMODR(MSR,R$) goes TOP-DOWN
C    ttx(MP$), ZP2(R$) are layer altitudes, BOTTOM-UP -> ttr(R$)
C
C    then use GEOS 4 SURFACE temps: TEMPG42D(NP$,3) are sfc, 1km, 2km;
C
C                                               use UNSMOOOOOOTHED TEMPS here - TEMPBX(NP$,MP$)

         do 211 ik=1,MP$
 211        ttx(ik) = TEMPBX(ij,ik)

            CALL LINTERP(zpp, mp$, ttx, ZP2, R$, 0, TTR)

         do 214 ik=1,R$
 214        TMODR(1,R$+1-ik) = TTR(ik)


C
C  interpolate TEMPG42D(NP$,3) (sfc, 1km, 2km - zzg4(3) ) to lowest 2 levels of radiation grid for current latitude
C  load in level 1 ONLY, ie CENTER of bottom layer, ZPL(R$) (usually 0.5 km), ttg4(3), zzg4(3), zzgr(2), ttgr(2)
C                                              eg, ttgr(2) is for 0.5075 km and 1.52250 km
         do 702 ik=1,3
 702        ttg4(ik) = TEMPG42D(ij,ik)

         CALL LINTERP(zzg4, 3, ttg4, zzgr, 2, 0, ttgr)


ccpw36           TMODR(1,R$) = ttgr(1)

C  lowest radiation grid point is .5075 km, so take avg of NCEP at sfc (~0km) and lowest model level (1.015 km)
C
            TMODR(1,R$) = (TEMPG42D(ij,1) + TEMPBX(ij,1))/2.


cccpw32         TMODR(1,R$) = (TEMPG42D(ij,2) - TEMPG42D(ij,1)) * zpl(R$) 
cccpw32     >                                 + TEMPG42D(ij,1)



c  surface temperature (K) - use TD GEOS 4 sfc values here for current day (not bottom level of model)
C      do the same for the land or ocean surface temperature (K)  tg(msr,nsurf)
C      and the vegetation temperature ,  tv(msr,nsurf) ,           TNCEPSFC(NP$,2), TEMPG42D(NP$,3) -> tsurf(MSR)
C      also set fractional cover of sub-grid regions,  fs(msr,nsurf)=1.

         do 285 ik=1,nsurf
         do 285 ins=1,MSR
            tsurf(ins) = TEMPG42D(ij,1)     !TNCEPSFC(ij,1)     ! TMODR(ins,M$+1)
            tg(ins,ik) = TEMPG42D(ij,1)     !TNCEPSFC(ij,1)     ! TMODR(ins,M$+1)
            tv(ins,ik) = TEMPG42D(ij,1)     !TNCEPSFC(ij,1)     ! TMODR(ins,M$+1)
 285        fs(ins,ik) = 1.


c
C   get constituent profiles - CNCH(S$C,NP$,R$) is mixing ratio (PARTS/PART) for current latitude/day
C     YPP(NP$) is latitude ,  ZPL(R$) goes top-down ,  CNCH goes top-down
C     also convert water vapor and ozone to gm/gm, use Joan's conversions
C     put CO2, CH4, N2O, F11, F12, F22 IN PARTS/PART
C
C     also for RADIATION CODE ONLY, limit Ozone in POLAR mesosphere to 7 ppmv to avoid blowups 
C     this is ONLY necessary for the first month or so where the initial conditions have 
C     VERY LARGE OZONE (hundreds of ppmv) in polar night, otherwise, mesospheric ozone stays < ~3.5 ppmv
C
      do 500 ik=1,R$
         vmro3(ik) = CNCH(4,ij,ik)
         if (zpl(ik) .ge. 60.  .and.  ABS(ypp(ij)) .ge. 60.) then
              if (vmro3(ik) .ge. 7.e-6) vmro3(ik) = 7.e-6
         endif

         oa(1,ik) = vmro3(ik)*1.657
         wa(1,ik) = CNCH(15,ij,ik)*0.622

         co2vmr(1,ik) = CNCH(20,ij,ik)
         vmrch4(1,ik) = CNCH(18,ij,ik)
         vmrn2o(1,ik) = CNCH(11,ij,ik)
         vmrcfc11(1,ik) = CNCH(34,ij,ik)
         vmrcfc12(1,ik) = CNCH(35,ij,ik)
         vmrcfc22(1,ik) = CNCH(52,ij,ik)
 500   CONTINUE


c
c...compute solar heating, need to march through diurnal cycle for call to SORAD - do every 1 hour
C   NOTE: I tried this with 1/2 hour integration, differences were small around the equinoxes (<~.5K), 
C       larger at the solstices - w/ 1 hour integration, temps were up to ~3K warmer in the summer mesosphere
C         1-2K colder in the winter mesosphere, zero line right at the equator, 
C         only a few tenths K difference in the stratosphere.
C
C         solar heating just needs COSZ = COS(Solar zenith angle) - (adapted from SETUPA)
C
C  also, for SOLAR ZENITH ANGLE - problem with directly overhead sun and COSZ=1, it may actually be 
C     numerically slightly greater than 1 (eg, 1.0000001), so that ACOS(COSZ) = NaN, so reset COSZ = 1.
C                                         SOLDEC is the SOLAR DECLINATION in DEGREES
       sind = SIN(soldec*dtr)
       cosd = COS(soldec*dtr)

       sinl = SIN(ypp(ij)*DTR)
       cosl = COS(ypp(ij)*DTR) + 1.e-20

                                          ! initialize diurnal avg array, flxbda(msr,R$+1,15)
       do 504 ik=1,R$+1
 504      flxda(1,ik) = 0.0

       do 514 ib=1,15
       do 514 ik=1,R$+1
 514      flxbda(1,ik,ib) = 0.0

C                                             the 15 factor converts hours-degrees, DTR converts degrees-rads
       DO 700 ICLOCK=1,24
C                                                                             COSZ = COS(SZA)
         COSZ(1) = sind*sinl + cosd*cosl*COS(15.*(iclock + 12.)*DTR)
         if (COSZ(1) .gt. 1.) COSZ(1) = 1.
         if (COSZ(1) .lt. -1.) COSZ(1) = -1.
                                                            ! use flx(msr,R$+1) - all sky flux, flxb(msr,R$+1,15)
         if (COSZ(1) .gt. 0.) then 

            CALL SORAD (msr, R$, cosz, prmode, tmodr, wa, oa, co2vmr,
     *         overcast,cwc,fcld,ict,icb,reff,
     *         naero, saerosols, raero, rhaero,
     *         rsuvbm,rsuvdf,rsirbm,rsirdf, hjfluxday,
     *        flx,flc,fdiruv,fdifuv,fdirpar,fdifpar,fdirir,fdifir,flxb)

          else
C                                               set solar heating = 0.0 if sun is below the horizon
              do 502 ik=1,R$+1
 502             flx(1,ik) = 0.0

              do 503 ib=1,15
              do 503 ik=1,R$+1
 503             flxb(1,ik,ib) = 0.0

          endif

C  sum up diurnal avg, and apply COS(SZA) and TOA flux (1367 W/m2), flxbda(msr,R$+1,15), flx(1,k)=flxb(1,k,14)
C
C   substitute  1367 W/m2 here with TOA value from J. Lean's data for proper month
C                                                        use  sjtoa defined above
 
             do 517 ik=1,R$+1
 517             flxda(1,ik) = flxda(1,ik) + flx(1,ik)*sjtoa*COSZ(1)/24.

              do 577 ib=1,15
              do 577 ik=1,R$+1
 577   flxbda(1,ik,ib) = flxbda(1,ik,ib)+flxb(1,ik,ib)*sjtoa*COSZ(1)/24.

 700    CONTINUE
C                     ! flxda(msr,R$+1), flxir(msr,R$+1) - all sky flux (same as clear sky flux w/ NO clouds)


c...compute longwave cooling, 1X per day, now include cloud parameters

           CALL IRRAD (msr, R$, prmode, tmodr, wa, oa, tsurf, co2vmr,
     *        high, trace, vmrn2o, vmrch4, vmrcfc11, vmrcfc12, vmrcfc22,
     *                  overcast,cwc,fcld,ict,icb,reff,
     *                  nsurf, fs, tg, eg, tv, ev, rv,
     *                  naero, saerosols, raero, rhaero,
!     *                  aerosol,na,taual,ssaal,asyal,
     *                  flxir, flcir, flair, dfdts, sfcem, taudiag)



C  for output, load in flxir(msr,R$+1), flcir(msr,R$+1), flair(msr,R$+1) -> irflux(NP$,R$+1,3) in W/m^2

          do 757 ik=1,R$+1
            irflux(ij,ik,1) = flxir(1,ik)
            irflux(ij,ik,2) = flcir(1,ik)
 757        irflux(ij,ik,3) = flair(1,ik)


C      
C  fluxes are across grid box edges, compute heating and cooling rates at centers of boxes (grid points)
C      flxda(msr,R$+1), flxir(msr,R$+1) are in W/m^2, prmode(1,R$+1) -> HEATT(R$), COOLT(R$)
C      to convert to K/sec:  multiply by g/Cp, use pressures in Pascals
C
          do 711 ik=1,R$
             heatt(ik) = -9.81/1004.*(flxda(1,ik+1) -  flxda(1,ik))/
     >                        (100.*(prmode(1,ik+1) - prmode(1,ik)))

             coolt(ik) = 9.81/1004.*(flxir(1,ik+1)  -  flxir(1,ik))/
     >                        (100.*(prmode(1,ik+1) - prmode(1,ik)))
 711      CONTINUE

C
C  also do for heating bands separately:  flxbda(msr,R$+1,15),  heatuv(NP$,R$,15) is in K/sec, reverse to go bottom up 
C                                          and for bands 5-8 set to zero above 60 km (remove the junk), use ZP2(R$)
          do 721 ib=1,15
          do 721 ik=1,R$
             hhh = -9.81/1004.*(flxbda(1,ik+1,ib) - flxbda(1,ik,ib))/
     >                        (100.*(prmode(1,ik+1) - prmode(1,ik)))
             heatuv(ij,R$+1-ik,ib) = hhh
 721      CONTINUE

          do 722 ib=5,8
          do 722 ik=1,R$
             if (ZP2(ik) .ge. 60.) heatuv(ij,ik,ib) = 0.
 722      CONTINUE


C  CO2 contribution (ib=10) gets very large above 80 km (maybe not realistic)
C    so reset to constant gradient, i80 is 80 km index from above, but DON'T go negative
C    also need to recompute total DF - heatuv(11) = 9+10:
 
          delco2 = heatuv(ij,i80,10) - heatuv(ij,i80-1,10)
          do 723 ik=i80+1,R$
             heat10 = heatuv(ij,ik,10)                                ! save unadjusted value

             heatuv(ij,ik,10) = heatuv(ij,ik-1,10) + delco2
             if (heatuv(ij,ik,10) .le. 0.) heatuv(ij,ik,10) = 0.

         heatuv(ij,ik,11) = heatuv(ij,ik,11) - heat10 + heatuv(ij,ik,10)
 723      CONTINUE


C  also sum up UV bands 1-4 and IR bands 6-8 for diagnostics. Also sum up bands 5-10 (that which is NOT computed from
C  the photolysis rates) and these will be added into the total heating rates - ALL in K/sec;  heatuv(NP$,R$,15)
C
C   w/ the adjustments above, also re-total heatuv(14), which is the total from the flx array in SORAD
C   this is close to what I get when totaling the individual bands together, with some small 
C                                                     differences in the upper mesosphere (no big whoop)

          do 725 ik=1,R$
             heatuv(ij,ik,12) = heatuv(ij,ik,1) + heatuv(ij,ik,2)
     >                        + heatuv(ij,ik,3) + heatuv(ij,ik,4)

             heatuv(ij,ik,13) = heatuv(ij,ik,6) + heatuv(ij,ik,7)
     >                        + heatuv(ij,ik,8)

             heatuv(ij,ik,14) = heatuv(ij,ik,11) + heatuv(ij,ik,12)
     >                        + heatuv(ij,ik,13) + heatuv(ij,ik,5)

             heatuv(ij,ik,15) = heatuv(ij,ik,5) + heatuv(ij,ik,6)
     >                        + heatuv(ij,ik,7) + heatuv(ij,ik,8) 
     >                        + heatuv(ij,ik,9) + heatuv(ij,ik,10)
 725      CONTINUE


c  
C  NOW reverse heating/cooling to go bottom-up -> HEATG5E(NP$,R$), COOLG5E(NP$,R$) are in K/sec
C    w/ photolysis calc, just include heatuv(15) here (IR bands ONLY), which is ALREADY bottom-up
C                                          ;  heatuv(NP$,R$,15)
          do 717 ik=1,R$
cccheatph           HEATG5E(ij,R$+1-ik) = heatt(ik) - !old

             HEATG5E(ij,ik) = heatuv(ij,ik,15)
             COOLG5E(ij,R$+1-ik) = coolt(ik)
 717      CONTINUE



C  ABOVE 85 km:  damp out Cooling rate ONLY - USE NEWTONIAN COOLING as done in NEWRAD7; VNC(M$) is 1/sec
C    use ZP2(R$).  (don't need to damp out heating w/ phot calc.)
C         use TTOTHR(R$), VNCR(R$), THGR(R$) are BOTTOM-UP; TMODR(1,R$) is TOP-DOWN

       ZP0X = 85.

       do 701 ik=2,R$
          if (zp2(ik) .GT. ZP0X) then

cccheatph             damp = EXP(-(zp2(ik)-ZP0X)/4.)
cccheatph             HEATG5E(ij,ik) = HEATG5E(ij,ik)*damp

             WGT = EXP(-(zp2(ik)-ZP0X)/4.)
             COOLG5E(ij,ik) = COOLG5E(ij,ik)*WGT +
     >  (1.-WGT)*VNCR(ik-1)*(TMODR(1,R$+1-ik+1)-THGR(ik-1)/TTOTHR(ik-1))

          endif
 701   CONTINUE



C  Calculate cooling in bottom 2 LAYERS (.5-1.5 km) using relaxation to TD GEOS 4 - TEMPG42D(NP$,3) are sfc, 1km, 2km 
C
C  temperature field (or enhanced temperature field).
C    TMODR(1,R$) is the model calculated temperature field. TMODR(MSR,R$) goes TOP-DOWN
C  Relaxation timescale is 1 day (get closer to OBS compared to timescale of 2 days per Joan)
C
C    do this for the bottom 2 layers will cover the bottom level of the EXTENDED DYNAMICS GRID (MP$)
C       when COOLG5E(NP$,R$) is interpolated to coupled model grid, BOTTOM-UP
C
C   (THIS WAS TAKEN FROM TROPHT in PRMLIB.f), RELAX SRFC. TEMP. to GEOS 4 Time dependent temps (theta)
C
C   ttgr(2) is TEMPG42D(NP$,3) interpolated to lowest 2 levels of radiation grid for current latitude
C                         NOTE: TMODR(1,R$)=ttgr(1) defined above, so the RELAXATION = 0 at ik=1 here
C
C    keep this in even when forcing bottom boundary TEMPS to NCEP/GEOS 4 values - for consistency in streamfunction
C
         SRFNMC = 1.0000

         do 703 ik=1,2
            COOLG5E(ij,ik) = COOLG5E(ij,ik) +
     >          SRFNMC * (TMODR(1,R$+1-ik) - TTGR(ik))/(1.*DAYL)
 703     CONTINUE


       IF (SRFNMC.LE.0.0) THEN 
         WRITE(6,'(5X,"**** NO SURFACE COOL TO OBS. !!!!")' )
       ENDIF

C                                         !  load in NETG5E(NP$,R$) in K/sec for diagostics ONLY - NOT USED
cccheatph       do 705 ik=1,R$
cccheatph 705       netg5e(ij,ik) = heatg5e(ij,ik) - coolg5e(ij,ik)


 2000     CONTINUE
c                       ...end of lat loop,  IJ = 1, NP$



C   interpolate to coupled model INNER grid (N$,M$)
C       HEATG5E(NP$,R$), COOLG5E(NP$,R$) are in K/sec, and are BOTTOM-UP - use ZP2(R$) goes BOTTOM-UP
C   --> HEATXX(N$,M$), COOLXX(N$,M$)         ypp(NP$) ->  YP(N$), zp(M$)
C

           CALL BINTERP(ypp, NP$, zp2, R$, HEATG5E, 
     >                  yp, N$, zp, M$, 0, 0, HEATXX)

           CALL BINTERP(ypp, NP$, zp2, R$, COOLG5E, 
     >                  yp, N$, zp, M$, 0, 0, COOLXX)


C
C  add in latent heating, G4LHDAYC(N$,M$) in K/day, convert to K/sec
C  this is TD Latent heating from WACCM (seasonal cycle) + GEOS 4 trend/CO2 sensitivity
C
C  add TOTAL diurnally avged heating rate from photolysis calc. for UV + VIS
c      PHEATC(N$,M$) goes bottom-up  (K/sec)
C 
C   load in arrays HEATG5(N$,M$), COOLG5(N$,M$) for return to coupled model dynamics (TWODS)
C
       DO K=1,M$
       DO J=1,N$ 
          HEATG5(J,K) = HEATXX(J,K) + PHEATC(J,K) + G4LHDAYC(J,K)/86400.
          COOLG5(J,K) = COOLXX(J,K)
       END DO
       END DO



C
C interpolate heating and cooling rates to chemistry latitude/altitude grid for output 
C   - use BINTERP - must have ALTITUDE arrays INCREASING in value!!!!!!!
C
C  use HEATG5E(NP$,R$), COOLG5E(NP$,R$) are in K/sec, and are BOTTOM-UP - use ZP2(R$) goes BOTTOM-UP, YPP(NP$)
C
C  NOTE: if CONSTIT grid goes BEYOND heating rate grid, values are set to last grid point available
C
C   ychem(NS$), zchem(MS$), HEATCHM(NS$,MS$), COOLCHM(NS$,MS$)
C
C                                    - NOTE: HEATG5E and HEATCHM here are ONLY for the IR bands

           CALL BINTERP(ypp, NP$, zp2, R$, HEATG5E, 
     >                  YCHEM, NS$, ZCHEM, MS$, 0, 0, HEATCHM)


           CALL BINTERP(ypp, NP$, zp2, R$, COOLG5E, 
     >                  YCHEM, NS$, ZCHEM, MS$, 0, 0, COOLCHM)



C  interpolate the 15 bands from SORAD.f to chemistry grid for output in MAIN;  cnc1(NS$,MS$), cnc2(NP$,R$)
C  heatuv(NP$,R$,15) -> SORADHEAT(15,NS$,MS$), in K/sec; ZP2(R$) goes bottom-up ; zchem(MS$)
C

        do 488 ib=1,15

           do 489 ik=1,R$
           do 489 ij=1,NP$
 489	       cnc2(ij,ik) = heatuv(ij,ik,ib)


           CALL BINTERP(ypp, NP$, zp2, R$, CNC2,
     >                  ychem, NS$, zchem, MS$, 0, 0, CNC1)


           do 491 ik=1,MS$
           do 491 ij=1,NS$
 491         SORADHEAT(ib,ij,ik) = cnc1(ij,ik)

 488	 CONTINUE



C  write out IRFLUX(NP$,R$+1,3) array every mid-month to fort.223 - in W/m^2
C    also YPP(NP$), ZPE(R$+1=57), PRMODE(MSR=1,R$+1),   ZPL(R$) goes TOP-DOWN, ZP2(R$) goes BOTTOM-UP
C    for reference,  also write out COOLG5E(NP$,R$) and HEATG5E(NP$,R$) (which is ONLY for the IR bands)
C
          if (mod(idoy360-15,30) .eq. 0.0) then 
              write (223) NP$, R$
              write (223) ypp, zpe, prmode, zpl, zp2
              write (223) irflux
              write (223) heatg5e
              write (223) coolg5e
          endif


      RETURN
      END




      SUBROUTINE WACL(YPP,NP$, ZPL, ZP2, R$, IDOY360, CLMR, CLFR, CREFF)

C
C  routine to get WACCM cloud parameters for CURRENT DAY
C        interpolate to current radiation grid, TOP-DOWN
C

C
C  COMMON for WACCM cloud parameters:  WWLIQ, WWICE are in kg/kg;  WWCLD is fraction
C        WWRLIQ, WWRICE are in microns  (from TEMPIN);  timesw(14)

      COMMON/CWCLOUDS/wwliq(46,33,14), wwice(46,33,14), wwcld(46,33,14),
     >                wwrliq(46,33,14), wwrice(46,33,14),
     >                wwlat(46), wwpr(33), timesw(14)



      INTEGER NP$, R$, IDOY360
      REAL YPP(NP$), ZPL(R$), ZP2(R$), ttime(1), ttout(1), ttw14(14)
      REAL wwz(33), wwin(46,33), wwout(NP$,R$)
      REAL clmr(NP$,R$,3), clfr(NP$,R$), creff(NP$,R$,3)


        do 117 ik=1,33
 117       wwz(ik) = 7.*ALOG(1000./wwpr(ik))


        ttime(1) = IDOY360*1.


C  LIQUID MIXING RATIO:   wwliq(46,33,14) -> clmr(NP$,R$,2)

        do 150 ik=1,33
        do 150 ij=1,46

           do 151 iii=1,14
 151          ttw14(iii) = wwliq(ij,ik,iii)

           CALL LINTERP(timesw, 14, ttw14, ttime, 1, 0, ttout)

           wwin(ij,ik) = ttout(1)
 150    CONTINUE


C  interpolate wwin(46,33) to RADIATION GRID (NP$,R$), first BOTTOM-UP
C    then reverse to go TOP-DOWN :   ZP2(R$) goes BOTTOM-UP,  ZPL(R$) goes TOP-DOWN

        CALL BINTERP(wwlat, 46, wwz, 33, WWIN,
     >               ypp, NP$, zp2, R$, 0, 0, WWOUT)


C  load into cloud water mixing ratio (gm/gm, or kg/kg) ->  CLMR(NP$,R$,3) goes TOP-DOWN
C                                1=ice particles;  2=liquid drops;  3=rain drops (not used)
        do 175 ik=1,R$
        do 175 ij=1,NP$
 175       clmr(ij,R$+1-ik,2) = wwout(ij,ik)




C  ICE MIXING RATIO:  wwice(46,33,14) -> clmr(NP$,R$,1)

        do 160 ik=1,33
        do 160 ij=1,46

           do 161 iii=1,14
 161          ttw14(iii) = wwice(ij,ik,iii)

           CALL LINTERP(timesw, 14, ttw14, ttime, 1, 0, ttout)

           wwin(ij,ik) = ttout(1)
 160    CONTINUE


        CALL BINTERP(wwlat, 46, wwz, 33, WWIN,
     >               ypp, NP$, zp2, R$, 0, 0, WWOUT)


        do 176 ik=1,R$
        do 176 ij=1,NP$
           clmr(ij,R$+1-ik,1) = wwout(ij,ik)
 176       clmr(ij,R$+1-ik,3) = 0.0




C  CLOUD FRACTION:   wwcld(46,33,14) -> clfr(NP$,R$)

        do 170 ik=1,33
        do 170 ij=1,46

           do 171 iii=1,14
 171          ttw14(iii) = wwcld(ij,ik,iii)

           CALL LINTERP(timesw, 14, ttw14, ttime, 1, 0, ttout)

           wwin(ij,ik) = ttout(1)
 170    CONTINUE


        CALL BINTERP(wwlat, 46, wwz, 33, WWIN,
     >               ypp, NP$, zp2, R$, 0, 0, WWOUT)


        do 177 ik=1,R$
        do 177 ij=1,NP$
 177       clfr(ij,R$+1-ik) = wwout(ij,ik)




C  LIQUID EFFECTIVE RADIUS:  wwrliq(46,33,14) -> creff(NP$,R$,2), limit to 8-14 um

        do 180 ik=1,33
        do 180 ij=1,46

           do 181 iii=1,14
 181          ttw14(iii) = wwrliq(ij,ik,iii)

           CALL LINTERP(timesw, 14, ttw14, ttime, 1, 0, ttout)

           if (ttout(1) .le. 8.) ttout(1) = 8.
           if (ttout(1) .ge. 14.) ttout(1) = 14.

           wwin(ij,ik) = ttout(1)
 180    CONTINUE


        CALL BINTERP(wwlat, 46, wwz, 33, WWIN,
     >               ypp, NP$, zp2, R$, 0, 0, WWOUT)


C  this should be set to 14um in upper troposphere, so at each latitude
C  find highest level = 14, and set everything above to 14, then limit to 8-14 um below
C                                  wwout(NP$,R$) is BOTTOM-UP
        do 275 ij=1,NP$
          ikww=R$
          do 276 ik=1,R$
 276         if (wwout(ij,ik) .eq. 14.) ikww=ik

          do 277 ik=ikww,R$
 277         wwout(ij,ik) = 14.
 275    CONTINUE


        do 188 ik=1,R$
        do 188 ij=1,NP$

           if (wwout(ij,ik) .le. 8.) wwout(ij,ik) = 8.
           if (wwout(ij,ik) .ge. 14.) wwout(ij,ik) = 14.

           creff(ij,R$+1-ik,2) = wwout(ij,ik)
 188    CONTINUE




C  ICE EFFECTIVE RADIUS:  wwrice(46,33,14) -> creff(NP$,R$,1), limit to 5-300 um

        do 190 ik=1,33
        do 190 ij=1,46

           do 191 iii=1,14
 191          ttw14(iii) = wwrice(ij,ik,iii)

           CALL LINTERP(timesw, 14, ttw14, ttime, 1, 0, ttout)

           if (ttout(1) .le. 5.) ttout(1) = 5.
           if (ttout(1) .ge. 300.) ttout(1) = 300.

           wwin(ij,ik) = ttout(1)
 190    CONTINUE


        CALL BINTERP(wwlat, 46, wwz, 33, WWIN,
     >               ypp, NP$, zp2, R$, 0, 0, WWOUT)


        do 198 ik=1,R$
        do 198 ij=1,NP$

           if (wwout(ij,ik) .le. 5.) wwout(ij,ik) = 5.
           if (wwout(ij,ik) .ge. 300.) wwout(ij,ik) = 300.

           creff(ij,R$+1-ik,1) = wwout(ij,ik)
           creff(ij,R$+1-ik,3) = 0.0
 198    CONTINUE


ccpw40          write (696) NP$, R$, idoy360
ccpw40          write (696) ypp, zpl, zp2, wwlat, wwpr, wwz
ccpw40          write (696) wwliq, wwice, wwcld, wwrliq, wwrice
ccpw40          write (696) clmr, clfr, creff


      RETURN
      END




      SUBROUTINE G5CL(YPP, NP$, ZPL, ZP2, R$, IDOY360, CO2DAY, G5CLDAY)

C
C  routine to get GEOS 5 cloud parameters for CURRENT DAY
C        interpolate to current radiation grid, TOP-DOWN
C

C  COMMON for GEOS 5 TD cloud parameters from TEMPIN
C       1 = mass_fraction_of_cloud_liquid_water (kg/kg)
C       2 = mass_fraction_of_cloud_ice_water    (kg/kg)
C       3 = cloud_fraction
C       4 = liquid_cloud_particle_effective_radius (meters)
C       5 = ice_cloud_particle_effective_radius (meters)
C
C  2nd index:  1-14 are 14 month seasonal cycle; 
C              15 = CO2 sensitivity (in cloud parameter/ppmv)

      COMMON/CG5CL/g5clin(5,15,37,30), yg5cl(37), zg5cl(30), timecl(14)


      INTEGER NP$, R$, IDOY360
      REAL YPP(NP$), ZPL(R$), ZP2(R$), ttime(1), ttout(1), ttw14(14)
      REAL gcin(37,30), gcout(NP$,R$), clfac(5), cllim(5)

      REAL g5clday(5,NP$,R$)


      DATA clfac / 1., 1., 1., 1.e6, 1.e6/
      DATA cllim / 1.e-25, 1.e-25, 3.e-11, 5.e-6, 2.e-5/


      ttime(1) = IDOY360*1.


C  interpolate each parameter to current day,   add CO2 contribution (468.35 ppmv is the 1950-2100 avg)
C
C    then interpolate to RADIATION GRID (NP$,R$), first BOTTOM-UP
C      and reverse to go TOP-DOWN :   ZP2(R$) goes BOTTOM-UP,  ZPL(R$) goes TOP-DOWN


        do 222 iis=1,5

        do 1150 ik=1,30
        do 1150 ij=1,37

           do 1151 iii=1,14
 1151          ttw14(iii) = g5clin(iis,iii,ij,ik)

           CALL LINTERP(timecl, 14, ttw14, ttime, 1, 0, ttout)

           gcin(ij,ik) = ttout(1) + (co2day-468.35)*g5clin(iis,15,ij,ik)
 1150    CONTINUE


C  interpolate gcin(37,30) to RADIATION GRID (NP$,R$), first BOTTOM-UP

        CALL BINTERP(yg5cl, 37, zg5cl, 30, GCIN,
     >               ypp, NP$, zp2, R$, 0, 0, GCOUT)


C  load into G5CLDAY(5,NP$,R$) goes TOP-DOWN,  convert radii to um

          do 475 ik=1,R$
          do 475 ij=1,NP$
 475         g5clday(iis,ij,R$+1-ik) = gcout(ij,ik)*clfac(iis)

C                            set lower limits (ensures non-negatives)
          do 477 ik=1,R$
          do 477 ij=1,NP$
             if (g5clday(iis,ij,ik) .le. cllim(iis)) 
     >                  g5clday(iis,ij,ik) = cllim(iis)
 477      CONTINUE

 222    CONTINUE


ccccccc          write (696) NP$, R$, idoy360
ccccccc          write (696) ypp, zpl, zp2, co2day
ccccccc          write (696) g5clin, yg5cl, zg5cl, timecl, g5clday


      RETURN
      END




      SUBROUTINE RADRED9(YPP, ij, PRLAY, msr,lprint,idoy360,cwc1,cloudi)

c...this routine prepares profiles for radiation calculation, just get cloud quantities if needed

      INCLUDE 'PARAM.INC'


      INTEGER msr, R$

      PARAMETER (R$=56)


c start here, load in proper latitude, etc.

      real latcld(72), camt(360,3,72), cpres(3,72), ctau(360,3,72)
      real cwcr(360,72,10,2), clfrac(360,72,10), zzcld(10)
      dimension cloudj(R$), taucld(R$), omgcld(3,R$), gcld(3,R$)
      real cwc1(R$,3), cloudi(R$)

      logical cldflg

      COMMON/CLOUDS1/latcld, camt, cpres, ctau, cwcr, clfrac, zzcld
      COMMON/CLOUDS2/cldflg, taucld, omgcld, gcld, cloudj


      logical lprint
 
      real YPP(NP$), PRLAY(R$)


c
c cloud processing, hardwired for 3 cloud layers
C    also intialize CWC1 array;  1=ICE crystals;  2=liquid droplets;  3=rain particles (not used)
c                                and initialize cloud fraction array cloudi(R$)
      do n=1,R$
        cloudi(n)=0.0
        cloudj(n)=0.0
	taucld(n)=0.0
          do it=1,3
	    omgcld(it,n)=0.0
	    gcld(it,n)=0.0
	    cwc1(n,it)=0.0
	  enddo  
      enddo 


      IF (cldflg) THEN

c first find lat index
c then initialize amount, and optical parameters for the 3 near-ir bands
c ic=1 and 2 are 2 lowest level clouds, assumed to be water
c ic=3 is top level cloud, assumed ice

c omgcld and gcld uses Chou et al. parameterization (J. Climate, 11, 202-214,
c 1998) with the effective radii given by ISCCP (10 um for water, 
c 30 um for ice);   camt(360,3,72), ctau(360,3,72) now on 360 day-year (EF, NOV. 2007) - read in TEMPIN
C
C   NOTE that the cloud optical thickness (taucld) here is for the VISIBLE
C           need to convert for IR, as in IRRAD (R&S,p.2273)
C                                                    convert these to the LAYER PRESSURES (NOT the EDGES)
          ilat = nearidx9(ypp(ij),latcld,72)
	  if (ilat .gt. 72) ilat=72

C  this is used for loading in cloud optical thickness, etc, per Joan
C
          do ic=1,3
            ip = nearidx9(cpres(ic,ilat), prlay, R$)
            cloudj(ip) = camt(idoy360,ic,ilat)
	    if(cloudj(ip).gt.0.) then 
	      taucld(ip)=ctau(idoy360,ic,ilat)
	      if(ic.eq.1 .or. ic.eq.2) then
	         omgcld(1,ip)=1.0
	         omgcld(2,ip)=0.992
	         omgcld(3,ip)=0.846
		 gcld(1,ip)=0.854
		 gcld(2,ip)=0.844
		 gcld(3,ip)=0.866
              endif
	      if(ic.eq.3) then
	         omgcld(1,ip)=1.0
		 omgcld(2,ip)=0.977
		 omgcld(3,ip)=0.830
		 gcld(1,ip)=0.782
		 gcld(2,ip)=0.800
		 gcld(3,ip)=0.875
	      endif
	    endif
          enddo

C  and load in CWC and cloud fraction arrays:    cwc1(R$,3), cloudi(R$)
C  use cwcr(360,72,10,2) - index 1=ICE crystals;  2=liquid droplets; clfrac(360,72,10), zzcld(10),  
C  find proper altitude index, use PRLAY(R$) which goes TOP-DOWN, and convert to altitude
C  NOTE: this works fine even though zzlay starts at 93 km, 
C  everything above level 10 is set to ikc=10, which is ALWAYS 0.0, 
C                                           and limits of 1 and 10 are already set in NEARIDX9
          do 777 ik=1,R$
              zzlay = 7.*alog(1000./prlay(ik))
              ikc = nearidx9(zzlay, zzcld, 10)

              cwc1(ik,1) = cwcr(idoy360,ilat,ikc,1)
              cwc1(ik,2) = cwcr(idoy360,ilat,ikc,2)
              cloudi(ik) = clfrac(idoy360,ilat,ikc)
 777      CONTINUE

      END IF
C              end CLDFLG loop


C
c * * *
      if(.not.lprint) go to 425
      write(6,*) 'printout from radred9'
      write(6,*) 'cloud quantities, cloudi, taucld,omgcld,gcld'
      do n=1,R$
        write(6,*) n,cloudj(n),taucld(n),omgcld(3,n),gcld(3,n)
      enddo
  425 continue

  410 format(' ',i3,2x,f10.2,f9.2,2x,f10.1,2x,e10.4,2x,e10.4)

      RETURN
      END


      integer function nearidx9(x,a,dim)

c Finds the index of the nearest value in an array
c T.Atwater 7/25/91

c Input:
c    x    -- real value to figure index for
c    a    -- real array to find index in; must be montonically increasing
c    dim  -- dimension of a
c Output:
c    returns index of array value closest to x
c        returns 1 if x is less than a(1)
c        returns dim if x is greater than a(dim)

      real x
      integer dim
      real a(dim)

      nearidx9 = 1
      if (x.le.a(1)) return
      do 10 i=1,dim-1
         if (x.le.a(i+1)) then
            if (abs(x-a(i)).lt.abs(x-a(i+1))) then
               nearidx9=i
            else
               nearidx9=i+1
            end if
            return
         end if
10    continue
      nearidx9 = dim

      RETURN
      END


C
C   **************************   SMOOTHING  SUBROUTINES ****************************************
C
C           these are the same as SMOOTH5/SMOOTH5C in GWAVE, except in REAL*4
C

        SUBROUTINE SMOOTH54(u1, u2, m, n)
 
c
c  subroutine to smooth 2D array,  u1(m,n),  and  return  smooooothed array u2(m,n)
C     smooth in the identical manner as smooth5.pro, adapted from smooth5.pro  (EF, 7/27/95)
C
C  *****  weights are:  .05, .075, .5
C     since, eg  x(-1,-1) is sqrt(2) further from x(0,0) than is x(-1,0)
C

        integer m, n
        REAL*4 u1(m,n), u2(m,n)


C  FIRST   SMOOTH all INTERIOR points

       do 200 ik=2,n-1
       do 200 ij=2,m-1

       u2(ij,ik) = .05*(u1(ij-1,ik-1) + u1(ij-1,ik+1) + u1(ij+1,ik-1) +
     c     u1(ij+1,ik+1)) + .075*(u1(ij-1,ik) + u1(ij,ik+1) +
     c     u1(ij,ik-1) + u1(ij+1,ik)) + .5*u1(ij,ik)
 200       continue


C  *****Compute smoooooooothed boundaries  -- new  (EF, 6/26/95)
C     use more weight to surrounding points for more smoothing, here and for corner points
  
      do 300 ik=2,n-1
          u2(1,ik) = .1*(u1(2,ik-1) + u1(2,ik+1)) + .35*u1(1,ik) +
     c              .15*(u1(1,ik-1) + u1(1,ik+1) + u1(2,ik))
 300  continue

      do 305 ik=2,n-1
          u2(m,ik) = .1*(u1(m-1,ik-1) + u1(m-1,ik+1)) + .35*u1(m,ik) +
     c              .15*(u1(m,ik-1) + u1(m,ik+1) + u1(m-1,ik))
 305  continue


      do 320 ij=2,m-1
          u2(ij,1) = .1*(u1(ij-1,2) + u1(ij+1,2)) + .35*u1(ij,1) +
     c              .15*(u1(ij-1,1) + u1(ij+1,1) + u1(ij,2))
 320  continue

      do 325 ij=2,m-1
          u2(ij,n) = .1*(u1(ij-1,n-1) + u1(ij+1,n-1)) + .35*u1(ij,n) +
     c              .15*(u1(ij-1,n) + u1(ij+1,n) + u1(ij,n-1))
 325  continue


c  *****Compute smoooooooooothed corner points

      u2(1,1) = .16*u1(2,2) + .24*(u1(1,2) + u1(2,1)) + .36*u1(1,1)
      u2(m,1) = .16*u1(m-1,2) + .24*(u1(m,2) + u1(m-1,1)) + .36*u1(m,1)

      u2(1,n) = .16*u1(2,n-1) + .24*(u1(2,n) + u1(1,n-1)) + .36*u1(1,n)
      u2(m,n)= .16*u1(m-1,n-1) +.24*(u1(m-1,n) + u1(m,n-1)) +.36*u1(m,n)
  

        RETURN
      END


C
C   **************************   SMOOTHING  SUBROUTINE - for COUPLED MODEL    *******************************
C                                                         less smoothing than SMOOTH54 above

        SUBROUTINE SMOOTH54C(u1, u2, m, n)
 
c
c  subroutine to smooth 2D array,  u1(m,n),  and  return  smooooothed array u2(m,n)
C
C  *****  weights are:  .03, .045, .7
C     since, eg  x(-1,-1) is sqrt(2) further from x(0,0) than is x(-1,0)
C

        integer m, n
        REAL*4 u1(m,n), u2(m,n)


C  FIRST   SMOOTH all INTERIOR points

       do 200 ik=2,n-1
       do 200 ij=2,m-1

       u2(ij,ik) = .03*(u1(ij-1,ik-1) + u1(ij-1,ik+1) + u1(ij+1,ik-1) +
     c     u1(ij+1,ik+1)) + .045*(u1(ij-1,ik) + u1(ij,ik+1) +
     c     u1(ij,ik-1) + u1(ij+1,ik)) + .7*u1(ij,ik)
 200       continue


C  *****Compute smoooooooothed boundaries  -- new  (EF, 6/26/95)
C     use more weight to surrounding points for more smoothing, here and for corner points
  
      do 300 ik=2,n-1
          u2(1,ik) = .06*(u1(2,ik-1) + u1(2,ik+1)) + .61*u1(1,ik) +
     c               .09*(u1(1,ik-1) + u1(1,ik+1) + u1(2,ik))
 300  continue

      do 305 ik=2,n-1
          u2(m,ik) = .06*(u1(m-1,ik-1) + u1(m-1,ik+1)) + .61*u1(m,ik) +
     c               .09*(u1(m,ik-1) + u1(m,ik+1) + u1(m-1,ik))
 305  continue


      do 320 ij=2,m-1
          u2(ij,1) = .06*(u1(ij-1,2) + u1(ij+1,2)) + .61*u1(ij,1) +
     c               .09*(u1(ij-1,1) + u1(ij+1,1) + u1(ij,2))
 320  continue

      do 325 ij=2,m-1
          u2(ij,n) = .06*(u1(ij-1,n-1) + u1(ij+1,n-1)) + .61*u1(ij,n) +
     c               .09*(u1(ij-1,n) + u1(ij+1,n) + u1(ij,n-1))
 325  continue


c  *****Compute smoooooooooothed corner points

      u2(1,1) = .1*u1(2,2)    + .15*(u1(1,2) + u1(2,1)) +     .6*u1(1,1)
      u2(m,1) = .1*u1(m-1,2)  + .15*(u1(m,2) + u1(m-1,1)) +   .6*u1(m,1)

      u2(1,n) = .1*u1(2,n-1)  + .15*(u1(2,n) + u1(1,n-1)) +   .6*u1(1,n)
      u2(m,n) = .1*u1(m-1,n-1) +.15*(u1(m-1,n) + u1(m,n-1)) + .6*u1(m,n)
  

        RETURN
      END
