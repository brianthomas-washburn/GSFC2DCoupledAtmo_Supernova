
       SUBROUTINE TEMPIN


       include "com2d.h"

C                                                                      common from clouds.inc for NEWRAD9
       real latcld(72), camt(360,3,72), cpres(3,72), ctau(360,3,72)
       real cwcr(360,72,10,2), clfrac(360,72,10), zzcld(10)
       COMMON/CLOUDS1/latcld, camt, cpres, ctau, cwcr, clfrac, zzcld


       REAL*4 ZZIN(37,59), ZZIN18(18,58), ZZINL(18), ZZIN37L(37)
       REAL*4 ZZ7(91,117), ZZ8(91), ZZ9(73,235), ZZ91(L$S,Z$S)
       REAL*4 ZZOUT(L$,Z$X), ZZOUTS(L$S,Z$S), ZZOUTL(L$), ZZOUTSL(L$S)
       REAL*4 ZZOUT1(L$+1,Z$X), llstr(91), zzstr(117), zx1(L$,10)
       REAL*4 otint(19,49), otout(L$+1,Z$), zz10(10), z4521(45,21)

       REAL*4 ozbcin(360, 45, 21), latbcin(45), zzbcin(21), zz6(45)
       REAL*4 depin(14,45,7), latdep(45), emisin(14,45,3), ohout(L$,Z$)
       REAL*4 depr8in(10,14,91), ddlat(91), ddz(16), yy91(91)
       REAL*4 scavin(20,14,91,16), yysc(91,16), yylz(L$,Z$)
       REAL*4 dukzr(14,45,60), latdu(45), zzdu(60), yz45(45,60)
       REAL*4 ohr(14,180,60), latoh(180), zzoh(60), ohint(180,60)
       REAL*4 WHACK(91,117,360)

       REAL*8 sjflux, sjfluxa, hjflux, hjfluxa


cccwconv48       REAL*4 qvt0(46,57,14), qwt0(46,57,14)


C   COMMON for SPARC OH field (#/cm3)
C
       COMMON/CCOH/ohm(14,L$,Z$), ohday(L$,Z$)


C  COMMON FOR Carbon-14 SETVEL is the Settling velocity for Strontium-90
       COMMON/CARB14/C14IN(18,46), SR90(18,46), SETVEL(18,46)

C  
C  COMMONs for COUPLED MODEL X* BBC, fixed model Kzz, Kyy in cm2/sec,   and latent heating from GEOS4 and WACCM

       COMMON/CCHI/BBCOUP(91,21,360), LATBBC(91), ZBBC(21), 
     >             KZZ9FR(45,89,360), latkzz(45), zzkzz(89), 
     >             KYY9FR(46,88,360), latkyy(46), zzkyy(88)


       COMMON/CFKZZ/ylatf(45), zaltf(89), timefz(14)


C                                                            ! WLH, WHACK are in K/sec
       COMMON/CCLH/WLH(91,117,360)
ccwconv175 - OLD       COMMON/CCLH/LHC(91,117,360), WLH(91,117,360), WHACK(91,117,360)

       COMMON/CCHACK/WHACK72(91,117,72)

C                                                   COMMON for coupled model geop hgt wave amps:
       COMMON/CGPBC/ GPBC(6,2,73,30,360), ypgp(73)

C                                                             common for KYYOT on coupled model dynamics grid
       COMMON/CKYYOT/ XKYYOTC(36,46,360), ypot(36), zpot(46)


C                                                       common for fixed model clim TEMP, EHF (K/sec)

       COMMON/CEHFX/tempfx(91,117,74), ehfx(91,117,74),  
     >              latfx(91),zzfx(117)


c   solar cycle variation used in UV heating in COUPLED MODEL

       COMMON/SCYCLE/ SOLCYCHR(5,50,360)


c   solar flux data from Judith Lean - ALL ARRAYS are REAL*8 ;  SJFLUXA are the long term avgs (1954-2008)
C   index 1 = time;  2-40 are the 39 model wavelength bins (in ph/cm^2/sec); 41=TOA flux (W/m2)
C
C   HJFLUX(12,1524) are for the heating rate intervals; HJFLUXA(12) are the 1954-2008 avgs (in mW/m2)
C
       COMMON/CSJFLUX/ SJFLUX(41,1524), SJFLUXA(41),
     >                 HJFLUX(12,1524), HJFLUXA(12)       !  , itj1, itj2



c  common for UBAR climatology, used in both fixed and COUPLED models
C                   and UBAR-CO2 Sensitivity from GEOS5 (in m/sec/ppmv)

       COMMON/CUBAR/ UBARCL(91,117,74)

       COMMON/CUBARCO2/ UBARCO2(91,117)


C  common for fixed model dynamical fields from MERRA data, clim avg for 1979-2010

       COMMON/CLMERRA/ xkyyclm(91,117,74), epclm(91,117,74), 
     >           tempclm(91,117,74), ubarclm(91,117,74),
     >           ehfclm(91,117,74), ehfzclm(91,117,74), qyclm(91,117,74)



c  common for TD GEOS 4 surface temperatures used in COUPLED model, 1950-2100, 
c          g4sens - sensitivity to CO2, w/ seas cycle

       COMMON/CTEMPG4/ TEMPG4(1800,91,3), latg4(91), timeg4(1800), 
     >                 g4sens(15,91,3), timesf(14)



c  common for TD GEOS 4 latent heating used in COUPLED model, 1950-2100, 
c          xlhsens - sensitivity to CO2, w/ seas cycle, and GEOS5 latent heating

       COMMON/CG4LH/ xlhsens(15,91,117), G4LHDAY(91,117)
       COMMON/CG5LH/ xlhsens5(15,91,117)


c  common for TD GEOS 4 tropospheric H2O 1950-2100 in ppmv, seas cycle + sensitivity to CO2

       COMMON/CG4H2O/ g4h2o(16,91,76), latg4h(91), zh76(76), timesfh(14)



c      common for coupled model temperature - NCEP offset (lat-hgt-season) from 1985-2005:

       COMMON/CTOFFS/ xlatoffs(45), zzoffs(76)
cccccccccccc       COMMON/CTOFFS/ TEMPOFFS(360,45,76), xlatoffs(45), zzoffs(76)


c      common for coupled model sfc albedo from TOMS (season-lat) 

       COMMON/CTALBEDO/ TALBEDO(15,180), ALBDAY(180), ALBLAT(180)


c      common for coupled model TBAR/UBAR - NCEP differences (lat-hgt-season) from 1979-2006:

       COMMON/CUTDIFF/ TDIFFR(14,45,76), UDIFFR(14,45,76), xlatf(45), 
     >                 zzff(76)


c      common for NCEP DELF for baroclinic waves ONLY - for coupled model - NOT Currently used

ccbaro       COMMON/CDFBARO/ DFBARO(73,33,360), BALAT(73), BZZ(33)


C
C  COMMON for WACCM cloud parameters:  WWLIQ, WWICE are in kg/kg;  WWCLD is fraction
C        WWRLIQ, WWRICE are in microns  

      COMMON/CWCLOUDS/wwliq(46,33,14), wwice(46,33,14), wwcld(46,33,14),
     >                wwrliq(46,33,14), wwrice(46,33,14),
     >                wwlat(46), wwpr(33), timesw(14)


C   COMMON for WACCM oro Gwaves DELF (m/sec^2)

       COMMON/CWOGW/wogw(46,57,14), platg(46), zzwg(57), tog14(14)


C  COMMON for GEOS 5 TD cloud parameters

       COMMON/CG5CL/g5clin(5,15,37,30), yg5cl(37), zg5cl(30), timecl(14)



C   common for GMI sfc deposition, emissions, interpolated to L$ latitude grid:
C          DEP in cm/sec;   EMIS in #/cm2/sec - need to divide by DELTAZ
C
       COMMON/CDEP/tdep(14), DEP0(14,L$,7), EMIS0(14,L$,3), 
     >             DEPDAY(L$,7), EMISDAY(L$,3)


C   common for CGCM-COMBO T1R8 sfc deposition (dry depos is the TOTAL DEPOSITION)
C          DEPR8 in cm/sec;  - divide by DELTAZ to get 1/sec
C
C  plus the wet scavenging profile from T1R8, SCAVR8(20,14,L$,Z$) is in 1/sec
C
       COMMON/CDEPR8/DEPR8(10,14,L$), SCAVR8(20,14,L$,Z$),
     >               SCAVDAY(20,L$,Z$)


C
C  common for the vertical gradients of the DU/km ozone climatology, for Kzz adjustment
C
       COMMON/CDUKZ/dukz(14,45,89)


C
C  common for HALOE/MLS tropical tropopause H2O (ppmv) for 1991-2010
C
       COMMON/CH2OFM/h2ofm(240,45,76), t240(240)



C
C  COMMON for WACCM convective mass flux for Kzz
C
ccc46       COMMON/CWCONKZ/cmfw(46,25,14), plat(46), zzw(25), tg14(14)


C
C  COMMON for WACCM combined eddy heating rate corrections (theta heating) in K/sec
C
cccwconv48       COMMON/CQWT/qwtt(46,57,14), platw(46), zz57(57), tq14(14)


C
C  COMMON for ERA-40, w'Th'
C
cerw      COMMON/CERWT/erwtr(73,49,74), eelat(73), eez(49), e74(74)


c   common for NCEP v'T', u'v' for baroclinic waves ONLY - for coupled model
C                   vt is in K/sec;  uv is in m/sec^2
ccvtbaro
ccvtbaro       COMMON/CVTBARO/ vt37(73,33,37), uv37(73,33,37),
ccvtbaro     >                 BALAT1(73), BZZ1(33), time37(37)



       SAVE


C
C  time series of tropical tropopause H2O (ppmv) from HALOE/MLS for 1992-2010
C    h2ofm(240,45,76), t240(240)
C        - this was written out in IDL on Linux, so use "little_endian"
C
       OPEN (78, FILE = 'hal_mls_h2o_ttrop.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='stream', STATUS = 'OLD')
           READ (78) h2ofm
           READ (78) t240
        CLOSE (78)



C
C  CLIMATOLOGICAL Geop hgt wave amps (waves 1-5) for coupled model (in METERS) -  ypgp(73) are the latitudes
C  1st index: 1=zonal avg, 2-6 are waves 1-5;                2nd index:  1=AMP,  2=Phase (degree longitude);
C  3rd index: latitudes (90S-90N, by 2.5 deg - NCEP grid);   4th index:  1 - 30 km, by 1 km
C  for 360 days - this is the climatological avg for 1979-2007
C                                   - this was written out in IDL on Linux, so use "little_endian"
C
        OPEN (78, FILE = 'zzbc_2dcoup.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='stream', STATUS = 'OLD')
           READ (78) GPBC
           READ (78) YPGP
        CLOSE (78)


c  read in coupled model TBAR/UBAR - NCEP differences (lat-hgt-season) avgd over 1979-2006:
C     TDIFFR(14,45,76), UDIFFR(14,45,76), xlatf(45), zzff(76) in COMMON above (in K, m/sec)
C                                       - this was written out in IDL on Linux, so use "little_endian"
C
        OPEN (78, FILE = 'utdiff.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='stream', STATUS = 'OLD')
           READ (78) TDIFFR
           READ (78) UDIFFR
           READ (78) XLATF
           READ (78) ZZFF
        CLOSE (78)



c READ IN  NCEP DELF for baroclinic waves ONLY - for coupled model (in m/sec^2) - NOT currently used
C   DFBARO(73,33,360), BALAT(73), BZZ(33) are in COMMON above
C                                      - this was written out in IDL on Linux, so use "little_endian"

ccbaro        OPEN (78, FILE = '/misc/kah03/chj25/nmc07/nmc_barocl_2d.xdr',
ccbaro     >        convert="little_endian", FORM = 'unformatted', 
ccbaro     >        ACCESS='SEQUENTIAL', STATUS = 'OLD')
ccbaro           READ (78) DFBARO
ccbaro           READ (78) BALAT
ccbaro           READ (78) BZZ
ccbaro        CLOSE (78)




c READ IN  NCEP v'T', u'v' for baroclinic waves ONLY - for coupled model
C   vt37(73,33,37), uv37(73,33,37), BALAT1(73), BZZ1(33), time37(37)  are in COMMON above
C                                      - this was written out in IDL on Linux, so use "little_endian"
C   vt is in K/sec;  uv is in m/sec^2
C
ccvtbaro
ccvtbaro        OPEN (78, FILE = '/misc/kah03/chj25/nmc07/nmc_barocl_vt_2d.xdr',
ccvtbaro     >        convert="little_endian", FORM = 'unformatted', 
ccvtbaro     >        ACCESS='SEQUENTIAL', STATUS = 'OLD')
ccvtbaro           READ (78) vt37
ccvtbaro           READ (78) uv37
ccvtbaro           READ (78) BALAT1
ccvtbaro           READ (78) BZZ1
ccvtbaro           READ (78) time37
ccvtbaro       CLOSE (78)



C
C  Read in  latent heating  from  Newell, 1974  10-DAY MEANS;  LHIN(37,7,36)
c
        OPEN (78, FILE = 'latheat_base8d.xdr', 
     >    convert="big_endian", 
     >    FORM = 'unformatted', ACCESS='stream',STATUS = 'OLD')
           READ (78) LHIN
        CLOSE (78)

C  LH10(L$S,Z$S,36) is for STREAMF, Z$S levels, ZZOUTS(L$S,Z$S)

        do 330 iim=1,36

          do 255 ij=1,37
            do 257 ik=1,7  
 257	       zzin(ij,ik) = lhin(ij,ik,iim)
            do 259 ik=8,59 
 259	       zzin(ij,ik) = 0.0
 255	  CONTINUE

            CALL BINTERP(xlat5, 37, zz59, 59, zzin, 
     >                   LATST4, L$S, ZSTR4, Z$S, 0, 0, ZZOUTS)

            DO 535 ik=1,Z$S
            DO 535 ij=1,L$S
 535	      lh10(ij,ik,iim) = zzouts(ij,ik)
 330	 CONTINUE


C
C  read in LATENT HEATING from GCM output, LHGCM(L$S,Z$S,72) in K/day,  in COMMON on STREAMF grid 
C                                                                
        OPEN (78, FILE = 'gcmlh_9ba.xdr', 
     >        convert="big_endian", FORM = 'unformatted', 
     >        ACCESS='stream', STATUS = 'OLD')
           READ (78) LHGCM
        CLOSE (78)

C
C  also read in daily values for coupled model - LHC(91,117,360) in K/day, in COMMON ABOVE on STREAMF grid
C                                           - this was written out in IDL on Linux, so use "little_endian"
C
ccwconv175 - NO LONGER USED
ccwconv175
ccwconv175        OPEN (78, FILE = '/misc/kah03/chj25/gmi/gcmlh_9ba_day.xdr', 
ccwconv175     >        convert="little_endian", FORM = 'unformatted', 
ccwconv175     >        ACCESS='SEQUENTIAL', STATUS = 'OLD')
ccwconv175           READ (78) LHC
ccwconv175        CLOSE (78)



C
C  read in daily values of WACCM latent heating for coupled model - 
C    WLH, WHACK(91,117,360) in K/sec, in COMMON ABOVE on STREAMF grid
C            - this was written out in IDL on Linux, so use "little_endian"
C
C   this was cleaned up for WCONV147, with WLH and HACK combined into one array (offline)
C                         WLH is in K/sec
C
      OPEN (78, FILE = 'waccm_lh_coup.xdr', 
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='stream', STATUS = 'OLD')
           READ (78) WLH
ccwconv147           READ (78) WHACK
        CLOSE (78)



C
C  for MERRA IAV run, add in HACK latent heating to get a bit faster circ/younger age
C  so read in old 2008 file to get WHACK(91,117,360) in K/sec (the 2nd array)
C     WHACK is in COMMON ABOVE on STREAMF grid, create 5-day avgs for use in STREAMF
C            - this was written out in IDL on Linux, so use "little_endian"
C
      OPEN(78,FILE='waccm_lh_coup_2008.xdr'
     >      , convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='stream', STATUS = 'OLD')
           READ (78) WHACK
           READ (78) WHACK
        CLOSE (78)



C   WHACK72(91,117,72) are 5-day averages in K/sec

        do 7710 ik=1,117
        do 7710 ij=1,91

          do 7715 imm=1,72
            ib = imm*5 - 4
            WHACK72(ij,ik,imm) = (whack(ij,ik,ib) + whack(ij,ik,ib+1) + 
     >     whack(ij,ik,ib+2) + whack(ij,ik,ib+3) + whack(ij,ik,ib+4))/5.
 7715     CONTINUE

 7710   CONTINUE



C
C   read in daily CLIMATOLOGICAL values of X* BBC for COUPLED MODEL - written out in IDL on Linux
C      LATBBC(91), ZBBC(21), BBCOUP(91,21,360) in m2/sec, all in COMMON ABOVE on STREAMF grid
C

       OPEN (78, FILE = 'chi_coupled_9gb.xdr', 
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
           READ (78) BBCOUP
           READ (78) LATBBC
           READ (78) ZBBC
       CLOSE (78)


C
C   read in daily CLIMATOLOGICAL Kzz from fixed model 9FR run for COUPLED MODEL - written out in IDL on Linux
C      latkzz(45), zzkzz(89), KZZ9FR(45,89,360) in cm2/sec, all in COMMON ABOVE
C

       OPEN (78, FILE = 'kzz_9fr_coup.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
           READ (78) KZZ9FR
           READ (78) LATKZZ
           READ (78) ZZKZZ
       CLOSE (78)



C   load latkzz(45), zzkzz(89) into arrays ylatf(45), zaltf(89) 
C      for use in SETDAILY for fixed model Kzz tropspheric adjustment
C     (to be consistent w/ coupled model)

           do 888 ij=1,45
 888          ylatf(ij) = latkzz(ij)

           do 889 ik=1,89
 889          zaltf(ik) = zzkzz(ik)



C
C   read in daily CLIMATOLOGICAL Kyy from fixed model 9FR run for COUPLED MODEL - written out in IDL on Linux
C       latkyy(46), zzkyy(88), KYY9FR(46,88,360) in cm2/sec, all in COMMON ABOVE
C

       OPEN (78, FILE = 'kyy_9fr_coup.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
           READ (78) KYY9FR
           READ (78) LATKYY
           READ (78) ZZKYY
       CLOSE (78)



C
C   read in climatological dynamical fields from fixed model 9FR run, on STREAMF grid
C     generated in IDL on Linux, so use "little_endian"
C
C       tempfx(91,117,74), ehfx(91,117,74) is in K/sec,  latfx(91), zzfx(117) all in COMMON above
C       UBARCL(91,117,74), QYCL(91,117,72) (1.e11 1/m/sec), EPCL(91,117,72) (m/sec/day) in COMMON
C 
        OPEN (78, FILE = 'ubarclim_9fr.xdr', 
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
           READ (78) tempfx
           READ (78) ehfx
           READ (78) UBARCL
           READ (78) EPCL
           READ (78) QYCL
           READ (78) latfx
           READ (78) zzfx
        CLOSE (78)


C
C   read in climatological dynamical fields from MERRA data (9GT run), on STREAMF grid
C     generated in IDL on Linux, so use "little_endian"
C
C     xkyyclm(91,117,74), epclm(91,117,74), tempclm(91,117,74), ubarclm(91,117,74)
C     ehfclm(91,117,74), ehfzclm(91,117,74), qyclm(91,117,74) all in COMMON above
C
C     QYCLM in  1.e11 1/m/sec;  EPCLM in m/sec/day;  EHF, EHFZ in K/day
C 
       OPEN (78, FILE='kyy_merra_clim.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
           READ (78) xkyyclm
           READ (78) epclm
           READ (78) tempclm
           READ (78) ubarclm
           READ (78) ehfclm
           READ (78) ehfzclm
           READ (78) qyclm
        CLOSE (78)



C
C  read in GEOS4/NCEP surface temperatures used for COUPLED model radiation scheme (GEOS 5 scheme);
C   TEMPG4(1800,91,3), latg4(91), timeg4(1800)  - this was written out in IDL on Linux, so use "little_endian"
C
        OPEN (78, FILE = 'g4temp_2dmod.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
           READ (78) tempg4
           READ (78) latg4
           READ (78) timeg4
        CLOSE (78)


C
C  read in GEOS4/NCEP surface temperature sensitivity to CO2  used for COUPLED model radiation scheme (GEOS 5 scheme);
C   G4SENS(15,91,3), latg4(91)  - this was written out in IDL on Linux, so use "little_endian"
C
        OPEN (78, FILE = 'g4temp_sens.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
           READ (78) g4sens
           READ (78) latg4
        CLOSE (78)



C
C  read in GEOS4 latent heating sensitivity to CO2  used for COUPLED model
C     xlhsens(15,91,117)  - this was written out in IDL on Linux, so use "little_endian"
C          this was cleaned up for wconv147: 
C
        OPEN (78, FILE = 'g4_lh_2dmod.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
           READ (78) xlhsens
        CLOSE (78)


C
C  read in GEOS5 latent heating sensitivity to CO2  used for COUPLED model - in K/day
C     xlhsens5(15,91,117)  - this was written out in IDL on Linux, so use "little_endian"
C       THIS has also been cleaned up as in WCONV147
C
        OPEN (78, FILE = 'geos5_lh_2dmod.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
           READ (78) xlhsens5
        CLOSE (78)


C
C  read in GEOS4 H2O (ppmv) + sensitivity to CO2    
C     g4h2o(16,91,76), latg4h(91), zh76(76)  - this was written out in IDL on Linux, so use "little_endian"
C
        OPEN (78, FILE = 'g4_h2o_2dmod.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
           READ (78) g4h2o
           READ (78) latg4h
           READ (78) zh76
        CLOSE (78)


C                         define zzoffs(76)
       do 221 ik=1,76
          zzoffs(ik) = zh76(ik)
 221   CONTINUE


C
C    also define time index for 14 months:  timesf(14) in COMMON,  for G4SENS(15,91,3), and XLHSENS(15,91,117)
C
        do 700 itj=1,14
           timesw(itj)  = itj*30. - 15. - 30.
ccc46           tg14(itj)  = timesw(itj)
ccc48           tq14(itj)  = timesw(itj)

           timecl(itj)  = timesw(itj)
           timefz(itj)  = timesw(itj)

           tog14(itj)  = timesw(itj)
           timesf(itj)  = timesw(itj)
           tdep(itj) = timesw(itj)
 700       timesfh(itj) = timesw(itj)


C
C  read in TOMS monthly mean surface Albedo UV-VIS used for COUPLED model
C    talbedo(15,180), 1-14 are months, 15=latitude  - this was written out in IDL on Linux, so use "little_endian"
C
        OPEN (78, FILE = 'srf_rflec.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
           READ (78) talbedo
        CLOSE (78)



C
C  read in coupled model temperature - NCEP offset (lat-hgt-season) from 1985-2005: for use in coupled model TPROB
C     TEMPOFFS(360,45,76), xlatoffs(45), zzoffs(76) in COMMON 
C       - this was written out in IDL on Linux, so use "little_endian"
C
ccccc        OPEN (78, FILE = '/misc/kah05/fleming/nmc07k/tempoffs.xdr',
ccccc     >        convert="little_endian", FORM = 'unformatted', 
ccccc     >        ACCESS='SEQUENTIAL', STATUS = 'OLD')
ccccc           READ (78) tempoffs
ccccc           READ (78) xlatoffs
ccccc           READ (78) zzoffs
ccccc        CLOSE (78)

C
C
C  Read in  new NMC tropopause heights, which have been interpolated to daily values
C        array  TROPHTIN(18,360) is in COMMON  (these are tropopause pressures in mbar, range: 333.128-98 mb)
c 
         OPEN (74, FILE = 'trop2d_base8e.xdr', 
     >     convert="big_endian",
     >     form = 'unformatted', access = 'stream', status = 'old')
         READ (74) TROPHTIN
         CLOSE (74)

C  TROPHTF(L$,360),  LAT4(L$) , interpolate to new grid logarithmically

        do 340 iid=1,360
          do 265 ij=1,18
 265	       zzinl(ij) = trophtin(ij,iid)

            CALL LINTERP(LATIN, 18, ZZINL, LAT4, L$, 1, ZZOUTL)

            DO 545 ij=1,L$
 545	      trophtf(ij,iid) = zzoutl(ij)
 340	 CONTINUE


C
C  READ IN daily tropospheric H2O values from Oort, 1983 -- TROPWVIN(18,11,360) defined in COMMON
C
        OPEN (560, FILE = 'troph2oday_2d.dat', convert="big_endian",
     >      FORM = 'UNFORMATTED', ACCESS='STREAM', STATUS = 'old')
        READ (560) TROPWVIN
        CLOSE (560)

C
C  TROPWV(L$,Z$X,360), ZALT(Z$X), ZZOUT(L$,Z$X), ZZIN18(18,58), defined at chemistry model grid points
Ctrop      don't use this to save memory - OOrt clim now contained in HALOE clim read in below
Ctrop   
Ctrop        do 350 iim=1,360
Ctrop
Ctrop          do 275 ij=1,18
Ctrop            do 276 ik=1,11 
Ctrop 276	       zzin18(ij,ik) = TROPWVIN(ij,ik,iim)
Ctrop            do 277 ik=12,58
Ctrop 277	       zzin18(ij,ik) = 0.0
Ctrop 275	  CONTINUE
Ctrop
Ctrop            CALL BINTERP(LATIN, 18, zz58, 58, zzin18, 
Ctrop     >                   LAT4, L$, ZALT, Z$X, 0, 0, ZZOUT)
Ctrop
Ctrop            DO 555 ik=1,Z$X
Ctrop            DO 555 ij=1,L$
Ctrop 555	      TROPWV(ij,ik,iim) = zzout(ij,ik)
Ctrop 350	 CONTINUE


C
C   READ in file of bottom boundary scaling factor for gravity wave model 
C     (based on Wu and Waters MLS results at 33 km), GWBCIN(13,37,72) defined in COMMON
C
       OPEN (63,FILE = 'gwbc_base9ba.xdr', 
     >    convert="big_endian",
     >    FORM = 'UNFORMATTED', ACCESS='STREAM', STATUS = 'OLD')
         READ (63) GWBCIN
       CLOSE (63)

c                             ZZIN37L(37); ZZOUTSL(L$S) , GWBC(13,L$S,72)
        do 360 iid=1,72
        do 360 iix=1,13
          do 295 ij=1,37
 295	       zzin37l(ij) = GWBCIN(iix,ij,iid)

            CALL LINTERP(XLAT5, 37, ZZIN37L, LATST4, L$S, 0, ZZOUTSL)

            DO 445 ij=1,L$S
 445	      GWBC(iix,ij,iid) = zzoutsl(ij)
 360	 CONTINUE


C
C  Read in bottom boundary conditions NOW FOR LEVEL 2 (762 mb) for stream function;
C        array  BBCADJ(37,36) is REAL*8 and is in COMMON  - this read-in table is NO LONGER USED
c 
CEF       OPEN (75, FILE = '/misc/kah02/newmod/base8x/chibbc2_base8xf.xdr',
CEF     *     form = 'system', access = 'sequential', status = 'old')
CEF         READ (75) BBCADJ
CEF         CLOSE (75)
CEF

C
C  Read in bottom boundary conditions for CO2; array CO2BCIN(18,101,360) is in COMMON
C    daily values for 1950-2050 (101 years)
c 
       OPEN (79, FILE = 'co2bc_day.xdr', 
     >     convert="big_endian",
     >     form = 'unformatted', access = 'stream', status = 'old')
         READ (79) CO2BCIN
       CLOSE (79)

!-- BThomas 2018 - changing the CO2 BCs for paleo (pre-industrial, ish) conditions.  
!   Basically, original data (co2bc_day.xdr) increases over the 1950-2050 range.  
!   I want it to be constant.  So, taking first year (360 days, 18 lat points), 
!    and copying into every year.
       do iid=1,360
          do iiy=2,101
             do ij=1,18
                co2bcin(ij,iiy,iid) = co2bcin(ij,1,iid)
             end do
          end do
       end do
!-- end BThomas edit --

C  CO2BC(L$,101,360),  LAT4(L$) ,  ZZINL(18), ZZOUTL(L$),

        do 380 iid=1,360
        do 380 iiy=1,101
          do 665 ij=1,18
 665	       zzinl(ij) = co2bcin(ij,iiy,iid)

            CALL LINTERP(LATIN, 18, ZZINL, LAT4, L$, 0, ZZOUTL)

            DO 945 ij=1,L$
 945	      co2bc(ij,iiy,iid) = zzoutl(ij)
 380	 CONTINUE


C
C  Read in daily tropical tropopause boundary conditions for H2O from HALOE, both yearly and 
C      climatological values ;  arrays H2OBC(18,360,9), H2OBCCL(18,360) are in COMMON
c 
       OPEN (79, FILE = 'h2obc.xdr', 
     >       convert="big_endian",
     >       form = 'unformatted', access = 'stream', 
     >       status = 'old')
         READ (79) H2OBC 
         READ (79) H2OBCCL
       CLOSE (79)


C
C  read in daily ERA40/HALOE High res H2O climatology - used for TROPICAL TROPOPAUSE
C     halh2or(45,92,360), lathal(45), zzhal(92), are in COMMON
C     to save memory, this is now interpolated to proper grid each day in RAINOUT
C                                - this was written out in IDL on Linux, so use "little_endian"

       OPEN (79, FILE = 'h2o_uars_era40.xdr',
     >       convert="little_endian",
     >       form = 'unformatted', access = 'stream', 
     >       status = 'old')
         READ (79) halh2or
         READ (79) lathal
         READ (79) zzhal
       CLOSE (79)



c  
C  READ IN initial conditions for Carbon-14 and Sr-90, and load in
C     C14IN(18,46), SR90(18,46), SETVEL(18,46) in COMMON block CARB14 ABOVE
C
        OPEN (61, FILE = 'c14bc.xdr', 
     >        convert="big_endian",
     >        FORM = 'UNFORMATTED', ACCESS='STREAM', STATUS = 'old')
          READ (61) C14IN
          READ (61) SR90
          READ (61) SETVEL
        CLOSE (61)

C
C  now INTERPOLATE C14IN(18,46) to NEW GRID,  LATIN(18), ZZ46(46), all defined in COMMON and are REAL*4
C       LAT4(L$),  ZALT90(Z$) are both REAL*4,   C14(L$,Z$) is in COMMON (REAL*4)
C
C  C14 is Mix rat, so do LINEAR INTERPOLATION in both Lat and pres/alt

          CALL BINTERP(latin, 18, zz46, 46, C14IN, 
     >                 LAT4, L$, ZALT90, Z$, 0, 0, C14)


C
C  read in HARVARD Ox dry deposition velocities, zonally averaged for each model day, in m/sec
C     HARVARD_DDOX(91,360), LAT_DDOX(91),  in COMMON
C
        OPEN (78, 
     >    FILE = 'harvard_ddox.xdr', 
     >    convert="big_endian",
     >    FORM = 'unformatted', ACCESS='stream', STATUS = 'OLD')
           READ (78) HARVARD_DDOX
           READ (78) LAT_DDOX
        CLOSE (78)

C  interpolate to model latitudes - OXDEP(L$,360) in COMMON ,   ZZ8(91),   ZZOUTL(L$)

        do 440 iim=1,360

          do 577 ij=1,91
 577	       zz8(ij) = HARVARD_DDOX(ij,iim)

            CALL LINTERP(LAT_DDOX, 91, ZZ8, LAT4, L$, 0, ZZOUTL)

            DO 735 ij=1,L$
 735	      OXDEP(ij,iim) = ZZOUTL(ij)
 440	 CONTINUE



C
C  read in ozone climatology (ppmv) for each model day,  INTERPOLATE to OZBC(360,L$,10) in COMMON
C    
C    OZBCIN(360,45,21), latbcin(45) - every 4 degrees (88S-88N), zzbcin(21) 0.5-20.5 km
C
        OPEN (477,  FILE = 'ozbc.xdr', 
     >    convert="big_endian",
     >    FORM = 'unformatted', ACCESS='STREAM', STATUS = 'OLD')
           READ (477) OZBCIN
           READ (477) latbcin
           READ (477) zzbcin
        CLOSE (477)

C
C  interpolate to model grid for lowest 10 levels (whatever they are)
C   zz10(10), z4521(45,21) -> zx1(L$,10), LAT4(L$);  OZBC(360,L$,10) in COMMON (ppmv)

        do 1333 ik=1,10
 1333      zz10(ik) = zalt(ik)


        do 1440 iim=1,360

          do 1577 ik=1,21
          do 1577 ij=1,45
 1577	       z4521(ij,ik) = OZBCIN(iim,ij,ik)

             CALL BINTERP(LATBCIN, 45, ZZBCIN, 21, z4521,
     >                    LAT4, L$, ZZ10, 10, 0, 0, ZX1)

            DO 1735 ik=1,10
            DO 1735 ij=1,L$
 1735	      OZBC(iim,ij,ik) = zx1(ij,ik)

 1440	 CONTINUE


C
C   read in OH climatology (#/cm3), INTERPOLATE to model grid
C        written out in LINUX IDL, so use "little_endian" 
C   ohr(14,180,60), latoh(180), zzoh(60) -> ohm(14,L$,Z$)
C
        OPEN (477,  FILE = 'oh2d.xdr',
     >         convert="little_endian", FORM = 'unformatted', 
     >         ACCESS='STREAM', STATUS = 'OLD')
           READ (477) ohr
           READ (477) latoh
           READ (477) zzoh
        CLOSE (477)

c             ohint(180,60); ohout(L$,Z$); LAT4(L$); ZALT90(Z$)

        do 2330 iim=1,14
          do 2257 ik=1,60
          do 2257 ij=1,180
 2257	       ohint(ij,ik) = ohr(iim,ij,ik)

            CALL BINTERP(latoh, 180, zzoh, 60, ohint,
     >                   LAT4, L$, ZALT90, Z$, 0, 0, ohout)

C  set to nominal value above 60 km (it's not used)

           DO 2537 ik=1,Z$
           DO 2537 ij=1,L$
       	      ohm(iim,ij,ik) = ohout(ij,ik)
              if (zalt90(ik) .ge. 60.) ohm(iim,ij,ik) = 1.e6
 2537	   CONTINUE
 2330   CONTINUE



C interpolate to model streamfunction grid,  CHI_REAN(L$S,Z$S,72,24) IN COMMON for STREAMF, in m^2/sec 
C    L$S, Z$S levels, ZZOUTS(L$S,Z$S), ZZ9(73,235) , latitudes go 90S-90N, poles = 0.0
cc
cc        do 230 iiy=1,24
cc        do 230 iim=1,72
cc 
cc           do 777 ik=1,235
cc           do 777 ij=1,73
cc  777	       zz9(ij,ik) = CHI7(ij,ik,iim,iiy)
cc 
cc             CALL BINTERP(LAT7, 73, ZZZ7, 235, zz9, 
cc      >                   LATST4, L$S, ZSTR4, Z$S, 0, 0, ZZOUTS)
cc 
cc            DO 733 ik=1,Z$S
cc            DO 733 ij=1,L$S
cc  733	      CHI_REAN(ij,ik,iim,iiy) = zzouts(ij,ik)
cc  230	 CONTINUE


C
C  Read in daily solar cycle variation for 1950-2006 (57 years) as a function of 
C    the 39 wavelengths of the model (IL$=39);   SOLCYCR(IL$,360,57) is in COMMON
c 
       OPEN (177, 
     >   FILE = 'solcyc_1950-2006.xdr', 
     >       convert="big_endian",
     >       form = 'unformatted', access = 'stream', 
     >       status = 'old')
         READ (177) solcycr
       CLOSE (177)



c
c read in solar cycle UV variation from off-line calculations FOR COUPLED MODEL HEATING RATES, based on F10.7
c   solcyc(5,50,360) is in COMMON ABOVE, where 50 is year, 1957-2006  for 360 model days
C   and the 1st index is the UV band:  1=SR Bands;   2=Herzberg;  3=Hartley;   4=Huggins-1;    5=Huggins-2
C
C                                     written out in LINUX IDL, so use "little_endian" (EF, NOv. 2007)
       OPEN (523, 
     >   file='solcyc_var-57-06_360.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
          READ (523) solcychr
       CLOSE (523)



C
C  Read in monthly solar flux variations for 1882-2008 from Judith Lean
C     SJFLUX(41,1524), SJFLUXA(41)  -  written out in LINUX IDL, so use "little_endian"

       OPEN (177, FILE = 'phot_flux.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
         READ (177) sjflux
         READ (177) sjfluxa
       CLOSE (177)


C
C  Read in monthly solar flux variations for 1882-2008 from Judith Lean - for heating rate grid
C  HJFLUX(12,1524), HJFLUXA(12) - values are in mW/m2;  written out in LINUX IDL, so use "little_endian"

       OPEN (177, FILE = 'tsi_heat.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
         READ (177) hjflux
         READ (177) hjfluxa
       CLOSE (177)



C  for SORCE steady state runs: define ITJ1, ITJ2 for Judith Lean's data:
C    use May 2003 (<310), April 2004(>310) for MAX, July 2008 for MIN - ITS1/ITS2
C
cccc9hb        itj1 = (2003-1882)*12 + 5   ! (2008-1882)*12 + 7
cccc9hb        itj2 = (2004-1882)*12 + 4   ! (2008-1882)*12 + 7

cccc9hb        its1 = 5                    ! 67   (July 2008)
cccc9hb        its2 = 16                   ! 67



c 
c  read in JOan's cloud data parameters for COUPLED MODEL HEATING rates (NEWRAD9), interpolated to 360 days
C       camt(360,3,72), ctau(360,3,72),  cpres(3,72),  latcld(72) all in COMMON ABOVE - Nov. 2007 (EF)
C                                                            written out in LINUX IDL, so use "little_endian"
C
        OPEN (677, FILE = 'cloud_360.dat',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
           READ (677) latcld
           READ (677) camt
           READ (677) cpres
           READ (677) ctau
        CLOSE (677)


c 
c  also read in cloud water content derived from CTAU above, and cloud fraction for COUPLED MODEL HEATING rates for 360 days
C    cwcr(360,72,10,2), clfrac(360,72,10), zzcld(10) are the altitudes, in COMMON ABOVE - Nov. 2008 (EF)
C                                                            written out in LINUX IDL, so use "little_endian"

        OPEN (677, FILE = 'cwc_360.xdr',
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
           READ (677) cwcr
           READ (677) clfrac
           READ (677) zzcld
        CLOSE (677)




C
C  read in WACCM cloud parameters:
C    wwliq(46,33,14), wwice(46,33,14), wwcld(46,33,14), wwrliq(46,33,14), wwrice(46,33,14)
C    WWLIQ, WWICE are in kg/kg;  WWCLD is fraction ; wwlat(46), wwpr(33) all in COMMON
C    WWRLIQ, WWRICE are in microns                  written out in LINUX IDL, so use "little_endian"

      OPEN(677,FILE='waccm_clouds.xdr',
     >         convert="little_endian", FORM = 'unformatted', 
     >         ACCESS='STREAM', STATUS = 'OLD')
           READ (677) wwliq
           READ (677) wwice
           READ (677) wwcld
           READ (677) wwrliq
           READ (677) wwrice

           READ (677) wwlat
           READ (677) wwpr
      CLOSE (677)


C
C  read in WACCM convective mass flux for Kzz
C    cmfw(46,25,14) is in kg/m2/sec,  plat(46), zzw(25), all in COMMON above
C
C     written out in LINUX IDL, so use "little_endian"
C
ccc46      OPEN(677,FILE='/misc/kah04/fleming/waccm/gwave/waccm_cmflx.xdr',
ccc46     >         convert="little_endian", FORM = 'unformatted', 
ccc46     >         ACCESS='SEQUENTIAL', STATUS = 'OLD')
ccc46           READ (677) cmfw
ccc46           READ (677) plat
ccc46           READ (677) zzw
ccc46      CLOSE (677)



C  read in WACCM orographic gravity wave DELF (m/sec^2)
C    wogw(46,57,14), platg(46), zzwg(57) all in COMMON above
C
C     written out in LINUX IDL, so use "little_endian"
C
      OPEN (677,file='waccm_gwaves.xdr',
     >         convert="little_endian", FORM = 'unformatted', 
     >         ACCESS='STREAM', STATUS = 'OLD')
           READ (677) wogw
           READ (677) platg
           READ (677) zzwg
      CLOSE (677)


C
C  read in WACCM  eddy heating rate (theta heating)
C  qvt0(46,57,14), qwt0(46,57,14) definded above ,  platw(46), zz57(57) in COMMON
C  NOTE: the 1st qvt0 is epw (u'w' - not used), so read qvt0 2X to save memory
C
C    w(46,57,14), QWT(46,57,14), QVT(46,57,14),plat(46), zzw(57) all in COMMON above
C
C    QWT, QVT (K/sec) -  written out in LINUX IDL, so use "little_endian"
C
cccwconv48      OPEN(677,FILE='/misc/kah04/fleming/waccm/gwave/waccm_conv.xdr',
cccwconv48     >         convert="little_endian", FORM = 'unformatted', 
cccwconv48     >         ACCESS='SEQUENTIAL', STATUS = 'OLD')
cccwconv48           READ (677) qvt0
cccwconv48           READ (677) qvt0
cccwconv48           READ (677) qwt0
cccwconv48           READ (677) platw
cccwconv48           READ (677) zz57
cccwconv48      CLOSE (677)


C   just use vertical eddy heating (theta heating) in K/sec (qwt0)
C                                   qwtt(46,57,14) is in COMMON
cccwconv48            do 717 im1=1,14
cccwconv48            do 717 ik=1,57
cccwconv48            do 717 ij=1,46
cccwconv48 717           qwtt(ij,ik,im1) = qwt0(ij,ik,im1)   ! + qvt0(ij,ik,im1) 



C
C  read in ERA-40, theta heating rate, w'th': erwtr(73,49,74) in K/sec (theta heating)
C     eelat(73), eez(49), e74(74) all in COMMON above, define time index below
C
C              -  written out in LINUX IDL, so use "little_endian"
C
cerw      OPEN(677,FILE='/misc/kah03/chj26/era04/era_wt.xdr',
cerw     >         convert="little_endian", FORM = 'unformatted', 
cerw   >         ACCESS='SEQUENTIAL', STATUS = 'OLD')
cerw           READ (677) erwtr
cerw           READ (677) eelat
cerw           READ (677) eez
cerw      CLOSE (677)
cerw
cerw
cerw        do 7117 itj=1,74
cerw 7117       e74(itj)  = itj*5. - 7.


C  read in UBAR - CO2 sensitivity factor from GEOS 5 - UBARCO2(91,117) in m/sec/ppmv
c                            written out in LINUX IDL, so use "little_endian"
C
      OPEN(677,FILE='geos5_ubarsens_2dmod.xdr',
     >         convert="little_endian", FORM = 'unformatted', 
     >         ACCESS='STREAM', STATUS = 'OLD')
           READ (677) UBARCO2
      CLOSE (677)

 

C
C  read in GEOS5 TD cloud parameters:   written out in LINUX IDL, so use "little_endian"
C
C    g5clin(5,15,37,30), yg5cl(37), zg5cl(30) - all in COMMON
C       1 = mass_fraction_of_cloud_liquid_water (kg/kg)
C       2 = mass_fraction_of_cloud_ice_water    (kg/kg)
C       3 = cloud_fraction
C       4 = liquid_cloud_particle_effective_radius (meters)
C       5 = ice_cloud_particle_effective_radius (meters)
C
C  2nd index:  1-14 are 14 month seasonal cycle; 
C              15 = CO2 sensitivity (in cloud parameter/ppmv)
C
      OPEN(677,FILE='g5_clouds_2dmod.xdr',
     >         convert="little_endian", FORM = 'unformatted', 
     >         ACCESS='STREAM', STATUS = 'OLD')
           READ (677) g5clin
           READ (677) yg5cl
           READ (677) zg5cl
      CLOSE (677)


C
C
C  *********************************************************************************
C
C
C   read in UT/LS Kyys from OT (PM,1987), interpolate to model grid
C      kyyotr(19,49,72), latot(19), zzot(49), kyyot(L$+1,Z$,72) is in cm2/sec - all are in COMMON
C
       OPEN (275, FILE = 'kyyot.xdr', 
     >     convert="big_endian",
     >     form = 'unformatted', access = 'stream', status = 'old')
         READ (275) kyyotr
         READ (275) latot
         READ (275) zzot
       CLOSE (275)

c                                      ;  otint(19,49);  otout(L$+1,Z$); LATEG4(L$+1); ZALT90(Z$)
        do 1330 iim=1,72
          do 1257 ik=1,49  
          do 1257 ij=1,19
 1257	       otint(ij,ik) = kyyotr(ij,ik,iim)

            CALL BINTERP(latot, 19, zzot, 49, otint,
     >                   LATEG4, L$+1, ZALT90, Z$, 0, 0, otout)

            DO 1535 ik=1,Z$
            DO 1535 ij=1,L$+1
 1535	      kyyot(ij,ik,iim) = otout(ij,ik)
 1330	 CONTINUE



C
C   read in UT/LS Kyys from OT (PM,1987), interpolated to COUPLED model DYNAMICS grid - daily values
C      xkyyotc(36,46,360) in cm2/sec, ypot(36), zpot(46)  in COMMON above
C                                            - - this was written out in IDL on Linux, so use "little_endian"
C
       OPEN (275, FILE = 'kyyot-2dcoup.xdr', 
     >        convert="little_endian", FORM = 'unformatted', 
     >        ACCESS='STREAM', STATUS = 'OLD')
         READ (275) xkyyotc
         READ (275) ypot
         READ (275) zpot
       CLOSE (275)




C
C   read in GMI surface deposition, emissions (from 2004-2007 avg of AURA run)
C   interpolate to current latitude grid - written out in IDL on Linux, so use "little_endian"
C               depin(14,45,7), latdep(45)
C
       OPEN (177, FILE = 'dep_2d_2011.xdr',
     >       convert="little_endian", form = 'unformatted', 
     >       ACCESS = 'STREAM', STATUS = 'OLD')
         READ (177) depin
       CLOSE (177)

C                         define latdep, xlatoffs(45)
       do 121 ij=1,45
          latdep(ij) = (ij-1)*4. - 88.
          xlatoffs(ij) = latdep(ij)
 121   CONTINUE




C  interpolate depositions to model latitudes - DEP0(14,L$,7) in COMMON in cm/sec
C                           -->   zz6(45), LAT4(L$) ==>  ZZOUTL(L$)

        do 2440 ik=1,7
        do 2440 iim=1,14

          do 1677 ij=1,45
 1677	       zz6(ij) = DEPIN(iim,ij,ik)

            CALL LINTERP(LATDEP, 45, ZZ6, LAT4, L$, 0, ZZOUTL)

            DO 1835 ij=1,L$
 1835	      DEP0(iim,ij,ik) = ZZOUTL(ij)

 2440	 CONTINUE


C
C   GMI surface EMISSIONS - written out in IDL on Linux, so use "little_endian"
C       emisin(14,45,3), latdep(45) -> EMIS0(14,L$,3) in COMMON in #/cm2/sec
C
       OPEN (177, FILE = 'emis_2d.xdr',
     >       convert="little_endian", form = 'unformatted', 
     >       ACCESS = 'STREAM', STATUS = 'OLD')
         READ (177) emisin
       CLOSE (177)


        do 3440 ik=1,3
        do 3440 iim=1,14

          do 1777 ij=1,45
 1777	       zz6(ij) = EMISIN(iim,ij,ik)

            CALL LINTERP(LATDEP, 45, ZZ6, LAT4, L$, 0, ZZOUTL)

            DO 1935 ij=1,L$
 1935	      EMIS0(iim,ij,ik) = ZZOUTL(ij)

 3440	 CONTINUE


C
C   read in surface DRY (TOTAL) DEPOSITION from CGCM T1R8 (combo), avgd over 2001-2005
C   interpolate to current latitude grid - written out in IDL on Linux, so use "little_endian"
C      depr8in(10,14,91) is in cm/sec ;    ddlat(91)
C
       OPEN (177, FILE = 'dep_2d_t1r8.xdr',
     >       convert="little_endian", form = 'unformatted', 
     >       ACCESS = 'STREAM', STATUS = 'OLD')
         READ (177) depr8in
         READ (177) ddlat
       CLOSE (177)



C  interpolate T1R8 surface deposition to model latitudes - DEPR8(10,14,L$) in COMMON in cm/sec
C                   -->   yy91(91), LAT4(L$) ==>  ZZOUTL(L$)

        do 4440 iim=1,14
        do 4440 iii=1,10

          do 7677 ij=1,91
 7677	       yy91(ij) = DEPR8IN(iii,iim,ij)

            CALL LINTERP(DDLAT, 91, YY91, LAT4, L$, 0, ZZOUTL)

            DO 1836 ij=1,L$
 1836	      DEPR8(iii,iim,ij) = ZZOUTL(ij)

 4440	 CONTINUE



C
C   read in wet scavenging from CGCM T1R8 (combo), avgd over 2001-2005
C   interpolate to current model grid - written out in IDL on Linux, so use "little_endian"
C      SCAVIN(20,14,91,16) is in 1/sec ;    ddlat(91), ddz(16)
C
       OPEN (177, FILE = 'scav_2d_t1r8.xdr',
     >       convert="little_endian", form = 'unformatted', 
     >       ACCESS = 'STREAM', STATUS = 'OLD')
         READ (177) scavin
         READ (177) ddlat
         READ (177) ddz
       CLOSE (177)


C  interpolate T1R8 wet scavenging to current model grid; 
C    yysc(91,16);  LAT4(L$), ZALT90(Z$) ==>  YYLZ(L$,Z$)
C    THERE are 0's so must do LINEAR interpolation in ALTITUDE
C    set to ZERO above 16 km;   SCAVR8(20,14,L$,Z$) in COMMON in 1/sec

        do 4444 iim=1,14
        do 4444 iii=1,20

          do 7577 ik=1,16
          do 7577 ij=1,91
 7577	       yysc(ij,ik) = SCAVIN(iii,iim,ij,ik)

          CALL BINTERP(ddlat, 91, ddz, 16, YYSC, 
     >                 LAT4, L$, ZALT90, Z$, 0, 0, YYLZ)

          do 7777 ik=1,Z$
          do 7777 ij=1,L$
               scavr8(iii,iim,ij,ik) = YYLZ(ij,ik)
               if (zalt90(ik) .ge. 16.) scavr8(iii,iim,ij,ik) = 0.
 7777     CONTINUE

 4444	 CONTINUE



C
C  read in vertical gradients of the DU/km ozone climatology, for Kzz adjustment
C     dukzr(14,45,60), latdu(45), zzdu(60)
C
C  this is ALREADY on the 45x89 Kzz grid used in XCOUPLED, (off by 0.5 km - close enough)
C     so don't need to interpolate here, but expand to 89 levels array -> dukz(14,45,89)
C     - written out in IDL on Linux, so use "little_endian"
C
       OPEN (177, FILE = 'ozclim_dukz.xdr',
     >       convert="little_endian", form = 'unformatted', 
     >       ACCESS = 'STREAM', STATUS = 'OLD')
         READ (177) dukzr
         READ (177) latdu
         READ (177) zzdu
       CLOSE (177)


       do 3300 ik=1,60
       do 3300 ij=1,45
       do 3300 iim=1,14
 3300       dukz(iim,ij,ik) = dukzr(iim,ij,ik)


C
C
cc ***************************************************************************************
c
C  also need to load in initial TPROB45(211,45,76,72) and LPROB45(181,45,76,72) w/ climatology here
C     this then gets interpolated to the proper day and grid in SETDAILY
C     (where ITP1(L$,Z$) and ITP2(L$,Z$) are also defined) which is then
C     used for the gas phase reaction rates (REACTION/REACTIONG), called right after SETDAILY
C     (this is all done BEFORE the "do 100" loop)
C

       CALL GETPROB


       RETURN
       END
