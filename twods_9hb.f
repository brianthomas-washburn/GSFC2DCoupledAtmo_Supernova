c      PROGRAM TWOD

      SUBROUTINE TWODS

!      USE degree_trig

      include "comcfg.h"
      include "timer.inc"


C     TWO-D DYNAMICAL MODEL BY M.R. SCHOEBERL
C                VERS 2.2

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
      INCLUDE 'COMMONR.INC'  
      INCLUDE 'COMMONW.INC'
      
      logical flag1
      
ccelf      include 'test_value.inc'
ccelf      logical nantest

      DIMENSION gflux(n$,m$)
      DIMENSION W1V(N$), W2V(N$), HEAT0(N$,M$), COOL0(N$,M$)
      DIMENSION UBST(N$,M$,2),THST(N$,M$,2)

C EF ubar output

      DIMENSION UB1(N$,M$), UB2(N$,M$), UB3(N$,M$), 
     >          ub1out(ns$,ms$), ub2out(ns$,ms$), ub3out(ns$,ms$)

C     DUM1X already declared in COMMOND.INC
C      DIMENSION DUM1X(NP$,MP$),ISWTMP(NSWTCHS),ISW2TMP(ISW$)
      DIMENSION ISWTMP(NSWTCHS),ISW2TMP(ISW$)
      
      COMMON/CORIOL/ CFA(N$,M$),CFB(N$,M$)

      COMMON /SLAKO/ AY(N$,M$),CY(N$,M$),B(N$,M$),AZ(N$,M$),CZ(N$,M$)
     1,F(N$,M$)

      LOGICAL REST
      
      logical co2var
      common/co2mr/co2vmr(45),co2lw(48),co2var

      COMMON/CBBC1/BBC1(N$,4)

      COMMON/CKZZ1/XKZZFX(NP$,MP$)

      COMMON/CKYY1/ BGKYY0(N$,M$), BGKYY(N$,M$), EPOTX(NP$,MP$),
     >              XKYYFX(N$,M$)

      COMMON/CLATH/GLHDAYC(N$,M$), WLHDAYC(N$,M$), HACKDAYC(N$,M$),
     >             G4LHDAYC(N$,M$), LTNHT(N$,M$),  TDIFFDAY(N$,M$),
     >             G4LHDAYCX(NP$,MP$)

                                                                        ! COMMON for CONSTITS -mixing ratio
      COMMON/CCHEM/CNC(S$C,NS$,MS$), HEATCHM(NS$,MS$), COOLCHM(NS$,MS$),
     >             HEATCHML(NS$,MS$), DRAGCHM(7,NS$,MS$)

C                                                  COMMON for coupled model geop hgt wave amps, from TEMPIN:
      COMMON/CGPBC/ GPBC(6,2,73,30,360), ypgp(73)
      COMMON/CGPBD/ GPBC1(6,2,N$)

      COMMON/CTNCEP/TNCEP(N$,M$), TNCEPSFC(NP$,2), TEMPG42D(NP$,3), 
     >              UT1KM(NP$,2), U3KM(NP$,2)

      COMMON/CTNCEPE/TNCEPE(NP$,MP$), XN2E(NP$,MP$), dtedz(NP$,MP$)


      COMMON/CUNCEP/UNCEP(NP$,MP$), UDIFFDAY(NP$,MP$)

      COMMON/CEPDD/ EPDD(N$,M$), EPDDX(NP$,MP$), EPBARO(NP$,MP$)


c   common for NCEP v'T', u'v' for baroclinic waves ONLY - for coupled model
C                   vt is in K/sec;  uv is in m/sec^2 (from TEMPIN)

ccvtbaro       COMMON/CVTBARO/ vt37(73,33,37), uv37(73,33,37),  
ccvtbaro     >                 BALAT1(73), BZZ1(33), time37(37)


       COMMON/CVTUV/ VTC(NP$,MP$), UVC(NP$,MP$)


c  common for coupled model temperature - NCEP offset (lat-hgt-season) from 1985-2005:  from TEMPIN

       COMMON/CTOFFS/ xlatoffs(45), zzoffs(76)
ccccccccccccccc       COMMON/CTOFFS/ TEMPOFFS(360,45,76), xlatoffs(45), zzoffs(76)


c            common for TD GEOS 4 latent heating,   G4LHDAY(91,117) is for current day in K/day

       COMMON/CG4LH/ xlhsens(15,91,117), G4LHDAY(91,117)


C   COMMON for WACCM oro Gwaves DELF (m/sec^2) from TEMPIN

       COMMON/CWOGW/wogw(46,57,14), platg(46), zzwg(57), tog14(14)

       COMMON/CWOGWD/ OGWD(NP$,MP$)


C   common for trop heating correction (diff w/ NCEP temps), store daily values in TTDIFFYR

       COMMON/CTDIFF/ TTDIFFYR(N$,M$,360)


C  common for year counter, steady state

       COMMON/CIYRSS/IYRCT, ISTST
ccccccc       COMMON/CDELF/ IYRCT0


C
C  COMMON for WACCM combined eddy heating rate (theta heating) in K/sec
C
cccwconv48       COMMON/CQWT/qwtt(46,57,14), platw(46), zz57(57), tq14(14)
cccwconv48
cccwconv48       COMMON/CQWTD/qwttd(NP$,MP$)




C  COMMON for WACCM convective mass flux for Kzz
C
ccc46       COMMON/CWCONKZ/cmfw(46,25,14), plat(46), zzw(25), tg14(14)
ccc46
ccc46       COMMON/CWCONKZD/cmfd(NP$,M$)


C
C  COMMON for WACCM convective u'w' DELF, convective (theta) heating corrections (from TEMPIN)
C    and interpolated to current day/dynamics grid
C    EPWD(NP$,MP$) is m/sec^2  ;  QWTD(NP$,MP$), QVTD(NP$,MP$) are K/sec (theta heating)

cccc      COMMON/CWCON/epw(46,57,14), qwt(46,57,14), qvt(46,57,14), 
cccc     >             plat(46), zzw(57), tg14(14)
cccc
cccc      COMMON/CWCOND/epwd(NP$,MP$), qwtd(NP$,MP$), qvtd(NP$,MP$)

C
C  COMMON for ERA-40, w'Th' (from TEMPIN) - in K/sec (theta heating)
C
ccert      COMMON/CERWT/erwtr(73,49,74), eelat(73), eez(49), e74(74)



C     NSTEP initialized to 0 is currently necessary for CALL to TRCINIT
C     from SUBROUTINE TWODS.

      DATA NSTEP/0/,IPRINT/0/

cjer added h2oin 4/14/97 to use model h2o in radiation
      dimension h2oin(ns$,ms$),h2otmp1(ns$,ms$)
cjer added co2in 1/6/99 to use model co2 in radiation
      dimension co2in(ns$,ms$),co2tmp1(ns$,ms$), o3in(ns$,ms$)

c     Array trr declaration added 10/7/94 by PEM to avoid compilation
c     error under IRIX 5.2.
      real kyyout,kyzout,kzzout,kzztemp
      dimension tout(ns$,ms$),wout(ns$,ms$),
     1          kyyout(ns$,ms$),kyzout(ns$,ms$),kzzout(ns$,ms$),
     2          yin(ns$),zin(ms$),kzztemp(n$,m$)
      real trr(ms$), LTNHT

      REAL BBCDAY(91,21), LATBBC(91), ZBBC(21), zp10(4), EHFDAYC(N$,M$)
      REAL KZZTH(45,89), xlatkzz(45), zzkzz(89), latlh(91), zzlh(117)
      REAL KYYDAY(46,88), xlatkyy(46), zzkyy(88)

      REAL TDIFFD(45,76),UDIFFD(45,76), XLATF(45), ZZFF(76), qwt2(46,57)
      REAL BARODAY(73,33), BALAT(73), BZZ(33), erwtd(73,49), cmf1(46,25)
      REAL epw1(46,57), qwt1(46,57), qvt1(46,57), t14(14), t74(74)

      REAL LHDAY(91,117), WLHDAY(91,117),HACKDAY(91,117),UBARDAY(91,117)
      REAL EHFDAY(91,117), TEMPDAY(91,117), TEMPDAYSF(91), xtemp(NP$)
      REAL HEAT9(N$,M$), COOL9(N$,M$), HEATG5(N$,M$), COOLG5(N$,M$)

      REAL GP73(73), GPN(N$), GPAP(6,2,N$), wogw1(46,57)
      REAL BDRAG0(NP$,MP$), BDRAG1(NS$,MS$), TEMPG4DAY(91,3)
      REAL zzz3(3), zzz1(1), UTOUT(NP$,1)
      REAL ttime(1), ttout(1), tt37(37), vtday(73,33), uvday(73,33)

      REAL t4576(45,76), TMODOFF(NP$,MP$), glon(72), gxx(72), rtime(N$)

C                                              arrays for Lin-Rood transport
      real*8 xcold(NP$,MP$), xcnew(NP$,MP$)
      real*8 xxs0(NP$,MP$), xxs1(NP$,MP$)

      REAL*8 TNCEPE8(NP$,MP$), dtedz8(NP$,MP$), dze8


C fixed model DELF for current day (m/sec/day), in PARAM.INC, NS$=L$=45; MS$=Z$=76)
C                                                and Z$X=Z$+12,  use YIN(NS$) for latitudes
      REAL EPDAY(NS$,MS$+12), ZALTC(MS$+12)


C  output back to XCOUPLED:

      REAL ypf(N$), yppf(NP$), zpf(M$), zppf(MP$)
      REAL tempb(N$,M$), kyyb(N$,M$), wsb(N$,M$), ttdiff(N$,M$)
      REAL ub3x(NP$,MP$), kzzb(NP$,MP$), tempbxo(NP$,MP$)


C     ms0, ns0 added 7/13/95 by PEM to avoid conflict between ms$, ns$
C     declared in PARAM.INC and as parameters to 'entry twodr', below.
C     IRIX 6.0 Fortran viewed this as an illegal double declration of the
C     identifiers.

      !INTEGER ms0, ns0, NS$, MS$, MSX
      INTEGER ms0, ns0, MSX
C
C   NOTE:  ns0, ns$,    ms0, ms$  should be equal


C
C      SET UP CONSTANTS AND INITIALIZE ARRAYS
C      NSTEP IS STEP COUNTER
C      IPRINT IS PRINT DUMP COUNTER

C      SET UP SWITCHES

      CALL SWIT(REST)


C     SET UP OTHER CONSTANTS

      NDBUG= ISW2(40)
      NSTEPS=ISW(1)
      IDUMP0=ISW(2)
      IHEATS=ISW(3)
      IGWS=ISW(4)
      IMODUL=ABS(ISW(5))


      IDELT=ISW(6)
      IKYY=ISW(8)
      IKYZ=ISW(9)
      ITHERMW=ISW(15)
      IZRES=ISW(17)


      DO NN=1,NSWTCHS
         ISWTMP(NN)=ISW(NN)
      END DO
      DO NN=1,ISW$
         ISW2TMP(NN)=ISW2(NN)
      END DO


      NDBUG=0
      IF (IDUMP .LT. 2) NDBUG=1


c     VERTICAL RESOLUTION
     
      ASOR=1.00     
      ZRES=0.01*IZRES

C     TIME AVERAGING WEIGHTS

      w1=0.01*isw(11)
      w2=1.-w1
 


c      PRINT 7788,W1
7788  FORMAT(' W1 = ',F10.2)

C                               TEMPBX(NP$,MP$) in COMMOND.INC initialized in SETCON
      CALL SETCON(IHEATS)

      CALL START(REST,IHEATS)


C     The assignments to NDAY, NYEAR, below are not really used, but are
C     necessary for calls to OUTA/OUTP in SUBROUTINE THE_TR and
C     TRCINIT.

      IF (.NOT.REST) then

        IFILE=0
        NDAY=0
        NYEAR=1
C        PRINT*, NSTEP, NDAY, NYEAR
C        PRINT*, DECL, DAYL, DT, MONTH
C        PRINT*, IDAYX, ISW(10)
C        CALL TIMEX (NSTEP,NDAY,NYEAR,DECL,DAYL,DT,MONTH,IDAYX,ISW(10))

C********************************************************
C--   SET UP THETA AS TRACER FOR PRATHER SCHEME ---------

        CALL THE_TR

        DO K=1,MP$
        DO J=1,NP$
           XC00(J,K,1)=XC(J,K,1)  ! INITIAL POT. TEMP. 
           XC00(J,K,2)=XC(J,K,2)  ! INITIAL ZONAL MOM. 
        END DO
        END DO
 
C******************************************************** 

      CALL TRCINIT(2)

      endif


      IF( IDUMP0 .GT. 0) then 
       IDUMP= IDUMP0 ! USE NEW IDUMP IF >0
      ENDIF


        TMAX=0.
        do 2223 k=1,m$
        do 2223 j=1,n$
        TMAX=AMAX1(TMAX,th(J,K,1))
 2223   UMAX=AMAX1(UMAX,ub(J,K,1))
c        write(*,*) ' bbb ', umax,tmax
C     RESTORE THE SWITCHES IF RESTARTED

      ISW(1)=NSTEPS
      ISW(2)=IDUMP
      ISW(3)=IHEATS
      ISW(4)=IGWS
      ISW(5)=IMODUL
      ISW(6)=IDELT
      ISW(8)=IKYY
      ISW(9)=IKYZ


      DO NN=1,NSWTCHS
         ISW(NN)=ISWTMP(NN)
      END DO
c      DO NN=1,ISW$
c         ISW2(NN)=ISW2TMP(NN)
c      END DO

      NRDLS=0


C  load in background Kyy as in fixed model: 
c    stored in COMMON/CKYY1/BGKYY0(N$,M$),BGKYY(N$,M$),EPOTX(NP$,MP$),XKYYFX(N$,M$)

        CALL KYYBACK(yp, zp)


C
C  define Fixed model STREAMF latitudes and altitudes for latent heating:   latlh(91), zzlh(117)
C                                    and zp10(4) for X* BBC;   zzz3(3) for sfc temps, UBAR
        do 701 ij=1,91
 701       latlh(ij) = (ij-1)*2.-90.

        do 702 ik=1,117
 702       zzlh(ik) = ik-1.

        do 703 ik=1,4
 703       zp10(ik) = zp(ik)

        do 704 ik=1,3
 704       zzz3(ik) = ik - 1.


C  load in latitudinal Rayleigh friction profile for use in GETMOM: 
C      rfric(NP$) is in COMMONC.INC, key on YPP(NP$), based roughly on OROGRAPHY
C
        do 707 ij=1,NP$
          rfric(ij) = 1.
          if (ypp(ij) .le. -65.) rfric(ij) = .4
          if (ypp(ij) .gt. -65. .and. ypp(ij) .le. -50.) rfric(ij) = .4
          if (ypp(ij) .gt. -50. .and. ypp(ij) .le. -30.) rfric(ij) = .4
          if (ypp(ij) .gt. -30. .and. ypp(ij) .le. 0.)   rfric(ij) = .4

ccrfric          if (ypp(ij) .gt. -30. .and. ypp(ij) .le. 0.)
ccrfric     >         rfric(ij) = .5 + .25*( (SIND(ypp(ij)*3.))**2 )

cccwconv48       if (ypp(ij) .gt. 0. .and. ypp(ij) .le. 30.)
cccwconv48     >         rfric(ij) = .4 + .35*( (SIND(ypp(ij)*3.))**2 )
cccwconv48
cccwconv48          if (ypp(ij) .ge. 30. .and. ypp(ij) .le. 80.) rfric(ij) = .75
cccwconv48          if (ypp(ij) .gt. 80.) rfric(ij) = .75


          if (ypp(ij) .gt. 0.) rfric(ij)= .4 + .3*( (SIND(ypp(ij)))**4)

cccwconv49          rfric(ij) = .4
 707    CONTINUE

ccrfric          print *, 'IN TWODS:  RFRIC = ', rfric
          print *, 'IN TWODS:  RFRIC = ', rfric


        return

c******************************************************************************
c
c entry point from chem code, inside main chem time loop
c
c******************************************************************************



      ENTRY TWODR(DTIN, NS0, MS0, YIN, ZIN, 
     >      LHDAY, WLHDAY, HACKDAY, TEMPDAY, UBARDAY, TEMPG4DAY, EHFDAY,
     >      BBCDAY, LATBBC, ZBBC, KZZTH, XLATKZZ, ZZKZZ,
     >                            KYYDAY, XLATKYY, ZZKYY, EPDAY, ZALTC,
     >      TDIFFD, UDIFFD, XLATF, ZZFF, BARODAY, BALAT, BZZ,
     >      ypf, zpf, yppf, zppf, tempbxo, ub3x, kyyb, kzzb, wsb,ttdiff)


C                    CNC(S$C,NS$,MS$) (mixing ratio, PARTS/PART) in COMMON, co2in, o3in, h2oin(ns$,ms$)
      do k=1,ms$
      do j=1,ns$

        co2in(j,k) = cnc(20,j,k)
        o3in (j,k) = cnc( 4,j,k) * 1.0E+06
        h2oin(j,k) = cnc(15,j,k) * 1.0e+06 

        if(co2in(j,k) .lt. 10.e-6) then
	  print *,'in twods: co2in lt 10.e-6 after entry at j,k = ',
     *	    j,k,co2in(j,k)
	  stop
	endif
	if(co2in(j,k) .gt. 10000.e-6) then
	  print *,'in twods: co2in gt 10000.e-6 after entry at j,k = ',
     *	    j,k,co2in(j,k)   
          stop
        endif      
      enddo
      enddo

      ITHERMW=ISW(15)
      IZRES=ISW(17)

      ndbug= isw2(40)
      nsteps=isw(1)
      idump0=isw(2)
      iheats=isw(3)
      igws=isw(4)
      imodul=abs(isw(5))
      idelt=isw(6)
      ikyy=isw(8)
      ikyz=isw(9)
      ithermw=isw(15)
      izres=isw(17)


ccelf                  also initialize ASOR here for proper call to GETMIX below
      ASOR=1.00     


C     Here we assign values from the timer.inc COMMON blocks to the
C     dynamics timekeeping variables.

      NDAY  = IDOY365
      NYEAR = IYOR
      MONTH = IMOY365
      IDAYX = IDOM365
      DECL  = SOLDEC

C      nsteps=dtin/dt
      NSTEPS = DYPRCH


C      DT0=DAYL/ISW(6)
C      FRACD=DT0/DAYL
C      DELT=2.*DT0
C      UMAX=0.    
C      DT=DT0
C      NCHKTS=0
C      ICHKAG=IDELT

      DT0    = DAYL/(DYSTPD)
      FRACD  = DT0/DAYL
      DELT   = 2.0*DT0
      UMAX   = 0.0
      DT     = DT0



C   interpolate eddy heat flux for current day: (K/sec) to coupled model grid
C                                                              EHFDAY(91,117) -> EHFDAYC(N$,M$)

        CALL BINTERP(latlh, 91, zzlh, 117, EHFDAY,
     >               yp, N$, zp, M$, 0, 0, EHFDAYC)



C   interpolate GCM latent heating to coupled model dynamics grid for current day:
C                       LHDAY(91,117) in K/day => GLHDAYC(N$,M$),  yp(N$), zp(M$)
C
        CALL BINTERP(latlh, 91, zzlh, 117, LHDAY, 
     >               yp, N$, zp, M$, 0, 0, GLHDAYC)



C  also interpolate WACCM latent heating to coupled model dynamics grid for current day:
C                             WLHDAY(91,117) in K/day => WLHDAYC(N$,M$),  yp(N$), zp(M$)
C
        CALL BINTERP(latlh, 91, zzlh, 117, WLHDAY, 
     >               yp, N$, zp, M$, 0, 0, WLHDAYC)


C                                                       HACKDAY(91,117) in K/day -> HACKDAYC(N$,M$)
        CALL BINTERP(latlh, 91, zzlh, 117, HACKDAY, 
     >               yp, N$, zp, M$, 0, 0, HACKDAYC)



c   TD GEOS 4 latent heating - G4LHDAY(91,117) is for current day (K/day)
C                           -> G4LHDAYC(N$,M$)
C
        CALL BINTERP(latlh, 91, zzlh, 117, G4LHDAY,
     >               yp, N$, zp, M$, 0, 0, G4LHDAYC)


c   TD GEOS 4 latent heating - G4LHDAY(91,117) is for current day (K/day)
C                           -> G4LHDAYCX(NP$,MP$) for extended grid
C
        CALL BINTERP(latlh, 91, zzlh, 117, G4LHDAY,
     >               ypp, NP$, zpp, MP$, 0, 0, G4LHDAYCX)


c   TD GEOS 4 latent heating - G4LHDAY(91,117) is for current day (K/day)
C        -> also interpolate to CHEMISTRY grid for output - HEATCHML(NS$,MS$) is in K/day
C
        CALL BINTERP(latlh, 91, zzlh, 117, G4LHDAY,
     >               YIN, NS$, ZIN, MS$, 0, 0, HEATCHML)



C   interpolate UBAR to extended coupled model grid - UBARDAY(91,117) in m/sec -> UNCEP(NP$,MP$)
                                
        CALL BINTERP(latlh, 91, zzlh, 117, UBARDAY,
     >               ypp, NP$, zpp, MP$, 0, 0, UNCEP)




C  TEMPDAY(91,117) -> TNCEP(N$,M$) - clim NCEP temps for current day,
C  also interpolate 0-1 km temps to model lats - TEMPDAYSF(91) -> TNCEPSFC(NP$,2),  xtemp(NP$)
C                                                                     
        CALL BINTERP(latlh, 91, zzlh, 117, TEMPDAY, 
     >               yp, N$, zp, M$, 0, 0, TNCEP)


        DO 488 ik=1,2 
            do 489 ij=1,91
 489           tempdaysf(ij) = tempday(ij,ik)

            CALL LINTERP(latlh, 91, TEMPDAYSF, ypp, NP$, 0, xtemp)

            do 491 ij=1,NP$
 491           TNCEPSFC(ij,ik) = xtemp(ij)
 488    CONTINUE



C    TEMPG4DAY(91,3) ->  GEOS 4 sfc, 1km, 2km temps for current day, interpolate to model latitudes
C                - TEMPDAYSF(91) -> TEMPG42D(NP$,3),  xtemp(NP$)
C                                                                     

        DO 1488 ik=1,3
            do 1489 ij=1,91
 1489           tempdaysf(ij) = TEMPG4DAY(ij,ik)

            CALL LINTERP(latlh, 91, TEMPDAYSF, ypp, NP$, 0, xtemp)

            do 1491 ij=1,NP$
 1491           TEMPG42D(ij,ik) = xtemp(ij)
 1488    CONTINUE


C  interpolate TEMPG42D(91,3) (sfc, 1km, 2km temps for current day) 
C     to extended latitudes at lowest grid point - zpp(1)
C     zzz1(1) -> UTOUT(NP$,1) ;  load into  UT1KM(NP$,2),  1=TEMP;  2=UBAR (m/sec)

        zzz1(1) = zpp(1)


        CALL BINTERP(latlh, 91, zzz3, 3, TEMPG4DAY,
     >               ypp, NP$, zzz1, 1, 0, 0, UTOUT)

         do 1588 ij=1,NP$
 1588       UT1KM(ij,1) = utout(ij,1)



C  also interpolate NCEP UBAR to lowest extended grid point, UBARDAY(91,117) in m/sec
                                
        CALL BINTERP(latlh, 91, zzlh, 117, UBARDAY,
     >               ypp, NP$, zzz1, 1, 0, 0, UTOUT)

         do 1591 ij=1,NP$
 1591       UT1KM(ij,2) = utout(ij,1)



C  also interpolate NCEP UBAR to 2-3rd extended grid point (3-5km), UBARDAY(91,117) in m/sec
C                             zzz1(1) -> UTOUT(NP$,1);   U3KM(NP$,2)
        zzz1(1) = zpp(2)
                                
        CALL BINTERP(latlh, 91, zzlh, 117, UBARDAY,
     >               ypp, NP$, zzz1, 1, 0, 0, UTOUT)

         do 1593 ij=1,NP$
 1593       U3KM(ij,1) = utout(ij,1)


        zzz1(1) = zpp(3)
                                
        CALL BINTERP(latlh, 91, zzlh, 117, UBARDAY,
     >               ypp, NP$, zzz1, 1, 0, 0, UTOUT)

         do 1594 ij=1,NP$
 1594       U3KM(ij,2) = utout(ij,1)




C  TEMPDAY(91,117) -> TNCEPE(NP$,MP$) - clim NCEP temps for current day on EXTENDED COUPLED MODEL GRID
C    convert to REAL*8 TNCEPE8(NP$,MP$) for DERV4, and compute BV freq - XN2E(NP$,MP$)
C                                                            dtedz8(NP$,MP$)

        CALL BINTERP(latlh, 91, zzlh, 117, TEMPDAY, 
     >               ypp, NP$, zpp, MP$, 0, 0, TNCEPE)


        xkx = 2./7.
        hhx = 7000.
        dze8 = zres*hhx

        DO 567 K=1,MP$
        DO 567 J=1,NP$
 567       TNCEPE8(J,K) = TNCEPE(J,K)*(EXP(zpp(K)/7.))**xkx    ! this is THETA


        CALL DERV4(1, tncepe8, dtedz8, dze8, NP$, MP$, 2)


        DO 767 K=1,MP$
        DO 767 J=1,NP$
 767       DTEDZ(J,K) = DTEDZ8(J,K)*1000.     ! this is dtheta/dz in REAL*4


c  now get XN2 (in 1/sec-day)

        DO 577 K=1,MP$
        DO 577 J=1,NP$
 577       TNCEPE8(J,K) = TNCEPE(J,K)

        CALL DERV4(1, tncepe8, dtedz8, dze8, NP$, MP$, 2)

        DO 777 K=1,MP$
        DO 777 J=1,NP$
 777     XN2E(J,K) = 287./hhx*(DTEDZ8(J,K) + xkx*TNCEPE(J,K)/hhx)*86400.



c
C   interpolate X* BBC from NCEP (fixed model) to coupled model dynamics grid for current day:
C      LATBBC(91), ZBBC(21), zp10(4);  BBCDAY(91,21) -> BBC1(N$,4) are in m2/sec  (in COMMON above)
C
C   USING X* BBC from NCEP in the COUPLED MODEL DOES NOT WORK!!!!!
C
        CALL BINTERP(latbbc, 91, zbbc, 21, BBCDAY, 
     >               yp, N$, zp10, 4, 0, 0, BBC1)



C  interpolate fixed model Kzzs to coupled model dynamics grid for current day:   do LOG INTERPOLATION in ALT
C      xlatkzz(45), zzkzz(89)   KZZTH(45,89) -> XKZZFX(NP$,MP$) in cm2/sec  (in COMMON above)
C
        CALL BINTERP(xlatkzz, 45, zzkzz, 89, KZZTH, 
     >               ypp, NP$, zpp, MP$, 0, 1, XKZZFX)



C  interpolate fixed model Kyys to coupled model dynamics grid for current day: LINEAR INTERPOLATION
C  xlatkyy(46), zzkyy(88)   KYYDAY(46,88) -> XKYYFX(N$,M$) in cm2/sec  (in COMMON above)
C
        CALL BINTERP(xlatkyy, 46, zzkyy, 88, KYYDAY, 
     >               yp, N$, zp, M$, 0, 0, XKYYFX)



C  interpolate fixed model DELF to coupled model dynamics grid for current day:
C      YIN(NS$), ZALTC(MS$+12) ; EPDAY(NS$,MS$+12) -> EPDD(N$,M$) in m/sec/day in COMMON ABOVE
C                                           also interpolate to extended grid - EPDDX(NP$,MP$)
        MSX = MS$+12
        CALL BINTERP(yin, NS$, zaltc, MSX, EPDAY, 
     >               yp, N$, zp, M$, 0, 0, EPDD)


        CALL BINTERP(yin, NS$, zaltc, MSX, EPDAY, 
     >               ypp, NP$, zpp, MP$, 0, 0, EPDDX)




C  interpolate NCEP DELF from BAROCLINIC EDDIES only to coupled model dynamics grid for current day:
C  BARODAY(73,33), BALAT(73), BZZ(33)    for momentum source ONLY (NO effect on KYY)
C                                   => EPBARO(NP$,MP$) is in m/sec^2
C
        CALL BINTERP(balat, 73, bzz, 33, BARODAY,
     >               ypp, NP$, zpp, MP$, 0, 0, EPBARO)



        ttime(1) = IDOY360*1.


C
C  INTERPOLATE WACCM orographic gravity wave DELF,to current day
C   and then to dynamics grid
C
C   wogw(46,57,14), platg(46), zzwg(57), tog14(14), t14(14)
C      --> wogw1(46,57),  ttout(1) -> OGWD(NP$,MP$) in COMMON above 

         do 7710 ik=1,57
         do 7710 ij=1,46 

            do 7711 iiw=1,14
 7711          t14(iiw) = wogw(ij,ik,iiw)

            CALL LINTERP(tog14, 14, t14, ttime, 1, 0, ttout)

            wogw1(ij,ik) = ttout(1)

 7710    CONTINUE


         CALL BINTERP(platg, 46, zzwg, 57, WOGW1,
     >                ypp, NP$, zpp, MP$, 0, 0, OGWD)



C
C  INTERPOLATE WACCM total eddy heating rate (theta heating, K/sec) to current day
C     and then to dynamics grid  , ttout(1)
C
C   qwtt(46,57,14), platw(46), zz57(57), tq14(14), t14(14) --> qwt2(46,57)
C                                QWTTD(NP$,MP$) is in COMMON above 
cccwconv48         do 5077 ik=1,57
cccwconv48         do 5077 ij=1,46 
cccwconv48
cccwconv48            do 5072 iiw=1,14
cccwconv48 5072          t14(iiw) = qwtt(ij,ik,iiw)
cccwconv48
cccwconv48            CALL LINTERP(tq14, 14, t14, ttime, 1, 0, ttout)
cccwconv48
cccwconv48            qwt2(ij,ik) = ttout(1)
cccwconv48 5077    CONTINUE
cccwconv48
cccwconv48
cccwconv48         CALL BINTERP(platw, 46, zz57, 57, QWT2,
cccwconv48     >                ypp, NP$, zpp, MP$, 0, 0, QWTTD)
cccwconv48
cccwconv48
cccwconv48C                             ramp down QWTTD to zero above 25 km, max of 20%
cccwconv48        do 5075 ik=1,MP$
cccwconv48            zwt = (23. - zpp(ik))/25.
cccwconv48            if (zwt .le. 0.) zwt = 0.
cccwconv48            if (zwt .ge. .2) zwt = .2
cccwconv48      
cccwconv48            do 5076 ij=1,NP$
cccwconv48 5076            QWTTD(ij,ik) = QWTTD(ij,ik)*zwt
cccwconv48 5075   CONTINUE





C  interpolate NCEP v'T'/dy, u'v'/dy from BAROCLINIC EDDIES only to coupled model dynamics grid for current day:
C  vt37(73,33,37), uv37(73,33,37), BALAT1(73), BZZ1(33), time37(37) is 1-365;  vt is in K/sec;  uv is in m/sec^2
C                                              tt37(37), vtday(73,33), uvday(73,33) -> VTC(NP$,MP$), UVC(NP$,MP$)
ccvtbaro
ccvtbaro         do 1010 ik=1,33
ccvtbaro         do 1010 ij=1,73
ccvtbaro
ccvtbaro            do 1011 iiw=1,37
ccvtbaro 1011          tt37(iiw) = vt37(ij,ik,iiw)
ccvtbaro
ccvtbaro            CALL LINTERP(time37, 37, tt37, ttime, 1, 0, ttout)
ccvtbaro
ccvtbaro            vtday(ij,ik) = ttout(1)
ccvtbaro
ccvtbaro
ccvtbaro            do 1012 iiw=1,37
ccvtbaro 1012          tt37(iiw) = uv37(ij,ik,iiw)
ccvtbaro
ccvtbaro            CALL LINTERP(time37, 37, tt37, ttime, 1, 0, ttout)
ccvtbaro
ccvtbaro            uvday(ij,ik) = ttout(1)
ccvtbaro 1010    CONTINUE
ccvtbaro
ccvtbaro
ccvtbaroC                                                            
ccvtbaro        CALL BINTERP(balat1, 73, bzz1, 33, VTDAY,
ccvtbaro     >               ypp, NP$, zpp, MP$, 0, 0, VTC)
ccvtbaro
ccvtbaro                                                            
ccvtbaro        CALL BINTERP(balat1, 73, bzz1, 33, UVDAY,
ccvtbaro     >               ypp, NP$, zpp, MP$, 0, 0, UVC)
ccvtbaro

ccvtbaro for output:
ccvtbaro        CALL BINTERP(balat1, 73, bzz1, 33, VTDAY,
ccvtbaro     >               YIN, NS$, ZIN, MS$, 0, 0, HEATCHML)


C
C  INTERPOLATE WACCM mass flux Kzz to current day, and then to dynamics grid
C       also need to set minimum for cmf1 (so it's non-zero)
C   cmfw(46,25,14) -> cmf1(46,25) -> cmfd(NP$,M$),  plat(46), zzw(25), tg14(14)
C                                                   tg14(14), t14(14)
ccc46         do 7010 ik=1,25
ccc46         do 7010 ij=1,46 
ccc46
ccc46            do 7011 iiw=1,14
ccc467011          t14(iiw) = cmfw(ij,ik,iiw)
ccc46
ccc46            CALL LINTERP(tg14, 14, t14, ttime, 1, 0, ttout)
ccc46
ccc46            cmf1(ij,ik) = ttout(1)
ccc46            if (cmf1(ij,ik) .le. 1.e-6) cmf1(ij,ik) = 1.e-6
ccc467010    CONTINUE
ccc46
ccc46
ccc46        CALL BINTERP(plat, 46, zzw, 25, CMF1,
ccc46     >               ypp, NP$, zp, M$, 0, 1, CMFD)



C
C  INTERPOLATE WACCM orographic gravity wave DELF, convective heating rate to current day
C   and then to dynamics grid
C
C   epw(46,57,14), qwt(46,57,14), qvt(46,57,14), plat(46), zzw(57), tg14(14), t14(14)
C                      --> epw1(46,57), qwt1(46,57), qvt1(46,57)  , ttout(1)
cccc         do 3010 ik=1,57
cccc         do 3010 ij=1,46 
cccc
cccc            do 3011 iiw=1,14
cccc3011          t14(iiw) = epw(ij,ik,iiw)
cccc
cccc            CALL LINTERP(tg14, 14, t14, ttime, 1, 0, ttout)
cccc
cccc            epw1(ij,ik) = ttout(1)
cccc
cccc
cccc            do 3012 iiw=1,14
cccc 3012          t14(iiw) = qwt(ij,ik,iiw)
cccc
cccc            CALL LINTERP(tg14, 14, t14, ttime, 1, 0, ttout)
cccc
cccc            qwt1(ij,ik) = ttout(1)
cccc
cccc
cccc            do 3017 iiw=1,14
cccc 3017          t14(iiw) = qvt(ij,ik,iiw)
cccc
cccc            CALL LINTERP(tg14, 14, t14, ttime, 1, 0, ttout)
cccc
cccc            qvt1(ij,ik) = ttout(1)
cccc 3010    CONTINUE
cccc
cccc
C   EPWD(NP$,MP$), QWTD(NP$,MP$), QVTD(NP$,MP$) in COMMON above 
cccc
cccc        CALL BINTERP(plat, 46, zzw, 57, EPW1,
cccc     >               ypp, NP$, zpp, MP$, 0, 0, EPWD)
cccc
cccc
cccc        CALL BINTERP(plat, 46, zzw, 57, QWT1,
cccc     >               ypp, NP$, zpp, MP$, 0, 0, QWTD)
cccc
cccc        CALL BINTERP(plat, 46, zzw, 57, QVT1,
cccc     >               ypp, NP$, zpp, MP$, 0, 0, QVTD)
cccc
cccc
C                                 - reduce by 1/2, and ramp down QWTD to zero above 25 km
cccc        do 3030 ik=1,MP$
cccc            zwt = (25. - zpp(ik))/10.
cccc            if (zwt .le. 0.) zwt = 0.
cccc            if (zwt .ge. .5) zwt = .5
cccc      
Ccc16            do 3031 ij=1,NP$
Ccc16 3031            QWTD(ij,ik) = QWTD(ij,ik)*zwt
cccc 3030   CONTINUE
cccc
C     - reduce by 3/4, and ramp down EPWD to zero above 15 km
C       - also limit negatives to -.5 m/sec/day,    EPWD is in m/sec^2 - ZERO-out
cccc
cccc        ewmin = -.5/86400.
cccc
cccc        do 6030 ik=1,MP$
cccc            zwt = (20. - zpp(ik))/10.
cccc            if (zwt .le. 0.) zwt = 0.
cccc            if (zwt .ge. .25) zwt = .25
cccc      
cccc            do 6031 ij=1,NP$
cccc                 EPWD(ij,ik) = 0.  ! EPWD(ij,ik)*zwt
cccc                 if (EPWD(ij,ik) .le. ewmin) EPWD(ij,ik) = ewmin
cccc 6031       CONTINUE
cccc 6030   CONTINUE


C
C  INTERPOLATE ERA40 w'th' heating rate to current day, then to dynamics grid
C
C   erwtr(73,49,74), eelat(73), eez(49), e74(74) ->  t74(74), erwtd(73,49)
C
ccert         do 4010 ik=1,49
ccert         do 4010 ij=1,73 
ccert
ccert            do 4011 iiw=1,74
ccert4011          t74(iiw) = erwtr(ij,ik,iiw)
ccert
ccert            CALL LINTERP(e74, 74, t74, ttime, 1, 0, ttout)
ccert
ccert            erwtd(ij,ik) = ttout(1)
ccert 4010    CONTINUE
ccert
C                                                   --> ERWT(NP$,MP$) 
ccert        CALL BINTERP(eelat, 73, eez, 49, ERWTD,
ccert     >               ypp, NP$, zpp, MP$, 0, 0, ERWT)




C  interpolate coupled model TBAR/UBAR - NCEP differences to coupled model dynamics grid for current day:
C    TDIFFD(45,76), UDIFFD(45,76), XLATF(45), ZZFF(76)
C
C    temperature diffs goes into Latent heating (in K) => TDIFFDAY(N$,M$),  yp(N$), zp(M$)
C

        CALL BINTERP(xlatf, 45, zzff, 76, TDIFFD, 
     >               yp, N$, zp, M$, 0, 0, TDIFFDAY)

C                                                             => UDIFFDAY(NP$,MP$)
        CALL BINTERP(xlatf, 45, zzff, 76, UDIFFD, 
     >               ypp, NP$, zpp, MP$, 0, 0, UDIFFDAY)




C  interpolate NCEP REANALYSIS CLIMATOLOGICAL GEOP hgts to coupled model latitude grid for current day
C  need to interpolate AMPS and PHASES, then convert back to Aks and Bks, INTERPOLATING Aks,Bks is NOT RIGHT!
C      ypgp(73), GPBC(6,2,73,30,360) -> GPAP(6,2,N$) (in COMMON above),   IDOY360, IYEARC are in timer.inc
C
C   1st index: 1=zonal avg, 2-6 are waves 1-5;   2nd index:  1=AMP (METERS),  2=PHASE (degrees longitude)
C   73 latitudes (90S-90N, by 2.5 deg - NCEP grid);   4th index:  altitude, 1 - 30 km, by 1 km
C
C
Cccbc   NOTE: tried using individual years, eg, 1958-2006, got some interannual variability in model resuls
Cccbc                                                        but not enough to warrant using individual years
ccbc        iybc = 50
ccbc        if (iyearc .ge. 23  .and.  iyearc .le. 71) iybc = iyearc - 22
ccbc
ccbcC                           ! set BC to clim. avg       ! GP73(73), GPN(N$)
ccbc        iybc = 50


        do 401 iih=1,2
        do 401 iiw=1,6
C                                                          just use 2 km value from NCEP geop hgts
          do 402 iij=1,73
 402	       gp73(iij) = GPBC(iiw,iih,iij,2,IDOY360)


Cwconv33 - fix NH Pwave phases to that in the SH with 180 day phase shift
CCWCONV156 - KEEP THIS - THIS IS NEEDED in NH!!

          if (iih .eq. 2  .and. iiw .ge. 2) then 
                if (IDOY360 .ge. 181.) then
                   do 405 ij=38,73
 405                  gp73(ij) = GPBC(iiw,iih,74-ij,2,IDOY360-180)
                endif

                if (IDOY360 .le. 180.) then
                   do 406 ij=38,73
 406                  gp73(ij) = GPBC(iiw,iih,74-ij,2,IDOY360+180)
                endif
          endif


           CALL LINTERP(ypgp, 73, GP73, yp, N$, 0, GPN)

          do 403 iij=1,N$
 403	       GPAP(iiw,iih,iij) = GPN(iij)

 401    CONTINUE


C  WCONV175 - DECREASE WAVES 1-4 AMPs by 10% at 30N-90N ONLY, ramp up a bit

          do 410 iij=1,N$
             amfac = 1.
             if (yp(iij) .ge. 30.) amfac = .95
             if (yp(iij) .ge. 35.) amfac = .9

             GPAP(2,1,iij) = GPAP(2,1,iij)*amfac
             GPAP(3,1,iij) = GPAP(3,1,iij)*amfac
             GPAP(4,1,iij) = GPAP(4,1,iij)*amfac
             GPAP(5,1,iij) = GPAP(5,1,iij)*amfac
 410	  CONTINUE



ccwonc25c    NOW CONVERT to 1=Ak and 2=Bk  (2nd index)  from AMP and PHASE:  GPAP(6,2,N$) -> GPBC1(6,2,N$)
ccwonc25                                    - THIS IS NOT CORRECT.....
ccwonc25        do 407 iij=1,N$
ccwonc25        do 417 iiw=2,6
ccwonc25            harr = iiw-1.
ccwonc25            ckk  = GPAP(iiw,1,iij)
ccwonc25            ckk2 = ckk**2
ccwonc25            tph2 = ( TAND( GPAP(iiw,2,iij)*harr) )**2
ccwonc25
ccwonc25            akk = SQRT( ckk2/(1. + tph2) )
ccwonc25            bkk = SQRT( ABS(ckk2 - akk**2) )                    ! do ABS value here just to ensure non-NANs
ccwonc25
ccwonc25           if (akk .ge. ckk) then 
ccwonc25               akk = ckk
ccwonc25               bkk = 0.0
ccwonc25            endif
ccwonc25            
ccwonc25            GPBC1(iiw,1,iij) = akk
ccwonc25            GPBC1(iiw,2,iij) = bkk
ccwonc25 417     CONTINUE
C ccwonc25                                                    ; also load in zonal avg
ccwonc25            GPBC1(1,1,iij) = GPAP(1,1,iij)
ccwonc25 407     CONTINUE


c    NOW CONVERT to 1=Ak and 2=Bk  (2nd index)  from AMP and PHASE:  GPAP(6,2,N$) -> GPBC1(6,2,N$)
C     glon(72), gxx(72),  do by expanding in longitude, then getting Ak's and Bk's
C                                    - THIS IS THE CORRECT WAY of DOING IT!!!!!
        do 801 ix=1,72
 801       glon(ix) = (ix-1)*5.


        do 807 iij=1,N$
        do 817 iiw=2,6
            harr = iiw-1.

            do 804 ix=1,72
               phlon = glon(ix) - gpap(iiw,2,iij)
 804           gxx(ix) = gpap(iiw,1,iij)*COSD(360./360.*harr*phlon)

            akk = 0.
            bkk = 0.

            do 710 ix=1,72
 710           akk = akk + gxx(ix)*SIND(harr*glon(ix))/36.

            do 711 ix=1,72
 711           bkk = bkk + gxx(ix)*COSD(harr*glon(ix))/36.
            
            GPBC1(iiw,1,iij) = akk
            GPBC1(iiw,2,iij) = bkk
 817     CONTINUE
C                                                     ; also load in zonal avg
            GPBC1(1,1,iij) = GPAP(1,1,iij)
 807     CONTINUE



C   
C  and interpolate NCEP - model temperture offset to NP$,MP$ grid, for current day - TMODOFF(NP$,MP$)
C     TEMPOFFS(360,45,76), xlatoffs(45), zzoffs(76) -> t4576(45,76)
C
ccccccc         do 617 ik=1,76
ccccccc         do 617 ij=1,45
ccccccc 617        t4576(ij,ik) = TEMPOFFS(IDOY360,ij,ik)
ccccccc
ccccccc         CALL BINTERP(xlatoffs, 45, zzoffs, 76, T4576, 
ccccccc     >                ypp, NP$, zpp, MP$, 0, 0, TMODOFF)
ccccccc


C  initialize TTDIFFYR(N$,M$,360) array (in COMMON)
C
       IF (IYRCT .eq. 1  .and.  IDOY360 .le. 2) then
           do 770 iit=1,360
           do 770 ik=1,M$
           do 770 ij=1,N$
 770          TTDIFFYR(ij,ik,iit) = 0.
       ENDIF



      NDBUG=0
      IF (IDUMP .LT. 2) NDBUG=1

 
c      PRINT 515,NSTEP,NDAY,NYEAR,DECL,DT
515   FORMAT(' START TIME STEP ',I5,/,' DAY,YEAR = ',2I5,/,
     1'  DECLINATION ANGLE = ',F10.2,/,' DT =',f10.0,' sec')


C ------------------  MAIN TIME LOOP -----------------------------

      DO 1000 LL=1, NSTEPS


C     Correct NSTEP calculation for new timing (timer.f).
C      NSTEP=NSTEP+1
      NSTEP = CHSTPD*(CHSTEP - 1) + LL


        CALL TABINIT  ! TIME DEP. TRACER INITIALIZATION
     

        TMAX=0.
        THMAX=0.
        do 2233 k=1,m$
        do 2233 j=1,n$
        TMAX=AMAX1(TMAX,T(J,K))
 2233   THMAX=AMAX1(THMAX,th(J,K,1))



C     FIRST GET GEOPOTENTIAL HEIGHT AND TEMPERATURE

c getalt in twodlib.f
       CALL GETALT2


ccpw32  also get temperatures on extended model grid for RADIATION - TEMPBX(NP$,MP$) in COMMOND.INC
ccpw32      compute as in GETALT2 -  THGX(MP$), TTOTHX(MP$) are in COMMONC.INC
ccpw32      THX(NP$,MP$) is in COMMONT.INC, updated in UTFRMX on previous step
ccpw32
ccpw32   THIS IS DONE on previous step in 517 loop
ccpw32
ccpw32         DO 527 K=1,MP$
ccpw32         DO 527 J=1,NP$
ccpw32 527        TEMPBX(J,K) = (THX(J,K) + THGX(K))/TTOTHX(K)
ccpw32


        TMAX=0.
        THMAX=0.
        do 2231 k=1,m$
        do 2231 j=1,n$
        TMAX=AMAX1(TMAX,T(J,K))
 2231   THMAX=AMAX1(THMAX,th(J,K,1))


C     OUTPUT THE TEMPERATURE

      CALL OUTA(9,NSTEP,NDAY,NYEAR,YP,ZP,T,N$,M$,0)

C     GET THE OZONE DISTRIBUTION

c Switch from chem o3 grid 'o3in' to dynm o3 grid 'o3m'
        call regrid(o3m,yp,zp,n$,m$,o3in,yin,zin,ns0,ms0,0)

c Limit o3 in polar night to 1.0 ppmv as input to radiation code only
c Make sure o3 is always >0
        do 977 k=1,m$
        do 977 j=1,n$
                if (zp(k).gt.60.0) o3m(j,k)=amin1(o3m(j,k),1.0)
977     o3m(j,k)=amax1(o3m(j,k),.01)

cjer to use model h2o in radiation, switch from chem h2o grid 'h2oin' to
c  dynm h2o grid 'h2om'

      do k=1,ms0
      do j=1,ns0
        h2otmp1(j,k)=alog(h2oin(j,k))
      end do
      end do

        call regrid(h2om,yp,zp,n$,m$,h2otmp1,yin,zin,ns0,ms0,0)

      do k=1,m$
      do j=1,n$
        h2om(j,k)=exp(h2om(j,k))
        if(h2om(j,k) .lt. 1.) h2om(j,k)=1.0
      end do
      end do

cjer to use model co2 in radiation, switch from chem co2 grid 'co2in' to
c  dynm co2 grid 'co2m'  (stored in common/o3/, which is in COMMOND.INC)

      
      if(co2var) call regrid(co2m,yp,zp,n$,m$,co2in,yin,zin,ns0,ms0,0)

      do j=1,ns$
      do k=1,ms$
        if(co2in(j,k) .lt. 10.e-6) then
	  print *,'in twods: co2in lt 10.e-6 after regrid at j,k = ',
     *	    j,k,co2in(j,k)
	  stop
	endif
	if(co2in(j,k) .gt. 10000.e-6) then
	  print *,'in twods: co2in gt 10000.e-6 after regrid at j,k = ',
     *	    j,k,co2in(j,k)   
          stop
        endif      
      enddo
      enddo
      
      CALL GETO3


C     Modified 9/17/96 by PEM.  Radiation calculation is performed at least
C     once per day.

      IDELH=ISW(22)
      IF (IDELH .GT. IDELT) IDELH = IDELT

C     End Modified 9/17/96.


C     Since NSTEP now counts the dynamics timesteps of each day, and IDELH
C     is at most the same as the number of dynamics timesteps per day, then
C     this should still work correctly:
C     ie, NSTEP = 1->12, IDELH=12, so this gets called 1x per day on 1st step (LL=1) - EF, 5/09

      IF ( MOD(NSTEP,IDELH) .EQ. 1 .OR. IDELH .EQ. 1 ) THEN 
              
C--      NEWTONIAN COOLING:
         IF(ABS(IHEATS).EQ.1) THEN
            CALL GETHET
            CALL GETCOL
            ENDIF


C--      ROSENFIELD'S HEATING & COOLING: - THIS IS WHAT IS USED for NEWRAD7 - NOV. 2007 (EF)
C                                          ie, ISW(3) = 0 ;   and ISW(28)=1
C
C   then in TROPHT, WACCM latent heating is used
C
         IF(ABS(IHEATS).EQ.0) THEN
            CALL GETRAD (YIN, ZIN)
c            CALL GETRAD(flag1)


C  interpolate HEAT, COOL(N$,M$) arrays to chemistry grid for output, first convert back to real temp heating
c                           HEATCHM(NS$,MS$), COOLCHM(NS$,MS$) in COMMON above, goes to XCOUPLED for output
C     interpolate latent heating LTNHT(N$,M$) -> HEATCHML(NS$,MS$) (K/sec) ALREADY DONE ABOVE
C
C  Also, calculate cooling in lowest model layer in terms of relaxation to NMC - TNCEP(N$,M$)
C  temperature field (or enhanced temperature field). T(N$,M$) is the model calculated temperature field.
C  Relaxation timescale is 1 day.   (THIS WAS TAKEN FROM TROPHT in PRMLIB.f)
C

            DO 222 K=1,M$
            DO 222 J=1,N$
               HEAT(J,K) = HEAT(J,K)/TTOTH(K)
 222           COOL(J,K) = COOL(J,K)/TTOTH(K)


            SRFNMC = 1.0000
                                                             ! RELAX SRFC. TEMP. to NCEP observed t (theta)
            do 3737 J=1,N$
                COOL(J,1) = COOL(J,1) +
     >              SRFNMC * (T(J,1) - TNCEP(J,1))/(1.*DAYL)
 3737       CONTINUE

            IF (SRFNMC.LE.0.0) THEN 
              WRITE(6,'(5X,"**** NO SURFACE COOL TO OBS. !!!!")' )
            ENDIF



            CALL BINTERP(yp, N$, zp, M$, heat, 
     >                   YIN, NS$, ZIN, MS$, 0, 0, HEATCHM)

            CALL BINTERP(yp, N$, zp, M$, cool, 
     >                   YIN, NS$, ZIN, MS$, 0, 0, COOLCHM)

            DO 227 K=1,M$
            DO 227 J=1,N$
               HEAT(J,K) = HEAT(J,K)*TTOTH(K)
 227           COOL(J,K) = COOL(J,K)*TTOTH(K)


            IF(ISW(28).EQ.1) CALL TROPHT

ccpw32            CALL BINTERP(yp, N$, zp, M$, LTNHT, 
ccpw32     >                   YIN, NS$, ZIN, MS$, 0, 0, HEATCHML)

         ENDIF



C--     NEWRAD9 (adopted from OFFLINE CODE) - ROSENFIELD'S HEATING & COOLING: - NOV. 2007 (EF)
C                                          ie, ISW(3) = 9 ;   and ISW(28)=1
C
C  pass the temperature field T(N$,M$) on the dynamics grid,
C    also pass the solar declination (in DEGREES), IDOY360, IYEARC (15=1950), updated in TIMETICK
C                                                      also pass dayl, TTOTH(M$), VNC(M$), THG(M$)
C
         IF(ABS(IHEATS).EQ.9) THEN
ccrad9            CALL RADIATE9 (yp, zp, T, YIN, ZIN, SOLDEC, IDOY360, IYEARC, 
ccrad9     >                     dayl, TTOTH, VNC, THG, HEAT9, COOL9)
ccrad9
ccrad9
C  Calculate cooling in lowest model layer in terms of relaxation to NMC - TNCEP(N$,M$)
C  temperature field (or enhanced temperature field). T(N$,M$) is the model calculated temperature field.
C  Relaxation timescale is 1 day.   (THIS WAS TAKEN FROM TROPHT in PRMLIB.f),   heat9, cool9(N$,M$)

            SRFNMC = 1.0000
                                                             ! RELAX SRFC. TEMP. to NCEP observed t (theta)
            do 2727 J=1,N$
                COOL9(J,1) = COOL9(J,1) +
     >              SRFNMC * (T(J,1) - TNCEP(J,1))/(1.*DAYL)
 2727       CONTINUE

ccrad9            IF (SRFNMC.LE.0.0) THEN 
ccrad9              WRITE(6,'(5X,"**** NO SURFACE COOL TO OBS. !!!!")' )
ccrad9            ENDIF

c                                                                  get HEATCHM, COOLCHM for output
            CALL BINTERP(yp, N$, zp, M$, heat9, 
     >                   YIN, NS$, ZIN, MS$, 0, 0, HEATCHM)

            CALL BINTERP(yp, N$, zp, M$, cool9, 
     >                   YIN, NS$, ZIN, MS$, 0, 0, COOLCHM)

C                    heat9, cool9(N$,M$) are returned from RADIATE9, reload into HEAT, COOL(N$,M$) arrays
C                                                         ! multiply by TTOTH as in NEWRAD7
            DO 2772 K=1,M$
            DO 2772 J=1,N$
               HEAT(J,K) = HEAT9(J,K)*TTOTH(K)
               COOL(J,K) = COOL9(J,K)*TTOTH(K)
 2772       CONTINUE
C                                                WACCM latent heating is now used in TROPHT
            IF(ISW(28).EQ.1) CALL TROPHT

ccpw32                                                   INTERPOLATION TO HEATCHML ALREADY DONE ABOVE
ccpw32            CALL BINTERP(yp, N$, zp, M$, LTNHT, 
ccpw32     >                   YIN, NS$, ZIN, MS$, 0, 0, HEATCHML)
         ENDIF



C
C   GEOS5  HEATING & COOLING: - JAN. 2008 (EF),  ie, ISW(3) = 7 ;   and ISW(28)=1
C
C  pass the temperature field T(N$,M$) on the dynamics grid,
C    NOW use TEMPBX(NP$,MP$) (COMMMOD.INC) temps on extended grid, 
C    computed ABOVE from previous time step/day - May 2009
C
C    also pass the solar declination (in DEGREES), IDOY360, IYEARC (15=1950), updated in TIMETICK
C                                                  also pass dayl, TTOTHX(MP$), VNC(M$), THGX(MP$) (in COMMON)
C
        IF(ABS(IHEATS).EQ.7) THEN
         CALL RADIATE_G5 (yp, zp, ypp, zpp, T, TEMPBX, YIN, ZIN, SOLDEC,
     >         IDOY360, IYEARC, dayl, TTOTHX, VNC, THGX, HEATG5, COOLG5)

C
C  HEATG5, COOLG5(N$,M$) (K/sec) are returned from RADIATE_G5, reload into HEAT, COOL(N$,M$) arrays 
C      multiply TTOTH(M$) factor as in NEWRAD7/NEWRAD9 to convert to THETA heating
C
C      if needed, add in EHFDAYC(N$,M$) (K/sec) ;  also CALL TROPHT to load any unnessessary stuff
C

            DO 2777 K=1,M$
            DO 2777 J=1,N$
               HEAT(J,K) = HEATG5(J,K)*TTOTH(K)    ! + EHFDAYC(J,K)*TTOTH(K) - better NOT to use
               COOL(J,K) = COOLG5(J,K)*TTOTH(K)
 2777       CONTINUE
C                                             GEOS4 latent heating now used in RADIATE 
C                                             interpolate to HEATCHML(NS$,MS$) is done ABOVE
            IF(ISW(28).EQ.1) CALL TROPHT

ccpw32          CALL BINTERP(yp, N$, zp, M$, LTNHT, 
ccpw32     >                   YIN, NS$, ZIN, MS$, 0, 0, HEATCHML)


C  load HEATG5, COOLG5(NP$,MP$) into HEATG5X(NP$,MP$), COOLG5X(NP$,MP$) (COMMOND.INC)
C  these are in K/sec in THETA heating, for use in GETSOR2 to compute THETA (instead of N$,M$ grid)
C   NOT USED - THIS DOESN'T WORK WELL - NEED SMOOOOOOTHED HEATING....

ccpw33         DO 5677 K=1,MP$
ccpw33         DO 5677 J=1,NP$
ccpw33            HEATG5X(J,K) = HEATG5(J,K)
ccpw33 5677       COOLG5X(J,K) = COOLG5(J,K)


C WCONV232 - add in extra tropospheric heat source as the NEGATIVE of the difference from NCEP temps
C   do for all tropospheric M$ levels (2.03 km and up, wherever latent heating is > 0.)
C   use updated model temps -  TEMPB(N$,M$) and TNCEP(N$,M$)
C   apply where latent heating is gt 0., G4LHDAYC(N$,M$) is in K/day
C   ttdiff(N$,M$) is in K/sec, add into HEAT(N$,M$) array as theta heating
C   load into daily array TTDIFFYR(N$,M$,360) for first 3 years ONLY if TD run, ALWAYS for SSTATE
C
C  WCONV237 - 30 day timescale for 90S-Eq, ramp up to 5-day time scale in polar NH - rtime(N$)

         iiadj = 0
       IF (ISTST .eq. 0 .and. IYRCT .le. 3  .or. ISTST .eq. 1) iiadj = 1


         do 5767 J=1,N$
            rtime(J) = 30.
            if (yp(J) .gt. 0.) rtime(J) = 5.*(1.+ 5.*((COSD(yp(j)))**6))
 5767    CONTINUE


         DO 5677 K=1,M$
         DO 5677 J=1,N$
            ttadj = -(TEMPB(J,K) - TNCEP(J,K))/(rtime(J)*86400.)

            if (G4LHDAYC(J,K) .gt. 0.  .and.  iiadj .eq. 1)
     >          TTDIFFYR(J,K,idoy360) = ttadj

            ttdiff(J,K) = TTDIFFYR(J,K,idoy360)
            HEAT(J,K) = HEAT(J,K) + ttdiff(J,K)*TTOTH(K)
 5677    CONTINUE

          
       ENDIF




C--     ZHU'S NLTE COOLING:
         IF(ABS(IHEATS).EQ.2) THEN
            CALL RADNLTE(LL)
            IF(ISW(28).EQ.1) CALL TROPHT
            ENDIF

C--      UNHOLY AND UNNATURAL MARRIAGE OF ZHU'S 
C--      AND ROSENFIELD'S RADIATION MODELS:
         IF(ABS(IHEATS).EQ.4) THEN
            CALL WED_RAD(LL)
            IF(ISW(28).EQ.1) CALL TROPHT
            ENDIF


         DO 2232 K=1,M$
            DO 2232 J=1,N$
               HEAT0(J,K)=HEAT(J,K)/TTOTH(K) ! CHANGE BACK TO
               COOL0(J,K)=COOL(J,K)/TTOTH(K) ! REAL TEMP. HEATING for output
2232           CONTINUE



C        NRDLS is no longer used.      
         NRDLS=NSTEP

         ENDIF


C     This also should still work correctly:, ie, convert back,
C     ie, HEAT0, COOL0 are ONLY used for OUTA below;  ie, NSTEP = 1->12, IDELH=12, 
ccpw32 
ccpw32  - DON'T NEED THIS - HEAT, COOL arrays are defined above
ccpw32
ccpw32      IF ( MOD(NSTEP,IDELH) .NE. 1 .AND. IDELH .NE. 1 ) THEN 
ccpw32
ccpw32         DO  K=1,M$
ccpw32            DO  J=1,N$
ccpw32               HEAT(J,K)=HEAT0(J,K)*TTOTH(K) 
ccpw32               COOL(J,K)=COOL0(J,K)*TTOTH(K) 
ccpw32               END DO
ccpw32            END DO
ccpw32
ccpw32         ENDIF



c write out real temperature heating
      CALL OUTA(1,NSTEP,NDAY,NYEAR,YP,ZP,HEAT0,N$,M$,0)
      CALL OUTA(2,NSTEP,NDAY,NYEAR,YP,ZP,COOL0,N$,M$,0)


C     GET MOMENTUM SOURCES

c      print *,'before getmomx'
      CALL GETMOMX
c      print *,'after getmomx'

      CALL PRMINMAX(PSIC(1,1,1),'PWV_AMP-MN',N$, M$,1)



C     CHECK FOR INERTIAL INSTABILITIES, MODIFY MERIDIONAL MOMENTUM FLUX,
C     (UBY ) IF NECESSARY

c      CALL ZINST2
       
cjer get stream function and velocities

      CALL XSTRM

      CALL OUTA(5,NSTEP,NDAY,NYEAR,YP,ZP,PSI,N$,M$,0)
      CALL OUTA(17 ,NSTEP,NDAY,NYEAR,YP,ZP,WS,N$,M$,0)
      CALL OUTA(18 ,NSTEP,NDAY,NYEAR,YP,ZP,VS,N$,M$,0)


      CALL OUTA(8 ,NSTEP,NDAY,NYEAR,YP,ZP,UB(1,1,IT1),N$,M$,0)

c      CALL OUTA(19,NSTEP,NDAY,NYEAR,YP,ZP,TH(1,1,IT1),N$,M$,0)


c     ***************************************************************       
C     ADVANCE THE CONSTITUENT FIELDS
c     ***************************************************************

C     If there are no transported dynamical tracers (ISW(7) < 1), then
C     jump over the tracer advection stuff (label 501).

      IF (ISW(7).LT.1) GOTO 501

C     SET UP VELOCITIES ON NEW GRID FOR TRANSPORT (WRBR, VRBR)
      call gettv
      
      CALL OUTA(6 ,NSTEP,NDAY,NYEAR,YPP,ZPP,VRBR,NP$,MP$,0)
      CALL OUTA(7 ,NSTEP,NDAY,NYEAR,YPP,ZPP,WRBR,NP$,MP$,0)



C --------  DO TRACER ADVECTION AND SOURCE AND SINK CALCULATION
C --------  TO EACH TRACER -----------------------------------

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          Note:                                     c
c                X is mass distribution              c
c               XC is mixing ratio                   c
c               X0 is initial distribution (mass)    c
c          These distributions are defined in start  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Loop is over the number of transported dynamical tracers,  ISW(7) = 2 (theta, ang momentum)


      DO NN=1,ISW(7)
    
C                          -- STORE OFF OLD MIXING RATIO,  xcold(NP$,MP$) (REAL*8)
          do k=1,mp$
          do j=1,np$
              xct(j,k,nn)=xc(j,k,nn)
              xcold(j,k) = x(j,k,nn)/RHO0(J,K)
          end do
          end do
    
C         NSTAR is no longer used.
ccelf          NSTAR=N

C                          -- RESET BOX MASSES TO RHO0
          CALL SETMASS



C  -- ADVECT  THETA, Angular momentum using the Lin and Rood scheme  - EF  (March 2008)
C     (modified from the constituent transport from the fixed model) 

       
          CALL TPCORE_DYN(np$, mp$, zres, dt, xcold, wrbr, xcnew)

C
C   after advection, reset bottom boundary (1.015 km) here to NCEP/GEOS 4 TEMPs
C   XCNEW(NP$,MP$) is THETA;   use UT1KM(NP$,2),  1=TEMP;  2=UBAR (m/sec) for current day
C                                    TTOTHX(MP$) is in COMMONC.INC
          if (NN .eq. 1) then
             do j=1,np$
                xcnew(j,1) = ut1km(j,1)*TTOTHX(1)
             end do
          endif


C  also reset angular momentum at sfc w/ NCEP UBAR: use XC00(NP$,MP$,NCON=2)
C                                             CST(NP$), convert UT1KM(NP$,2) to cm/sec
C
c  TURN ON FOR WCONV75 (was turned OFF for PW40, WCONV85), also do for levels 2 (WCONV82) - U3KM(NP$,2)
C
ccwconv85          if (NN .eq. 2) then
ccwconv85             do j=1,np$
ccwconv85                xcnew(j,1) = 100.*ut1km(J,2)*CST(J) + xc00(J,1,2)
ccwconv85                xcnew(j,2) = 100.*u3km(J,1)*CST(J)  + xc00(J,2,2)
ccwconv82                xcnew(j,3) = 100.*u3km(J,2)*CST(J)  + xc00(J,3,2)
ccwconv85             end do
ccwconv85          endif	



C                     -- UPDATE ADVECTED MIXING RATIO - xcnew(NP$,MP$),  REAL*8
C                        load into mass array X for proper initialization in GETMIX below
          do k=1,mp$
          do j=1,np$
              x(j,k,nn) = xcnew(j,k)*RHO0(J,K)
          end do
          end do



C                 -- ADVECT THE CONSTITUENT WITH
c                 -- Modified PRATHER'S SCHEME - NO LONGER USED (EF, March 2008)
cjer in trnsp.f
ccelfdyn          CALL DYNMY2_1(NN)
ccelfdyn          CALL DYNMX2_1(NN)
c...updated "mass" now in X array
C
              
C                          --  COMPUTE DIFFUSION AND CHEMISTRY BELOW

        IF ((NN .EQ. 1) .OR. (NN .EQ. 2)) THEN
           DO 5055 K=1,MP$
           DO 5055 J=1,NP$

c            nantest=test_nan( xc(j,k,nn) )
c            if(nantest) then
c              PRINT*, "NAN before GETMIX:"
c              PRINT*, "      J = ", J
c              PRINT*, "      K = ", K
c              PRINT*, "      N = ", N
c              PRINT*, "      X = ", X(J, K, NN)
c              PRINT*, "   RHO0 = ", RHO0(J, K)
c              STOP
c              ENDIF

           IF (XC(J, K, NN) .GE. 1.0E+15) THEN
              PRINT*, "XC blowup before GETMIX:"
              PRINT*, "XC before GETMIX"
              PRINT*, "      J = ", J
              PRINT*, "      K = ", K
              PRINT*, "      N = ", N
              PRINT*, "      X = ", X(J, K, NN)
              PRINT*, "   RHO0 = ", RHO0(J, K)
              STOP
              ENDIF

5055      CONTINUE
          ENDIF
	  
c getmix is in trcrlib.f
          CALL GETMIX(ll,NN,ASOR)

      END DO  ! LOOP OVER TRACER INDEX NN


501   CONTINUE


      do k=1,mp$
      do j=1,np$
         dum1x(j,k)=SM(j,k)/RHO0(j,k)
      end do
      end do
cc      CALL OUTA(18,NSTEP,NDAY,NYEAR,YPP,ZPP,dum1x ,NP$,MP$,0)


C  apply NCEP temperature adjustment here  -  TMODOFF is MODEL-NCEP so SUBTRACT this from XC array to adjust
C     use TMODOFF(NP$,MP$) from above, convert to THETA using ZPP(MP$) which is in KM
C     operate on the XC(NP$,MP$,NCON) - this was last updated in GETMIX  -  theta is NCON=1
C     use NSTEPS = DYSTPD = 12, so DT = DT0 = 7180.333 seconds,  DAYL = 86164 seconds
C     use change 1x per day on final time step, ie when LL = NSTEPS - NOT USED
cadj
cadj        if (LL .eq. NSTEPS) then
cadj           do 620 ik=1,MP$
cadj           do 620 ij=1,NP$
cadj               thetadj = TMODOFF(ij,ik) * exp(zpp(ik)*.288/7.)
cadj               XC(ij,ik,1) = XC(ij,ik,1) - thetadj/2.     ! *DT/(DAYL*1.)*8.
cadj 620       CONTINUE
cadj        ENDIF


ccc       type *, '  '
ccc       type *, ' IN TWODS: '
ccc       type *, 'LL, DT0, DT, DAYL, DYSTPD = ', LL, DT0, dt, DAYL, DYSTPD, IDOY360
ccc       type *, 'LL, DT, DYSTPD, IDOY360 = ', LL, dt, DYSTPD, IDOY360
ccc       type *, '  '

ccc       write (524) np$, mp$, M$
ccc       write (524) ttoth, zz, ypp, zpp
ccc       write (524) xc, xc00


C pw82 -  SMOOTH Fields on FINAL DIURNAL TIME STEP, first convert to TEMP, UBAR -  U3KM(NP$,2)
C    TTOTHX(MP$), xcnew(NP$,MP$), xxs0(NP$,MP$), xxs1(NP$,MP$) are REAL*8  - XC(NP$,MP$,NCON=2)
C
       IF (LL .eq. NSTEPS) then
            do k=1,mp$
            do j=1,np$
               xxs0(j,k) = xc(j,k,1)/TTOTHX(K)
            end do
            end do
    
            CALL SMOOTH5(xxs0, xxs1, NP$, MP$)

            do k=1,mp$
            do j=1,np$
               xc(j,k,1) = xxs1(j,k)*TTOTHX(K)
            end do
            end do

            do j=1,np$
               XC(J,1,1) = ut1km(J,1)*TTOTHX(1)
            end do


            do k=1,mp$
            do j=1,np$
               xxs0(j,k) = (xc(j,k,2) - xc00(j,k,2))/CST(J)
            end do
            end do
      
            CALL SMOOTH5(xxs0, xxs1, NP$, MP$)

            do k=1,mp$
            do j=1,np$
               xc(j,k,2) = xxs1(j,k)*CST(J) + xc00(j,k,2)
            end do
            end do

ccwconv85            do j=1,np$
ccwconv85               xc(j,1,2) = 100.*ut1km(J,2)*CST(J) + xc00(J,1,2)
ccwconv85               xc(j,2,2) = 100.*u3km(J,1)*CST(J)  + xc00(J,2,2)
ccwconv82               xc(j,3,2) = 100.*u3km(J,2)*CST(J)  + xc00(J,3,2)
ccwconv85            end do
       ENDIF



C     GET U BAR, THETA Prime, & THERMAL WIND 

c utfrmx in twodlib.f
      CALL UTFRMX     !(XC,XC00,UB,TH)

      CALL JGETUB
      CALL THERMW


c      CALL OUTA(20,NSTEP,NDAY,NYEAR,YP ,ZP ,UTHRM ,N$, M$ ,0)
      CALL OUTA(13,NSTEP,NDAY,NYEAR,YPP,ZPP,EPFLX ,NP$,MP$,0)


C     CYCLE THE INDICIES  -   EF: load in UBAR1,2,3 for output DIMENSION UB1(N$,M$), UB2(N$,M$), UB3(N$,M$)
C                                                             tempb(N$,M$), kyyb(N$,M$), wsb(N$,M$)
C  NOTE:  
C
C   TH(it1=3) and UB(it1=3) are updated in UTFRMX  - index IT1=3, IT3=1 are in COMMON, defined in START
C   TH(it3=1) is then used in GETALT2 to compute temperature array T at beginning of DO 1000 NSTEP loop above
C     but UB(it3=1) doesn't seem to be used
C
C   NOTE the TH array is THETA minus the global mean!!!
C
C     index IT2=2, but TH(2) doesn't seem to be used. UB(2) is used in THERMW
C
C                   TEMB, KYYB, WSB, KZZB are returned to the chemistry, compute temp from theta as in GETALT2
C
         DO 507 K=1,M$
         DO 507 J=1,N$
            TH(J,K,1) = TH(J,K,3)     ! - make these TH(3), UB(3) from UTFRMX -  previously TH(J,K,2), UB(J,K,2)
            UB(J,K,1) = UB(J,K,3)

            TH(J,K,2)=TH(J,K,3)
            UB(J,K,2)=UB(J,K,3)
C                                                      update temperature for chemistry from theta TH(1) as in GETALT2
            TEMPB(J,K) = (TH(J,K,1) + THG(K))/TTOTH(K)

            KYYB(J,K) = KYY(J,K)
            WSB(J,K) = WS(J,K)
C
C                                    these are ONLY used for REGRID below
            UB1(J,K)=UB(J,K,1)
            UB2(J,K)=UB(J,K,2)
            UB3(J,K)=UB(J,K,3)

C                                                these are NOT used
ccc            THST(J,K,1)=THST(J,K,2)
ccc            UBST(J,K,1)=UBST(J,K,2)
ccc            THST(J,K,2)=TH(J,K,1)
ccc            UBST(J,K,2)=UB(J,K,1)
507     CONTINUE


C  for output, use THX(NP$,MP$), UBX(NP$,MP$) arrays which are in extended grid (COMMONT.INC)
C   directly from theta, angular momentum computation (no interpolation)
C     -> tempbxo(NP$,MP$), ub3x(NP$,MP$), kzzb(NP$,MP$),  TGX(MP$), THGX(MP$), TTOTHX(MP$) are in COMMONC.INC
C                         also load in TEMPBX(NP$,MP$) in COMMOND.INC for use in RADIATION on next step
         DO 517 K=1,MP$
         DO 517 J=1,NP$
             TEMPBX(J,K) = (THX(J,K) + THGX(K))/TTOTHX(K)
            TEMPBXO(J,K) = (THX(J,K) + THGX(K))/TTOTHX(K)
            UB3X(J,K) = UBX(J,K)
            KZZB(J,K) = KZZ(J,K)
 517     CONTINUE



C   This problem should be eliminated with the new timing package:
c      (the following comment is from AD)
c      there is a timing problem if nsteps = 1 (same day time step in 
c      chem code as in dynamical model; we call timex 1 time above.
c      here we should call it nsteps-1 times, remember ll goes from 1
c      to nsteps

c      CALL OUTA(999,NSTEP,NDAY,NYEAR,YP(1),ZP(1),DUM1(1,1),1,1,0)


      CALL COMFORT  !PRINTS OUT MAX,MIN'S DURING RUN


C     CONDITIONAL RESTART
C     IMODUL contains the frequency (in days, 360-day calender) at which
C     the dynamics restart file is written.  It is set from ISW(5).  The
C     NSTEP.EQ.1 condition assures that the restart file is only written
C     ONCE on the appropriate day.
     
      IF (IMODUL.GT.0) THEN
         IF ((MOD(IDOY360, IMODUL).EQ.0).AND.(NSTEP.EQ.1))
     2      CALL DUMPA(.FALSE.)
         ENDIF


1000  CONTINUE




C  now interpolate DRAG terms from final time step to chemistry grid for output,  bdrag0(NP$,MP$)
C     dragcoup(7,NP$,MP$) is in COMMOND.INC,  YPP(NP$),ZPP(MP$) in COMMONC.INC
C       ->  DRAGCHM(7,NS$,MS$), yin(ns$), zin(ms$), BDRAG1(NS$,MS$)  units are in cm/sec^2
C

       do 875 iig = 2,6

           do 877 ik=1,MP$
           do 877 ij=1,NP$
 877          bdrag0(ij,ik) = dragcoup(iig,ij,ik)


           CALL BINTERP(ypp, NP$, zpp, MP$, BDRAG0,
     >                  yin, NS$, zin, MS$, 0, 0, BDRAG1)

          
           do 878 ik=1,MS$
           do 878 ij=1,NS$
 878          dragchm(iig,ij,ik) = bdrag1(ij,ik)

 875    CONTINUE



C  and interpolate QBARY from final time step to chemistry grid for output
C   QBRYW(N$,M$) is in COMMONW.INC,  YP(N$), ZP(M$) in COMMONC.INC
C   yin(ns$), zin(ms$), BDRAG1(NS$,MS$) here is QBARY, units are in 1/(cm-sec), load into DRAGCHM(7,NS$,MS$)


           CALL BINTERP(yp, N$, zp, M$, QBRYW,
     >                  yin, NS$, zin, MS$, 0, 0, BDRAG1)


           do 1777 ik=1,MS$
           do 1777 ij=1,NS$
 1777         dragchm(1,ij,ik) = bdrag1(ij,ik)




ccelf         write (778) N$, M$
ccelf         write (778) yp, zp, ypp, zpp
ccelf         write (778) gwmom
ccelf         write (778) epflx
ccelf         write (778) qbar
ccelf         write (778) qbary
ccelf         write (778) gwkzz
ccelf         write (778) kyy
ccelf         write (778) uthrm
ccelf         write (778) ub
ccelf         write (778) t
ccelf         write (778) psi
ccelf         write (778) ws




c Translate data to chem grids for output

        do 932 k=1,ms0
932     trr(k)=zin(k)+(zin(2)-zin(1))*.5

c put in a quick fix to get kzz fields right:

        do 1234 k=1,n$
           do 1234 j=1,m$
              if( kzz(k,j) .lt. 0.0) then
                print *,' negative kzz before interp: k,j,kzz = ',
     *              k,j,kzz(k,j)
              end if
              kzztemp(k,j)=alog( kzz(k,j) )
 1234   continue


C  EF:
        call regrid(ub1out,yin,zin,ns0,ms0,ub1,yp,zp,n$,m$,0)
        call regrid(ub2out,yin,zin,ns0,ms0,ub2,yp,zp,n$,m$,0)
        call regrid(ub3out,yin,zin,ns0,ms0,ub3,yp,zp,n$,m$,0)

        call regrid(tout,yin,zin,ns0,ms0,t,yp,zp,n$,m$,0)
        call regrid(wout,yin,zin,ns0,ms0,ws,yp,zp,n$,m$,0)
        call regrid(kyyout,yin,zin,ns0,ms0,kyy,yp,zp,n$,m$,0)
        
        call regrid(kzzout,yin,zin,ns0,ms0,kzztemp,yp,zp,n$,m$,0)
        do k=1,ns$
        do j=1,ms$
        kzzout(k,j)=exp( kzzout(k,j) )
        if( kzzout(k,j) .lt. 0.0) then
          print *,' negative kzz after interp: k,j,kzzout = ',
     *        k,j,kzzout(k,j)
        end if
        enddo
        enddo


c Limits on kyyout

c also limit kzz to be less than 5e4 above troposphere
        do k=1,ms0
         do j=1,ns0
                kyyout(j,k)=amin1(amax1(kyyout(j,k),3.0e9),1.8e10)
		if(zin(k).gt.15.) kzzout(j,k)=amin1(kzzout(j,k),5.e4)
         enddo
        enddo
	
        
         do ij=1,N$
           ypf(ij) = yp(ij)
         enddo

         do ij=1,NP$
           yppf(ij) = ypp(ij)
         enddo

         do ij=1,M$
           zpf(ij) = zp(ij)
         enddo

         do ij=1,MP$
           zppf(ij) = zpp(ij)
         enddo


      RETURN
      END
