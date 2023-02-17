C
C  /misc/kah02/9g1> em main_9ga.f &
C  [3] 7127
C
C
	PROGRAM MAIN
C       MAIN COMMON ASSEMBLY FOR SPEC TRANS 2-D MODEL  -- MAIN_9GA.F 
C            has new PPM transport (TPCORE2D) call
C       also writes out final 360 days of mix rat at top level
C

         include "com2d.h"
         include "com_aerd.h"
         include "com_aerg.h"
         include "com_aers.h"


! BThomas Nov2017 - some stuff for outputting data I want:
      INTEGER outS,outX,outnum
      INTEGER nsp 
      PARAMETER(nsp=7)   !number of consists in ncolgd to be output
        !CAREFUL! nsp must equal NSPC in colden subroutine
      PARAMETER(outS=80) !65) !number of constituents to output; 
        ! S$ is larger than this, see 2Dmodel_ConstituentList_constit15_9ua.pdf
      PARAMETER(outX=6)  !number of extra values to output; such as temperature
      PARAMETER(outnum=outS+outX)
      real*4 outdatarr(outnum,L$,Z$)


C  Common for mapping into GSFC code - NOTE: CAN't have any L$, Z$ here in these COMMONS

      COMMON/CAER1/ISPMAP(NFSP)
      COMMON/CAER2/DTIMEA(NTIME), DTIMEB(L$,NTIME), 
     >                PMTIME(NTIME), DLITE(NHT)
      COMMON /AERHET1/KHETA(NKR,18,46), GAMAER(NKR,18,46)

      REAL KHETA, GAMAER, PMTIME, DLITE
      DOUBLE PRECISION DTIMEA, DTIMEB

C  JAVAER are the diurnally avg J-coeffs   ,   JAER are the diurnally varying Js
C  CNDCA = CN array for 18+3 diurnal time steps for up to 1 daily iterations, CNDCA is in Num dens., NFSP=40
C
ccc         REAL*8 JAVAER(NJR+3,L$,Z$), CNDCA(S$,ntime+3,L$,Z$)
ccc         REAL*8 JAER(NJR+3,NTIME,L$,Z$), OXG(L$,Z$), ch4out(L$,Z$)
C
Cj  common for diurnally avged J-coeffs - ONLY USE if dayime avg J's are needed
cj      COMMON/JAER1/JAVAER(NJR,18,46), AJOUT(2,ntime,18), 
cj     >             AJAVDAY(3,NJR,18,46)
cj         REAL AJOUT(2,ntime,18), AJAVDAY(3,NJR,18,46)
cj
cj      COMMON/JAER1/AJOUT(2,ntime,18)


      COMMON/CHIDSFAC/CHIDA(18,L$),SFACA(18,L$)

C INSERT SUBSONIC/SUPERSONIC 12/24/92   
      COMMON/EMISSION/H2OEM(L$,Z$), XNOXEM(L$,Z$), HCEM(L$,Z$), 
     >      COEM(L$,Z$), AGECN(3,L$,Z$)

C CNANOUT array to output C array (mix rat) for 15 days before NaN's
      COMMON/COMNAN/cnanout(15,S$,L$,Z$) 


C      common block for the coupled model dynamics: dimensions are as in SETDAILY

      COMMON/TRCOUP/tcoup(L$,Z$X), ubard(L$,Z$X), kyycoup(L$+1,Z$X), 
     >              kzzcoup(L$,Z$X+1), wcoup(L$,Z$X+1), DRAGOUT(7,L$,Z$)


c  common of diurnally avged heating rates computed from the photolysis in PHOTHEAT (K/sec)
C  SORADHEAT(15,L$,Z$) are the individual bands from SORAD.f
c 
         COMMON/CPHEAT/PHEAT(5,L$,Z$), SORADHEAT(15,L$,Z$)

         REAL H2OJET,CH4JET,NOXJET,COJET

c	common/bcnox/ibc(l$,z$),bcclono2(l$,z$),bcbrono2(l$,z$)
c        common/noxc/bcclo24(l$,z$),bcn2o5(l$,z$),cno2(l$,z$),
c     c  bcbro24(l$,z$),gamloc(l$,z$),zzz(20,l$,z$),yy(4,l$,z$),
c     c  xxz(4,l$,z$),alpha(l$,z$),betaz(l$,z$),isw2(l$,z$),
c     c  isw(l$,z$),icount(l$,z$)

        real stime, extime, cputime(2),cnan(S$,L$,Z$X),cnan58(S$,L$,Z$X)
        REAL*4 ZZOUT(L$,Z$)
        REAL*4 HEATOUT(2,L$,Z$), COOLOUT(L$,Z$)

        dimension conv(l$,z$,t$),xsave(l$,z$),OUTO3(L$,Z$)
        DIMENSION NEWSP(30,18,Z$),OUTSP(18,Z$),JNEW(30),
     >    PHOTJOUT(20,L$,Z$),WAVELOW(39),WAVEHIGH(39),COLO3(L$,360)    ! ,COLO3(2,L$,360)

C AYelland 02Nov2020 - The following are used for the ionization and lighting modifers to the NOy production
C                      They are declared in com2d_9hb_ions.h and are used in main_9hb.f & solv_9hb.f

      OPEN(UNIT=1585, FILE='ionsflag.txt')
        READ(1585,*) ionsflag
      CLOSE(UNIT=1585)
      print *, "MAIN: ionsflag = ", ionsflag
      ionsflag_solv = 0

      OPEN(UNIT=1588, FILE='lightnflag.txt')
        READ(1588,*) lightnflag
      CLOSE(UNIT=1588)
      print *, "MAIN: lightnflag = ", lightnflag
      lightnflag_solv = 0

      OPEN(UNIT=1586, FILE='ionStart.txt')
        READ(1586,*) ionStart
      CLOSE(UNIT=1586)
      print *, "MAIN: ionStart, lightnStart = ", ionStart

      OPEN(UNIT=1587, FILE='runNum.txt')
        READ(1587,*) runNum
      CLOSE(UNIT=1587)
      print *, "MAIN: runNum = ", runNum


C CEF =  code that is turned off for testing purposes, but should be turned on in final model
!BThomas Nov.2017 - I have checked all CEF lines below and none should be turned on for my purposes

!BThomas Nov.2017 - Open output files to be written to daily. 
! 'alt-lat_daily' has altitude-latitude-time dependent stuff, including chemical
!                 constituent number densities and some temp/dynamics stuff
!                 along with the latitude and altitude arrays
	open(unit=5130,file='BToutput_alt-lat_daily.out',
     >       form='formatted',status='new')
	open(unit=5131,file='BToutput_coldens_daily.out',
     >       form='formatted',status='new')
	open(unit=5132,file='BToutput_coldenhno3_daily.out',
     >       form='formatted',status='new')
! 'alts-O3-NO2_TUV_daily' has some of the same data as 'alt-lat_daily'
!                         but I'm generating it separately for ease of use with TUV
	open(unit=5133,file='BToutput_alts-O3-NO2_TUV_daily.out',
     >       form='formatted',status='new')
	outP = 0 !initialize
!	
c initialize output array:
	 do ik=1,Z$
	    do ij=1,L$
	       do is=1,outnum
		  outdatarr(is,ij,ik) = 0.0
	       end do
	    end do
	 end do
!


C  SWITCH FOR ZONAL MEAN or PARCEL model temperature, latitude PDFs:
C                                                       ITPDF/ILPDF = 0 gives zonal average values
C                                                       ITPDF/ILPDF = 1 gives parcel PDFs
        ITPDF = 1
        ILPDF = 1

        IZONAVGT = 0
        IF (ITPDF .EQ. 0) IZONAVGT = 1
C                                          parameter IZONAVGT for using zonal avg temps in GAS/HET reactions
C                                          reset = 1 below if using coupled model

c  get inputs

	CALL CONTROL

ccc        print *,' CONTROL '

C    for coupled model, use ADJUSTED PDFs for GAS/HET/PSC reactions and INSOLATION correction
C                                   for SETDAILY, ICOUP is read in CONTROL
ccc        IF (ICOUP .EQ. 1) THEN 
ccc             IZONAVGT = 1
ccc             ILPDF = 0
ccc        ENDIF

C
C  initialize year counter for current run (for coupled model and chemical relaxations)
C  define ISTST for steady state (IYRCT, ISTST in COMMON), LSTSTATE defined in CONTROL.f

         IYRCT = 0

         ISTST = 0
         if (LSTSTATE) ISTST = 1


C  define INYEARS = number of years in run (INYEARS is in COMMON)

         INYEARS = INT(INDAYS/360)



C  set up some contants by calling input, also SETS up P and PRESS arrays (in COMMON) -- BASE8A

	CALL INPUT

        print *,' INPUT '


c  galactic cosmic rays
      !other ionization source input (e.g. ionization by gamma-rays, cosmic rays, or SEPs)
      !is read-in within time loop - AYelland 02Nov2020

	CALL GCRS



c  aerosol data - NOTE  read both distributions (appropriate for 6 mos)

	CALL AEROSOLS


C READ IN TD AEROSOL SURFACE AREAS (FROM D. CONSIDINE) - 8/20/96

        CALL AEROSOL_DBC


c reaction rates

           CALL RCREAD


c  fixed ozone?
CEF	YSO3FIX=.false.
CEF	if(YSO3FIX)CALL O3SBUV
	

c  aircraft fleets?
c	YSHSCT=.true.
c	if(YSHSCT)CALL EMISSREAD   

c  get boundary conditions
	CALL BCIN

        print *,' BCIN '


C   READ INITIAL PROFILES BY CALLING UNDUMP.

        CALL UNDUMP

        print *,' UNDUMP '



        CALL CHINITAER

        print *,' CHINITAER'

C	write(33)c
C
C  Initialize ozone column density array
        DO 650 ID=1,360
                DO 650 IJ=1,L$
                  COLO3(IJ,ID)=0.0E0

cc79                COLO3(1,IJ,ID)=0.0E0
cc79                COLO3(2,IJ,ID)=0.0E0
 650            CONTINUE



C  SET IYRLIFE FOR LIFETIME CALCULATION
	IYRLIFE=IYR
c                                           Note: DAY360 now goes 1.5-360.5, so IDAY360 => 1-360
   	DAY=(365./360.)*DAY360
   	DAYCHEK=DAY360
   	IDAY360=INT(DAY360)
   
   	DT=DT0*365./360.
                                       ! also do DT in REAL*8
        DT8=86400.d0*365.d0/360.d0
                                       ! get starting year (+ fraction of year) for fort.33 below
        yrst = (DAY360 - .5)/360. + iyr


C  GET PHOTOCHEMICAL INPUT DATA

	CALL SOLFLIN

	CALL CROSSIN
c        print *,' after call crossin'


	CALL CSIN_NO
C Read in NO photodissociation information from Minschwaner & Siskind (1993)


C  READ IN look-up tables.  TEMPS, DYNAMICS, HEATING RATES now read in READSTR below (for 1979-2010 run)

        CALL TEMPIN

        print *,' TEMPIN '
C
C

C  READ IN TABLE OF REDUCED SOLAR FLUXES GENERAGED BY RANDY KAWA - 1/24/95

        CALL PHOTIN


C
C  SET UP DYNAMICS ARRAYS  for  all  times  (the 72 5-day means), stored in COMMON
C       TROPKZZ sets up tropospheric Kzz's 
C     - do for Jan. 1979 up through Dec. 2010,  LDYNOUT is the logical for writing out dynamics
C

        LDYNOUT = .TRUE.
cc        LDYNOUT = .FALSE.
        if (icoup .EQ. 1) LDYNOUT = .FALSE.



C  set ICLIM (in COMMON) based on ICOUP from control.dat/control.f:
C
C     ICLIM = 1 (default) when:
C                ICOUP = 0  (fixed model, clim dynamics/tprob/lprob)
C                ICOUP = 1  (coupled model)
C
C     Exceptions:
C
C       ICOUP = -1 (fixed model, interannual dynamics/tprob/lprob)  => ICLIM = 0
C
C       ICOUP = -3 (fixed model testing, just use 1st time period
C                                            dynamics in SETDAILY)  => ICLIM = 3

        iclim = 1
        if (icoup .EQ. -1) iclim = 0

        if (icoup .EQ. -3) iclim = 3



C  when using already-calculcated transport fields, run STREAMF ONCE 
C     to initialize needed variables, and define transport 1X

                 ! in10 = # of 5-day dynamics periods per year to run,   inyr = number of years
       inyr = 1    ! 32    ! 2
       in10 = 1    ! 72

       itdyn = in10*inyr            ! total number of dynamics outputs written to fort.37, fort.38


C  READ IN yearly look-up tables of Kyys, TEMPS, and HEATING RATES in READSTR (for 1979-2010 run)

          do 445 iiyr = 1,inyr 

                 CALL READSTR(iiyr)

          do 447 im10 = 1,in10

                 im324 = (iiyr-1)*72 + im10

           CALL STREAMF(im10, iiyr, im324)

           CALL TROPKZZ(im10, iiyr, im324)

           IF (LDYNOUT) CALL DYNOUT(im10, iiyr, im324)
 447   CONTINUE

 445   CONTINUE

cc      print *,' STREAMF '



C  for both FIXED and COUPLED MODELS, read in clim transport fields for intializing in SETDAILY
C                                                and UBARCG(L$,Z$X,74)
        CALL GETDYN
        CALL UBARDIAG

C
C  for coupled model, need to initialize arrays in COMMON block TRCOUP, using clim transport fields
C   which were just read in GETDYN (use Jan), these then get loaded into chemistry arrays in SETDAILY (next)
C
C     tcoup(L$,Z$X), ubard(L$,Z$X), kyycoup(L$+1,Z$X), kzzcoup(L$,Z$X+1), wcoup(L$,Z$X+1),DRAGOUT(7,L$,Z$) 
C
C     TEMPALL(L$,Z$X,74), KYYALL(L$+1,Z$X,74), KZZTOT(L$,Z$X1,74), WALL(L$,Z$X1,74) are in COMMON
C
C
c     ALSO  set up coupled model inputs

       IF (ICOUP .EQ. 1) THEN 

           CALL XCOUPLED_IN(IYR, IDAY360)

           DO 7505 IK=1,Z$X
           DO 7505 IJ=1,L$
              TCOUP(IJ,IK) = TEMPALL(ij,ik,2)
 7505         UBARD(IJ,IK) =  UBARCG(ij,ik,2)

           DO 7707 IK=1,Z$
           DO 7707 IJ=1,L$
              DRAGOUT(1,IJ,IK) = QYCG(ij,ik,2)
 7707         DRAGOUT(2,IJ,IK) = EPCG(ij,ik,2)

           DO 7506 IK=1,Z$X
           DO 7506 IJ=1,L$+1
 7506         KYYCOUP(IJ,IK) = KYYALL(ij,ik,2)

           DO 7507 IK=1,Z$X1
           DO 7507 IJ=1,L$
              WCOUP(IJ,IK) = WALL(ij,ik,2)
              KZZCOUP(IJ,IK) = KZZTOT(ij,ik,2)
 7507  CONTINUE

       ENDIF


C
C
C  Call SETDAILY -- interpolates dynamics to daily values, and sets up aerosol and major species, ZKM, etc.
C    SETDAILY keys on IDAY360


         CALL SETDAILY

ccc        print *,' SETDAILY  '



	  CALL REACTION


		CALL COLDEN

C Initialize column ozone density into array on first day and store
                DO 1550 IJ=1,L$
                   ICOL=IDAY360
		   IF(ICOL.LT.1)GO TO 1550
                   IF(ICOL.GT.360)THEN
                      WRITE(6,2100)ICOL,IYR,IDAY360
                      GO TO 1550
                      ENDIF

                  COLO3(IJ,ICOL)=NCOLGD(2,IJ,1)

cc79                COLO3(1,IJ,ICOL)=NCOLGD(2,IJ,1)
cc79                COLO3(2,IJ,ICOL)=NCOLGD(3,IJ,1)
 1550            CONTINUE

 991             dum=1


C
C    initialize CNANOUT(15,S$,L$,Z$) array for writeout when NaN's occur (IN MIXING RATIO)
C
	do 915 ik=1,Z$
	do 915 ij=1,L$
	do 915 is=1,S$
 915       cnanout(1,is,ij,ik) = c(is,ij,ik)/m(ij,ik)


C  Reset special water vapor constituents at tropopause (90.3 mbar), 25S-25N 
C    to interannual HALOE values (PPMV), cn(75);  and climatological HALOE, cn(76) 
C
C    H2OBC(18,360,9), for 1992-April 2000;  H2OBCCL(18,360)
C    use 1992 values for 1960-1992
C
CEF      do 4000 ij=7,12
CEF        IF (IYR .le. 57) cn(75,ij,9) = h2obc(ij,iday360,1)*1.E-6*M(IJ,9) 
CEF        IF (IYR .ge. 58) 
CEF     >             cn(75,ij,9) = h2obc(ij,iday360,iyr-56)*1.E-6*M(IJ,9) 
CEF
CEF        c(75,ij,9) = cn(75,ij,9)
CEF
CEF        cn(76,ij,9) = h2obccl(ij,iday360)*1.E-6*M(IJ,9) 
CEF         c(76,ij,9) = cn(76,ij,9)
CEF 4000    CONTINUE
CEF

c  Initialize simple time increasing age of air to proper year everywhere - CN(78)
cc                              don't do with new initialization
cc        do 4500 ik=1,Z$
cc        do 4500 ij=1,L$
cc           c(78,ij,ik) = (1935.+iyr)*M(IJ,IK)
cc 4500      cn(78,ij,ik) = (1935.+iyr)*M(IJ,IK)


C
c  INitialize C-14, convert from mixing ratio to #/cm3,   , C14(L$,Z$) is in COMMON (REAL*4)  - CN(21)
C    then will re-intialize below, since Carbon-14 initial conditions are for OCTOBER 1963

          do 4504 ik=1,Z$
          do 4504 ij=1,L$
             c(21,ij,ik) = c14(ij,ik)*4.82E-18*M(IJ,IK) 
 4504        cn(21,ij,ik) = c14(ij,ik)*4.82E-18*M(IJ,IK) 


C
C    Initialize new extended level arrays of transported species only with initial profiles,
C          keep constant mixing ratio above level Z$ to start 
C
        do 2150 ik=1,Z$
        do 2150 ij=1,L$
        do 2150 ist=1,ITRANS
           c58(ist,ij,ik) = c(inttra(ist),ij,ik)
 2150      cn58(ist,ij,ik) = c(inttra(ist),ij,ik)

        do 2151 ik=Z$+1,Z$X
        do 2151 ij=1,L$
        do 2151 ist=1,ITRANS
           c58(ist,ij,ik) = c(inttra(ist),ij,Z$)/m(ij,Z$)*m(ij,ik)
 2151      cn58(ist,ij,ik) = c(inttra(ist),ij,Z$)/m(ij,Z$)*m(ij,ik)


C  also for new high resolution model, define indicies for 35 km, 40 km, and 58 km for HETCHEM and NATICE
C     ZALT(Z$X), ik35, ik58, ik40  in COMMON

          do 778 ik=1,Z$
             if (zalt(ik) .ge. 35.) then
                ik35 = ik-1
                goto 9999 
             endif
 778      CONTINUE

 9999     do 779 ik=1,Z$
             if (zalt(ik) .ge. 58.) then
                ik58 = ik-1
                goto 9988 
             endif
 779      CONTINUE

 9988     dum=1  

C                              also define index IKMES for 60 km to omit mesospheric chemistry in AER model
          do 777 ik=1,Z$
             if (zalt(ik) .ge. 60.) then
                ikmes = ik-1
                goto 9991 
             endif
 777      CONTINUE

 9991     dum=1  


          do 877 ik=1,Z$
 877         if (zalt(ik) .lt. 40.) ik40 = ik


C  also for new high resolution model, define indicies for 35N and 55N for SF6 input in SOLVER,  
C                                                                  ij35, ij55 in COMMON,  LAT(L$)
C  also define SP and NP indicies for NCB, ijsp, ijnp in COMMON
C    just take the first latitude equal to or equatorward of 65S, 65N

          do 5788 ij=L$,1,-1
             if (LAT(ij) .ge. 32.5) ij35 = ij
             if (LAT(ij) .ge. -65.) ijsp = ij
 5788     CONTINUE

          do 479 ij=1,L$
             if (LAT(ij) .le. 57.5) ij55 = ij
             if (LAT(ij) .le. 65.) ijnp = ij
 479      CONTINUE


ccccccc      print *,' ik35, ik58, ij35, ij55, ik40 = ',ik35,ik58,ij35,ij55,ik40


C  set Solid HNO3 and Solid H2O = 0.0 above 35 km
        do 4704 ik=ik35+1,Z$
        do 4704 ij=1,L$
           c(66,ij,ik) = 1.D-12
           c(67,ij,ik) = 1.D-12
           cn(66,ij,ik) = 1.D-12
           cn(67,ij,ik) = 1.D-12
 4704   CONTINUE

        do 4705 ik=ik35+1,Z$X
        do 4705 ij=1,L$
           c58(26,ij,ik) = 1.D-12
           c58(27,ij,ik) = 1.D-12
           cn58(26,ij,ik) = 1.D-12
           cn58(27,ij,ik) = 1.D-12
 4705   CONTINUE



C    and load in CNOUT(S$+8,L$,Z$) (REAL*4 defined in COMMON) for output of initial conditions
C                                             TEMP(L$,Z$X), ubard(L$,Z$X) in m/sec
C
C  also initialize HEATOUT(2,L$,Z$), COOLOUT(L$,Z$), DRAGOUT(7,L$,Z$) arrays as necessay for output
C                                    in K/day, m/sec/day, 1.e11 (1/m-s)
C
C                               initialize PHEAT(5,L$,Z$), SORADHEAT(15,L$,Z$) for output
C
	do 462 ik=1,Z$
	do 462 ij=1,L$

           heatout(1,ij,ik) = 0.0
           heatout(2,ij,ik) = 0.0
           coolout(ij,ik) = 0.0

           dragout(3,ij,ik) = 0.0
           dragout(4,ij,ik) = 0.0
           dragout(5,ij,ik) = 0.0
           dragout(6,ij,ik) = 0.0


           do 475 ihh=1,5
 475          pheat(ihh,ij,ik) = 0.

           do 476 ihh=1,15
 476          soradheat(ihh,ij,ik) = 0.
           

           do 452 is=1,S$
 452         cnout(is,ij,ik) = cn(is,ij,ik)

           cnout(S$+1,ij,ik) = TEMP(ij,ik)
           cnout(S$+2,ij,ik) = UBARD(ij,ik)
           cnout(S$+3,ij,ik) = HEATOUT(1,ij,ik)
           cnout(S$+4,ij,ik) = COOLOUT(ij,ik)

           cnout(S$+5,ij,ik) = DRAGOUT(1,ij,ik)
           cnout(S$+6,ij,ik) = DRAGOUT(2,ij,ik)
           cnout(S$+7,ij,ik) = DRAGOUT(3,ij,ik)
           cnout(S$+8,ij,ik) = DRAGOUT(4,ij,ik)
                                                 ! load Ray fric into CN(32), Fixed MOD DELF into CN(58)
           cnout(32,ij,ik) = DRAGOUT(5,ij,ik)
           cnout(58,ij,ik) = DRAGOUT(6,ij,ik)
462	continue


C  load in latent heating HEATOUT(2,ij,ik) into cnout(59) - washout freq - use levels 21-46
C
	do 463 ik=1,26
	do 463 ij=1,L$
 463       cnout(59,ij,ik+20) = HEATOUT(2,ij,ik)



C OUTPUT initial values, first output gridding and starting year; latitude grid for f14 total ozone file
C       LAT4(L$), LATEG4(L$+1), ZALT90(Z$), PRES90(Z$), ZALT(Z$X), ZALTE(Z$X1) are REAL*4
C
C    also w/ coupled model, now write out dynamics on (33)
C
           write (33) L$, Z$, yrst
           write (33) LAT4, LATEG4, zalt90, pres90, zalt, zalte

           write (33) cnout
ccccc           write (33) cnout, cdn

C                             ! W(L$,Z$X1), V(L$+1,Z$X) are REAL*8
C                             ! EKYY(L$+1,Z$X), EKYZ(L$+1,Z$X1), EKZZ(L$,Z$X1) are REAL*4
           write (33) W
           write (33) V
           write (33) EKYY
           write (33) EKYZ
           write (33) EKZZ




           write (17) L$, Z$, ntime
           write (17) LAT4, zalt90, pres90, zalte
ccccc           write (17) noonrat


           write (34) L$, Z$, Z$X
           write (34) LAT4, zalt90, pres90, zalte
           write (34) cn58


cc           write (547) L$, Z$
cc           write (547) LAT4, zalt90, pres90, zalte

cc           write (19) L$, Z$
cc           write (19) LAT4, zalt90, pres90, zalte

cc           write (22) L$, Z$
cc           write (22) LAT4, zalt90, pres90, zalte


           WRITE(14) L$, yrst
           WRITE(14) LAT4


C
C  daily output for heating rates on fort.222; PHEAT(5,L$,Z$), SORADHEAT(15,L$,Z$), HEATOUT(2,L$,Z$),COOLOUT(L$,Z$)
C
           write (222) L$, Z$
           write (222) LAT4, zalt90, pres90

           write (222) SORADHEAT
           write (222) PHEAT
           write (222) HEATOUT
           write (222) COOLOUT


C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c  TIME LOOP

ccj	ict = 0


c2ef
c2ef        print *, '   '
c2ef        print *, 'BEFORE  DO 100 '
c2ef        print *, 'ClONO2  =   ', cn(30,10,13)
c2ef        print *, 'OZONE  =   ', cn(4,10,13)
c2ef        print *, '   '
c2ef

!Some outputs (BT Feb.2018):
!Pressure array:
!	   open(5222,file="BT_pressure-array.out")
!	   do ik=1,Z$
!	      write(5222,'(1x,1i10,1f20.6)') ik,press(ik)
!	   end do
!	   close(5222)
!
!Pressure array into main output file:
	   do ik=1,Z$                                                                              
              write(5130,'(1x,1f20.6)') press(ik)                                           
           end do  
!Lat and alt arrays, into main output file:
	   do ij=1,L$
	      write(5130,'(1x,1e16.8)') lat4(ij)
	   end do
	   do ik = 1,Z$
	      write(5130,'(1x,1e16.8)') zalt90(ik)
	   end do
	   do ik=1,Z$
	      do ij=1,L$
		 write(5130,'(1x,1e16.8)') zkm(ij,ik)
	      end do
	   end do
!-- end outputs

        stime = etime(cputime)/3600.

C
C  with New AER fast chemistry, we'll use this order:
C      1) fast chemistry
C      2) transport
C      3) slow chemistry (SOLVER)


c  load up CH4 for output

cc        do 4514 ik=1,Z$
cc        do 4514 ij=1,L$
cc 4514      ch4out(ij,ik) = cn(18,ij,ik)
C

        print *,' CHEMISTRY = AERCHEM  '

	DO 100 I=1,INDAYS !Start of time loop here

c           DO 100 I=1,1
c  we are assuming we are using a 1 day timestep for chemistry, also in itday parameter

             ITDAY = I
	DAY360=DAY360+1.
                                         ! update IYRCT
	IF(DAY360 .GT. (YRL+1.))THEN
		DAY360=DAY360-YRL
c                    WRITE(14)COLO3
		IYR=IYR+1
                iyrct = iyrct + 1
	ENDIF

	DAY=(365./360.)*DAY360
	DYT360=DAY360-0.5+IYR*360.
	IDAY360=INT(DAY360)


cccc        print *,' icoup, iday360, iyr = ', ICOUP, IDAY360, IYR
      if (iday360.eq.1) print *,'ICOUP, iday360, iyr= ',
     >                           ICOUP,IDAY360,IYR
cccc      if (mod(iday360-1,30) .eq. 0.) print *,'iday360, iyr= ',IDAY360,IYR

c       print *,' timestep, day360, day, iyr, ict=',
c     *  i,day360,day,iyr,ict


C  interpolate to current day for use in fixed & coupled models (SAD, solar cycle)

       CALL DAILYINT

!Put time values into main output file (BT Feb.2018)
	write(5130,'(1x,1i10)') i
	write(5130,'(1x,1i10)') iday360	
!--end output

C ----------------------------------------------------------------------

C AYelland 02Nov2020 - Reading in cosmic ray ionization data to modify NOy production 
C and develop lightning multiplier

    !Read-in for first run at start day "ionStart"
      if ((runNum .eq. 1) .and. (I .eq. ionStart)) then
        
        if (ionsflag .eq. 1) then
          CALL ionSourceIn
            print *, "ionSource is Read-in ; during run 1"
C             print *, ionSource  !Checking Read-in
          ionsflag_solv = 1

          if (lightnflag .eq. 1) then

            CALL GCRionsIn
              print *, "GCRions is read-in ; during run 1"
C               print *, GCRions  !Checking Read-in

C           Calculate ionsRatio for lightning effect for all ionization read-in
            CALL ionsRatio_calc
              print *, "ionsRatio is calculated ; during run 1"
C               print *, ionsRatio  !Checking Read-in
            lightnflag_solv = 1
          end if
        endif
      endif

    !Read-in for all following runs at start of run
      if ((runNum .ne. 1) .and. (I .eq. 1)) then

        if (ionsflag .eq. 1) then
          CALL ionSourceIn
            print *, "ionSource is read-in ; after run 1"
C             print *, ionSource  !Checking Read-in
          ionsflag_solv = 1

          if (lightnflag .eq. 1) then

            CALL GCRionsIn
              print *, "GCRions is read-in ; after run 1"
C               print *, GCRions  !Checking Read-in

C           Calculate ionsRatio for lightning effect for all ionization read-in
            CALL ionsRatio_calc
              print *, "ionsRatio is calculated ; after run 1"
C               print *, ionsRatio  !Checking Read-in
            lightnflag_solv = 1
          end if
        endif
      endif
C end ions stuff
	
C ----------------------------------------------------------------------
c  if running coupled model, compute heating rates, dynamics

       IF (ICOUP .EQ. 1) CALL XCOUPLED (DT, IYR, IDAY360, IYRCT, Z$X,
     >                   LAT4, LATEG4, ZALT90, ZALT, ZALTE, EPCG, CN, M,
     > tcoup, ubard, kyycoup, kzzcoup, wcoup, HEATOUT, COOLOUT, DRAGOUT)


C
C  Now re-initialize Carbon-14 on Oct. 15, 1992 for INTERANNUAL transport,
C     and a 2nd time on Oct. 15, 2012 for CLIMATOLOGICAL transport
C     since the initial conditions are for OCTOBER 1963, -  convert from mixing ratio to #/cm3,
C     C14(L$,Z$) is in COMMON (REAL*4)  - CN(21), also initialize month counter, TMONC14

      if (iyr .eq. 57  .or.  iyr .eq. 77) then
        if (iday360 .eq. 285) then

           do 4505 ik=1,Z$
           do 4505 ij=1,L$
             c(21,ij,ik) = c14(ij,ik)*4.82E-18*M(IJ,IK) 
 4505        cn(21,ij,ik) = c14(ij,ik)*4.82E-18*M(IJ,IK) 
C
C  Also re-initialize new extended level arrays of transported species only with initial profiles,
C      keep constant mixing ratio above level Z$ to start 
C
           do 2250 ik=1,Z$
           do 2250 ij=1,L$
             c58(35,ij,ik) = c(21,ij,ik)
 2250        cn58(35,ij,ik) = c(21,ij,ik)

           do 2251 ik=Z$+1,Z$X
           do 2251 ij=1,L$
             c58(35,ij,ik) = c(21,ij,Z$)/m(ij,Z$)*m(ij,ik)
 2251        cn58(35,ij,ik) = c(21,ij,Z$)/m(ij,Z$)*m(ij,ik)

           TMONC14 = 0.0
        ENDIF
      ENDIF



c  set up boundary conditions if this is a time dependent run

ccccccc 	IF(.NOT.LSTSTATE) CALL BOUNDC
          CALL BOUNDC


C  Call SETDAILY -- interpolates dynamics to daily values, and sets up aerosol and major species, ZKM, etc.

          CALL SETDAILY 

ccc       print *,' SETDAILY  '


C  Now  CALL REACTION  every time step to define reaction rates and het. gammas every
C     day with changes in temperature and water vapor

 	  CALL REACTION

ccc        print *,' REACTION  '


	  CALL COLDEN
cccc        print *,' COLDEN  '

C Put column ozone density into array each day and store
                DO 550 IJ=1,L$
                   ICOL=IDAY360+1
                   ICOL=IDAY360
                   IF(ICOL.LT.1)THEN
                      ICOL=360
                   ENDIF
                   IF(ICOL.GT.360)THEN
                      WRITE(6,2100)ICOL,IYR,IDAY360
 2100                 FORMAT(' ICOL=',I5,' IYR=',I3,' IDAY360=',I5)
                      GO TO 550
		   ENDIF
c                   if(ij.eq.1)print *,' iday360(550)=',iday360,
c     *  ' icol=',icol,' day360=',day360

                COLO3(IJ,ICOL)=NCOLGD(2,IJ,1)

cc79                COLO3(1,IJ,ICOL)=NCOLGD(2,IJ,1)
cc79                COLO3(2,IJ,ICOL)=NCOLGD(3,IJ,1)
 550            CONTINUE


ccccc        print *,' COLO3 '


C  CALL AER diurnal CHemisry  (GSFC daytime avg chemistry NO LONGER USED)

           CALL AERCHEM


C
c  Now call separate SOLVER for extended levels using skeleton chemistry (no meat), no diurnal chemistry
C                                                      
 	do 920 ik=Z$+1,Z$X
 	do 920 ij=1,L$
  		CALL SOLVER58(IJ,IK)
 920	continue
 

         CALL RAINOUT




        IF(IYR.NE.IYRLIFE)YSLFPRNT=.TRUE.

 	CALL LIFETIME
c        print *,' after call lifetime'

        IF(IYR.NE.IYRLIFE)THEN
          IYRLIFE=IYR
          YSLFPRNT=.FALSE.
       ENDIF

c  update constituents, check for NaN's using ISNAN function (for Linux) 
C    also, Set limit to 1.E-12 number density
c    as was done in SOLVER, since HETCHEM may give zero's for Solid HNO3 and ice
C    (m at 115 km ~1.e12, and ~3.e19 at the ground, so the minimum mixing ratio will be ~3.E-31)
C    BUT DON'T ADJUST RAINOUT PARAMETERS
C      also load in REAL*4 CNOUT(S$+8,L$,Z$) array (in COMMON) for output

	do 950 ik=1,Z$
	do 950 ij=1,L$
	do 950 is=1,S$
          if (is .ne. 58  .and.  is .ne. 59) then
                IF (cn(is,ij,ik) .LT. 1.D-12) cn(is,ij,ik) = 1.D-12
          endif

C  Note: w/ new AER chemistry, the C and CN arrays are the same, except for the 
C    transported species (CN is updated in SOLVER), so update the C array here - just do for all

		c(is,ij,ik) = cn(is,ij,ik)
		cnout(is,ij,ik) = cn(is,ij,ik)

c  Now load in CNANOUT array, for i ge 15, move values back 1 day and load current day into index 15

        if (i .le. 14) cnanout(i+1,is,ij,ik) = c(is,ij,ik)/m(ij,ik)

        if (i .ge. 15) then 
              do 986 idayk=1,14 
 986              cnanout(idayk,is,ij,ik) = cnanout(idayk+1,is,ij,ik) 

                  cnanout(15,is,ij,ik) = c(is,ij,ik)/m(ij,ik)
        endif


          cnan(is,ij,ik) = 0.

        
C   w/ Linux INTEL FORTRAN, use ISNAN function  -  if c(is,ij,ik) in NAN, then ISNAN is true

          if (ISNAN( c(is,ij,ik) ) ) then
              cnan(is,ij,ik) = 1.
              write(6,866) is, ij, ik, c(is,ij,ik), i, iday360
866     format(1x,' is ij ik =',3i6,4x,'c = ',1PD12.3,1x,'i(day) =', I5,
     >              5x, 'iday360 = ', I3)
          endif

       if (ij .ge. 2) then
       if (cnan(is,ij,ik) .eq. 1.  .and.  cnan(is,ij-1,ik) .eq. 1.) then 
           write (32) L$, Z$
           write (32) LAT4, zalt90, pres90, zalte
           write (32) cnanout
           goto 8888
       endif
       endif

950	CONTINUE

C
C         also load in current TBAR, UBAR, HEAT, COOL, DRAG (coupled only) - CNOUT(S$+8,L$,Z$)
C                        in m/sec, K/day, m/sec/day,  1.e11 (1/m-s)
	do 762 ik=1,Z$
	do 762 ij=1,L$
           cnout(S$+1,ij,ik) = TEMP(ij,ik)
           cnout(S$+2,ij,ik) = UBARD(ij,ik)
           cnout(S$+3,ij,ik) = HEATOUT(1,ij,ik)
           cnout(S$+4,ij,ik) = COOLOUT(ij,ik)

           cnout(S$+5,ij,ik) = DRAGOUT(1,ij,ik)
           cnout(S$+6,ij,ik) = DRAGOUT(2,ij,ik)
           cnout(S$+7,ij,ik) = DRAGOUT(3,ij,ik)
           cnout(S$+8,ij,ik) = DRAGOUT(4,ij,ik)
!           ! load Ray fric into CN(32), Fixed MOD DELF into CN(58)
           cnout(32,ij,ik) = DRAGOUT(5,ij,ik)
           cnout(58,ij,ik) = DRAGOUT(6,ij,ik)
762	continue

C  load in latent heating HEATOUT(2,ij,ik) into cnout(59) - washout freq - use levels 21-46

	do 763 ik=1,26
	do 763 ij=1,L$
 763       cnout(59,ij,ik+20) = HEATOUT(2,ij,ik)


c  also update constituents for new c58 above and below 90 km and check for NaN's, load in 
c     c58 array below level Z$ just for cosmetic purposes and for output

	do 970 ik=1,Z$X
	do 970 ij=1,L$
	do 970 is=1,ITRANS

           IF (cn58(is,ij,ik) .LT. 1.D-12) cn58(is,ij,ik) = 1.D-12

           if (ik .le. Z$) then
               c58(is,ij,ik) = cn(inttra(is),ij,ik)
               cn58(is,ij,ik) = cn(inttra(is),ij,ik)
           endif
             if (ik .ge. Z$+1) c58(is,ij,ik) = cn58(is,ij,ik) 

           cnan58(is,ij,ik) = 0.

           if (ISNAN( c58(is,ij,ik) ) ) then
              cnan58(is,ij,ik) = 1.
              if (ik .ge. Z$+1)                                         ! only write out upper levels here
     >           write(6,866) is, ij, ik, c58(is,ij,ik), i, iday360
           endif

           if (ij .ge. 2) then
          if (cnan58(is,ij,ik) .eq. 1. .and. cnan58(is,ij-1,ik) .eq. 1.) 
     >             goto 8888
           endif

c                                                                    load in mixing ratio array for output
           mr58(is,ij,ik) = cn58(is,ij,ik)/m(ij,ik)


c   load in CH4 tracer for output to fort.34, mr58o(L$,Z$X)
           
           mr58o(ij,ik) = cn58(41,ij,ik)
970	continue


C  load in SF6 age and Simple age cn(77,78) and C-14 (CN(21)) into agecn array for output, AGECN(3,L$,Z$)
C    write out MIXING RATIO to fort.47 every 10 days for entire time series

	do 971 ik=1,Z$
	do 971 ij=1,L$
	  agecn(1,ij,ik) = cn(77,ij,ik)/m(ij,ik)
          agecn(2,ij,ik) = cn(78,ij,ik)/m(ij,ik)
 971      agecn(3,ij,ik) = cn(21,ij,ik)/m(ij,ik)/4.82E-18

        if (mod(iday360-5,10) .eq. 0.0) write(547) agecn

C
C   write out to fort.33 every 4 months through 1991 (iyr = 56), and then 
C      every 10 days at day 5,15, 25, 35, ... starting on year 57 (1992)
C        and write out total ozone for the year on the final day of each year 
C      NOW WRITING OUT REAL*4 CNOUT array,   CNOUT(S$+8,L$,Z$)
C
         if (iyr .le. 56) then
            if (mod(iday360-75,90) .eq. 0.0) then
ccccc                 write(33) cnout

cc                  write(19) ctprod
cc                  write(19) ctloss
cc                  write(19) chmpox
cc                  write(19) chmlox
cc                  write(19) rain1 
cc                  write(19) ice1 
cc
cc                  write(22) diffy
cc                  write(22) diffz 
cc                  write(22) advy
cc                  write(22) advz 
           end if
        end if

         if (iyr .ge. 57) then
            if (mod(iday360-5,10) .eq. 0.0) then
ccccccc                 write(33) cnout
cc                 write(19) ctprod
cc                 write(19) ctloss
cc                 write(19) chmpox
cc                 write(19) chmlox
cc                 write(19) rain1 
cc                 write(19) ice1 
cc
cc                 write(22) diffy
cc                 write(22) diffz 
cc                 write(22) advy
cc                 write(22) advz 
           end if
         end if



cccccc             write(33) cnout


          if (mod(iday360-15,30) .eq. 0.0) then 
                 write (33) cnout
                 write (33) W
                 write (33) V
                 write (33) EKYY
                 write (33) EKYZ
                 write (33) EKZZ

                 write (222) SORADHEAT
                 write (222) PHEAT
                 write (222) HEATOUT
                 write (222) COOLOUT
           endif


c write out CH4 tracer every 10 days

          if (mod(iday360-5,10) .eq. 0.0) write (34) mr58o, m 


ccccc         if (mod(iday360-15,30) .eq. 0.0) write(33) cnout, cdn


c            output every 10 days on 9th and 10th year (years 43-44) to unit 34

cccc       if (iyr .eq. 43  .or.  iyr .eq. 44) then
cccc         if (mod(iday360-5,10) .eq. 0.0) write(34) cnout
cccc       end if



         if (iday360 .eq. 360) WRITE(14) COLO3


C   And write out diffusive and advective flux terms to UNIT 22;  DIFFY, DIFFZ, ADVY, ADVZ(T$,L$,Z$X)
C        every 10 days at day 5,15, 25, 35, ... starting on year 57 (1992)

cef          if (iyr .ge. 57) then
cef               if (mod(iday360-5,10) .eq. 0.0) then
cef                    write(22) diffy
cef                    write(22) diffz 
cef                    write(22) advy
cef                    write(22) advz 
cef                end if
cef            end if



ccccc        if (iday360 .eq. 360) then
ccccc        if (mod(iyr-5,10) .eq. 0.0) then

C                                            write out to fort.17 every month for 1980 and 1995 ONLY
        if (iyr .eq. 45  .or.  iyr .eq. 60) then
ccccc        if (iyr .eq. 35  .or.  iyr .eq. 38) then
        if (mod(iday360-15,30) .eq. 0.0) then
cc             write (17) dtimeb
cc             write (17) cndc, m
cc             write (17) j

ccc             write (17) aerosol, nataer, aerice
ccc             write (17) kh
                                        ! khg(21,L$,Z$), ggsfc(21,L$,Z$), KHETA(NKR,18,46), GAMAER(NKR,18,46)
ccccc             write (17) gfaerout
ccccc             write (17) khaer
                                            ! hno3p1(7,L$,Z$), hno3l1(3,L$,Z$) 
cccccc             write (17) hno3p1, hno3l1

cccccccc             write(33) cnout

ccccc             write (17) noonrat
ccccc              write (17) j, jdc
        endif
        endif


c  OUTPUT in DUMP on DAY 360 every 5 years, in case the run crashes 
C    (or computers go down, or something......)

ccc         if (iyr .ge. 25  .and.  iday360 .eq. 360) then
         if (iday360 .eq. 360) then
             if (mod(iyr-5,5) .eq. 0.0) CALL DUMP
         endif

! BThomas Nov.2017 - Generating daily output of everything that might be useful for me:
!  This output includes chemical species identified in 2Dmodel_ConstituentList_constit15_9ua.pdf
!       along with index values for date, lat, alt and actual date, lat and alt values
!       as well as some other data (such as temperature, see below).
	 print *, "Writing output to files at timestep i = ",i
	 do ik=1,Z$
	    write(5130,'(1x,1i10)') ik
	    do ij=1,L$
	       write(5130,'(1x,1i10)') ij
!	       print *,      
!     >               i,iday360,ij,lat4(ij),ik,zalt90(ik),zkm(ij,ik)
               !load in constituent (chemical species) data:
	       do is=1,outS
		  outdatarr(is,ij,ik) = cn(is,ij,ik)
	       end do
               !load in some other data:
	       ! BThomas Nov2017 - I have some uncertainty about the following.
               !  I THINK they are what I've labled them as...
	       outdatarr(outS+1,ij,ik) = TEMP(ij,ik) !temperature, K (confident)
	       outdatarr(outS+2,ij,ik) = UBARD(ij,ik) !mean zonal wind speed, m/s (confident) 
	       outdatarr(outS+3,ij,ik) = HEATOUT(1,ij,ik) !heating rate? K/day
	       outdatarr(outS+4,ij,ik) = COOLOUT(ij,ik) !cooling rate? K/day
	       outdatarr(outS+5,ij,ik) = W(ij,ik) !vertical velocity cm/s (pretty confident) 
	       outdatarr(outS+6,ij,ik) = V(ij,ik) !horizontal velocity (north-south?) cm/s (pretty confident)
	       !Write out to a file, all in one long column.  
               !  Will process into a more usable form using an external NCL routine.
	       do is=1,outnum
		  write(5130,'(1x,1e16.8)') outdatarr(is,ij,ik)
	       end do !is loop
!
	    end do !ij loop
	 end do !ik loop
!
!        !Write out column density values:
         ! in com2d.h, NCOLGD(10,L$,Z$); See "colden" subroutine for more info (BThomasNov2017)
         ! the number of items to write out in format statement is nsp, see above
	 !Also output here (to separate file) is ground-level column density
         ! of HNO3 rainout.
	 !Note that ik=1 here means ground level
	 do ij=1,L$
	    write(5131,'(1x,5e16.8)') 
     >            (ncolgd(is,ij,1),is=1,nsp)
	    !print *, i,ij, "coldenhno3(ij,1): ",coldenhno3(ij,1)
	    write(5132,'(1x,1e16.8)') coldenhno3(ij,1)
	 end do
	 !output just the O3 and NO2 concentrations (alt-lat dependent), with index/grid info
	 ! O3 is cn index 4, NO2 index 6
	 ! this file for use with TUV model
	 do ij=1,L$
	    do ik=1,Z$
	       write(5133,1960) i, iday360, ij, lat4(ij), zkm(ij,ik),
     *                      cn(4,ij,ik), cn(6,ij,ik)
 1960	       format(1i7,1i5,1i5,1f7.1,1f10.4,2e16.8)
	    end do
	 end do
	 
100	CONTINUE !End of time loop

	close(5130) !BThomas Nov.2017 - close my daily output files
	close(5131) 
	close(5132) 
	close(5133) 

C        print *,' before call dump'
            CALL DUMP
C        print *,' after call dump'

c  also write out total ozone for Jan-April 2000
ccccccccccccccccccccccccccccccccccccc       WRITE(14)COLO3


         IF(DAY360 .GT. 90.)THEN
C Get printout of lifetime at the end here only if the model has run
C for at least 90 days.
           YSLFPRNT=.TRUE.
ccc9HB           CALL LIFETIME
         ENDIF

 8888   extime = etime(cputime)/3600.
       write(6,6741) stime, cputime(1)/3600., cputime(2)/3600., extime
6741   format(1x, 'START-UP CPU =', F7.3, 7X, 'USER CPU =', F7.3, 4X,
     >            'SYS CPU =', F7.3, 10X, 'TOTAL CPU (hrs) =', F7.3)

	STOP
	END
