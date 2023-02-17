      SUBROUTINE GETMOMX

c     compute eddy momentum inputs due to gravity waves and
c     planetary waves. The kyy and kzz fields are also generated here
C
!       USE degree_trig

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONT.INC'
      INCLUDE 'COMMONW.INC'
      
ccelf      include 'test_value.inc'
ccelf      logical nantest1, nantest2
      
        REAL KYYMAX, KZZMAX, GWDRG(NP$,MP$), epexx(NP$,MP$)
        REAL*8 GWDRGS0(NP$,MP$), GWDRGS1(NP$,MP$)

        REAL fks(NP$,MP$), fkq(NP$,MP$), fkq2(NP$,MP$), frq(NP$,MP$)
        REAL fegwc(NP$,MP$), fktot(NP$,MP$)

        REAL pzz(MP$), EPXE(NP$,MP$), EPXW(NP$,MP$), OGG(NP$,MP$)


        COMMON/CKYY1/ BGKYY0(N$,M$), BGKYY(N$,M$), EPOTX(NP$,MP$),
     >                XKYYFX(N$,M$)
        

c common of NCEP UBAR in m/sec, and coupled model - NCEP UBAR difference for current day (m/sec)

        COMMON/CUNCEP/UNCEP(NP$,MP$), UDIFFDAY(NP$,MP$)


C  common for fixed model DELFs for current day (m/sec/day); EPBARO - from barclinic eddies only (m/sec^2)

        COMMON/CEPDD/ EPDD(N$,M$), EPDDX(NP$,MP$), EPBARO(NP$,MP$)


C   common for PW DRAG computed from NCEP UBAR for coupled model (from XCOUPLED)
C     also store total daily drag in array EPXT

cccc        COMMON/CDELF/ IYRCT0
        COMMON/CEPXT/ EPXT(NP$,MP$,360)


C  common for model compute eddy heating/EP flux correction: 
C      EHFC is in K/sec (theta heating); EPZZ is in cm/sec^2

       COMMON/CEHFC/ ehfc(NP$,MP$), epzz(NP$,MP$)


C  COMMON for WACCM combined eddy heating rate (theta heating) in K/sec
C
cccwconv48       COMMON/CQWTD/qwttd(NP$,MP$)


C   COMMON for WACCM oro Gwaves DELF for current day (m/sec^2) from TEMPIN/TWODS - WCONV97

       COMMON/CWOGWD/ OGWD(NP$,MP$)


C  common for year counter, steady state

       COMMON/CIYRSS/IYRCT, ISTST



C     Added 9/18/96 by PEM to allow access to new timing (timer.f).

      include "timer.inc"


      LOGICAL PWIN_FLAG
      DATA    PWIN_FLAG /.FALSE./

C     End Added 9/18/96.

        
C------------GRAVITY WAVE SECTION---------------------------

c       gravity wave scheme turned on after nswd DAYS

C       nswd=30
        nswd=isw(12)
C        nswitch=nswd*dayl/dt
        nswitch = nswd
        ralx=20.
        

C  GRAVITY WAVE PARAMETERIZATIONS

C Modified 9/18/96 by PEM to use new timing (timer.f).  Note that
C  since ISW(12), and hence NSWITCH, is currently -10 by default, and
C  IGRAVW is set to 0, this section is not used.  The gravity wave
C  turn-on day (ISW(12)) is now in terms of the 360-day idealized calendar.

 
C Modified 9/18/95 by PEM for new JTB mountain waves:
C
C
C   EF - OCt 2009, if computing OGW - DON'T CALL JMTNWV (just Julio's cubic drag law)
C      and initialize arrays GWDRG and GWKZZ below (these are loaded in JMNTWV)
C      then substitute in WACCM orographic GWAVE drag below
C

ccwconv97           CALL JMTNWV(GWDRG)               !Get Mountain Waves


           CALL GWDRAG(1)                   !Get Grav. Waves
           
           DO K=1, MP$                      !Add In Grav. Waves
           DO J=1, NP$

               GWDRG(J,K) = 0.0
               GWKZZ(J,K) = 0.0
           
           GWDRG(J,K) = GWDRG(J,K) + DUM1X(J,K)
           GWKZZ(J,K) = GWKZZ(J,K) + DUM2X(J,K)


c   nantest1=test_nan( gwdrg(j,k) )
c   nantest2=test_nan( gwkzz(j,k) )
c   if(nantest1 .or. nantest2) then
c      PRINT*, "NAN for GWDRG/GWKZZ in GETMOMX:"
c      PRINT*, "   J = ", J, ", K = ", K
c      PRINT*, "   GWDRG = ", GWDRG(J,K)
c      PRINT*, "   GWKZZ = ", GWKZZ(J,K)
c      STOP
c   ENDIF


C  load in REAL*8 array for SMOOTH5

           GWDRGS0(J,K) = GWDRG(J,K)

           END DO
           END DO 

c  WRITE(*,*) '  USED NEW JTB MTN WAVES '

C-----------------------------------------------------------------

C 
C  for 9GC run, smooth Julio's gravity wave drag array
C   REAL*4 GWDRG(NP$,MP$), use SMOOTH5;  REAL*8 GWDRGS0(NP$,MP$), GWDRGS1(NP$,MP$)


           CALL SMOOTH5(GWDRGS0, GWDRGS1, NP$, MP$)


           DO K=1, MP$
           DO J=1, NP$
              GWDRG(J,K) = GWDRGS1(J,K)
           END DO
           END DO
 



C  WACCM orog. gravity waves - OGWD(NP$,MP$) is in m/sec2, convert to cm/sec2 - OGG(NP$,MP$)
C    also check against model UBAR - UBX(NP$,MP$) in cm/sec 
C    since this drags UBAR back to 0 :
C
C    if UBAR > 0 m/sec, OGW must be < 0
C    if UBAR < 0 m/sec, OGW must be > 0
C
                                ! ulim in cm/sec;   gmax, glim in cm/sec^2
cogg           ulim = 0.*100.
           glim = -.5*100/86400.
           gmax = -2.*100/86400.

           DO K=1, MP$
           DO J=1, NP$
              OGG(J,K) = OGWD(J,K)*100.    !/2.

              if (UBX(J,K) .ge. 0. .and. OGG(J,K) .ge. 0.) OGG(J,K) = 0.
              if (UBX(J,K) .le. 0. .and. OGG(J,K) .le. 0.) OGG(J,K) = 0.

cogg    if (UBX(J,K) .ge. ulim) then
cogg       if (OGG(J,K) .gt. 0.) OGG(J,K) = 0.
cogg    endif

cogg    if (UBX(J,K) .le. -ulim) then
cogg       if (OGG(J,K) .lt. 0.) OGG(J,K) = 0.
cogg    endif

C                ! reduce OGG by 4X if < 0.
cogg    if (OGG(J,K) .lt. 0.) OGG(J,K) = OGG(J,K)/4.

C
C  set limits on OGG:  reduce by 4X if LT -.5 m/sec/day, but keep it at least -.5 m/sec/day
C   then limit to -2 m/sec/day
C
              if (OGG(J,K) .lt. glim) then
                  gg0 = OGG(J,K)/4.
                  gg1 = AMIN1(glim, gg0)
                  OGG(J,K) = AMAX1(gmax, gg1)
              endif

CC  wconv99 - zero out WACCM OGW

              OGG(J,K) = 0.0

           END DO
           END DO 



CC
C  ****************  START   EXTRA MOMENTUM  LOADING   *************
C
C
C  compute relaxation to NCEP UBAR for current day - IYRCT, IDOY360 (in timer.inc) (both in COMMON)
C  UNCEP(NP$,MP$) is in m/sec;    UBX(NP$,MP$) in cm/sec    - EPXT(NP$,MP$,360) is in cm/sec^2
C  for TD run compute ONLY for the 1st 3 years, then reuse the daily values after the 3rd year)
C  for SSTate run, compute for all years
C     in 1st year, IYRCT = 1;   also set up z mask, pzz(MP$) for westerly accel, initialize
C
       IF (IYRCT .eq. 1  .and.  IDOY360 .le. 2) then
           do 770 iit=1,360
           do 770 K=1,MP$
           do 770 J=1,NP$
 770          EPXT(J,K,iit) = 0.0

           do 771 K=1,MP$
              pzz(K) = (13.-zpp(K))/5.
              if (pzz(K) .le. 0.) pzz(K) = 0.
              if (pzz(K) .ge. 1.) pzz(K) = 1.
 771       CONTINUE
       ENDIF


C  set latitude dependent timescale, pyy is in days (equat=60, poles=5-10 max)
C  larger in NH, seasonal dependence in SH (only for OCt 15-Dec)
C    pyy = 60. - yhem*( (TANH(ypp(j)*cvt*2.))**6 ) - MAX at poles

       cvt = 3.1415927/180.

       IF (ISTST .eq. 0  .and.  IYRCT .le. 3  .or.  ISTST .eq. 1) then
          DO 775 K=1,MP$
          DO 775 J=1,NP$
cpw86              yhem = 59.   ! 57.   ! 55.
cpw86
cpw86              if (ypp(j) .lt. 0.) then
cpw86                 yysh = 0.
cpw86
cpw86                   if (idoy360 .ge. 285) then
cpw86                      yysh = (idoy360 - 285.)/30.
cpw86                      if (idoy360 .ge. 345) yysh = (360. - idoy360)/15.
cpw86                      if (yysh .ge. 1.) yysh = 1.
cpw86                      if (yysh .le. 0.) yysh = 0.
cpw86                   endif
cpw86
cpw86                 yhem = -50. + yysh*109.
cpw86              endif
C                                      gives max at 60S,N
cc60max              pyy = 60. - yhem*(.5-(COSD((ypp(j)*3.))/2.))

c      gives max at 70S,N, small at low lats

cpw86              pyf = (.5-(COSD(ypp(j)*2.5))/2.)**2
cpw86              pyy = 60. - yhem*pyf
cpw86
cpw86              if (pyy .le. 1.) pyy = 1.
cpw86              if (pyy .ge. 120.) pyy = 120.


C  if westerly accel, only use below 13 km and set to 30-day time scale everywhere
cpw86
cpw86              if (exmom .gt. 0.) then
cpw86                   exmom = -(UBX(J,K) - UNCEP(J,K)*100.)/(DAYL*30.)
cpw86                    EPXT(J,K,idoy360) = exmom*pzz(K)
cpw86              endif



C  EPXT is in cm/sec^2, reverse sign (it's a drag term)
c  WCONV189 - set THIS to 0 EVERYWHERE


cccc    EPXT(J,K,idoy360) = 0.  ! -(UBX(J,K) - UNCEP(J,K)*100.)/(DAYL*360.)


C
cpw86   if (exmom .le. 0.  .or.  epncep .lt. 0.) then
cpw86          exmom0 = exmom
cpw86          exmom = AMIN1(exmom0, epncep)   ! take larger of exmom, 25% NCEP DELF
cpw86
cpw86          exmom1 = exmom
cpw86          ulim = -20. *100.               ! ulim in cm/sec
cpw86          if (UBX(J,K) .lt. 0.) then
cpw86            exmom = exmom1/(-ulim)*( UBX(J,K) - ulim)
cpw86            if (UBX(J,K) .lt. ulim) exmom = 0.
cpw86          endif
cpw86
cpw86   endif

 775     CONTINUE
       ENDIF

C  reset to 0 for the first 30 days to allow for spin up

cpw86       if (IYRCT .eq. 1  .and.  IDOY360 .le. 30) then
cpw86           do 776 K=1,MP$
cpw86           do 776 J=1,NP$
cpw86 776          EPXT(J,K,idoy360) = 0.0
cpw86       endif


C  load in strong lower trop relaxation (< 3.5 km) to NCEP UBAR here - time scale of 24 hours
C  need to do this EVERY DAY (not just in first 3 years) 
C  maybe need to limit to <.009 cm/sec^2 (7.75476 m/sec/day) ?
C        UBX(NP$,MP$) in cm/sec ;  UNCEP(NP$,MP$)
C
CCWCONV91 - do relaxation at 1-3 km easterly accel only in tropics:
CCWCONV91   add additional accel for 40-80N, reduce > 80N

       do 757 K=1,MP$
          if (ZPP(K) .le. 3.5) then

          do 776 J=1,NP$
              sfcr = -(UBX(J,K) - UNCEP(J,K)*100.)/(DAYL*1.)
              if (ypp(j) .ge. 40.  .and.  ypp(J) .le. 80.) then
                  sfcr1 = sfcr
                  sfcr = sfcr1*(1. + 3.*((COSD(ypp(j)*4.9-90.))**2))
              endif
              if (ypp(j) .gt. 80.) sfcr = sfcr/2.

              EPXT(J,K,idoy360) = sfcr
       if (sfcr .lt. 0. .and. ABS(ypp(j)) .gt. 25.) EPXT(J,K,idoy360)=0.
 776     CONTINUE

         endif
 757    CONTINUE



C  load in extra momentum for current day here, separate easterly/westerly
C     EPXE(NP$,MP$), EPXW(NP$,MP$) are in cm/sec^2  -  first initialize
C
         DO 777 K=1,MP$
         DO 777 J=1,NP$
             EPXE(J,K) = 0.0
             EPXW(J,K) = 0.0

             relaxu = EPXT(J,K,idoy360)
c
c  Use up to 35% NCEP DELFS at high lats (current day), ramp down to zero above 60 km and below 15 km
c    EPDDX(NP$,MP$) are the NCEP DELFs for current day in m/sec/day - Convert to cm/sec^2
C    weight by SIN(lat)^6 , zero-out equatorward of 30 deg and IN SH AND IN NH
C
             pwmx = .35
             pwmn = .00   ! .05

             zzpw = pwmx - (ZPP(K)-50.)/35.
             if (zzpw .lt. pwmn) zzpw = pwmn

             if (ZPP(K) .le. 23.) then
               zzpw = pwmx - (23.-ZPP(K))/25.
               if (zzpw .lt. 0.) zzpw = 0.
             endif

             if (zzpw .gt. pwmx) zzpw = pwmx
CC                                              ! zero out in tropics and in SH - wconv34
             pwy = (SIND(ypp(j)*1.23))**6
             if (ypp(j) .le. 30.) pwy = 0.
cccccccc             if (ABS(ypp(j)) .le. 30.) pwy = 0.

             obdelf = EPDDX(J,K)*100./86400.
             epncep = 0.  ! zzpw*pwy*obdelf  ! *0.

c  do full NCEP for easterlies, ramp down below 23 km, limit to -.5 m/s/d everywhere (not above 70 km)
c  even though this keys on model UBAR, do it time dependently since it's accounting 
c  for momentum sources in easterlies not computed in model, and include in Kyy
C
C  WCONV219- remove ALL NCEP DELF
C
             if (UNCEP(J,K) .le. 0.  .and.  UBX(J,K) .le. 0.) then
                zzmm = 1. - (23.-ZPP(K))/7.5
                if (zzmm .lt. 0.) zzmm = 0.
                if (zzmm .gt. 1.) zzmm = 1.

                pwy = (SIND(ypp(j)*1.8))**6   ! 1. - (COSD(ypp(j)))**4
                if (ABS(ypp(j)) .ge. 50.) pwy = 1. 
                epncep = 0.   ! zzmm*pwy*obdelf

                dfx = -.4*100./86400.
CWCONV219  if (epncep .le. dfx) epncep = dfx

ccwconv35      dfx = -4.*100./86400.
ccwconv35      if (zpp(k) .gt. 70. .and.  epncep .le. dfx) epncep = dfx
             endif

C
C  limit TOTAL DELF to Magnitude LE NCEP DELF - EPDDX(NP$,MP$) (m/sec/day)
C
C  use EPFLX(NP$,MP$) (cm/sec^2) from previous time step, which is OK since EPXT array is 
C  from the final diurnal time step, and the changes for 2 hour time step are very small 

             tdelf = EPFLX(J,K) + epncep

             if (tdelf .lt. obdelf) epncep = obdelf - EPFLX(J,K)
             if (EPFLX(J,K) .lt. obdelf) epncep = 0.         ! do this last


c  load into EPEXX(NP$,MP$), EPXE(NP$,MP$), EPXW(NP$,MP$) for Kyy/output
C
              EPEXX(J,K) = epncep

              if (relaxu .le. 0) EPXE(J,K) = relaxu
              if (relaxu .gt. 0) EPXW(J,K) = relaxu
 777    CONTINUE


c       if (IYRCT .eq. 1  .and.  idoy360 .eq. 1) then 
c       if (mod(idoy360-15, 30) .eq. 0.0) then 
c          write (798) N$, M$, iyrct, idoy360
c          write (798) yp, zp, ypp, zpp, pzz
c          write (798) ubx, uncep, EPXE, EPXW, EPXT
c       endif

CC  ***************   END  EXTRA MOMENTUM  LOADING   *************



C
C------------- P L A N E T A R Y * W A V E S ---------------------



C   Modified 9/18/96 by PEM to be consistent with new timing (timer.f).

c read in from P_wave_number.dat

C        IF (NSTEP.EQ.1) THEN
        IF (.NOT.PWIN_FLAG) THEN        !  Read in array containing
            READ(99,9900) MWVNS         !  actual physical zonal 
            CLOSE(UNIT=99)              !  planetary wave # and 
            write(6,9900) MWVNS         !  other stuff
            PWIN_FLAG = .TRUE.
            ENDIF
 9900   FORMAT( 4(I5,1X) )

C  End Modified 9/18/96.

C isw(46)=1 ,   ISW(60)=1, so   nswitch2 = 2
        nswd   =isw(46)
        NOPWDRG=ISW(60)
C        nswitch = nswd*dayl/dt          !  Turn P-waves on at day nswd
C        nswitch2=(nswd+nopwdrg)*dayl/dt !  Begin to affect Flow NOPWDRG days later

        nswitch  = nswd
        nswitch2 = nswd + nopwdrg


C Modified 9/18/96 by PEM to use new timing (timer.f).  The planetary
C wave calculation turn-on day (NSWITCH) and the turn-on day for
C wave dissipation to begin affected the flow (NSWITCH2) are now in
C terms of the 360-day idealized calendar.


C IF ( (nstep.ge.nswitch) .and. (nswitch .ge. 0.) ) THEN


        IF ((IDOR360.GE.NSWITCH) .AND. (NSWITCH.GE.0)) THEN

C End Modified 9/18/96.


C  include day360 day number and IYR (IYEARC in timer.f) here - EF, Nov. 2007
C  also include EPEXX(NP$,MP$) - NCEP DELF for Kyy
C
         CALL JPWAVE(IDOY360, IYEARC, EPEXX, EPXE)


        ENDIF



4321  FORMAT( ' SRF.DRAG TIME(DAYS) (U and PSI) ',I3,
     2        ' ADDED EP FLUX DIV TO DRAG ',I3 )


                            ! NEPDV and ISW(32) = -1
      NSRFD=ISW(29)*ISW(35)
      NEPDV=ISW(32)
      
      
C  Modified 9/18/96 by PEM to use new timing (timer.f).

C      IF(NSTEP .LT. NSWITCH2) NEPDV=0
C      IF((NSTEP .EQ. NSWITCH2).AND.(NEPDV.NE.0)) THEN
C        WRITE(*,*)'  *********************************** '
C        WRITE(*,*)'  --------  TURNING EP DIVERGENCE -- '
C        WRITE(*,*)'  --------  ON IN UBAR EQ. --------- '
C        WRITE(*,*)'  *********************************** '
C      ENDIF

C IDOR360 is Day-of-run, Initialized to 1 but not reset at the start of a new model year, nswitch2= 2
C
       IF (IDOR360 .LT. NSWITCH2) NEPDV = 0
       
C  End Modified 9/18/96.


c  WRITE(6,4321) NSRFD, NEPDV


C
C   get tropical momentum sources  - EF (May 2009), units are in cm/sec^2

cckw       CALL KWAVE_COUP(NP$, MP$, yp, ypp, zp, zpp, t, ubx,
cckw     >                 FKS, FKQ, FKQ2, FRQ, FEGWC)



      DO K=1,MP$
      DO J=1,NP$


C  sum up tropical wave sources here, apply factors - fktot(NP$,MP$), in cm/sec^2
C     if EQ GWs NOT USED, fegwc=0 ;   don't use fkq2??

cckw         fktot(j,k) = fks(j,k) + fkq(j,k) + fkq2(j,k) 
cckw     >              + frq(j,k) + fegwc(j,k)

           fktot(j,k) = 0.0



C  Set up Rayleigh friction - RFRIC(NP$) (defined in TWODS, COMMONC.INC)
C  reduce:  below 30 km if westerly, below 60 km if easterly, UBX(NP$,MP$)

         if (UBX(j,k) .ge. 0.) then
            rfad = 1. - (30. - zpp(K))/20.
            if (rfad .lt. .2) rfad = .2
            if (rfad .gt. 1.) rfad = 1.

            RALCO = RALP(K)*RFRIC(J)*rfad
         endif

         if (UBX(j,k) .lt. 0.) then
            rfad = 1. - (60. - zpp(K))/20.
            if (rfad .lt. .2) rfad = .2
            if (rfad .gt. 1.) rfad = 1.

            RALCO = RALP(K)*rfad
         endif


Ccexm add in winter/spring easterly GW momentum source, weight by sin(lat)^6, max at high lats
CCcexm    zero-out w/ adding fixed model DELFs, and zero-out < 40 degrees
CCcexm 
Ccexm         seasrf = 0.
Ccexm         ampm = 0.
Ccexm
Ccexm         if (YPP(J) .lt. 0.) then
Ccexm           ampm = 0.   ! .25   ! .5
Ccexm
Ccexm           if (idoy360 .ge. 270.  .and.  idoy360 .le. 300.)
Ccexm     >         seasrf = SIND((idoy360-270.)*3.)
Ccexm
Ccexm           if (idoy360 .ge. 300.  .and.  idoy360 .le. 330.) seasrf = 1.
Ccexm
Ccexm           if (idoy360 .ge. 330.)    !  .and.  idoy360 .le. 90.)
Ccexm     >         seasrf = SIND((idoy360-60.)*3.)
Ccexm         endif
Ccexm
Ccexm
Ccexm         if (YPP(J) .ge. 0.) then
Ccexm             ampm =  0.   ! .75   ! 1.
Ccexm
Ccexm             if (idoy360 .ge. 315.) seasrf = SIND((idoy360-315.)*2.)
Ccexm
Ccexm             if (idoy360 .le. 105.) seasrf = 1.
Ccexm 
Ccexm             if (idoy360 .ge. 105.  .and.  idoy360 .le. 135.)
Ccexm     >          seasrf = SIND((idoy360-75.)*3.)
Ccexm          endif
Ccexm
Ccexm          seasrf = seasrf * ((SIND(ypp(j)))**6)        ! ((SIND(ypp(j)))**2)
Ccexm          if (ABS(ypp(j)) .le. 40.) seasrf = 0.
Ccexm
CcexmC  ZRFAC  - alt weighting - ZPP(MP$)
CcexmC
Ccexm         zrfac = (ZPP(k) - 15.)/7.5
Ccexm         if (zrfac .le. 0.) zrfac = 0.
Ccexm         if (zrfac .ge. 1.) zrfac = 1.
Ccexm

         ABSLAT = ABS(YPP(J))

C Added 10/26/95 by PEM to zero-out equat rayl damping:
C    IF (ABSLAT.LE.15.0) THEN
C         RALCO = 0.0
C       ELSEIF (ABSLAT.LE.30.0) THEN
C         RALCO = ((ABSLAT - 15.0)/(15.0))*RALP(K)
C       ELSE
C         RALCO=RALP(K)
C    ENDIF

C    IF (ABSLAT.LE.25) THEN
C         RALCO = EXP(-1.0*((ABSLAT - 25.0)/12.5)**2.0)*RALP(K)
C      ELSE
C         RALCO = RALP(K)
C    ENDIF
C        End added 10/25/95, PEM.

C
C  replace  EPFLX(NP$,MP$) w/ EPOTX(NP$,MP$) if updated w/ KYYOT where applies
c  DON'T DO THIS - get anomlous EASTERLIES in tropical upper troposphere - EF
c
c   if (zpp(k) .ge. 7.  .and.  zpp(k) .le. 19.) then
c    if (ABS(ypp(j)) .le. 50.) then
c     if (EPOTX(j,k) .lt. EPFLX(j,k)) EPFLX(j,k) = EPOTX(j,k)
c     endif
c    endif


C
C  do WEAK relaxation of model UBAR to NCEP climatology for current day to get better results
C    use 30 day timescale;  UBX(NP$,MP$) in cm/sec;  UNCEP(NP$,MP$) is m/sec
C
C  also add in extra easterly momentum source of AMPM (m/sec/day), convert to cm/sec^2 
C   weight by seasrf and altitude, ZRFAC - this contribution SHOULD BE POSITIVE
C   NO CONTRIBUTION when using FIXED MODEL DELF
C
C
Cccpw4 - use  pre-determined NCEP-model diffs here to avoid long term contribution
C        UDIFFDAY(NP$,MP$) is in m/sec convert to cm/sec - use 30 day timescale
C
         relaxu = 0.   ! UDIFFDAY(J,K)*100./(DAYL*30.)
Cccpw7     >          + ampm*100./(DAYL*1.) * seasrf * zrfac

Cccpw4         relaxu = (UBX(J,K) - UNCEP(J,K)*100.)/(DAYL*30.)
Cccpw4     >          + ampm*100./(DAYL*1.) * seasrf * zrfac


CCpw23
CCpw23
CCpw23  Include contribution from NCEP DELFs when model UBAR is easterly:
CCpw23    UBX(NP$,MP$) in cm/sec;  UNCEP(NP$,MP$) is m/sec;   
ccpw23    EPDDX(NP$,MP$) are the fixed model DELFs in m/sec/day - Convert to cm/sec^2
CCpw23       don't do in extratropical troposphere

CCpw23        if (UBX(J,K) .lt. 0.  .and.  UNCEP(J,K) .lt. 0.) 
CCpw23         epu = 0.
CCpw23         if (UNCEP(J,K) .lt. 0.) epu = .3*EPDDX(J,K)*100./86400.
CCpw23         IF (ABS(YPP(J)) .ge. 40. .and.  ZPP(K) .LE. 10.) epu = 0.



         DRAGX(J,K)=-RALCO*UBX(J,K)*CST(J)     !(XC(J,K,2)-XC00(J,K,2))
     1              -NEPDV*EPFLX(J,K)*CST(J)   !PWAVE DRIVING
     2              -GWDRG(J,K)*CST(J)         !GWAVE DRIVING
     3              -relaxu*CST(J)             ! relaxu = 0.0 currently
     4             + fktot(J,K)*CST(J)         ! tropical wave driving
     5             + EPXE(J,K)*CST(J)          ! EPXE(NP$,MP$) easterly relax to UNCEP, cm/sec^2
     6             + EPXW(J,K)*CST(J)          ! EPXW(NP$,MP$) westerly relax to UNCEP, cm/sec^2
     7             + EPEXX(J,K)*CST(J)         ! EPEXX(NP$,MP$) NCEP DELF, cm/sec^2
     7             + OGG(J,K)*CST(J)           ! OGG(NP$,MP$) WACCM OGWwave, cm/sec^2


C  set SURFACE Drag as ONLY relaxation to NCEP UBAR (everything else should be 0)

CCWCONV66        if (ZPP(K) .le. 3.5) DRAGX(J,K) = (EPXE(J,K) + EPXW(J,K))*CST(J)

         if (K .eq. 1) DRAGX(J,1) = (EPXE(J,1) + EPXW(J,1))*CST(J)




c  nantest1=test_nan( dragx(j,k) )
c   if(nantest1) then
c    PRINT*, "NAN for DRAGX in GETMOMX:"
c    PRINT*, "   J = ", J, ", K = ", K
c    print*, 'ralco,ubx,cst = ',ralco,ubx(j,k),cst(j)
c    print*, 'nepdv,epflx = ',nepdv,epflx(j,k)
c    print*, 'gwdrg = ',gwdrg(j,k)
c    STOP
c   ENDIF


c  also load in DRAG contributions here for output - use final time step for output, NO COSINE!
C  dragcoup(7,NP$,MP$) in COMMOND.INC; EPEXX(NP$,MP$) is NCEP DELF;  EHFC(NP$,MP$), TTOTHX(MP$)
C
          dragcoup(5,J,K) = -RALCO*UBX(J,K)
          dragcoup(2,J,K) = -NEPDV*EPFLX(J,K) + EPEXX(J,K)
          dragcoup(3,J,K) = -GWDRG(J,K)
          dragcoup(4,J,K) = EPXE(J,K) + EPXW(J,K)
          dragcoup(6,J,K) = EPEXX(J,K)



C  OLD - surface drag to zero (time scale of 12 hours; changed from 1 day), UBX(NP$,MP$)
C  gives a bit more westerly mom at surface (relax to NCEP UBAR didn't work as well)
C   load into dragcoup(3) for output;  UNCEP(NP$,MP$) is m/sec

ccpw107   IF(NSRFD.GT.0) then
ccpw107     DRAGX(J,1) = -( UBX(J,1)*CST(J) )/( NSRFD*DAYL/2.)
ccpw107     DRAGX(J,1)= -( (UBX(J,1)- UNCEP(J,1)*100.) )/(DAYL/4.)*CST(J)

ccpw107     dragcoup(3,J,1) = DRAGX(J,1)/CST(J)
ccpw107   ENDIF

c     DRAGX(J,2)=-( UBX(J,2)*CST(J) )/(10.0*DAYL)

      END DO
      END DO



      DO K=1,MP$
      DO J=1,NP$
         DUM1X(J,K)=  DRAGX(J,K)  
cc     1             -THERMCR*ISW(15)*THRMWX(J,K)*CST(J) 
c  nantest1=test_nan( dum1x(j,k) )
c  if(nantest1) then
c    PRINT*, "NAN for DUM1X in GETMOMX:"
c    PRINT*, "   J = ", J, ", K = ", K
c    print*, "dragx,dum1x",dragx(j,k),dum1x(j,k)
c    STOP
c  ENDIF
      END DO
      END DO


C
C   FINISHED CALCULATING MOMENTUM SOURCES AND/OR SINKS
C   NOW INTERPOLATE TO RADIATION GRID 
C   ---------------------------------
C
c
c  IF UB is tracer then GWMOM ~ DRAGX
c  IF M  is tracer then GWMOM ~ DRAGX/C(j); ISW(33) = 0  currently (Feb 2008)
c 
      DO K=1,M$
      DO J=1,N$
         GWMOM(J,K)= 
     >               ( DUM1X(J,K) + DUM1X(J+1,K)
     >               + DUM1X(J,K+1)+DUM1X(J+1,K+1) )
     >                     /(4.*C(J))

     >             + ( XMIX(J,K,2)  + XMIX(J+1,K,2)
     >               + XMIX(J,K+1,2)+ XMIX(J+1,K+1,2) )
     >             *ISW(33)/(4.*C(J))


c  nantest1=test_nan( gwmom(j,k) )
c  if(nantest1) then
c     PRINT*, "NAN #1 for GWMOM in GETMOMX:"
c     PRINT*, "   J = ", J, ", K = ", K
c     print*, "dum1x(j,k),dum1x(j+1,k),dum1x(j,k+1),dum1x(j+1,k+1)",
c   *     dum1x(j,k),dum1x(j+1,k),dum1x(j,k+1),dum1x(j+1,k+1)
c     print*, "c(j) = ",c(j)
c     print*, "xmix(j,k ,xmix(j+1,k ,xmix(j,k+1 ,xmix(j+1,k+1 = ",
c   *       XMIX(J,K,2),XMIX(J+1,K,2),XMIX(J,K+1,2),XMIX(J+1,K+1,2) 
c   STOP
c  ENDIF


      END DO
      END DO 

        CALL OUTA(13,NSTEP,NDAY,NYEAR,YPP,ZPP,EPFLX,NP$,MP$,0)
        CALL OUTA(22,NSTEP,NDAY,NYEAR,YP ,ZP ,GWMOM,N$ ,M$ ,0)

       
      DO K=1,M$
      DO J=1,N$
         GWMOM(J,K)=  GWMOM(J,K)
     1             -THERMCR*ISW(15)*THRMW(J,K) 

c  nantest1=test_nan( gwmom(j,k) )
c  if(nantest1) then
c    PRINT*, "NAN #2 for GWMOM in GETMOMX:"
c     PRINT*, "   J = ", J, ", K = ", K
c    STOP
c  ENDIF

      END DO
      END DO


C      CALL OUTA(22,NSTEP,NDAY,NYEAR,YPP,ZPP,DRAGX,Np$,Mp$,0)
c      CALL OUTA(50,NSTEP,NDAY,NYEAR,YP,ZP,THRMW,N$,M$,0)

      RETURN


      END



	SUBROUTINE JPWAVE(IDOY360, IYEARC, EPEXX, EPXE)

C	PROGRAM TO SOLVE THE PLANETARY WAVE EQUATION ONCE GIVEN
C	THE MEAN ZONAL WIND (UBAR) AND THE ZONALLY DECOMPOSED OROGRAPHIC HEIGHTS

	INCLUDE 'PARAM.INC'
	INCLUDE 'COMMONC.INC'
	INCLUDE 'COMMOND.INC'
	INCLUDE 'COMMONW.INC'
	INCLUDE 'COMMONT.INC'

celf      include 'test_value.inc'
celf      logical nantest

C   common for CLIMATOLOGICAL geop hgt wave (1-5)  coeffs (in METERS), read in TEMPIN
C   1st index: 1=zonal avg, 2-6 are waves 1-5;   2nd index:  1=Ak,  2=Bk;  
C   these have been PROPERLY interpolated to current dynamics grid in TWODS

        COMMON/CGPBD/ GPBC1(6,2,N$)


C  common for fixed model DELFs for current day (m/sec/day); 
C    EPBARO - from barclinic eddies only (m/sec^2)

        COMMON/CEPDD/ EPDD(N$,M$), EPDDX(NP$,MP$), EPBARO(NP$,MP$)

      
	common/slakwk3/ 
     2   phi(n$,m$),PHAZ(n$,M$),PAMP(n$,M$)
     3  ,bigq(n$,m$),XNU1(N$,m$),XNU0(N$,m$),
     4   PHICR(N$,M$),epfy(n$,m$),epfz(n$,m$),tempx1(n$,m$)

	COMPLEX AYC(n$,M$),AZC(n$,M$),cyc(n$,M$),CZC(n$,M$),
     1   phic(n$,m$),BC(n$,M$),FC(n$,M$)
     2  ,phi,bigq,ai,PHICR
     3  ,COE1,COE2,COEZ,ff,rr,ss,XNU0,XNU1,XNU0l,XNU1l
     4  ,coef,s4,s5,BCA(N$,M$),BCB(N$,M$),fg(n$)

        COMPLEX PSAO(N$,M$),PSBO(N$,M$),PSICO(N$,M$)
     2  , DUM1C(N$,MMAX$)


        REAL*8 EPFLXS0(NP$,MP$), EPFLXS1(NP$,MP$), EPFLXS2(NP$,MP$)

        real ubw1(N$,M$), epstr(N$,M$), epex(N$,M$), epexx(NP$,MP$)
        real pwyy(N$), EPXE(NP$,MP$), EPXE1(N$,M$)
        real thc(N$,M$), thce(NP$,MP$), delf10(N$,M$)
        real pwbf(n$,m$),xn(m$),deq(m$),dum1v(n$)
        REAL PSICR(N$,M$,MMAX$), PSICI(N$,M$,MMAX$)

c       Array XNSTMP added 10/7/94 by PEM to avoid compilation error
c       under IRIX 5.2.
        REAL HTPWR(N$,MMAX$), HTPWI(N$,MMAX$), MVEC(MMAX$)
        REAL    EPFY2(N$,M$,MMAX$),EPFZ2(N$,M$,MMAX$),UNMC(24,40)
        REAL    EPDV (N$,M$,MMAX$),XN2TMP(M$),XNSTMP(MP$),SURF(N$,M$)
        REAL YNMC(24),ZNMC(40)
        COMPLEX DPHIDX(N$,M$),DPHIDY(N$,M$),DPHIDZ(N$,M$) 
        COMPLEX HEIGHT_L(N$,MMAX$) 
        
	COMPLEX ONEC/(1.0, 0.0)/, EYEC/(0.0, -1.0)/, EYEP/(0.0, 1.0)/


c       Added 8/7/95 by PEM:
        INTEGER IWRINT, IDOY360, IYEARC
        DATA    IWRINT /5/

	LOGICAL WERROR,rest
	DATA WERROR/.false./,ai/(0.,1.)/

c	simple time dependent planetary wave model 

C       SET UP CONSTANTS AND ZERO OUT THE TOTAL PW FORCING (PWMOM)

c        CALL PRMINMAX(PSIC(1,1,1),'PWV_AMP-00',N$, M$,1)

        DMPMX   = 1./DAYL
        EQDMP0  = ISW(54)/100.
        AMP_KLG = ISW(55)/100.        ! ISW(55) currently = 0.0  (EF, Nov. 2007)

        
C   Removed 9/18/96 by PEM.  These are no longer used with the new timing (timer.f).
C        NDAY0   = NDAY
C        NSTEP0  = NSTEP

C   End Removed 9/18/96.


c    WRITE(*,*) 'EQ. DAMPING FACTOR  ---> ',EQDMP0
c    WRITE(*,*) 'AMPLITUDE REDUCTION ---> ',AMP_KLG

       

C ---STORE OFF N-SQUARED PROFILES SO THAT THEY 
C ---CAN BE HARMLESSLY FUDGED FOR THE PWAVE CODE
C
        DO J=1,M$
           XN2TMP(J)=XN2(J)
        END DO
        DO J=1,MP$
           XNSTMP(J)=XN2ST(J)
        END DO



CC        M=1



C  Modified 12/29/95 by PEM.  This section now performs the task of
C    zeroing out the wave foring height field poleward of the cutoff
C    latitude specified by switch #63 in the switches.dat file.
C    Forcing heights are extended to 90 degrees latitude during November in the southern hemisphere.

cjer        ZERLAT = ISW(63)
cjer heights come in here in cm

        DO MSS=1,MMAX$
           DO J  =1,N$

c              IF (YP(J).GE.ZERLAT) THEN
c                    HEIGHT_L(J,MSS) = 0.0
c                       DUM1C(J,MSS) = 0.0

c                 ELSEIF (((YP(J).LE.-ZERLAT).AND.
c     2                    ((MONTH.LT.10).OR.
c     3                     (MONTH.GT.11))).OR.
c     4                  (YP(J).LE.-80.0)) THEN
c                    HEIGHT_L(J,MSS) = 0.0
c                       DUM1C(J,MSS) = 0.0

c                 ELSE
                    HEIGHT_L(J,MSS) = HEIGHT(J,MSS)
                       DUM1C(J,MSS) = HEIGHT(J,MSS)
c                 ENDIF

              END DO 
           END DO 


C  End Modified 12/29/95, by PEM.

cjer or use my code for zeroing topog poleward of 60 for may thru
c october, 80 otherwise
c  if (nstep.eq.1) write(6,*) ' TIME VARYING CUT OFF LATITUDE'
c     DO MSS=1,MMAX$
c     DO J  =1,N$
cjer zero topog poleward of 60 for may through october, 80 otherwise
c       tempht=HEIGHT(J,MSS)

c       if(yp(j) .lt. -80.) tempht=0.
c       if(yp(j) .gt.  80.) tempht=0.
c       if(nday.ge.121 .and. nday.le.304 .and. yp(j).le.-60.) tempht=0.  

c       HEIGHT_L(J,MSS)=  tempht
c          DUM1C(J,MSS)=  tempht

c        HEIGHT_L(J,MSS)=  HEIGHT(J,MSS)
c          DUM1C(J,MSS)=  HEIGHT(J,MSS)
c     END DO 
c     END DO 


      IF ( ISW(43) .EQ. 1  ) THEN 
C   FIRST SMOOTH THE CRAP OUT OF THE TOPOGRAPHIC FORCING
C

        DO MSS=1,MMAX$
        DO J  =1,4
           HEIGHT_L(J,   MSS)= 0.
           HEIGHT_L(N$-J,MSS)= 0.
        END DO 
        END DO 

        DO MSS=1,MMAX$
        DO J  =4,N$-3
           DUM1C(J,MSS)=( HEIGHT_L(J-1,MSS) + HEIGHT_L(J+1,MSS) +
     2                  HEIGHT_L(J-2,MSS) + HEIGHT_L(J+2,MSS) +
     3                  HEIGHT_L(J-3,MSS) + HEIGHT_L(J+3,MSS) +
     4                        HEIGHT_L(J,MSS)  )/7.0 
        END DO 
        END DO 

        DO MSS=1,MMAX$
        DO J  =3,N$-2
           DUM1C(J,MSS)=( DUM1C(J-1,MSS) + DUM1C(J+1,MSS) +
     1                  DUM1C(J-2,MSS) + DUM1C(J+2,MSS) +
     2                        DUM1C(J,MSS)  )/5.0 
        END DO 
        END DO 

        DO ISMTH=1,3
        DO MSS=1,MMAX$
        DO J  =2,N$-1
           DUM1C(J,MSS)=( DUM1C(J-1,MSS) + DUM1C(J+1,MSS) +
     2                        DUM1C(J,MSS)  )/3.0 
        END DO 
        END DO 
        END DO 
      END IF  !!! FINISHED CRAP SMOOTHING

        DO MSS=1,MMAX$
        DO J  =1,N$
           IF(YP(J).LT.0.0) THEN 
             AMP_KLG=MWVNS(2,MSS)*(1.000/1000.)
             HEIGHT_L(J,MSS)= DUM1C(J,MSS)*AMP_KLG
           ENDIF
           IF(YP(J).GE.0.0) THEN 
             AMP_KLG=MWVNS(3,MSS)*(1.000/1000.)
             HEIGHT_L(J,MSS)= DUM1C(J,MSS)*AMP_KLG
           ENDIF
        END DO 
        END DO 

cjer looks like heights are now in decimeters
c    amp_klgs were roughly .1 - .3, dum1c in cm 

        DO MSS=1,MMAX$
        DO J  =1,N$
           HTPWR(J,MSS)=  REAL(HEIGHT_L(J,MSS))
           HTPWI(J,MSS)= AIMAG(HEIGHT_L(J,MSS))
        END DO 
        END DO 

c        CALL OUTP(66,NSTEP,NDAY,NYEAR,YP,ZP(1),HTPWR,N$,1,MMAX$,0)
c        CALL OUTP(67,NSTEP,NDAY,NYEAR,YP,ZP(1),HTPWI,N$,1,MMAX$,0)

C  CONSTANT HEIGHT FORCING FOR P-WAVES
c

        DT00=DT   !SAVE ORIGINAL TIME STEP


        CALL PRMINMAX(UBX,  'UBX-PWV---',NP$, MP$,1)



cc        CALL GETQB




        DT=DT00/(ISW(48))

        DAYFR =( DT/DAYL )



        IDAYFR=int( 1.00/DAYFR )

c      WRITE(*,*) ' REDUCED DT FROM',DT00,' TO',DT,' FOR PWAVE'
c      WRITE(*,*) ' PWAVE DT IS ',IDAYFR,' OF A DAY '
c      WRITE(*,*)' USING  IMPL. TRAP. J-P WAVES '

	dtp=dt*.5
	DRY=DY/A
	ss=.25/H-1./dz
	rr=.25/H+1./dz

	DO 105 K=1,M$
	xn(K)=sqrt(xn2(K))
	DO 105 J=1,n$
	pwmom(J,K)=0.
        kyy(j,k)=0.
105	CONTINUE


        DO K=1,M$
        DO J=1,N$
           UBW(J,K)=( UBX(J,  K) + UBX(J,  K+1) 
     >             + UBX(J+1,K) + UBX(J+1,K+1) )*0.25
        END DO
        END DO


        IF ( ISW(45) .EQ. 0 ) THEN
C*************** READ IN NMC ZONALLY AVG'D WIND FOR TESTING PURPOSES *********
        WRITE(*,*) ' .... READING NMC WIND FIELD  '

       open(unit=597, file='nmc_ubar.dat', status='old',
     >    convert="big_endian", form='unformatted')
         READ(597) UNMC
       close (597)

ccelf        READ(97) UNMC
   
        DO J=1,24
           YNMC(J)=7.2*J-90.
        END DO
        DO K=1,40
           ZNMC(K)=K*2.52
        END DO

        CALL Regrid(UBW,YP,ZP,N$,M$,UNMC,YNMC,ZNMC,24,40,0)

        DO K=1,M$
        DO J=1,N$
           UBW(J,K)=UBW(J,K)*100.  ! CONVERT TO CM/S
cc           UBW(J,K)=UBW(J,K)*0.+2000.  ! CONVERT TO CM/S
        END DO
        END DO
C***************************************
        END IF   


        IF ( ISW(45) .EQ. 2 ) THEN
C*** SET UP UNIFORM CONSTANT ******
C*** WIND FOR TESTING PURPOSES ****
        WRITE(*,*) ' .... MAKING UNIFORM WIND FIELD  '

        DO K=1,M$
           XN2(k)=1.e-4
           XN2ST(k)=1.e-4
        DO J=1,N$
           UBW(J,K)=2000.  ! CONVERT TO CM/S
        END DO
        END DO
        XN2ST(mP$)=1.e-4
C***************************************
        END IF   
  
        
        KSRF=ISW(44)
        DO K=1,KSRF
        DO J=1,N$
           UBW(J,K)=UBW(J,KSRF)
           XN2(K)  =XN2(KSRF)
           XN2ST(K)=XN2ST(KSRF)
        END DO
        END DO
        DO K=2,KSRF
        DO J=1,N$
           UBW(J,K)=(UBW(J,K-1)+UBW(J,K)+UBW(J,K+1))/3.0
           XN2(K)  =(XN2(K-1)+XN2(K)+XN2(K+1))/3.0
           XN2ST(K)=(XN2ST(K-1)+XN2ST(K)+XN2ST(K+1))/3.0
        END DO
        END DO



        CALL MRS2RBR(UBW,QUBX)

        CALL GETQB



        DO K=1,M$
        DO J=1,N$
           QBRYW(J,K)=( QBRYWX(J,  K) + QBRYWX(J+1,K+1) +
     >                  QBRYWX(J+1,K) + QBRYWX(J,  K+1)   
     >                                         )*0.25           

        END DO
        END DO





        DO K=1,M$
        DO J=1,N$
           UBW(J,K)=AMAX1(UBW(J,K), 200.  )  ! keep UBW slightly positive
c           IF(UBW(J,K) .LE. 200.) QBRYW(J,K)=0.0
        END DO
        END DO

        DO K=1,M$
        DO J=1,N$
           QBRYW(J,K)=AMAX1(QBRYW(J,K), 3.e-14)  !restrain qbary, changed from 1.e-14 => 3.e-14 (May 09)
        END DO
        END DO





C CALCULATE Phi BY LOOPING OVER ZONAL WAVENUMBERS


        DO NPWV=1,ISW(47)    !LOOP OVER P-WAVE SOLVER   




	DO 300 MSS=1,ISW(53)          ! MMAX$
C

        M=MWVNS(1,MSS)    ! MSS refers to INDEX in (*,*,MMAX$) arrays
c                         ! M refers to actual PHYSICAL wavenumber
c                         

        CALL BKGDMP
        CALL ENHDMP(MSS)


c	print 666,m
666	format(' m=',i5,' <<<<******* ZONAL**ZONAL**ZONAL P-WAVE # ')
C	print 851,(height(j,MSS),j=1,n$)
851	format(' height =',2g14.4)

C        IF(ISW(50).EQ.0)  WRITE(6,*)'  NO SPHERICITY IN PWAVE '
C        IF(ISW(50).EQ.1)  WRITE(6,*)' YES SPHERICITY IN PWAVE '

852     FORMAT( '*** SPHERICITY',I2,';  COTANGENT TERM',I2,
     >          ';  TANGENT TERM',I2,';  WAVENUMBERS',I2,' ***')

c        WRITE(6,852) ISW(50),ISW(51),ISW(52),ISW(53) 

	do 100 k=1,m$

cc     dvnc=dtp*VNC(K)*0.                         ! NEWT. COOLING (0)
cc       upper limit the raleigh friction
        rral=ral(k)*0.                             ! RAYL. DRAG    (0)
c        if (ral(k).gt.(1./(dayl*30.))) rral=1./(dayl*10.)
        deq(k)= amin1(amax1((rral/tomeg)**2,.06),.1)



        DVNC0=0.*dtp  !*PWVDMP         ! (1.00/(5.0*DAYL) )         
        DRAL0=0.*dtp  !*PWVDMP         ! (1.00/(5.0*DAYL) )         


	do 100 j=1,n$

           j0=20
           IF(ISW(50).EQ.1) J0=J

C       DRAL=DTP*rral

        S4=1.+ ai*ubw(j,k)*m*dtp/(A*C(J0))
        s5=1.- ai*ubw(j,k)*m*dtp/(A*C(J0))


        DRAL=DRAL0 + EQDMP(J,K)*DTP
        DVNC=DVNC0 + EQDMP(J,K)*DTP

        XNU0(j,k)= S4 + DRAL                       
        XNU1(j,k)= S4 + DVNC                        

        XNU0L=     s5 - DRAL                        
        XNU1L=     s5 - DVNC                        

        BIGQ(J,K)=ai*dtp*M*QBRYW(J0,K)/(A*C(J0))

        FC(j,k)= XNU1L*psa(j,k,MSS) + XNU0L*psb(j,k,MSS)
     >         -  bigq(j,k)*psic(j,k,MSS)

        PSAO(J,K)=PSA(J,K,MSS)
        PSBO(J,K)=PSB(J,K,MSS)
        PSICO(J,K)=PSIC(J,K,MSS)

100	continue

C	CALCULATE COEFICIENTS FOR CSLAK
	DO 200 k=1,M$
c        deq=.05
	DO 200 j=1,n$

           j0=20 !j   ! INDEX FOR LATITUDE DEP. CONSTANTS E.G. S,C,CF     
           IF(ISW(50).EQ.1) J0=J

        req= 1.000                                           

	COE1=  XNU0(j,k)                                      

	coe2= -XNU0(j,k)*(
     >    ISW(51)*2.*( C(J0)/S(J0) ) + ISW(52)*( S(J0)/C(J0) )    )  
c                     = cot(LAT) =              = tan(LAT) =
c
	COEZ=  XNU1(J,K)*cf(j0)*cf(j0)/xn2(k)

	  AZC(j,k) = COEZ/(DZ*dz)
	  CZC(j,k) = COEZ/(DZ*dz) 

	  AYC(j,k) = (coe1/(Dy*dy)) + coe2/(2.0*DY*A)
	  cyc(j,k) = (coe1/(DY*dy)) - coe2/(2.0*DY*A)

	  BCA(j,k) = -2.0*COEZ/(DZ*DZ)-COEZ*0.25/(H*H)
	  BCB(J,K) = -2.0*coe1/(Dy*dy)-M*M*coe1/(A*C(J0)*A*C(J0))

	  bc(j,k)=bca(j,k)+bcb(j,k)+bigq(j,k)

	  psic(j,k,MSS)=0.

c	set up upper b.c.
C             t is constant
c	if (k.eq.m$) then
c	bc(j,k)=bc(j,k)-azc(j,k)*ss/rr
c	endif

c	set up lower boundary condition

	  IF(k.EQ.1) THEN

C   FLOW OVER TOPOGRAPHY

c	fg(j) = -HEIGHT_L(j,MSS)*xn2(k)
ccc	need this ?bca(j,k) = bca(j,k)-rr*czc(j,k)/ss	
c	bc(j,k)=bc(j,k)-rr*czc(j,k)/ss	
c	fc(j,k)=fc(j,k)-czc(j,k)*fg(j)/ss


C  for using NCEP HGTS, also include amplitude Kluge here (SH/NH) and apply below

           IF(YP(J).LT.0.0) AMP_KLG = MWVNS(2,MSS)/1000.
           IF(YP(J).GE.0.0) AMP_KLG = MWVNS(3,MSS)/1000.


C  over write array COMPLEX HEIGHT_L(N$,MMAX$) w/ NCEP REANALYSIS GEOP hgts (in METERS)
C  1st index: 1=zonal avg, 2-6 are waves 1-5;   2nd index:  1=Ak,  2=Bk
C
C  GPBC1(6,2,N$) has been PROPERLY interpolated to current dynamics grid for current day/year (in TWODS)
C  ramp up wave 1 amps at SH high lats ONLY - DON'T USE
C
        ampinc = 1.
cpw        if (MSS .eq. 1  .and.  yp(j) .le. 0.) then
cpw           ampinc = 1.5*(SIND(yp(j)))**2
cpw           if (ampinc .le. 1.) ampinc = 1.
cpw        endif
                                                                    ! convert TO CENTIMETERS
        HEIGHT_L(j,MSS) = ONEC*GPBC1(MSS+1,2,j)*100.*AMP_KLG * ampinc
     >                  + EYEC*GPBC1(MSS+1,1,j)*100.*AMP_KLG * ampinc


C   FIXED HEIGHT BC
C       FG(J)=FTIME*height_L(j,MSS)*GZ
        FG(J)=1.000*height_L(j,MSS)*GZ
	fc(j,k)=fc(j,k)-czc(j,k)*FG(J)
	endif

200	CONTINUE



C ELLIPTIC OPERATOR INVERSION


cc        CALL ELLGUAR(AYC, N$, M$)
cc        CALL ELLGUAR(CYC, N$, M$)

        DO K=1,M$
        DO J=1,N$
           DUM1(J,K)= REAL( AYC(J,K) )
           DUM2(J,K)=AIMAG( AYC(J,K) )
        END DO
        END DO

        CALL PRMINMAX(DUM1,  'AYC---RL--',N$, M$,1)
        CALL PRMINMAX(DUM2,  'AYC---IM--',N$, M$,1)

        DO K=1,M$
        DO J=1,N$
           DUM1(J,K)= REAL( CYC(J,K) )
           DUM2(J,K)=AIMAG( CYC(J,K) )
        END DO
        END DO

        CALL PRMINMAX(DUM1,  'CYC---RL--',N$, M$,1)
        CALL PRMINMAX(DUM2,  'CYC---IM--',N$, M$,1)



        IF(ISW(62) .EQ. 0) THEN
C     1    CALL CSLAK(AYC,AZC,BC,cyc,CZC,FC,PSIC(1,1,MSS),n$,M$,
C     2    WERROR,ERMAX,ERRET)
           PRINT*, "CSLAK is not currently supported."
           PRINT*, "Routines in file slak.f must be modified before"
           PRINT*, "it can be compiled and linked into the model."
           PRINT*, "Program Terminated."
           PRINT*, " "
           STOP
           ENDIF

        IF(ISW(62) .EQ. 1) 
     1    CALL CJSOR2(AYC,AZC,BC,cyc,CZC,FC,PSIC(1,1,MSS),n$,M$,
     2    WERROR,ERMAX,ERRET,NSTEP)

c	print 4455,erret
C4455	format(' returned error from cslak = ',2g14.4)



c       now COMPUTE AND store psa and psb

	do 115 k=1,m$
	do 115 j=1,n$


         j0=20
         IF(ISW(50).EQ.1) J0=J
         
        
          coe2= -1.00*(
     >    ISW(51)*2.*( C(J0)/S(J0) ) + ISW(52)*( S(J0)/C(J0) )    )  
c                     = cot(LAT) =               = tan(LAT) =
c


          req=1.000                      ! (deq(k) + s(j0)*s(j0))

	  AZC(j,k) = 1./(DZ*dz)*cf(j0)*cf(j0)/xn2(k)
	  CZC(j,k) = azc(j,k)

	  AYC(j,k) = 1./(DY*DY) + coe2/(2.0*DY*A)
	  cyc(j,k) = 1./(DY*DY) - coe2/(2.0*DY*A)

	  BCA(j,k) = (-2.0/(DZ*DZ)-0.25/(H*H))*cf(j0)*cf(j0)/xn2(k)

	  BCB(J,K) = (-2.0/(Dy*dy) )
     >                -M*M/(A*C(J0)*A*C(J0))
115	continue

	do 110 k=1,m$
	do 111 j=2,n1$
           PSB(J,K,MSS)=AYC(J,K)*PSIC(J+1,K,MSS)+BCB(J,K)*PSIC(J,K,MSS)
     >                  +cyc(J,K)*PSIC(J-1,K,MSS)
C     >                  +cyc(J,K)*1PSIC(J-1,K,MSS)
C BThomas 2Nov2017 - the varbiable 1PSIC does not appear anywhere but here (and below)
C                    I think it's a typo, and from symmetry (j+1,j,j-1), I think 
C                    it should just be PSIC
111	CONTINUE
	J=1
	PSB(J,K,MSS)=AYC(J,K)*PSIC(J+1,K,MSS)+BCB(J,K)*PSIC(J,K,MSS)
	J=N$
	PSB(J,K,MSS)=BCB(J,K)*PSIC(J,K,MSS)+cyc(J,K)*PSIC(J-1,K,MSS)
110	CONTINUE
	do 210 j=1,n$
	do 211 k=2,m1$
           PSA(J,K,MSS)=AZC(J,K)*PSIC(J,K+1,MSS)+BCA(J,K)*PSIC(J,K,MSS)
     >               +CZC(J,K)*PSIC(J,K-1,MSS)
C     >               +CZC(J,K)*1PSIC(J,K-1,MSS)
C Ditto above note (BThomas 2Nov2017)
211	CONTINUE

	K=1
C           FLOW OVER TOPOGRAPHY
c	PSA(J,K,MSS)=AZC(J,K)*PSIC(J,K+1,MSS)+(BCA(J,K)-czc(j,k)*rr/ss)
c     1*PSIC(J,K,MSS)+fg(j)*czc(j,k)/SS

C           FIXED HEIGHT BC                                           - psa(n$,m$,mmax$)
	PSA(J,K,MSS)=AZC(J,K)*PSIC(J,K+1,MSS)+BCA(J,K)*PSIC(J,K,MSS)
     1+fg(j)*czc(j,k)

c	zero amp upper bc
	K=M$
	psa(j,k,MSS)=bca(j,k)*psic(j,k,MSS) + czc(j,k)*psic(j,k-1,MSS)
c	t=0 upper bc
c	PSA(J,K,MSS)=(BCA(J,K)-ss*azc(j,k)/rr)*PSIC(J,K,MSS)+CZC(J,K)*PSIC(J,K-1,MSS)
210	CONTINUE

C       TIME SMOOTHING
c        w1A=.5 !9   !/3.
c        w2A=1.-w1
        w1a=.999     ! w1
        w2a=1.-W1A   ! w2
        
        DO 705 K=1,M$
        DO 705 J=1,N$
        PSA(J,K,MSS)= W1a*PSA(J,K,MSS) + W2a*PSAO(J,K)
        PSB(J,K,MSS)= W1a*PSB(J,K,MSS) + W2a*PSBO(J,K)
        PSIC(J,K,MSS)= W1a*PSIC(J,K,MSS) + W2a*PSICO(J,K)
705     CONTINUE


cpw69 - interpolate equatorial zone here in COMPLEX PSIC(N$,M$,mmax$) (values are 0. otherwise)
Cpw69   do REAL and IMAGINARY separately, restore to proper sign - iinc is -1 (SH) or +1 (NH)
Cpw70      - BEST NOT TO DO THIS
cpw70
cpw70        DO 707 J=1,N$
cpw70           if (ABS(YP(J)) .le. 5.) then
cpw70
cpw70              iinc = INT(SIGN(1., YP(J)))
cpw70
cpw70              DO 717 K=1,M$
cpw70                aa1 = REAL(PSIC(j+iinc,k,mss))
cpw70                aa2 = REAL(PSIC(j+iinc+iinc,k,mss))
cpw70                if (aa2 .ne. 0.) aar = SIGN(aa1*aa1/aa2, aa1)
cpw70                if (aa2 .eq. 0.) aar = aa1  ! 0.
cpw70
cpw70                aa1 = AIMAG(PSIC(j+iinc,k,mss))
cpw70                aa2 = AIMAG(PSIC(j+iinc+iinc,k,mss))
cpw70                if (aa2 .ne. 0.) aai = SIGN(aa1*aa1/aa2, aa1)
cpw70                if (aa2 .eq. 0.) aai = aa1  ! 0.
cpw70
cpw70                PSIC(j,k,mss) = ONEC*aar + EYEP*aai
cpw70 717          CONTINUE
cpw70
cpw70           endif
cpw70 707    CONTINUE



C	NOW GET u', v', w', phi',  FOR EACH WAVENUMBER
C       PHI IS THE GEOPOTENTIAL

	DO 610 J=1,M$
	DO 610 I=1,n$
	  Phi(I,J)=ez2h(j)*PSIC(I,J,MSS)
610	CONTINUE

        CALL GARCIA(PHI,M,MSS)



C*****************
C
        IF ( ISW(56) .EQ. 0 ) THEN
  
           CALL EPD1( PHI,M,EPDV(1,1,MSS),EPFY2(1,1,MSS)
     2                                 ,EPFZ2(1,1,MSS),N$,M$)

C           WRITE(6,*) 'DID EP-FLUX,DIV IN SUB. EPD1 '
             
        END IF   !EP CALCULATION <--- ISW(56)= 0
C*****************
C
  
        IF ( ISW(56) .EQ. 1 ) THEN
  
           CALL EPD2( PHI,M,EPDV(1,1,MSS),EPFY2(1,1,MSS)
     2                                 ,EPFZ2(1,1,MSS),N$,M$)

c           WRITE(6,*) 'DID EP-FLUX,DIV IN SUB. EPD2 **2*2*2** '
             
        END IF   !EP CALCULATION <--- ISW(56)= 0
C******************
C  

cc        CALL MATSQP(PHI,M,N$,M$,QPRIME,QPRMY)
        CALL MATSQP(PHI,M,MSS,N$,M$) 
    
        DO K=1,M$
        DO J=1,N$
           DUM1(J,K)=CABS(QPRMY(J,K))
        END DO
        END DO

        IF(M.EQ.1) IQUA=52
        IF(M.GT.1) IQUA=52+10000+M*100
c        CALL OUTA(IQUA,NSTEP,NDAY,NYEAR,YP,ZP,DUM1, N$,M$,0)

        IF(M.EQ.1) THEN
c        CALL OUTA(11,NSTEP,NDAY,NYEAR,YP,ZP,QBRYW,N$,M$,0)
        ENDIF

        IF(ISW(47).NE.1) THEN
        WRITE(*,*)' TESTING P-WAVE PARAM STEP=',NPWV,' DAY=',NDAYP
cc          RETURN


C Removed 9/18/96 by PEM. This section is not used with the new timing (timer.f).
C       NSTEP=NSTEP0+NPWV
C       NDAY=INT( NPWV*DAYFR )+NDAY0
C       NYEAR=1
C End Removed 9/18/96.

        END IF


300     CONTINUE   ! LOOP OVER WAVENUMBER


        CALL PRMINMAX(PSIC(1,1,1),'PWV_AMP---',N$, M$,1)

        DO MSS=1,MMAX$
        DO K=1,M$
        DO J=1,N$
           PSICI(J,K,MSS)=AIMAG(PSIC(J,K,MSS))  *ez2h(K) 
           PSICR(J,K,MSS)= REAL(PSIC(J,K,MSS))  *ez2h(K)
        END DO
        END DO 
        END DO 
        

c        CALL OUTP(59,NSTEP,NDAY,NYEAR,YP,ZP,UBW(1,1) ,N$,M$,1,0)
c        CALL OUTP(61,NSTEP,NDAY,NYEAR,YP,ZP,EQDMP(1,1),N$,M$,1,0)
c        CALL OUTP(55,NSTEP,NDAY,NYEAR,YP,ZP,PSICR,N$,M$,MMAX$,0)
c        CALL OUTP(56,NSTEP,NDAY,NYEAR,YP,ZP,PSICI,N$,M$,MMAX$,0)
C
c        CALL OUTP(60,NSTEP,NDAY,NYEAR,YP,ZP,QBRYW(1,1) 
c     2                                           ,N$,M$,1,0)
c        CALL OUTP(64,NSTEP,NDAY,NYEAR,YPP,ZPP,QBRYWX(1,1) 
c     2                                           ,NP$,MP$,1,0)
c        CALL OUTP(57,NSTEP,NDAY,NYEAR,YP,ZP,XN2(1) 
c     2                                           ,1 ,M$,1,0)
c        CALL OUTP(62,NSTEP,NDAY,NYEAR,YP,ZP,EPFY2,N$,M$,MMAX$,0)
c        CALL OUTP(63,NSTEP,NDAY,NYEAR,YP,ZP,EPFZ2,N$,M$,MMAX$,0)
c        CALL OUTP(65,NSTEP,NDAY,NYEAR,YP,ZP,EPDV ,N$,M$,MMAX$,0)

C


        END DO  ! END OF P-WAVE TEST LOOP (DONE ISW(47) TIMES)


        DT=DT00

c  WRITE(*,*) ' RESTORED DT TO',DT,' AFTER PWAVE'


        
        DO K=1,M$
        DO J=1,N$
           DUM1(J,K)=0.000
        END DO
        END DO 

        DO M=1,ISW(53)
        DO K=1,M$
        DO J=1,N$
           DUM1(J,K)=DUM1(J,K)+MWVNS(4,M)*EPDV(J,K,M)

c  nantest=test_nan( dum1(j,k) )
c    if(nantest) then
c     PRINT*, "NAN for dum1 in GETMOMX:"
c     PRINT*, "   J = ", J, ", K = ", K, 'M = ',M
c     print*, 'dum1,mwvns,epdv = ',dum1(j,k),mwvns(4,m),epdv(j,k,m)
c     STOP
c    ENDIF

C Added 6/4/96 by PEM.
c  IF ((MONTH.GE.10) .AND. (MONTH.LE.11)) THEN
c    IF ((J.LT.5) .AND. (K.GT.5) .AND. (K.LT.11)) THEN
c      DUM1(J, K) = 10.0*DUM1(J, K)
c    ENDIF
c  ENDIF
C End Added 6/4/96.

        END DO
        END DO
        END DO


C
C  get current theta - THCE(NP$,MP$)
C
         DO K=1,MP$
         DO J=1,NP$
            thce(j,k) = xc(j,k,1)
         END DO
         END DO


c  interpolate UBX(NP$,MP$) to UBW1(N$,M$) (cm/sec), and thce(NP$,MP$) -> thc(N$,M$)

        CALL BINTERP(ypp, NP$, zpp, MP$, UBX, yp, N$, zp, M$, 0,0, UBW1)
        CALL BINTERP(ypp, NP$, zpp, MP$, THCE, yp, N$, zp, M$, 0,0, THC)



C  adjust 80S-85S EP flux so it's at most 20% larger than 75S-80S - wconv89

        DO K=1,M$
        DO J=2,1,-1
           xx1 =     DUM1(J,K)
           xx2 = 1.2*DUM1(J+1,K)
           if (xx1 .le. xx2) DUM1(J,K) = xx2
        END DO
        END DO


C  adjust 85N TOTAL EP flux to be NO MORE THAN 70% of 80N (EF, Sep 2010) - wconv80
C  to avoid getting strange easterlies poleward of 80N in mid-winter
C
        DO K=1,M$
           do J=2,N$
              if (yp(J) .ge. 81.) then
                 epmax = .7*DUM1(J-1,K)
                 if (DUM1(J,K) .lt. epmax) DUM1(J,K) = epmax
              endif
           end do
        END DO


C  WCONV220 - set min High Lat EPFLX to -.05 m/sec/day if westerly above 17 km - DUM1(N$,M$)
C
       eplim = -.05*100./86400.

       DO K=1,M$
           if (ZP(K) .ge. 17.) then

           do J=N$-1,1,-1
              if (yp(J) .le. -50.  .and.  UBW1(J,K) .gt. 0.) then
CC                epmin = DUM1(J+1,K)*C(J)/C(J+1)
                 if (DUM1(J,K) .gt. eplim) DUM1(J,K) = eplim   ! AMIN1(epmin,eplim)
              endif
           end do

           do J=2,N$
              if (yp(J) .ge. 50.  .and.  UBW1(J,K) .gt. 0.) then
CC                epmin = DUM1(J-1,K)*C(J)/C(J-1)
                 if (DUM1(J,K) .gt. eplim) DUM1(J,K) = eplim   ! AMIN1(epmin,eplim)
              endif
           end do

           endif
        END DO



C   compute EPflux divergence using linear-balanced winds here:  - NOT USED

ccccccc  CALL EPD10(PSICR, PSICI, UBW1, THC, IDOY360, DELF10)


c  add in extra PWAVE momemtum source here to DUM1(N$,M$), which is NEGATIVE (EASTERLY)
C  this will impact also KYY
cckyyf8
cckyyf8   - DON'T DO THIS w/ FIXED MODEL Kyys/DELFs
cckyyf8
cckyyf8
cckyyf8
cckyyf8        DO 515 K=1,M$
cckyyf8        DO 515 J=1,N$
cckyyf8
cckyyf8           zrfac = (ZP(k) - 15.)/7.5
cckyyf8           if (zrfac .le. 0.) zrfac = 0.
cckyyf8           if (zrfac .ge. 1.) zrfac = 1.
cckyyf8
cckyyf8           seasrf = 0.
cckyyf8           ampm = 0.
cckyyf8
cckyyf8         if (YP(J) .lt. 0.) then
cckyyf8           ampm = .15
cckyyf8
cckyyf8           if (idoy360 .ge. 285.  .and.  idoy360 .le. 345.)
cckyyf8     >         seasrf = SIND((idoy360-285.)*3.)
cckyyf8         endif
cckyyf8
cckyyf8         if (YP(J) .ge. 0.) then
cckyyf8           ampm = 1.
cckyyf8
cckyyf8             if (idoy360 .ge. 330.) seasrf = SIND((idoy360-330.)*3.)
cckyyf8
cckyyf8             if (idoy360 .le. 90.) seasrf = 1.
cckyyf8 
cckyyf8             if (idoy360 .ge. 90.  .and.  idoy360 .le. 135.)
cckyyf8     >          seasrf = SIND((idoy360-45.)*2.)
cckyyf8          endif
cckyyf8
cckyyf8         seasrf = seasrf * ((SIND(yp(j)))**4)
cckyyf8         if (ABS(yp(j)) .le. 40.) seasrf = 0.
cckyyf8
cckyyf8         extram = -ampm*100./(DAYL*1.) * seasrf * zrfac
cckyyf8
C   don't apply correction if easterly UBW1(N$,M$) (cm/sec)
cckyyf8
cckyyf8         IF (UBW1(J,K) .LT. 0.) extram = 0.
cckyyf8         IF (extram .GT. 0.) extram = 0.
cckyyf8
cckyyf8         DUM1(J,K) = DUM1(J,K) + extram
cckyyf8
cckyyf8 515   CONTINUE
cckyyf8

c  Add in up to 30% of FIxED model DELFS here at high latitudes (current day), 
C  ramp down in mesosphere, first store coupled model DELFs in EPSTR(N$,M$)
c  EPDD(N$,M$) are the fixed model DELFs in m/sec/day - Convert to cm/sec^2 - weight by SIN(lat)^4
c  load into EPEX(N$,M$)
C
C  for PW23,apply NCEP DELF ONLY where MODEL UBAR is EASTERLY, pwyy(N$)=1 in NH; .5 in SH, include tropics
C
C  for PW26, apply correction based ONLY on MODEL UBAR, empirically derived from NCEP UBAR/DELF
C   this is a MODEST correction, use in Both SH/NH, key on UBW1(N$,M$) (cm/sec)
C
ccpw30        DO J=1,N$
ccpw30           pwamp = 1.
ccpw30           if (yp(j) .lt. 0.) pwamp = .25
ccpw30           if (ABS(yp(j)) .le. 40.) pwamp = 0.
ccpw30           pwyy(J) = pwamp * ((SIND(yp(j)))**4)
ccpw30        END DO
ccpw30
ccpw30
ccpw30        pwmn = -3.*100./86400.
ccpw30        pwmx =  0.*100./86400.
ccpw30
ccpw30        DO K=1,M$
ccpw30
ccpw26           zzpw = pwmx - (ZP(K)-50.)/100.
ccpw26           if (zzpw .lt. pwmn) zzpw = pwmn
ccpw26           if (zzpw .gt. pwmx) zzpw = pwmx
ccpw26             EPEX(J,K) = 0.0   ! zzpw*pwyy(J)*EPDD(J,K)*100./86400.
c
C   ramp down to zero in troposphere
ccpw30          zzpw = 1.
ccpw30          if (ZP(K) .le. 17.) then
ccpw30             zzpw = 1. - (17.-ZP(K))/5.
ccpw30             if (zzpw .lt. 0.) zzpw = 0.
ccpw30          endif
ccpw30
ccpw30C              add extra momentum here, convert to cm/sec^2, uux is in m/sec
ccpw30          DO J=1,N$
ccpw30            uux = UBW1(J,K)/100.
ccpw30            EPEX(J,K) = -zzpw*pwyy(J)*(EXP((uux+10.)/50.)-1.)*100./86400.
ccpw30            if (EPEX(J,K) .lt. pwmn) EPEX(J,K) = pwmn
ccpw30            if (EPEX(J,K) .gt. pwmx) EPEX(J,K) = pwmx
ccpw30
ccpw30            EPSTR(J,K) = DUM1(J,K)
ccpw30            DUM1(J,K) = EPSTR(J,K) + EPEX(J,K)
ccpw30          END DO
ccpw30        END DO 
ccpw30


Ccwconv48    zero-out EP flux in lower trop below 3.5 km - DUM1(N$,M$)

        ZEPFD = 3.5

        DO K=1,M$
        DO J=1,N$
           IF(ZP(K) .LT. ZEPFD) DUM1(J,K) = 0.0
        END DO
        END DO 


        CALL MRS2RBR( DUM1 , DUM1X )

cccwconv48        ZEPFD=1.00*ISW(64)

        DO K=1,MP$
        DO J=1,NP$
           IF (ZPP(K) .LT. ZEPFD) DUM1X(J,K)=0.00
           EPFLX(J,K) = DUM1X(J,K)
           EPFLXS0(J,K) = EPFLX(J,K)

C  load in REAL*8 EPFLXS0 array for smoothing

cccwconv48    IF(ZPP(K) .LT. ZEPFD) EPFLX(J,K)=0.00
        END DO
        END DO 
       

C----RESTORE N-SQUARED PROFILES
C----
        DO J=1,M$
           XN2(J)=XN2TMP(J)
        END DO
        DO J=1,MP$
           XN2ST(J)=XNSTMP(J)
        END DO

cc        CALL OUTA(13,NSTEP,NDAY,NYEAR,YP,ZP,EPDV(1,1,1),N$,M$,0)
cc        CALL OUTA(13,NSTEP,NDAY,NYEAR,YPP,ZPP,EPFLX,NP$,MP$,0)

        IF (ISW(59) .eq. 1) THEN
           CALL SMTHR3( EPFLX, N$, M$ )
           WRITE(6,*) '  >>******** SMOOTHING PWAVE EP FLUX DIV'
           WRITE(6,*) '  >>******** BEFORE APPLYING TO MOM. EQ.'
        ENDIF


C 
C  WCONV117 - smooth EPFLX 2X - use SMOOTH5
C  REAL*8 EPFLXS0(NP$,MP$), EPFLXS1(NP$,MP$), EPFLXS2(NP$,MP$)
C
ccwconv230  CALL SMOOTH5(EPFLXS0, EPFLXS1, NP$, MP$)
ccwconv230  CALL SMOOTH5(EPFLXS1, EPFLXS2, NP$, MP$)


C  WCONV117 - reload smoothed EPFLUX values to EPFLX(NP$,MP$) array, set limits (cm/sec^2)

        DO K=1,MP$
        DO J=1,NP$

ccwconv230  EPFLX(J,K) = EPFLXS2(J,K)

CC
C  changed here - EF - zero out DELF at ONLY IN mid-high latitude troposphere/LS
CCef 
cef   IF (UBX(J,K) .LT. 200.) THEN
cckyyf12           IF (UBX(J,K) .LT. 000.) THEN

           IF (ABS(YPP(J)) .ge. 40. .and.  ZPP(K) .LE. 20.) THEN
              IF (UBX(J,K) .LT. 0.) then 
                 EPFLX(J,K) = 0.
                 EPEXX(J,K) = 0.
              ENDIF
           ENDIF

           IF (EPFLX(J,K) .GT. 0.00) THEN
               EPFLX(J,K)=0.000
ccpw30               EPEXX(J,K) = 0.
           ENDIF
        END DO
        END DO 



C  interpolate NCEP DELF EPEXX(NP$,MP$) (cm/sec^2) -> EPXE1(N$,M$) (used for Kyy)
C
      CALL BINTERP(ypp, NP$, zpp, MP$, EPEXX, yp, N$, zp, M$, 0,0,EPXE1)



        DO K=1,M$
        DO J=1,N$
           DUM2(J,K)=QBRYW(J,K)  ! NEEDED TO PASS QBRYW TO JPWVKYY
        END DO
        END DO 

c        IKYYSB=ISW(21)
        IKYYSB=2

        IF(IKYYSB.EQ.1) CALL JPWVKYY   
        IF(IKYYSB.EQ.2) CALL JPWVKY2(IDOY360, EPXE1)

c        CALL OUTP(68,NSTEP,NDAY,NYEAR,YPP,ZPP,PWKYY,NP$,MP$,1,0)

cc    REAL PSICR(N$,M$,MMAX$), PSICI(N$,M$,MMAX$), T(N$,M$), UBX(NP$,MP$)
cc         DUM1(N$,M$), DUM1X(NP$,MP$), EPFLX(NP$,MP$)

ccelf       if (mod(idoy360-15, 30) .eq. 0.0) then 
ccelf          write (797) N$, M$, MMAX$, idoy360
ccelf          write (797) yp, zp, ypp, zpp, dz
ccelf          write (797) PSICR, PSICI, T, UBX
ccelf          write (797) DUM1, DUM1X, EPFLX
ccelf
ccelf          write (797) ez2h, phi, psic
ccelf
ccelf          write (797) height_L
ccelf       endif
ccelf
ccelf     write (797) height
ccelf       write (797) height_L
ccelf       write (797) gpbc1

C
C  get longitudinal temperature distribution from pwave geop fields:  units are cm^2/sec^2; DZ is cm
C    this will be used in SETDAILY (just use values from last dynamics time step for current day)
C    also compute eddy heating correction terms


         CALL TEMPLON(PSICR, PSICI, THCE, UBX, IDOY360)


	RETURN
	END

C  END JPWAVE
C



      SUBROUTINE JPWVKY2(IDOY360, EPXE1)
C -----
C     SUBROUTINE TO CALCULATE PLANETARY HORIZONTAL DIFFUSION
C     FROM EP-FLUX DIVERGENCE AND D(QBAR)/DY
C REQ:
C--     DUM1 SHOULD CONTAIN EPDV  BEFORE SUBR. IS CALLED
C--     DUM2 SHOULD CONTAIN QBRYW BEFORE SUBR. IS CALLED
 

	INCLUDE 'PARAM.INC'
	INCLUDE 'COMMONC.INC'
	INCLUDE 'COMMOND.INC'
	INCLUDE 'COMMONW.INC'
	INCLUDE 'COMMONT.INC'

        COMMON/CKYY1/ BGKYY0(N$,M$), BGKYY(N$,M$), EPOTX(NP$,MP$),
     >                XKYYFX(N$,M$)


        REAL*8 KYY8(N$,M$), KYY8S(N$,M$), KYY8S2(N$,M$)
        REAL*8 KYY8S3(N$,M$), KYY8S4(N$,M$)


        INTEGER IDOY360

        real kyylow, kyyhigh, EPXE1(N$,M$)
        REAL      FSCALE(N$,M$)
        LOGICAL   FSINIT
        DATA      FSINIT /.FALSE./


        IF (.NOT.FSINIT) THEN
           DO K=1, M$
              DO J=1, N$
                 ABSLAT = ABS(YP(J))
                 IF (ABSLAT .LE. 80.0) THEN
                       FSCALE(J,K)= 1.0
                    ELSE
C                       FSCALE(J,K)= EXP(-1.0*(((ABSLAT-80.0)/5.0)**2.0))
C                       FSCALE(J,K) = 0.0
                       FSCALE(J, K) = 1.0
                    ENDIF
                 END DO
              END DO
           FSINIT = .TRUE.
           ENDIF


C  here, DUM1(N$,M$) is E-P flux div,  DUM2(N$,M$) is QBRYW(N$,M$)  in COMMOND.INC
C    add in EPXE1(N$,M$) is the NCEP DELF (easterly ONLY) in cm/sec^2
C
        DO K=1,M$
        DO J=1,N$

           KYY(J,K) =
     1          -((DUM1(J,K) + EPXE1(J,K))/DUM2(J,K))*FSCALE(J,K)

           PWKYY(J,K)=KYY(J,K)

        END DO
        END DO 


cpw125 - reduce Kyy at SPole, 10-20 km for better ozone - don't do for WCONV57
       DO K=1,M$
          J=1
ccwconv57       if (ZP(K) .ge. 10. .and. ZP(K) .le. 20.) KYY(J,K) = .1*KYY(J+1,K)
cpw125        J=N$
cpw125        KYY(J,K) =  0.5*KYY(J-1,K)
       END DO


C  LIMIT KYY VALUES - changed to 1.e8/1.e11 - EF (Nov. 2007) -   KYY(N$,M$)
C    BGKYY now used w/ 9GD (similar to fixed model Kyymins) - BGKYY(N$,M$)
C    EPFLX(NP$,MP$), QBARY(Np$,MP$)
C
C    also load in KYYOT, DON'T do any EPFL adjustment in SKYYOT
C
        CALL SKYYOT(yp, zp, IDOY360, DUM2)     
      

cjer         kyylow=1.0e9
        kyylow  = 1.e8
        kyyhigh = 1.e11


        DO K=1,M$
        DO J=1,N$

ccc  IF(  UBX(J,K) .LT. 200.)      KYY(J,K)=kyylow
ccc  IF( DUM2(J,K) .LT. 1.501E-14) KYY(J,K)=kyylow

           if (  kyy(j,k) .lt. bgkyy(j,k) )   kyy(j,k) = bgkyy(j,k)
           if (  kyy(j,k) .lt. kyylow)        kyy(j,k) = kyylow


C                                load in REAL*8 array for SMOOTH5
           KYY8(J,K) = KYY(J,K)

        END DO
        END DO


C 
C  for 9GL run, smooth Kyys 2X - REAL*4 KYY(N$,M$) - use SMOOTH5
C  REAL*8 KYY8(N$,M$), KYY8S, KYY8S2, KYY8S3, KYY8S4

           CALL SMOOTH5(KYY8, KYY8S, N$, M$)
           CALL SMOOTH5(KYY8S, KYY8S2, N$, M$)

cc           CALL SMOOTH5(KYY8S2, KYY8S3, N$, M$)
cc           CALL SMOOTH5(KYY8S3, KYY8S4, N$, M$)


C   reload smoothed values to KYY array except use UNSMOOTHED VALUES at 18-28 km in tropics

           DO K=1, M$
           DO J=1, N$

              if (zp(k) .ge. 18.  .and.  zp(k) .le. 28.) then
                 if (ABS(yp(j)) .le. 15.) KYY8S2(J,K) = KYY8(J,K)
              endif

              KYY(J,K) = KYY8S2(J,K)


CCkyyfx  reset here to a blend of fixed/coupled model Kyys EVERYWHERE - XKYYFX(N$,M$)
CCkyyfx    KYY(J,K) = .65*KYY8S2(J,K) + .35*XKYYFX(J,K)


C  reset max Kyy here AFTER SMOOTHING

              IF (  KYY(J,K) .GT. kyyhigh)  KYY(J,K) = kyyhigh


C   set max lower trop Kyy = 5.e10
cc9gk  if (zp(k) .le. 5. .and. kyy(j,k) .gt. 5.e10)kyy(j,k)=5.e10

           END DO
           END DO 


        RETURN
        END



      SUBROUTINE JPWVKYY
C -----
C     SUBROUTINE TO CALCULATE PLANETARY HORIZONTAL DIFFUSION
C     FROM EP-FLUX DIVERGENCE AND D(QBAR)/DY
C REQ:
C--     DUM2 SHOULD CONTAIN QBRYW BEFORE SUBR. IS CALLED
 

	INCLUDE 'PARAM.INC'
	INCLUDE 'COMMONC.INC'
	INCLUDE 'COMMOND.INC'
	INCLUDE 'COMMONW.INC'
	INCLUDE 'COMMONT.INC'



        DO K=1,MP$
        DO J=1,NP$
           PWKYY(J,K) =
     1          -EPFLX(J,K)/QBRYWX(J,K)
        END DO
        END DO 


        DO K=1,M$
        DO J=1,N$
           KYY(J,K)  =( PWKYY(J,  K) + PWKYY(J+1,K+1) +
     >                  PWKYY(J+1,K) + PWKYY(J,  K+1)   
     >                                         )*0.25           

           IF( DUM2(J,K) .LT. 1.501E-14) KYY(J,K)=0.000  ! RESTRAIN
           IF(  KYY(J,K) .LT. 0.00)      KYY(J,K)=0.000  ! KYY
           IF(  KYY(J,K) .GT. 5.E12)     KYY(J,K)=5.E12  ! KYY
           IF(  UBX(J,K) .LT. 200.)      KYY(J,K)=0.000

cjer add lower limit to tropospheric kyy
           if(zp(k).lt.10. .and. kyy(j,k).lt.1.5e10) kyy(j,k)=1.5e10

        END DO
        END DO

        RETURN
        END

        SUBROUTINE EPD1(PHI,M,EPDV,EPFY,EPFZ,NS,MS)

	INCLUDE 'PARAM.INC'
	INCLUDE 'COMMONC.INC'
	INCLUDE 'COMMOND.INC'
	INCLUDE 'COMMONW.INC'
	INCLUDE 'COMMONT.INC'

        REAL    EPFY2(N$,M$),EPFZ2(N$,M$) 
        REAL    EPFY (NS,MS),EPFZ (NS,MS) 
        REAL    EPDV (NS,MS)
        COMPLEX DPHIDX(N$,M$),DPHIDY(N$,M$),DPHIDZ(N$,M$) 
        COMPLEX PHI(NS,MS),AI
C	DATA    WERROR/.false./,ai/(0.,1.)/
        LOGICAL WERROR = .false.
	DATA    ai/(0.,1.)/


	CALL CDERV4(1,phi,DPHIDZ,DZ,n$,M$,2)

	DO 215 I=1,n$
	  DPHIDZ(i,m$) = 0. !phicr(i-1,m$)
	  DPHIDZ(i,1 ) = (PHI(I,2)-PHI(I,1))/(DZ)
215	CONTINUE

C	COMPUTE U', V'


	CALL CDERV4(1,PHI,DPHIDY,dy,n$,M$,1)


	DO 228 J=1,M$
	  DPHIDY(1,J) = 0.5*PHI(2,J)/DY
	  DPHIDY(n$,J) = -0.5*PHI(n1$,J)/DY
 228	CONTINUE

	DO 230 J=1,M$
	DO 230 I=1,n$

C********************************************************
C                                                       *
C         DPHIDX(i,j) =  ai*m*phi(i,j)/(A*C(I))         *
C               ---- OK FOR REALLY,REALLY QG SYSTEM,    *
C               ---- I.E. A BETA PLANE,                 *
C               ---- BUT WE REALLY WANT                 *
C               ---- D(PHI)/D(LAMBDA).                  *
C               ---- IN A SPHERICAL GEOMETRY.  SO, LET: *
C                                                       *
	  DPHIDX(i,j) =  ai*m*phi(i,j)/(A)        !     *     
C                                                       *
C               ---- KEEP `A' IN DENOMINATOR FOR        *
C               ---- DIMENSIONAL REASONS                *
C               ---- C.F. AHL, EQ. 5.2.5, p. 231        *
C                                                       *
C********************************************************

 230	CONTINUE

C -- EP FLUX CALCULATION --

C
C
C   -- ALMOST QG APPROXIMATION OF EP-FLUX
C   -- (SPHERICAL NOT PLANAR, DEPENDING ON
C   -- ACTUAL DEFINITION OF DPHIDX)


	CALL ZMEAN_SUB(DPHIDX,DPHIDY,epfy,n$,M$)
	call zmean_sub(DPHIDX,DPHIDZ,epfz,n$,m$)

c	note the additional EP terms are neglected for now

	DO J=1,M$
	DO I=1,n$
	   epfy2(i,j)= epfy(i,j)/(CF(I)*CF(I))   !*rho(j)
	   epfz2(i,j)= epfz(i,j)/( xn2(J) )
        END DO
        END DO




        DO K=1,M$
        DO J=1,N$
           EPFZ(J,K)=
     >        EPFZ2(J,K)*RHO(K)
           EPFY(J,K)= 
     >        EPFY2(J,K)*RHO(K)*C(J)
        END DO
        END DO

        DO K=2,M$-1
        DO J=2,N$-1
           EPDV(J,K)=
     >        ( EPFY(J+1,K)-EPFY(J-1,K) )/( 2.*DY*C(J) )
     >     +  ( EPFZ(J,K+1)-EPFZ(J,K-1) )/( 2.*DZ )
        END DO
        END DO
C                        set bottom and top = 0. (OK)
        DO J=1,N$
           EPDV(J,M$)=0.00
           EPDV(J, 1)=0.00
        END DO

        DO K=1,M$
           EPDV(1, K)=0.00
           EPDV(N$,K)=0.00
        END DO

        DO K=2,M$-1
        DO J=2,N$-1
           EPDV(J,K)=EPDV(J,K)/( C(J)*RHO(K) )
        END DO
        END DO

C                  reset 85S, 85N = 80S,N (instead of 0 - EF, June 2009)
        DO K=1,M$
           EPDV(1, K) = EPDV(2, K)
           EPDV(N$,K) = EPDV(N$-1,K)
        END DO

c        CALL PRMINMAX(EPDV,  'EPDV_EPD1-',N$, M$,1)
c        CALL PRMINMAX(PHI ,  'PHI__EPD1-',N$, M$,1)

      RETURN
      END


     
        SUBROUTINE EPD2(PHI,M,EPDV,EPFY,EPFZ,NS,MS)

	INCLUDE 'PARAM.INC'
	INCLUDE 'COMMONC.INC'
	INCLUDE 'COMMOND.INC'
	INCLUDE 'COMMONW.INC'
	INCLUDE 'COMMONT.INC'

        REAL    EPFY2(N$+1,M$),EPFZ2(N$,M$+1) 
        REAL    EPFY (NS+1,MS),EPFZ (NS,MS+1) 
        REAL    EPDV (NS,MS)
        COMPLEX DPHIDX(N$,M$),DPHIDY(NP$,M$),DPHIDZ(N$,MP$) 
        COMPLEX PHI(NS,MS),AI,  PHIY(NP$,M$),  PHIZ(N$,MP$)
C	DATA WERROR/.false./,ai/(0.,1.)/
        LOGICAL WERROR = .false.
        DATA    ai/(0.,1.)/

C
C            //         |         |         |         |
C            /X----&----X----&----X----&----X----&----X--
C            //         |         |         |         |
C            //         |         |         |         |               
C            //    #    O    #....O....#    O    #    O
C            //         |    :    |    :    |         |
C            //         |    :    |    :    |         |
C            /X----&----X----&----X----&----X----&----X---
C            //         |    :    |    :    |         |
C            //         |    :    |    :    |         |               
C            //    #    O    #....O....#    O    #    O
C            //         |         |         |         |
C            //         |         |         |         |
C            /X----&----X----&----X----&----X----&----X---
C            //         |         |         |         |
C            //         |         |         |         |               
C            //    #    O    #    O    #    O    #    O
C            //         |         |         |         |
C            //         |         |         |         |
C            /X/////////X/////////X/////////X/////////X////
C
C           (90S, OKM)
C
C
C                      X-points ----> PHI, DPHIDX, EPDV (MRS grid)
C           
C                      &-points ----> DPHIDY, EPFY
C           
C                      O-points ----> DPHIDZ, EPFZ
C
C                      #-points ----> UBX, etc. (RBR grid)
C


	DO 215 J=2,M$
	DO 215 I=1,N$
	  DPHIDZ(I,J ) = (PHI(I,J)-PHI(I,J-1))/(DZ)
215	CONTINUE

        DO I=1,N$
           DPHIDZ(I,1 ) =DPHIDZ(I,2 ) 
           DPHIDZ(I,MP$) =DPHIDZ(I,M$) 
        END DO  

	DO 228 J=1,M$
	DO 228 I=2,N$
	  DPHIDY(I,J ) = (PHI(I,J)-PHI(I-1,J))/(DY)
228	CONTINUE

        DO J=1,M$
           DPHIDY(1, J) =DPHIDY(2, J) 
           DPHIDY(NP$,J) =DPHIDY(N$,J) 
        END DO  



	DO 230 J=1,M$
	DO 230 I=1,n$

C********************************************************
C                                                       *
C         DPHIDX(i,j) =  ai*m*phi(i,j)/(A*C(I))         *
C               ---- OK FOR REALLY,REALLY QG SYSTEM,    *
C               ---- I.E. A BETA PLANE,                 *
C               ---- BUT WE REALLY WANT                 *
C               ---- D(PHI)/D(LAMBDA).                  *
C               ---- IN A SPHERICAL GEOMETRY.  SO, LET: *
C                                                       *
	  DPHIDX(i,j) =  ai*m*phi(i,j)/(A)        !     *     
C                                                       *
C               ---- KEEP `A' IN DENOMINATOR FOR        *
C               ---- DIMENSIONAL REASONS                *
C               ---- C.F. AHL, EQ. 5.2.5, p. 231        *
C                                                       *
C********************************************************

 230	CONTINUE

C -- EP FLUX CALCULATION --

C
C
C  -- ALMOST QG APPROXIMATION OF EP-FLUX
C  -- (SPHERICAL NOT PLANAR, DEPENDING ON
C  -- ACTUAL DEFINITION OF DPHIDX)


        
	DO J=1,M$
	DO I=2,n$
           PHIY(I,J)=( DPHIDX(I,J)+DPHIDX(I-1,J) )/2.0
        END DO
        END DO
	DO J=1,M$
           PHIY(NP$,J)=DPHIDX(N$,J)
           PHIY(1  ,J)=DPHIDX(1 ,J)
        END DO
	CALL ZMEAN_SUB(PHIY,DPHIDY,epfy,NP$,M$)

	DO J=2,M$
	DO I=1,n$
           PHIZ(I,J)=( DPHIDX(I,J)+DPHIDX(I,J-1) )/2.0
        END DO
        END DO

	DO I=1,N$     ! J=1,M$; - corrected, EF, 7/09 ;  PHIZ(N$,MP$); DPHIDX(N$,M$)
           PHIZ(I,MP$)=DPHIDX(I,M$)
           PHIZ(I,1  )=DPHIDX(I,1 )
        END DO

	call ZMEAN_SUB(PHIZ,DPHIDZ,epfz,N$,MP$)

c  note the additional EP terms are neglected for now

	DO J=1,M$
	DO I=1,nP$
           IF(CFST(I) .NE. 0.0 ) THEN
	     epfy2(i,j)= epfy(i,j)/(CFST(I)*CFST(I))   !*rho(j)
           ENDIF
           IF(CFST(I) .EQ. 0.0 ) THEN
	     epfy2(i,j)= 0.00
           ENDIF
        END DO
        END DO

	DO J=1,MP$
	DO I=1,n$
	   epfz2(i,j)= epfz(i,j)/( xn2ST(J) )
        END DO
        END DO




        DO K=1,M$
        DO J=1,NP$
           EPFY(J,K)= 
     >        EPFY2(J,K)*RHO(K)*CST(J)
        END DO
        END DO

        DO K=1,MP$
        DO J=1,N$
           EPFZ(J,K)=
     >        EPFZ2(J,K)*RHOST(K)
        END DO
        END DO

        DO K=1,M$
        DO J=1,N$
           EPDV(J,K)=
     >        ( EPFY(J+1,K)-EPFY(J,K) )/( 1.*DY*C(J) )
     >     +  ( EPFZ(J,K+1)-EPFZ(J,K) )/( 1.*DZ )
        END DO
        END DO

        DO J=1,N$
           K=1
           EPDV(J,K)=0.
        END DO

        DO K=1,M$
        DO J=1,N$
           EPDV(J,K)=EPDV(J,K)/( C(J)*RHO(K) )
        END DO
        END DO


c        CALL PRMINMAX(EPFY,  'EPFY_EPD2-',NP$, M$,1)
c        CALL PRMINMAX(PHIY,  'PHIY_EPD2-',NP$, M$,1)
c        CALL PRMINMAX(DPHIDY,'DPDY_EPD2-',NP$, M$,1)
c        CALL PRMINMAX(EPFZ,  'EPFZ_EPD2-',N$, MP$,1)
c        CALL PRMINMAX(PHIZ,  'PHIZ_EPD2-',N$, MP$,1)
c        CALL PRMINMAX(DPHIDZ,'DPPZ_EPD2-',N$, MP$,1)
c        CALL PRMINMAX(EPDV,  'EPDV_EPD2-',N$, M$,1)
c        CALL PRMINMAX(PHI ,  'PHI__EPD2-',N$, M$,1)

      RETURN
      END




        SUBROUTINE MATSQP(PHI,M,MSS,NS,MS)   !,QPRIME,QPRMY)


        INCLUDE 'PARAM.INC'
        INCLUDE 'COMMONC.INC'
        INCLUDE 'COMMOND.INC'
        INCLUDE 'COMMONW.INC'
        INCLUDE 'COMMONT.INC'
        COMPLEX PHI(NS,MS)     !,QPRIME(NS,MS),QPRMY(NS,MS)
        COMPLEX DPHIDX(N$,M$),DPHIDY(N$,M$),DPHIDZ(N$,M$) 
        COMPLEX AI,  PHIY(N$,M$),  PHIZ(N$,M$)



C****---------------------------------------------
C**** D()/DXX TERM -------------------------------

        DO K=1,MS
        DO J=1,NS
           DPHIDX(J,K)=( -M*M*PHI(J,K) )/( CF(J)*(A*C(J))**2  )
        END DO
        END DO
C---------------------------------------------****
C---------------------------------------------****



C****---------------------------------------------
C**** D()/DZZ TERM -------------------------------

	CALL CDERV4(1,phi,DPHIDZ,DZ,n$,M$,2)

	DO I=1,n$
	  DPHIDZ(i,m$) = 0. !phicr(i-1,m$)
	  DPHIDZ(i,1 ) = (PHI(I,2)-PHI(I,1))/(DZ)
        END DO
 
        DO K=1,M$
        DO J=1,N$
           PHIZ(J,K)= ( RHO(K)*DPHIDZ(J,K)/XN2(K) )
        END DO
        END DO

	CALL CDERV4(1,phiZ,DPHIDZ,DZ,n$,M$,2)

	DO I=1,n$
	  DPHIDZ(i,m$) = 0. 
	  DPHIDZ(i,1 ) = (PHIZ(I,2)-PHIZ(I,1))/(DZ)
        END DO

        DO K=1,M$
        DO J=1,N$
           DPHIDZ(J,K)=  CF(J)*DPHIDZ(J,K)/RHO(K) 
        END DO
        END DO
C-------------------------------------------------****


C****-------------------------------------------------
C*** D()/DYY TERM ------------------------------------

	CALL CDERV4(1,PHI,DPHIDY,dy,n$,M$,1)

	DO J=1,M$
	  DPHIDY(1,J) = 0.5*PHI(2,J)/DY
	  DPHIDY(n$,J) = -0.5*PHI(n1$,J)/DY
        END DO

        DO K=1,MS
        DO J=1,NS
           PHIY(J,K)=( DPHIDY(J,K)*C(J) )/( CF(J)**2  )
        END DO
        END DO

	CALL CDERV4(1,PHIY,DPHIDY,dy,n$,M$,1)

	DO J=1,M$
	  DPHIDY(1,J) = 0.5*PHIY(2,J)/DY
	  DPHIDY(n$,J) = -0.5*PHIY(n1$,J)/DY
        END DO

        DO K=1,MS
        DO J=1,NS
           DPHIDY(J,K)=( DPHIDY(J,K)*CF(J) )/( C(J) )
        END DO
        END DO
C-------------------------------------------------****
C-------------------------------------------------****
        

C****-------------------------------------------------
C*** ADD THREE SECOND DER. TERMS TO GIVE QPRIME ------
C*** AND TAKE Y-DERIVATIVE ---------------------------

        DO K=1,MS
        DO J=1,NS
           QPRIME(J,K) = DPHIDX(J,K) + DPHIDY(J,K) + DPHIDZ(J,K)
        END DO
        END DO

	CALL CDERV4(1,QPRIME,QPRMY,dy,n$,M$,1)

	DO J=1,M$
	  QPRMY(1,J) = 0.5*QPRIME(2,J)/DY
	  QPRMY(n$,J) = -0.5*QPRIME(n1$,J)/DY
        END DO
C-------------------------------------------------****


C****-------------------------------------------------
C*** MASK OUT GARCIA DISSIPATION WHERE QPRMY IS ------
C*** LESS THAN QBRYW, I.E. WHERE PWAVE IS NOT 
C*** SATURATED   --
C*** MODIFIED to be where QPRMY < 2*QBRYW, per RG94, eq. 19 (EF, 4/09)
C*** -----------------------------------------------------------------
        DO K=1,M$
        DO J=1,N$
           DUM1(J,K) = 0.00
           SAT_PARAM = ( CABS( QPRMY(J,K) )/QBRYW(J,K) ) - 2.00
           SAT_PARAM = AMAX1( SAT_PARAM, 0.00 )  ! SAT_PARAM --> GE 0.0
           SAT_PARAM = AMIN1( SAT_PARAM, 1.00 )  ! SAT_PARAM --> LE 1.0
           IF(SAT_PARAM.GT.0.00) DUM1(J,K)=GDISS(J,K,MSS)*SAT_PARAM
           IF(DUM1(J,K).LT.0.00) DUM1(J,K)=0.00
           GDISS(J,K,MSS)=DUM1(J,K)
        END DO
        END DO

        IF(M.EQ.1) IQUA=46
        IF(M.GT.1) IQUA=46+10000+M*100
c      CALL OUTA(IQUA,NSTEP,NDAY,NYEAR,YP,ZP,GDISS(1,1,mss), N$,M$,0)



      RETURN
      END


        SUBROUTINE GARCIA(PHI,M,MSS)

C*************************************************************
C*   SUBROUTINE TO CALCULATE DISSIPATION NECESSARY TO        *
C*   MAINTAIN PLANETARY WAVES AT SATURATION.                 *
C*   REFS. (GARCIA, R.R., JAS, 48, 1405;  GARCIA ET AL.,     *
C*   JGR, 97, 12967)                                         *
C*************************************************************

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONW.INC'
      INCLUDE 'COMMONT.INC'
      REAL     RK(N$,M$),RL(N$,M$),RM(N$,M$)
      COMPLEX  PHI(N$,M$),DPHIDY(N$,M$),DPHIDZ(N$,M$),PHI_X

C
C  ---        FIRST CALCULATE WAVE NUMBERS   
C


      DO K=1,M$
      DO J=1,N$
 
         RK(J,K)= M / ( A*C(J) )

      END DO
      END DO

      CALL CDERV4(1,PHI,DPHIDY,DY,n$,M$,1)

      CALL CDERV4(1,phi,DPHIDZ,DZ,n$,M$,2)


      DO K=1,M$
      DO J=1,N$
         PHI_X=PHI(J,K)                          !+CMPLX( 1.E-6, 1.E-6 )

c         PHI_X=CMPLX( AMAX1( REAL(PHI_X),1.E-6 ) ,
c     >                AMAX1(AIMAG(PHI_X),1.E-6 )   )

   
         IF( CABS(PHI_X).LT. 1.E-6 ) PHI_X=CMPLX( 1.E-6, 1.E-6 )
         
         RL(J,K)= AIMAG( DPHIDY(J,K)/PHI_X )    ! MERID. WAVE #
         RM(J,K)= AIMAG( DPHIDZ(J,K)/PHI_X )    ! VERTICAL WAVE #

      END DO
      END DO


        IF(M.EQ.1) IQUA=41
        IF(M.GT.1) IQUA=41+10000+M*100
c        CALL OUTA(IQUA,NSTEP,NDAY,NYEAR,YP,ZP,RL , N$,M$,0)
        IF(M.EQ.1) IQUA=42
        IF(M.GT.1) IQUA=42+10000+M*100
c        CALL OUTA(IQUA,NSTEP,NDAY,NYEAR,YP,ZP,RM , N$,M$,0)


C  ---   CALCULATE GROUP VELOCITIES IN Y AND Z
C


      DO K=1,M$
      DO J=1,N$

         EPS   =  CF(J)**2 / XN2(K) 
         RBIGK2=  RK(J,K)**2 + RL(J,K)**2 + EPS*
     >          ( RM(J,K)**2+1./(4*H**2) )  

         DUM1(J,K)= 2.0 * RK(J,K) * RL(J,K) * UBW(J,K) 
     >                 / RBIGK2                   ! ------- GY

         DUM2(J,K)= 2.0 * RK(J,K) * EPS * RM(J,K) * UBW(J,K) 
     >                 / RBIGK2                   ! ------- GZ

      END DO
      END DO


      CALL DERV4B(1,UB(1,1,1),UBY,DY,n$,M$,1)
      CALL DERV4B(1,UB(1,1,1),UBZ,DZ,n$,M$,2)

      DO K=1,M$
      DO J=1,N$
         GDISS(J,K,MSS)=
     >    DUM1(J,K)*( (-5./4.)*( UBY(J,K)/UBW(J,K) ) ) +
     >    DUM2(J,K)*( (-3./4.)*( UBZ(J,K)/UBW(J,K) ) + 1./(2*H) )
      END DO
      END DO

        IF(M.EQ.1) IQUA=45
        IF(M.GT.1) IQUA=45+10000+M*100
c      CALL OUTA(IQUA,NSTEP,NDAY,NYEAR,YP,ZP,GDISS(1,1,MSS), N$,M$,0)

        IF(M.EQ.1) IQUA=43
        IF(M.GT.1) IQUA=43+10000+M*100
c      CALL OUTA(IQUA,NSTEP,NDAY,NYEAR,YP,ZP,DUM1, N$,M$,0)
        IF(M.EQ.1) IQUA=44
        IF(M.GT.1) IQUA=44+10000+M*100
c      CALL OUTA(IQUA,NSTEP,NDAY,NYEAR,YP,ZP,DUM2, N$,M$,0)

      RETURN
      END



      SUBROUTINE BKGDMP

C*************************************************************
C*   SUBROUTINE TO CALCULATE BACKGROUND DISSIPATION          *
C*   FOR PLANETARY WAVES                                     *
C*************************************************************

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONW.INC'
      INCLUDE 'COMMONT.INC'





C***BACKGROUND DAMPING***
C---
        DMPMX   = 1./DAYL
        EQDMP0  = ISW(54)/100.

        PWVDMP=0.000 

        IF(ISW(49).GT.0) then 
          PWVDMP=( 1.00/(DAYL*ISW(49)) )
        endif
C************************


        DO K=1,M$
        DO J=1,N$
  
           EQDMP(J,K)=EQDMP0*COS( YP(J)*PI/60. ) / DAYL

           IF( ABS(YP(J)) .GE. 30. ) EQDMP(J,K)=0.         

        END DO
        END DO


        PWDMPZ=1.*ISW(57)
        DMPTHK=1.*ISW(58)

        DO K=1,M$
        DO J=1,N$
  
           IF(ZP(K) .GT. PWDMPZ) THEN
             EQDMP(J,K)=EQDMP(J,K)+ (ZP(K)-PWDMPZ)/(DMPTHK*DAYL)
           ENDIF

           EQDMP(J,K) = AMIN1( EQDMP(J,K), DMPMX ) 
           EQDMP(J,K) = AMAX1( EQDMP(J,K), PWVDMP ) 

        END DO
        END DO


       
        RETURN
        END


      SUBROUTINE ENHDMP(MSS)

C*************************************************************
C*   SUBROUTINE TO MAKE ENHANCED DISSIPATION FOR             *
C*   PLANETARY WAVES AT CL'S AND WAVE BREAKING ZONES         *
C*************************************************************

      INCLUDE 'PARAM.INC'
      INCLUDE 'COMMONC.INC'
      INCLUDE 'COMMOND.INC'
      INCLUDE 'COMMONW.INC'
      INCLUDE 'COMMONT.INC'


        DMPMX   = 1./DAYL
        EQDMP0  = ISW(54)/100.

 
C**************************************************
C       ENHANCE DISSIPATION NEAR CRITICAL LINES 
        IF ( ISW(61) .GE. 1 )  THEN
c        WRITE(6,*) '  >>>> \/\/\/ ENHANCED CL DAMPING '
        DO K=1,M$
        DO J=1,N$
           IF( UBW(J,K) .LT. 400. ) THEN 
               EQDMP(J,K)=DMPMX
           END IF
        END DO
        END DO
        END IF
c----
C       SMOOTH ENHANCED CL DISSIPATION 
        IF ( ISW(61) .GE. 2 )  THEN
           DO NSMTH=1,ISW(61)-1
           CALL SMTHR3( EQDMP, N$, M$ )
c           WRITE(6,*) '  >>>> \/\/\/ SMOOTHING CL DAMPING ',NSMTH
           END DO
        END IF
C***************************************************

C***************************************************
C ENHANCE DISSIPATION WHERE |D(q')/Dy| > D(QBAR)/Dy

        IF ( ISW(65) .GE. 1 )  THEN
c        WRITE(6,*) '  >>>> \/\/\/ ENHANCED |QY| DAMPING '
c----
C       SMOOTH ENHANCED QY DISSIPATION 
        IF ( ISW(65) .GE. 2 )  THEN
           DO NSMTH=1,ISW(65)-1
           CALL SMTHR3( GDISS(1,1,MSS), N$, M$ )
c           WRITE(6,*) '  >>>> \/\/\/ SMOOTHING |QY| DAMPING ',NSMTH
           END DO
        END IF

        DO K=1,M$
        DO J=1,N$
           EQDMP(J,K) = EQDMP(J,K) + ISW(66)*GDISS(J,K,MSS)
        END DO
        END DO

        IF(MWVNS(1,MSS).EQ.1) IQUA=53
        IF(MWVNS(1,MSS).GT.1) IQUA=53+10000+MWVNS(1,MSS)*100
c        CALL OUTA(IQUA,NSTEP,NDAY,NYEAR,YP,ZP,EQDMP, N$,M$,0)
c
        END IF

C***************************************************


        RETURN
        END
