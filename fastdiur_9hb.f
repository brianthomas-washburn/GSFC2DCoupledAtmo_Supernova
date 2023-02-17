C **********************************************************************
C
cccc       SUBROUTINE FASTCHAER(DTIMEB,TFD,L$,Z$,S$,ikmes,JAVAER,JAER,CNDCA)
       SUBROUTINE FASTCHAER(DTIMEB,TFD,L$,Z$,S$,ikmes,JAVAER,JAER,CNDCA
     >   ,iday360, iyr, gaerosol, gnataer, gaerice, PH$, IL$, JZ39, J39)

C
C **********************************************************************
C             @(#)fastdiur.f	1.57  05/08/00 
C        DIURNAL VERSION OF FAST CHEMISTRY
C

      include  'com_aerg.h'
      include  'com_aers.h'
      include  'com_aerd.h'


      INTEGER L$, Z$, S$, ikmes, PH$, IL$


C  Common for mapping into GSFC code
C      CNDCA = CN array for 18+3 diurnal time steps for up to 1 daily iterations, NFSP=41
C  
      COMMON/CAER1/ISPMAP(NFSP)

      REAL TFD(L$)
      REAL*8 JAVAER(NJR+3,L$,Z$), CNDCA(S$,ntime+3,L$,Z$)
      REAL*8 JAER(NJR+3,NTIME,L$,Z$), xhcl, xclono2, OXG(L$,Z$)

C  wavelength dependent JZ39 array (for diagnostics ONLY), used in TIMEAVJ
C            J39 is diurnally averaged output from TIMEAVJ

      REAL*8 JZ39(18,PH$,5,L$,Z$), J39(PH$,5,L$,Z$)


      REAL gaerosol(l$,z$), gnataer(l$,z$), gaerice(l$,z$)
      REAL gfaer, gfnat, gfice

cj
Cj  AJOUT is just a subsample array for output (AER's J1 and J2 at 30 km)
cj
cj      COMMON/JAER1/AJOUT(2,ntime,18)
cj     >             AJAVDAY(3,NJR,18,46)
cj

      REAL*8 JEX
      COMMON/JEXAER/JEX(3)

      DOUBLE PRECISION SPNOON
      COMMON/SP/SPNOON(NTSP)

      DOUBLE PRECISION ALLFSP, ALLSSP, ALLPLQ, ALLJR
      COMMON /SPEC/ALLFSP(NLT,NHT,NFSP), ALLSSP(NLT,NHT,NSSP),
     $     ALLPLQ(NLT,NHT,NPLQ), ALLJR(NLT,NHT,NJR)

      PARAMETER (MLT=(NLT)/2)
      PARAMETER (IFSSP=NFSP+1,ILSSP=NCSP)
C
      LOGICAL  LASTYEAR,LSTOUT,LPRT,LPLAT
      LOGICAL  DOUTFS,NOUTFS,NOUTSS,FOUT,DIURNL
      LOGICAL  LTDBC,LTDFL,LSEASBC,LINTP,LHSCT,LAERO,LTDIST,LPSC,LCONV
      INTEGER  LOLA(50),PRTFSP(NFSP),PRTSSP(NSSP)
C
      REAL COSLAT(NLT2),SINLAT(NLT2),TANLAT(NLT),EE(NHT2)
      REAL DLITE(NHT)

      COMMON /AIR/TMP(NLT,NHT),PR(NLT,NHT),AIRVD(NLT,NHT)
c      REAL TMP(NLT,NHT),PR(NLT,NHT),AIRVD(NLT,NHT)

C

C         ALLFSP CONTAINS ALL FAST SPECIES MIXING RATIOS
c      REAL ALLFSP(NLT,NHT,NFSP)

      REAL FSPDAY(NLT,NHT,NFSP)
      REAL FSPNIT(NLT,NHT,NFSP)
      REAL FSPAVG(NLT,NHT,NFSP)
      REAL FSPRIS(NLT,NHT,NFSP)
      REAL FSPSET(NLT,NHT,NFSP)
      REAL FSPMID(NLT,NHT,NFSP)

C         ALLSSP CONTAINS ALL SLOW SPECIES MIXING RATIOS
c      REAL ALLSSP(NLT,NHT,NSSP)

C         PHOTOCHEMICAL RATES:  PRODUCTION, LINEAR AND QUADRATIC LOSS
C         ALL IN ONE ARRAY FOR 12 SPECIES

c      REAL ALLPLQ(NLT,NHT,NPLQ)
C         J-RATES FOR ALL SPATIAL POINTS

c      REAL ALLJR(NLT,NHT,NJR)

      REAL ALLJRNOON(NLT,NHT,NJR),ALLJRDAY(NLT,NHT,NJR)
C         ARRAY OF TROPOPAUSE HEIGHTS
      INTEGER ITROP(NLT)
C          TIME INTERVAL USED BY NWT2DD.
      DOUBLE PRECISION DTIMEA(NTIME),DTIMEB(L$,NTIME)
c      DOUBLE PRECISION DTIME(NTIME)
      DIMENSION PMTIME(NTIME)
C         ARRAY OF ALL SPECIES USED BY NWT2DD, NUMBER DENSITY
      DOUBLE PRECISION ALLSP(NTSP)

C          ARRAY OF ALL FAST SPECIES AROUND THE CLOCK AT ONE SPATIAL POINT
      DOUBLE PRECISION FSP(NFSP,NTIME+3)
C     DOUBLE PRECISION TOT(NFSP,NTIME+3)
C         J-RATES AND K-RATES AT ONE SPACIAL POINT
      DOUBLE PRECISION JR(NJR,NTIME+3), JRA(NJR,NTIME), JEXR(3,NTIME)
      DOUBLE PRECISION JRMOM(NJR),KR(NKR)
C         ARRAY OF PRODUCTION AND LOSS RATES AT ONE SPATIAL POINT
      DOUBLE PRECISION PRLOC(NCSP),LOLOC(NCSP),QLLOC(NCSP)
C         NWT2DD ERROR TOLERANCE ARRAY
      DOUBLE PRECISION EI1,EIPER
C         LOSS TERM POWER  FOR SPECIE NUMBER I
C         LTPW(I)=1 IF LINEAR
C                =2 IF QUADRATIC
      DIMENSION LTPW(NRSP),INDBC(NRSP)
C         NWT2DD ERROR EXPLANATIONS
      DIMENSION NWERR(3)
      DIMENSION PU(NHT),ZU(NHT),TU(NHT),DU(NHT),SU(NHT),O3(NHT)
C
      CHARACTER*24 TXTH2O,NWERR
      CHARACTER*24 SPID,PLQID,QID,JRID,REID
C
      COMMON /CAPTN/MODNO,NODAY,NSTOP,IYEAR,IDOFY
      COMMON /ALTITUDES/ZSTAR(NHT),ZZALT(NLT,NHT),ALAT(NLT)
      COMMON /SULFATE/SURFACE(NLT,NHT),AMASS(NLT,NHT),AH2SO4(NLT,NHT)
      COMMON /PSC/PSC1(NLT,NHT),SPSC1(NLT,NHT),VSED1(NLT,NHT1),
     $     PSC2(NLT,NHT),SPSC2(NLT,NHT),VSED2(NLT,NHT1)

c      COMMON /AIR/TMP,PR,AIRVD

C         J-RATES AND K-RATES IDS
      COMMON /QID/QID(NJR),JRID(NJR)
      COMMON /REAC/REID(NKR)
C         SPECIES NAMES
      COMMON /SID/SPID(NTSP),PLQID(NPLQ)

C         MIXING RATIOS FOR ALL CALCULATED SPECIES
c      COMMON /SPEC/ALLFSP,ALLSSP,ALLPLQ,ALLJR

      COMMON /JRDIUR/ALLJRNOON,ALLJRDAY
      COMMON /DEBUG/LOLA,PRTFSP,PRTSSP
      COMMON /DIUR/FSPDAY,FSPNIT,FSPAVG,FSPRIS,FSPSET,FSPMID
      COMMON /DIPARM/MINLAT,MAXLAT,MINALT,MAXALT,
     *               DOUTFS,NOUTFS,NOUTSS,FOUT,DIURNL
C         COMMON BLOCKS TO INTERFACE JACOB
      COMMON/PROD/PRLOC/LOSS/LOLOC/QLOS/QLLOC
      COMMON/JRATE/JRMOM/KRATE/KR/SP/ALLSP
      COMMON/SRATES/SNOXPR(NLT,NHT),SPANPR(NLT,NHT),SCLXPR(NLT,NHT),
     * SHCLPR(NLT,NHT)
      COMMON/CLOSURE/IFAM(NFAM),ICF(NFAM-1,NFP),ISF(NFSP,NFAM-1)
C         MODEL DIMENSIONS AND PARAMETERS SET IN SUBROUTINE PARAMS
      COMMON /COORDS/DH,DK,COSLAT,SINLAT,TANLAT,EE
C         SLOW SPECIES LOSS TERM POWER, BOUNDARY CONDITIONS
      COMMON /TREAT/LTPW,INDBC
      COMMON /LOGIC/LTDBC,LTDFL,LSEASBC,LINTP,LHSCT,LAERO,LTDIST,LPSC,
     $     LCONV
C         USED TO INTERPOLATE SPECIES FROM FILE
      COMMON /SPINT/DSPINT(NLT,NHT,NSPINT),LSTPSP(NSPINT),MRECSP(NSPINT)
      COMMON/NOZMIN/ANOZMIN(NLT,NHT)
C         ATMOSPHERIC PROFILE FOR UCI J-RATE CODE
      COMMON/FLAGS/IFL(90)
C               NEWTON ITERATIONS EVERY TIME STEP:
      DATA ITFS/200/

C
C    NWT2DD ERROR TOLERANCE ARRAY    
C
C    NOTE: after MUCH TESTING, it was determined that an ERROR TOLERANCE of 1.E-5 was best 
C          this gave as much as 5% difference in ozone from 1.E-4 in the lower stratosphere
C          but an ERROR TOLERANCE of 1.E-5 was, at most, 0.5% different than 1.E-6
C          200 iterations are used (in NWT2DD), but it seems that it usually only takes < 40 in most places
C
      DATA EI1/1.E-5/
cccc      DATA EI1/1.E-4/


      DATA EIPER/.005/
C
C ......................................................................
C                    EXECUTION BEGINS
C
C
C        GET THE SLOW SPECIES BY INTERPOLATION IF REQUIRED
C
C
C
C           GET AEROSOL SURFACE AREA FOR HETEROGENEOUS REACTIONS
ccc      IF(.NOT. LAERO) CALL GETSURF(SURFACE,IDOFY,LOLA(8))
C
C           GET PSC SURFACE AREA AND EQUILIBRIUM HNO3, H2O
ccc      IF(LPSC) CALL CALCPSC
C
C
C
C
C ++++++++++++++ START THE MAIN FAST CHEMISTRY LOOP ++++++++++++++++++++
C
C          LOOP OVER LATITUDES, GOING FROM MLT TO 1, THEN MLT+1 TO NLT


cc        print *, ' '
cc        print *, ' IN FASTDIUR, L$, TFD = ', L$
cc        print *, tfd
cc        print *, ' '
cc        print *, ' '


      DO 900 I=1,NLT
c      DO 900 I=4,4
C
         do J=1,NHT
            PU(J)=PR(I,J)
            ZU(J)=ZZALT(I,J)
            TU(J)=TMP(I,J)
            DU(J)=AIRVD(I,J)

c            SU(J)=SURFACE(I,J)
c            O3(J)=ALLSSP(I,J,IO3-NFSP)

         end do

c         print *,' pu=',(pu(j),j=1,nht)
c         print *,' zu=',(zu(j),j=1,nht)
c         print *,' tu=',(tu(j),j=1,nht)
c         print *,' du=',(du(j),j=1,nht)
C
C
C+++++++++ LOOP OVER ALTITUDES, FROM GROUND TO TOP OF ATMOSPHERE
      DO 800 J=1,NHT
C
C         AIR VOLUME DENSITY
      AIRM=AIRVD(I,J)

cpp         if(I .eq. 9)then
cpp            if(J .eq. 11 .or. J .eq. 21 .or. J .eq. 31)then
cpp      print *,' i=',i,' j=',j,' airm=',airm
cpp            endif
cpp            endif

C   also load in the diurnally varying extra Js  into JEXR(3,NTIME), JEX(3), first initialize

       do 4400 iclock=1,ntime
       do 4400 ijr=1,3
 4400     JEXR(ijr,iclock) = 0.0D0


C    NOTE:  SUNSET=7 and SUNRIS=12 ALWAYS (at ALL latitudes, days of the year) so set J's = 0 at night
C                                     although this should already be done in PHOTHEATIN initialization
       DO 5500 ICLOCK=1,18   ! SUNSET

           CALL JMAP(I,J,ICLOCK)

           DO 6050 IJR=1,NJR
             JR(IJR,ICLOCK)=JRMOM(IJR)
             if (iclock .ge. 8 .and. iclock .le. 11) JR(IJR,ICLOCK)=0.d0
 6050      CONTINUE

            DO 7050 IJR=1,3
              JEXR(IJR,ICLOCK) = JEX(IJR)
           if (iclock .ge. 8 .and. iclock .le. 11) JEXR(IJR,ICLOCK)=0.d0
 7050      CONTINUE


C   also do for diagnostic array JZ39(18,PH$,5,L$,Z$)

          if (iclock .ge. 8 .and. iclock .le. 11) then
             do 7175 ikk0=1,5
             do 7175 ijj0=1,PH$
 7175           JZ39(iclock,ijj0,ikk0,i,j) = 0.d0
          endif

 5500   CONTINUE


cccol      DO IJR=1,NJR
cccol      DO ICLOCK=1,SUNSET
cccol         JCLOCK=NTIME-ICLOCK+1
cccol         JR(IJR,JCLOCK)=JR(IJR,ICLOCK)
cccol      END DO
cccol      END DO
cccol
cccol      DO 7051 IJR=1,3
cccol      DO 7052 ICLOCK=1,SUNSET
cccol         JCLOCK=NTIME-ICLOCK+1
cccol         JEXR(IJR,JCLOCK) = JEXR(IJR,ICLOCK)
cccol 7052 CONTINUE
cccol 7051 CONTINUE

C     ONLY FOR POLAR DAY IS THERE DAYLIGHT AT SUNSET/RIS TIMES
cccccccc      if(i.eq.9 .and. j.eq.1)print *,' TFD(fastch)=', TFD(I)

      IF(TFD(I).LT.1.0) THEN
      DO IJR=1,NJR
         JR(IJR,SUNSET)=0.D0
         JR(IJR,SUNRIS)=0.D0
      END DO

        DO IJR=1,3
           JEXR(IJR,SUNSET)=0.D0
           JEXR(IJR,SUNRIS)=0.D0
        END DO

C   also do for diagnostic array JZ39(18,PH$,5,L$,Z$)

        do 7176 ikk0=1,5
        do 7176 ijj0=1,PH$
           JZ39(SUNSET,ijj0,ikk0,i,j) = 0.d0
           JZ39(SUNRIS,ijj0,ikk0,i,j) = 0.d0
 7176   CONTINUE
      END IF


C  load in JR(NJR,NTIME+3) array into new array JRA(NJR,NTIME) (both REAL*8)

         DO 761 IIC = 1,NTIME
         DO 761 IJR=1,NJR
 761        JRA(ijr,iic) = JR(ijr,iic)


Cj  for output/testing, load in AER's J1 and J2 at 30 km:  AJOUT(2,ntime,L$), JR(NJR,NTIME+3), JEXR(3,NTIME)
cj
cj       IF (J .eq. 15) then
cj          DO 771 iic = 1,ntime
cj            ajout(1,iic,i) = JEXR(1,iic)
cj 771        ajout(2,iic,i) = JEXR(2,iic)
cc            ajout(1,iic,i) = JR(1,iic)
cc 771        ajout(2,iic,i) = JR(2,iic)
cj       END IF
cj
cc JAVAER(NJR,18,46), AJOUT(2,ntime,46), AJAVDAY(3,NJR,18,46)
cj
cj         if(I .eq. 9 .and.  J .eq. 15) then
cj       print *,'IN FASTDIUR-1, i,j, JR, JRA, AJOUT, JAVAER, AJAVDAY =',
cj     > JR(1,2), JRA(1,2), AJOUT(1,2,I), JAVAER(1,i,j), ajavday(3,1,i,j)
cj         endif
cj

C
C          FAST SPECIES:  CONVERT MIXING RATIO TO VOLUME DENSITY
C
C         PUT INITIAL GUESS INTO ALLSP ARRAY
      DO 650 IFSP=1,NFSP
      ALLSP(IFSP)=ALLFSP(I,J,IFSP)
      SPNOON(IFSP)=ALLFSP(I,J,IFSP)
cpp          if(I .eq. 9)then
cpp             if(J .eq. 11 .or. J .eq. 21 .or. J .eq. 31)then
cpp                if(ifsp .eq. 23)print *,' spnoon(no2)=',spnoon(ifsp)
cpp                if(ifsp .eq. 1)print *,' spnoon(o1d)=',spnoon(ifsp)
cpp                endif
cpp                endif
c      print *,' ifsp=',ifsp,' allsp(ifsp)=',allsp(ifsp)
  650 CONTINUE
c  650 ALLSP(IFSP)=ALLFSP(I,J,IFSP)*AIRM
C

c      print *,' allsp(1,21,23)=',allsp(1),allsp(21),allsp(23)
c      print *,' spnoon(1,21,23)=',spnoon(1),spnoon(21),spnoon(23)

C         CALCULATED SLOW SPECIES:  MIXING RATIOS
      DO 660 ISP=IFSSP,NTSP
      ALLSP(ISP)=ALLSSP(I,J,ISP-NFSP)
      SPNOON(ISP)=ALLSSP(I,J,ISP-NFSP)
cpp          if(I .eq. 9)then
cpp             if(J .eq. 11 .or. J .eq. 21 .or. J .eq. 31)then
cpp                if(isp .eq. 64)print *,' spnoon(o3)=',spnoon(ifsp)
cpp                endif
cpp                endif
c      print *,' isp=',isp,' allsp(isp)=',allsp(isp)
c      ALLSP(ISP)=ALLSSP(I,J,ISP-NFSP)*AIRM
  660 CONTINUE
      ALLSP(IO2)=0.21*AIRM
      ALLSP(IM)=AIRM

C
C   CALCULATE K-RATES - AER KR array, load in Gas phase rates from the GSFC K array in KMAP
C                                          cc  AER routine KRGAS not called
C        INITIALIZE RATE ARRAY
      DO K=1,NKR
         KR(K)=0.0
      END DO


cccccccccc         CALL KRGAS(KR,AIRVD(I,J),TMP(I,J),PR(I,J))

C
C  load in GSFC K-rates into the AER K-rate array KR(NKR=285), K(RB$=195,L$,Z$), KH(RH$=14,L$,Z$)
C

         CALL KMAP(KR,NKR,I,J)



CCCCC      REAL gaerosol(l$,z$), gnataer(l$,z$), gaerice(l$,z$)
CCCCC
         gfaer = gaerosol(I,J)
         gfnat = gnataer(I,J)
         gfice = gaerice(I,J)
CCCCC
CCCCC         SURFACE(I,J) = gfaer
CCCCC         PSCSURF(I,J) = gfnat
CCCCC         ICESURF(I,J) = gfice


C  Get HET rates, first define HCl and ClONO2 for call to HETCHEM,  KH(RH$,18,L$,Z$)

          xhcl = allsp(IHCL)
          xclono2 = allsp(ICLNO3)

          CALL HETCHEM(KR,nkr,xhcl,xclono2,i,j)


cccccccccccc         CALL KRCALC(KR,I,J,iday360, gfaer, gfnat, gfice)



C
      do 1100 ii=1,ntime
         dtimea(ii)=dtimeb(i,ii)
 1100    continue
c      print *,' dtimea(fastdiur)=',(dtimea(ii),ii=1,ntime-1)


C         RUN NEWTON TO GET FAST SPECIES CONCENTRATIONS
C             ITFD IS NUMBER OF DAYS TO TIME STEP

CMESO  
c                                      don't do fast chemistry in mesosphere w/ AER model
ccccccccccc      if (j .ge. ikmes) goto 9999 

cccc         CALL NWT2DD(FSP,JR,DTIMEA,EI1,EIPER,ITFS,ITFD,I,J,IERR,LLOLA
cccc     >              ,iday360)


ccccc           CALL NWT2DD(FSP,JR,DTIMEA,EI1,EIPER,ITFS,ITFD,I,J,IERR,LLOLA, 

           CALL NWT2DD(FSP,JR,DTIMEA,EI1,EIPER,ITFS,ITFD,I,J,IERR,LOLA,
     >  gfaer, gfnat, gfice, TMP(I,J), PR(I,J), AIRVD(I,J),iday360, iyr)

c      print *,' after call nwt2dd'



C  get diurnally/night/daytime avg constituents

      CALL TIMEAV(FSP,NFSP,DTIMEA,I,J,TFD,L$)

CMESO  
 9999     dum=1  


      
C  get diurnally avgd Js, JRA(NJR,NTIME), and extra Js JEXR(3,NTIME) and diurnal Js - JAER(NJR+3,NTIME,L$,Z$)
C    also do wavelength dependent JZ39(18,PH$,5,L$,Z$) (for diagnostics ONLY)

      CALL TIMEAVJ(JRA,JEXR,DTIMEA,I,J,TFD,L$,Z$,JAVAER,JAER,
     >             PH$, IL$, JZ39, J39)
 

cpp          if(I .eq. 9)then
cpp             if(J .eq. 11 .or. J .eq. 21 .or. J .eq. 31)then
cpp        print *,' O1D=(after timeav)',(fsp(1,iii),iii=1,ntime+3),
cpp      *  ' ntime=',ntime
cpp        print *,' NO2=(after timeav)',(fsp(23,iii),iii=1,ntime+3),
cpp      *  ' ntime=',ntime
cpp        print *,' midnit=',midnit,' ris90=',ris90,' set90=',set90,
cpp      *  ' dayavg=',dayavg,' nitavg=',nitavg,' diuravg=',diuravg
cpp        endif
cpp        endif

c      CALL TIMEAV(FSP,NFSP,DTIME,DLITE(J))

c      print *,' after call timeav'

c      if(1.gt.0)stop
C
C         STORE THE FAST SPECIES MIXING RATIOS
      DO 740 IFSP=1,NFSP
C            NOONTIME CONTAINED IN ALLSP, ALSO IN FSP(IFSP,1)
         ALLFSP(I,J,IFSP)=FSP(IFSP,NTIME)/AIRM
         FSPMID(I,J,IFSP)=FSP(IFSP,MIDNIT)/AIRM
         FSPRIS(I,J,IFSP)=FSP(IFSP,RIS90)/AIRM
         FSPSET(I,J,IFSP)=FSP(IFSP,SET90)/AIRM
         FSPDAY(I,J,IFSP)=FSP(IFSP,DAYAVG)/AIRM
         FSPNIT(I,J,IFSP)=FSP(IFSP,NITAVG)/AIRM
         FSPAVG(I,J,IFSP)=FSP(IFSP,DIURAVG)/AIRM

C
C  Also load in FSP(NFSP - 40 , NTIME+3 - 18+3), into GSFC model array, for each diurnal time step
C          CNDCA(S$,ntime+3,L$,Z$) is in num. density, FSP(NFSP,NTIME+3),    NTIME = 18
C               CNDCA and FSP are REAL*8
C
C   NOTE: in com_aerd.h =>  DAYAVG=NTIME+1, NITAVG=NTIME+2, DIURAVG=NTIME+3
C
           do 777 ijtt=1,ntime+3
               if (ISPMAP(ifsp) .ne. 0)
     >             CNDCA(ISPMAP(ifsp),ijtt,i,j) = fsp(ifsp,ijtt)
 777        CONTINUE

C  reload AER FAST species into GSFC array,  - NO LONGER NEEDED
C     ALLFSP(18,46,40),  FSPMID, FSPRIS, FSPSET, FSPDAY, FSPNIT, FSPAVG(18,46,40)
C     GFDAY, GFNIT, GFAVG, GFRIS, GFSET(100,18,46), GFMID(100,18,46)  all in COMMON
cc
cc           GFMID(ISPMAP(ifsp),i,j) = FSPMID(i,j,ifsp)
cc           GFRIS(ISPMAP(ifsp),i,j) = FSPRIS(i,j,ifsp)
cc           GFSET(ISPMAP(ifsp),i,j) = FSPSET(i,j,ifsp)
cc           GFDAY(ISPMAP(ifsp),i,j) = FSPDAY(i,j,ifsp)
cc           GFNIT(ISPMAP(ifsp),i,j) = FSPNIT(i,j,ifsp)
cc           GFAVG(ISPMAP(ifsp),i,j) = FSPAVG(i,j,ifsp)

 740     CONTINUE

c      print *,' O=(in FASTCHAER)',ALLFSP(I,J,2),FSPMID(I,J,2),
c     *  FSPRIS(I,J,2),FSPSET(I,J,2),FSPDAY(I,J,2),FSPNIT(I,J,2),
c     *  FSPAVG(I,J,2)
cpp          if(I .eq. 9)then
cpp             if(J .eq. 11 .or. J .eq. 21 .or. J .eq. 31)then
cpp        print *,' NO2=(in FASTCHAER)',ALLFSP(I,J,23),FSPMID(I,J,23),
cpp      *  FSPRIS(I,J,23),FSPSET(I,J,23),FSPDAY(I,J,23),FSPNIT(I,J,23),
cpp      *  FSPAVG(I,J,23)
cpp        print *,' O1D=(in FASTCHAER)',ALLFSP(I,J,1),FSPMID(I,J,1),
cpp      *  FSPRIS(I,J,1),FSPSET(I,J,1),FSPDAY(I,J,1),FSPNIT(I,J,1),
cpp      *  FSPAVG(I,J,1)
cpp cpp        endif
cpp        endif
C         STORE THE SLOW SPECIES MIXING RATIOS UPDATED FOR EQUILIBRIUM
      DO 745 ISP=1,NSSP
      IASP=NFSP+ISP
      ALLSSP(I,J,ISP)=ALLSP(IASP)/AIRM
  745 CONTINUE

cpp          if(I .eq. 9)then
cpp             if(J .eq. 11 .or. J .eq. 21 .or. J .eq. 31)then
cpp       print *,' O3=(in FASTCHAER)',ALLSSP(I,J,24)
cpp        endif
cpp        endif

C
C
C
  800 CONTINUE
  900 CONTINUE
c      if(1.gt.0)stop
C
      RETURN
      END
