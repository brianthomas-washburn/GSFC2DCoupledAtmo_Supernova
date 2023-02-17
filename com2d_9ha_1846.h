C-*-Fortran-*-  -  this is used for emacs/Fortran mode capability, or just do M-x fortran-mode
C
	IMPLICIT real*4 (A-H,K-Z)

	INTEGER L$,Z$,L$1,Z$58,Z$59,S$,T$,RB$,RW$,RH$,RHT$,PH$,INBIG
	INTEGER Z$1
        INTEGER NYR$,NMON$,L$D,Z$X,Z$X1,Z$S,L$S,TH$

        REAL*8 delz, delzc, dth, dthst


C need to set the number of latitudes and pressure levels here
C
C       dth = latitude resolution in degrees, so L$ = number of latitudes =180/dth, eg., 45 lats = 4 deg res
c	Z$ = number of pressure levels (variable intervals), MUST CONFORM TO VALUES IN CONTROL.DAT)
C        and set the vertical resolution for constituent grid, DELZC (in meters, REAL*8):
C  just use the MINIMUM of deltaz1, deltaz2, and deltaz3:   DELZC=AMIN1(deltaz1, deltaz2, deltaz3)*1000.
C
        PARAMETER(dth=10.d0)
	PARAMETER(L$=IDNINT(180./dth)) 
	PARAMETER(Z$=46)
        PARAMETER(delzc=1.9908d0*1000.)


c	Z$X = number of grid points in extended model up to ~115 km (set to Z$ + 12 levels at 2 km res.)
c	Z$X1 = number of grid box edges in extended model up to ~115 km for STREAMF (just Z$X+1)

	PARAMETER(Z$X=Z$+12)
	PARAMETER(Z$X1=Z$+12+1)
	PARAMETER(Z$1=Z$+1)

C
C  Grid for STREAMFUNCTION CALCULATIONS:
C  
C   latitude resolution is now fixed at 2 degrees (dthst), number of latitudes L$S is computed from dthst
C   vertical resolution, DELZ (in meters, REAL*8), just use equal spacing of 1 km (1000 m):
C      the number of levels is Z$S, computed from DELZ up to 116 km - the extended model
C   
        PARAMETER(dthst=2.d0) 
	PARAMETER(L$S=IDNINT(180./dthst)+1)

        PARAMETER(delz=1.0d0*1000.)
        PARAMETER(Z$S=IDNINT(116.*1000./delz)+1)


C  TH$ is the number of potential temp levels for Kyy, done for 7-107 km, defined as 2x the vertical 
C                                                      resolution of minimum deltaz used (delzc)
ccccc        PARAMETER(TH$=IDINT((107-10)/(delzc/2000.))+1)
ccccc        PARAMETER(TH$=IDINT((107-7)/(delzc/2000.))+1)
ccccc        PARAMETER(TH$=IDINT((25-10)/(delzc/2000.))+1)  ! or maybe do only for 7-25 km??? or 10-25 km????


C   no longer used????
c	If L$=18 (original grid), need to define L$D = 37 , i.e., the number of 
C       latitudes for streamfunction and gravity wave scheme, otherwise set it to L$

	PARAMETER(L$D=37)
	PARAMETER(L$1=37)


C   no longer used
c	Z$59 = number of pressure levels for streamfunction and transport scheme (set to 59)
c	Z$58 = number of pressure levels for chemistry in extended 58 level model(set to 58)

        PARAMETER(Z$58=58)
        PARAMETER(Z$59=59)



c set the parameters:
c	RB$ = size of base reaction rate arrays
c 	RW$ = size of working reaction rate arrays
c 	RH$ = size of het reaction rate arrays
c 	RHT$ = total number of het reactions:  RH$*3  (sulfate + NAT + ICE reactios)
c 	S$ = number of constituents
c 	T$ = number of transported constituents	
c	PH$ = number of photolysis rates
c 	IL$ = number of wavelength intervals
c 	NYR$ = number of years for interannual run
c 	NMON$ = number of months for interannual run, usually NYR$*72

	PARAMETER(S$=80+10)
	PARAMETER(T$=77)

	PARAMETER(RB$=284)
	PARAMETER(RW$=30)
	PARAMETER(RH$=14)
        PARAMETER(RHT$=RH$*3)

	PARAMETER(IL$=39)
	PARAMETER(PH$=81)

	PARAMETER(NYR$=46)
	PARAMETER(NMON$=3312)


	COMMON/DAT100/INTTRA(T$), ITPDF, ILPDF, IZONAVGT, IYRCT


cc       COMMON/DAT200/KW(RW$,L$,Z$), KH(RH$,L$,Z$),      
cc           >         RNDAY(S$,L$,Z$), MR58(T$,L$,Z$X)

        COMMON/DAT200/KW(RW$,L$,Z$), KHGF(RH$,L$,Z$),
     >          RNDAY(S$,L$,Z$), MR58(T$,L$,Z$X), mr58o(L$,Z$X)


        COMMON/DKRATE/K, J, KH
        REAL*8 K(RB$,L$,Z$), J(PH$,L$,Z$), KH(RH$,18,L$,Z$)      ! KH =(RH$,NTIME,L$,Z$);  NTIME=18


        COMMON/JBL1/JDC 
        COMMON/JBL2/CNDC 
        REAL*8 JDC(PH$,18,L$,Z$), CNDC(S$,18+3,L$,Z$)         ! CNDC, JDC=(PH$,NTIME,L$,Z$);  NTIME=18


	COMMON/DAT300/TFD(L$)
	COMMON /CCNTRL/ CONVG,ITMAX,ITWCH,IKMAX,IJMAX,ICTT,JCT(S$),
     *		JSEQ(S$)
	COMMON /CCONS/ RD,WTA,RE,RCGS,RINV,RCSINV(L$),TANGENT(L$)


	COMMON/BKGD/AREA(L$),COLMOD(10),DENCLX(Z$,L$,2,12),
     *		CH4MR(12,L$,Z$),N2OMR(12,L$,Z$),H2ODEN(12,L$,Z$),
     *		TEMPSAMS(12,L$,Z$), globe

	common/aerosol/aerosol(l$,z$),hoff(10),nataer(l$,z$),
     *                   aerice(l$,z$), isadcl

	common/extra/t0,dzstar   
	COMMON /CONST1/ DZcm(z$),DY,DT1,DY1,DY2,DZ1,
     *		DZ2,DELT,C00,GZ ,DYX,DZX,XNORM,DAYL,YP(l$),
     1		ZP1(z$),ZZ(z$),RHO(z$),RSTAR(z$),h25,
     2		YPP(l$),ZPP(z$),EZ2H(z$),Sl(l$),DTHDZ(z$),
     3		THG(z$),TG(z$),TTOTH(z$), 
     4		CF(l$),BETA2(l$),XMM(l$,z$)

c  commons for reaction rate data 
      common/kk/k0(rb$),e(rb$),khi(rb$),ehi(rb$),khno3(3),ehno3(3)
      common/kkw/k0w(rw$),ew(rw$),khiw(rw$),ehiw(rw$),ibdy(rb$),
     c	ibdyw(rw$)

        REAL*8 gs, gn, gi, khg, ggsfc
        COMMON/gamhet/gs(rh$), gn(rh$), gi(rh$), 
     >                khg(rht$,L$,Z$), ggsfc(rht$,L$,Z$)

        common/gamhetg/gsg(rh$),gng(rh$),gig(rh$)

	COMMON/EDDY/EKZZI,EKYYI,FZZ(10),IZZ(10),FYY(10),IYY(10)
	COMMON/TRNSPORT/ITRANS, ICOUP
	COMMON/CTRL/INDAYS,IOUTP

	LOGICAL YSMSCAT,YSLFPRNT
	LOGICAL YSHSCT,YSO3FIX
	LOGICAL  lcolumns,lprofiles
C	LOGICAL LPHOT1,LPHOT2,LDAV,YSCORSR
	LOGICAL LDAV,YSCORSR, LDYNOUT

	CHARACTER*37 PHOTCHAR
	CHARACTER*10 LIFECHAR

c  commons for boundary contitions

	LOGICAL LSTSTATE,LBCMRTD,LBCLATTD,LBCMRSS,LBCLATSS

	COMMON/BCINFO/IYEARSS,IYEARBCTD(300),IYEARBCSS(10),IYRBC(300),
     *  ISPBC(S$),FRMOLWT(S$),BCTDINPUT(S$,300),BVAL(S$,L$), coi, noyi,
     *  BCSSINPUT(S$,10),TDLAT(S$,L$),SSLAT(S$,L$),TDLATR(S$,18),
     *  BVALG(S$,L$,10), IGSPEC(19), DLIGHT(L$,Z$), HNO3W(L$,Z$)

c                        TDLATR(S$,18) is for reading in BCs on orig. lat grid

	COMMON/BCLOGICAL/LSTSTATE,LBCMRTD(S$),LBCLATTD(S$),
     *  LBCMRSS(S$),LBCLATSS(S$)

c  commons for photolysis calculations

	COMMON /CFLUX/ WVL(IL$),FLUX0(IL$),FLUX(IL$),RFLUX(IL$,L$,Z$)

	COMMON /CJCOEF/ PHOTCHAR(PH$) 
	COMMON /CXSECT/ XSECT(PH$,IL$),XO2(IL$),XO3(IL$), IPHOT(PH$),
     *   XSECTTD(PH$,IL$,201), xspd(3,201,157), prpd(157), aprpd(157)
        COMMON /CNO2XS/  NO2T(IL$), XNO2(IL$)

	COMMON /CSRB1/ SR1(9,17),SR2(6,17),SRNO1(9,2),SRNO2(5,2),
     .  INOL(2),INOU(2),NOSIG(2,L$,Z$)
	COMMON /CSRB2/SRL(2),PRA(4),PRAL(4),SIGNO(2,2),SIGNOUSE(4),
     .  TLIM,SIGLM,ESIGLM,CORSR(17),ISR(2),JBL(3),JBU(3),ILLB,
     .  INOUSE(4)
	COMMON /CSRB3/ XSCHRUN(20,L$,Z$),ILLSR,ILUSR
	COMMON /CFIT1/ WO1D(2),WN2O5(2),A0(3),A1(3),A2(3),A3(3),A4(3),
     *  B0(3),B1(3),B2(3),B3(3),B4(3)
	COMMON /CFIT2/ ILO1D(2),IN2O5(2)
	COMMON /CFIT3/ TS320,TS300,TS230,TS225,TS220,TS180,WCH2O
        COMMON /CFIT4/ AS01,AS02,AS03,AS11,AS12,AS13,AS21,AS22,AS23,
     .  AS31,AS32,AS33,AS41,AS42,AS43,AS51,AS52,AS53,AS61,AS62,AS63
        common /cfit5/ hobrlam(3),hobra(3),hobrb(3)
	COMMON /CSOLAR1/ RSC(5),R27(4),RLAMB(6)
	COMMON /CSOLAR2/ D1,D2,D3,D4,D5,D6,DEC1,DEC2,DEC3,DEC4,DEC5,
     .  DEC6,DEC7,DEC8
	COMMON /CDAV/ FAVG(12,9,10),DECA(9),TAUA(12),LATR(10),ZREF,
     .  IDA,ILATR,ITA,ITAUA(1001)
	COMMON /CLOGSR/ YSCORSR
	COMMON /CLOC/ ZP,LATP,TS,T,DAY0,DAY,CLOCK,PRP,DECD
	COMMON /CPHYS/ SINL,COSL,SIND,COSD,CHID(L$),OMEGAD,HBAR,
     .  ZBAR,R0,YRL,NM,DPHI,SFAC(L$),ZGRZ(L$)
	COMMON /CDEN/ NCOL(10),ISNCOL(10),JNCOL,JSPD,ISP,
     .  ISO2,ISO3,ISNM
	COMMON /BEDO/ ALBB(Z$),YSMSCAT
	COMMON/MSCA/SNDZ(3,Z$),XNTAU(3,Z$),CHIMS,DECMS,LATMS,SO2,
     *  SO3,PHLUX,LAMMS,SFACMS,IKMS
	COMMON /CMAT/ JLRT1,JURT1,LDAV
	COMMON /FLUXTAU/ FLUXMULT(IL$,L$,Z$),TAUMULT(IL$,L$,Z$),
     *  IDAYMN(12)

	COMMON /CREST/TSS,CHNGX,IJCX,IKCX

	COMMON /CINDX/ PI,DTR,RTD,RPI,BK,AMU,CP,CL,ISCO2,ISCO3
	COMMON /CYEAR/ IYR,DYT,DAYST,DAYIN
 	COMMON /YEAR360/ DAY360,DYT360, IDAY360
 
	COMMON /SOLARCYC/ SOLCYCSET,PIYRUSE,FMAXSET,FMINSET,
     >                    SOLCYCR(IL$,360,57), iysc, wvmid(IL$)
 	COMMON /NOXGCR/ GCR(L$,Z$), NOXAIR(L$,Z$)
 	COMMON /GCRMM/ GCRMIN(18,46),GCRMAX(18,46)

C  COMMONS FOR NCB
	common/bcnox/ibc(l$,z$),bcclono2(l$,z$),bcbrono2(l$,z$)
        common/noxc/bcclo24(l$,z$),bcn2o5(l$,z$),cno2(l$,z$),
     c  bcbro24(l$,z$),gamloc(l$,z$),zzz(10,l$,z$),yy(5,l$,z$),
     c  xxz(4,l$,z$),alpha(l$,z$),betaz(l$,z$),isw2(l$,z$),
     c  isw(l$,z$),icount(l$,z$)

C  COMMONS FOR LIFETIME CALCULATION
	COMMON/LIFE/CLOSS(40,L$),CLOSSST(40,L$),CLOSSTR(40,L$),
     *		COLUMN(40,L$),RLOSS(40,L$,Z$),LIFECHAR(40),
     *		YSLFPRNT


C  COMMON FOR INTERACTIVE WATER VAPOR, also new ITDAY param
C   NOTE: don't use TROPWV, PRECIPYR to save memory, trop H2O from OOrt now contained in the 
C         HALOE climatolgy, which is now used in RAIN
        COMMON/OORTH2O/TROPWVIN(18,11,360), ITDAY,
     *         RAINDAY(L$,360), RAINYR(L$)
ccccc     >         , TROPWV(L$,Z$X,360), PRECIPYR(L$,Z$,360)


C  COMMON FOR Model Chemistry extension into Lower Thermosphere (levels 47-58, 90-115 km), and write-out
C     also for write out of all photochemical prod/loss, and Ox prod/loss , 
C     and H2O Loss from RAINOUT and ICE FALLOUT from PSC's (in HETCHEM) ALL IN MIX RAT/SEC
c                           30                            60          72   
        COMMON/THERMOS/CTLOSS(T$,L$,Z$), CTPROD(T$,L$,Z$),
     >       CHMPOX(2,L$,Z$), CHMLOX(8,L$,Z$), RAIN1(L$,Z$), ICE1(L$,Z$)

C                AND advection, diffusion terms (Y & Z), ALSO IN MIX RAT/SEC

        COMMON/COMTRANS/DIFFY(T$,L$,Z$X), DIFFZ(T$,L$,Z$X),
     >              ADVY(T$,L$,Z$X), ADVZ(T$,L$,Z$X), DIFFYZ(T$,L$,Z$X)

C  COMMON FOR CO2 run with DAILY time dependent BCs for 1950-2050, and TMONC (month counter) 
        COMMON/COMCO2/CO2BCIN(18,101,360), CO2BC(L$,101,360), 
     >                TMONC, TDAY, DAYCH

c                           30                            60          72
C
C  COMMON for tropical tropopause BCs for H2O, 1992-2000, and climatology
C     and daily Oort/HALOE high resolution H2O climatology -> convert to proper grid in TEMPIN
C
        COMMON/CH2OBC/H2OBC(18,360,9), H2OBCCL(18,360), 
     >     halh2or(45,92,360), lathal(45), zzhal(92), halh2o(L$,Z$)

C
C 
C  COMMONs for the new high resolution grid and the dynamics/transport routines
C
	COMMON /CGRID/ P(Z$), PRESS(Z$X), PHI(L$), PHIH(L$+1), DT, DT0,
     >      inz1, inz2, inz3, ZALT(Z$X), LDYNOUT, ITDYN, ik40, 
     > ik35, ik58, ikmes, ij35, ij55, ijsp, ijnp, ZALT90(Z$), PRES90(Z$)

	COMMON /CGRID8/deltaz1, deltaz2, deltaz3, deltaz4, DT8 
        REAL*8 deltaz1, deltaz2, deltaz3, deltaz4, DT8

        COMMON/DATABC/GWBC(13,L$S,72)

        COMMON/DATA12/LHIN(37,7,36), SBUVMRO3(12,L$,Z$), 
     >     TROPHTIN(18,360), GWBCIN(13,37,72),
     >     LH10(L$S,Z$S,36), TROPHTF(L$,360), ITROP360(L$,360)
C
C    TROPHTF = new NMC tropopause hgts for each day at model lats
C    ITROP360 = model index of trop hgts for each day at model lats
C    GWBC = BC's for gravity waves read in TEMPIN

C
C  SEPARATE COMMON for STREAMF/TEMPIN/GWAVE/:
C
C  NOTE: TEMP5IN, TEMPALLS  are the temperatures for the STREAMF/TROPKZZ calculations
C        TEMPALL are the temperatures for the chemistry model grid points
C 
        COMMON/DYNSTRM/KYYIN(L$S,Z$S,72), EPIN(L$S,Z$S,72), 
     >                 TEMP5IN(L$S,Z$S,72), UBARIN(L$S,Z$S,72),
     >                 EHFIN(L$S,Z$S,72), EHFZIN(L$S,Z$S,72),
     >      QBARY(L$S,Z$S,72), HEATIN(L$S,Z$S,72), FEGWIN(L$S,Z$S,72),
     >                         QYCL(91,117,72), EPCL(91,117,72)


C  COMMON for orthogonal tracer Kyys:

      COMMON/CKYY/kyyotr(19,49,72),latot(19),zzot(49), kyyot(L$+1,Z$,72)



C  COMMON for extra quantities on model grid:
C      -  DON'T need to store these for every time period
C
        COMMON/DYNT/VALL(L$+1,Z$X), KZZGW(L$,Z$X1), KZZTROP(L$,Z$X1), 
     >              KYZALL(L$+1,Z$X1), THSLPALL(L$+1,Z$X1)



       COMMON/DYNMOD/WALL(L$,Z$X1,74), KYYALL(L$+1,Z$X,74),  
     >   KZZTOT(L$,Z$X1,74), TEMPALL(L$,Z$X,74), TEMPALLS(L$S,Z$S),
     >   UBARCG(L$,Z$X,74), QYCG(L$,Z$X,74), EPCG(L$,Z$X,74),

     >   W11(L$,Z$X1,2), KZZ11(L$,Z$X1,2), KYY11(L$+1,Z$X,2), 
     >   TEMP11(L$,Z$X,2), 

     >   W74(L$,Z$X1,2), KZZ74(L$,Z$X1,2), KYY74(L$+1,Z$X,2), 
     >   TEMP74(L$,Z$X,2)


C    ! Kyz computed at grid box corners, KZZ is defined at top/bottom BOX EDGES

C  these are for the Kyz calculation: 
c        DPTDYALL(L$,Z$X,NMON$), DPTDZALL(L$,Z$X,NMON$)

C
C    common for climatological transport values:
C
        COMMON/DYNCL/WCL(L$,Z$X1,74), KYYCL(L$+1,Z$X,74),  
     >               KZZCL(L$,Z$X1,74), TEMPCL(L$,Z$X,74)
C
c
c                           30                            60          72   
C
        COMMON /CATM/ZGP(L$,Z$X), ZTOP(L$), EKZH(Z$58),
     >               PRES58(Z$58), ZZ58(Z$58)

        COMMON /DATLAR/NCOLGD(10,L$,Z$),      !! TEMPCH(L$,Z$X,36,NYR$),
     >                 ZKM(L$,Z$X), NCOLD(2,18+3,L$,Z$)

	COMMON /CTRDAT/ EKYY(L$+1,Z$X), EKYZ(L$+1,Z$X1), EKZZ(L$,Z$X1)


c  ZALTE(Z$X1), PRESSE(Z$X1) are the altitudes, and pressures of the box edges
C    for the extended model (~115 km), and W1 is for the constituents, from SETDAILY

        COMMON/AGRID/ZALTE(Z$X1), PRESSE(Z$X1), LATEG4(L$+1)


        COMMON/VWMASSAD/WBAR, W1, W, V
        REAL*8 WBAR(2,Z$X1), W1(L$,Z$X1), W(L$,Z$X1), V(L$+1,Z$X) 

        COMMON/TPCORE1/ZK, DY0, DZ0, XMR               !NEW COMMON here for PPM transport (8/96)
        REAL*8 XMR(L$,Z$X,T$), ZK, DY0, DZ0

C  TEMP(L$,Z$X) is for chemistry ONLY, at grid points, load into array TEMPCH(L$,Z$X,36,NYR$) - REAL*4
C  define all arrays for NEWDIFY and NEWDIFZ (including M and cn, c) to be REAL*8
C    also define new CNOUT array in REAL*4 for output of the CN array,
C    CTEST is for output of test run of 6 transported species + M,  ZALT8(Z$X) is ZALT(Z$X) in REAL*8
c    delyc is the DELta-Y for the constituent grid,  LATEG, COSEG are the latitudes, cos of box edges (sides)
C
      COMMON/NEWDIF1/TEMP, Q, DELTAZ, ST, ST2, PRESSE8, DP, H
      COMMON/NEWDIF2/DIFNY, DIFNYZ, DIFNZ, DIFTH
      COMMON/NEWDIF3/CNOUT, CTEST, LAT, COSC, TCOS, LATEG, COSEG, 
     >               C14, TMONC14, THSLP
      COMMON/NEWDIF4/radd, dphid, phid, phihd, delyc
      COMMON/NEWDIF5/ZALT8, ZALTE8, PRESS8, LAT181
c
      REAL*8 Q(L$,Z$X), DELTAZ(L$,Z$X), LAT(L$), COSC(L$)
      REAL*8 ST(T$,L$,Z$X), H, LATEG(L$+1), COSEG(L$+1)
      REAL*8 difny(T$,L$,Z$X), difnz(T$,L$,Z$X), difnyz(T$,L$,Z$X)     ! ,difth(T$,L$,TH$)
      REAL*8 TEMP(L$,Z$X), DP(Z$X), TCOS 
      REAL*8 ZALT8(Z$X), ZALTE8(Z$X1), PRESS8(Z$X), PRESSE8(Z$X1)
      REAL*8 CTEST(7,L$,Z$X), ST2(T$,L$,Z$X), LAT181(181)
      REAL*8 radd, dphid, phid(L$), phihd(L$+1), delyc
      REAL*4 CNOUT(S$+8,L$,Z$), C14(L$,Z$), THSLP(L$+1,Z$X1)

C        w/ coupled model,CNOUT now includes TBAR, UBAR, HEAT, COOL, DRAG, QBARY

      COMMON/CBL1/C, N2
      REAL*8 C(S$,L$,Z$), N2(L$,Z$)

      COMMON/CBL2/CN, M
      REAL*8 CN(S$,L$,Z$), M(L$,Z$X)

      COMMON/CBL3/C58
      REAL*8 C58(T$,L$,Z$X)

      COMMON/CBL4/CN58
      REAL*8 CN58(T$,L$,Z$X)

C
C   common for input files
C
      COMMON/CINPUT/CCIN(S$,18,46), CCIN1(18,46), LATIN(18), LAT5IN(37), 
     >              PRES46(46), ZZ46(46), PRES59IN, ZZ30(30)

C                                                   also do LAT, LATST, ZSTR arrays as REAL*4
      COMMON/CINPUT5/pres59(59), zz59(59), xlat5(37), 
     >               LAT4(L$), LATST4(L$S), ZSTR4(Z$S)

C   common for look-up tables of prod/loss arrays
C
cef      COMMON/tracert/ch4lin(360,18,46), h2olossin(360,18,46), 
cef     >      h2oprodin(360,18,46), oxprodin(360,18,46), oxlin(360,18,46), 
cef     >      ch4l(360,L$,Z$), h2oloss(360,L$,Z$), h2oprod(360,L$,Z$), 
cef     >      oxprod(360,L$,Z$), oxl(360,L$,Z$)

C   SOME THINGS IN COMMON FOR STREAMF COMPUTATION

       COMMON/NEWSF1/XN2, UBAR1, DUDZ, XKZZT, FXTT, FK, FKQ, FKQ2, FRQ, 
     >  RAD, OM1, HH, RR1, KAP, DELY, DELY2, DELZ2, LATST, BBCADJ, 
     >  zstr, prst, rhost, grst, xpi, qp1, cose, tana, ff0, WALLA8

       REAL*8 XN2(L$S,Z$S), UBAR1(L$S,Z$S), DUDZ(L$S,Z$S),XKZZT(L$S,Z$S)
       REAL*8 FXTT(L$S,Z$S), FK(L$S,Z$S), FKQ(L$S,Z$S), FRQ(L$S,Z$S)
       REAL*8 RAD, OM1, HH, RR1, KAP, DELY, DELY2, DELZ2, xpi, qp1
       REAL*8 zstr(Z$S), prst(Z$S), rhost(Z$S), grst(Z$S),BBCADJ(L$S,72)
       REAL*8 latst(L$S), cose(L$S), tana(L$S), ff0(L$S),WALLA8(L$S,Z$S)
       REAL*8 FKQ2(L$S,Z$S)
       

C   for DYNOUT, KEEP AT REAL*4 --
C
       COMMON/DYNEXTRA/WALLA(L$S,Z$S), VALLA(L$S,Z$S), UBARALL(L$S,Z$S),
     >    HEATALL(L$S,Z$S), CHIALL(L$S,Z$S), FXTALL(L$S,Z$S),
     >    FKALL(L$S,Z$S), FKQALL(L$S,Z$S), FRQALL(L$S,Z$S),
     >    FKQ2ALL(L$S,Z$S), FEGWALL(L$S,Z$S),
     >    CFFALL(L$S,Z$S), QVALL(L$S,Z$S), KRUALL(L$S,Z$S),
     > EHFALL(L$S,Z$S),EHFZALL(L$S,Z$S),QTOTALL(L$S,Z$S),LHALL(L$S,Z$S),
     >    DDHALL(L$S,Z$S), DMOM(L$S,Z$S), GWH1(L$S,Z$S), GWH2(L$S,Z$S)

C  COMMON for new diurnal AER chemistry, NOONRAT array is the ratio of noontime radical/family
C
       COMMON/CDIUR8/NOONRAT
       REAL*8 NOONRAT(S$,L$,Z$)


C  COMMON FOR EXTRA VARIABLES, SINCE COMPILER DOESN'T SEEM TO LIKE ADDING ON EXISTING COMMON BLOCKS
C    PLRAT=ratio:photoequil members/family, 1=CHx, 2=HOx;  
C    CHXRAT=ratio of CHx members to total family - this is kept constant during polar night
C         iclim=1 => climatology is used, ie, im324=1,36 in SETDAILY


      COMMON/CEXTRA/INBC, iclim, PLRAT(2,L$,Z$), CHXRAT(5,L$,Z$), 
     >         hno3p1(7,L$,Z$), hno3l1(3,L$,Z$), 
     >         GFAEROUT(3,L$,Z$), CDN(2,S$,L$,Z$), IAGMAP(382)

C       ! IAGMAP maps AER-GSFC reaction rates

C    ! khaer is only used if you need to output diurnal cycle of reaction rates
ccccc      COMMON/CEXTRA1/khaer(207,18,L$,Z$)



ccccc        COMMON/C9AA/TEMP9AA(18,58,324), KYY9AA(18,58,324), 
ccccc     >              W9AA(18,59,324), KZZ9AA(18,59,324)


C  COMMON for GCM output of latent heating, ozone climatology, and Harvard Ox dry deposition

      COMMON/CGCM/LHGCM(L$S,Z$S,72), OZBC(360,L$,10),
     >            HARVARD_DDOX(91,360), LAT_DDOX(91), OXDEP(L$,360)

cccccccc     >   DTMPGCM(91,117,36), LATGCM(91), ZZGCM(117),



C  COMMON FOR EXTRA VARIABLES, SINCE COMPILER DOESN'T SEEM TO LIKE ADDING ON EXISTING COMMON BLOCKS
C    PLRAT=ratio:photoequil members/family, 1=CHx, 2=HOx;  
C    CHXRAT=ratio of CHx members to total family - this is kept constant during polar night
c    HOXRAT = ratio of HOx members to total family - this is kept constant during polar night

        COMMON/HOXCM/OHNEW(L$,Z$), HNEW(L$,Z$), HOXRAT(3,L$,Z$), 
     >                       HOXSTEP(25,4,L$,Z$), HOXIT(20,3,L$,Z$)

C
C   Common for PDFs: LAT181(181) (above), and ijind(181) are used in PDFINSOL
C      now includes COMMON for reading in PDFs, and the starting/ending indicies

      COMMON/CPDFR/TPROB45(125,45,76,72), LPROB45(70,45,76,72), 
     >   itr1(45,76,72), itr2(45,76,72), ilr1(45,76,72), ilr2(45,76,72)

      COMMON/CPDF/TPROBD(211,45,76), LPROBD(181,45,76), ijind(181),
     >            TPROB(211,L$,Z$), LPROB(181,L$,Z$), ijl(L$), ikz(Z$),
     >            ITP1(L$,Z$), ITP2(L$,Z$), ILP1(L$,Z$), ILP2(L$,Z$), 
     >            cdum(187,22), 
     >            TPROB_PSC(211,L$,Z$), ITPSC1(L$,Z$), ITPSC2(L$,Z$)
