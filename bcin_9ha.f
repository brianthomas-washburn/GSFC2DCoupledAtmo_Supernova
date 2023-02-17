	SUBROUTINE BCIN


        include "com2d.h"


c  common for GMI sfc constits, lightning

        COMMON/CGMI/ gmi2dbc(14,19,L$,10), timesb(14),
     >               gmilt(14,45,20), xlatlt(45), zzlt(20)


        REAL tdbc18int(18), tdbcint(L$), gmi2dbcr(14,19,45,10), zzgmi(4)
        REAL gmilat(45), bbin(45), bbout(L$)


	SAVE


C  this is from bcin_wmo02.f (for the WMO-2002 BCs) and has been adopted for 9aa (variable resolution)
C
C
C                         INBC is the number of species that have SS BCs read in in the BC table
C
C
c   TDLATR(S$,18) is for reading in BCs on orig. lat grid


        INBC = 21

C INITIALIZE SOME ARRAYS
	DO 100 IS=1,S$
	LBCMRTD(IS)=.FALSE.
	LBCLATTD(IS)=.FALSE.
	LBCMRSS(IS)=.FALSE.
	LBCLATSS(IS)=.FALSE.
	DO 100 IL=1,L$
	TDLAT(IS,IL)=1.0
	SSLAT(IS,IL)=1.0
	BVAL(IS,IL)=0.0
100	CONTINUE


	DO 101 IS=1,S$
	DO 101 IL=1,18
 101	   TDLATR(IS,IL)=1.0


C READ IN HEADER INFORMATION

	READ(20,1910)
	READ(20,1910)
	READ(20,1910)
1910	FORMAT(1X)

C READ IN MOLECULAR WEIGHTS OF SPECIES.  MOLECULAR WEIGHTS GM/MOLE.
C THESE ARE PRIMARILY USED IN FLUX BOUNDARY CONDITIONS.
	READ(20,1910)
	READ(20,1910)
	READ(20,1950)(FRMOLWT(II),II=1,10)
	READ(20,1910)
	READ(20,1910)
	READ(20,1950)(FRMOLWT(II+10),II=1,5)
1950	FORMAT(10F8.2)

C NOW FOR TIME-DEPENDENT CONDITIONS
	READ(20,1910)
	READ(20,1910)
	READ(20,1910)

C READ IN IBC, YEAR NUMBER (FOR USE IN 2D MODEL), OF BOUNDARY CONDITIONS
C INPUT
	READ(20,1915)IEND
        write(6,*)'in bcin..., iend = ',iend
1915	FORMAT(37X,I3)
	READ(20,1925)(IYRBC(II),II=1,IEND)
1925	FORMAT(16I5)
	READ(20,1910)

C READ IN INSTART (STARTING INDEX) AND INFINISH (FINISHING INDEX)
	READ(20,1960)INSTART,INFINISH
1960	FORMAT(42X,I3,15X,I3)
        write(6,*)'in bcin..., instart, infinish = ',instart,infinish

C READ IN YEARS OF INPUT FOR TIME-DEPENDENT BOUNDARY CONDITIONS
	READ(20,981)(IYEARBCTD(JJ),JJ=INSTART,INFINISH)
981	FORMAT(10I9) 
        write(6,*)'in bcin..., iyearbctd = ',(iyearbctd(jj),jj=instart,
     c infinish)
	READ(20,1910)

C NOW READ IN BOUNDARY CONDITIONS FOR THE 30 SPECIES THAT HAVE TIME-DEPENDENT B.C.'S from WMO-2010
c                             ! BCTDINPUT(S$,300) - includes HFCs (Sept. 2011)
	DO 1970 II=1,30
	READ(20,2000)IS,LBCMRTD(IS),LBCLATTD(IS)
        write(6,*)'in bcin ..., ii, is, lbcmrtd ,lbclattd = ',ii,is,
     c             lbcmrtd(is),lbclattd(is)
	ISPBC(II)=IS
2000	FORMAT(22X,I2,1X,L1,L1)

        if (II .le. 27) then 
            READ(20,1981)(BCTDINPUT(IS,JJ),JJ=INSTART,INFINISH)
1981	    FORMAT(10F9.5)
        endif
                                     ! slightly different format statement for N2O, CH4, CO2
        if (II .ge. 28) then 
            READ(20,1982)(BCTDINPUT(IS,JJ),JJ=INSTART,INFINISH)
1982	    FORMAT(10F9.2)
        endif

	DO 1972 JJ=1,INSTART-1
	BCTDINPUT(IS,JJ)=BCTDINPUT(IS,INSTART)
1972	CONTINUE
	DO 1974 JJ=INFINISH+1,IEND
	BCTDINPUT(IS,JJ)=BCTDINPUT(IS,INFINISH)
1974	CONTINUE

	if(is.eq.18)print *,' bctdinput=',(bctdinput(is,jj),jj=1,iend)

C ONLY ALLOW LATITUDE DEPENDENCE FOR CFC11 - THIS WILL BE USED FOR OTHER
C CFC CONSTITUENTS WHEN FLUX BOUNDARY CONDITIONS
	IF(LBCLATTD(IS).AND.II.EQ.1)
     *  READ(20,1980)(TDLAT(IS,JJ),JJ=1,L$)
1980	FORMAT(8F10.4)

	READ(20,1910)

1970	CONTINUE


C NOW FOR STEADY-STATE CONDITIONS
	READ(20,1910)
	READ(20,1910)

C READ IN YEARS OF INPUT FOR STEADY-STATE BOUNDARY CONDITIONS
	READ(20,980)(IYEARBCSS(JJ),JJ=1,7)
980	FORMAT(I6,7I10)

	READ(20,1910)

C NOW READ IN BOUNDARY CONDITIONS FOR THE 21 (INBC) SPECIES THAT HAVE STEADY-
C STATE B.C.'S FOR FIVE DIFFERENT YEARS (this has NOT been updated for WMO-2002)
	DO 2070 II=1,inbc
	READ(20,2000)IS,LBCMRSS(IS),LBCLATSS(IS)
C	print *,' is=',is,' lmr=',lbcmrss(is),' llat=',
C     *  lbclatss(is)
	READ(20,1980)(BCSSINPUT(IS,JJ),JJ=1,7)
	IF(LBCLATSS(IS))READ(20,1980)(TDLAT(IS,JJ),JJ=1,L$)
C	print *,' tdlat=',(tdlat(is,jj),jj=1,l$)
2070	CONTINUE


	READ(20,1910)
	READ(20,1910)
	READ(20,1910)

C DO OTHER TRANSPORTED GASES
	READ(20,2060)IGAS
2060	FORMAT(69X,I2)
c	print *,' igas=',igas
C NOW READ IN BOUNDARY CONDITIONS FOR THE OTHER TRANSPORTED GASES
	DO 3070 II=1,IGAS
	READ(20,2000)IS,LBCMRSS(IS),LBCLATSS(IS)
c	print *,' is=',is,' lmr=',lbcmrss(is),' llat=',
c     *  lbclatss(is)
	READ(20,1980)BCSSINPUT(IS,1)
c	print *,' bcssi=',bcssinput(is,1)
C
C - Need to keep to 18 latitudes for reading in here, then interpolate to new latitude grid - EF, 10/01
c    TDLATR(S$,18) is for reading in BCs on orig. lat grid
	IF(LBCLATSS(IS))READ(20,1980)(TDLATR(IS,JJ),JJ=1,18)
c	print *,' tdlat=',(tdlat(is,jj),jj=1,l$)

C  now linearly interpolate TDLATR to new latitude grid, tdbc18int(18), tdbcint(L$)

        do 504 ij=1,18
 504	   tdbc18int(ij) = TDLATR(IS,ij)

          CALL LINTERP(LATIN, 18, tdbc18int, LAT4, L$, 0, tdbcint)

        do 505 ij=1,L$
 505	   TDLAT(IS,ij) = tdbcint(ij)
C
C SET B.C.'S FOR OTHER TRANSPORTED GASES
	IF(LBCMRSS(IS))THEN
	DO 3080 IL=1,L$
	BVAL(IS,IL)=TDLAT(IS,IL)*BCSSINPUT(IS,1)*1.E-9
c	print *,' is=',is,' il=',il,' tdlat=',tdlat(is,il),
c     *  ' bval=',bval(is,il)
c	print *,' bcssi=',bcssinput(is,1)
3080	CONTINUE
	ENDIF

	IF(.NOT.LBCMRSS(IS))THEN
	DO 3085 IL=1,L$
	BVAL(IS,IL)=TDLAT(IS,IL)*BCSSINPUT(IS,1)
c	print *,' is=',is,' il=',il,' tdlat=',tdlat(is,il),
c     *  ' bval=',bval(is,il)
c	print *,' bcssi=',bcssinput(is,1)
3085	CONTINUE
	ENDIF

3070	CONTINUE


C THIS SECTION OF THE CODE IS FOR A STEADY-STATE B.C.
	IF(LSTSTATE)THEN
           print *,' enter bcin lststate .eq. true'
           print *,' iyearss',iyearss
	DO 5000 IC=1,7

	IF(IYEARSS .EQ. IYEARBCSS(IC))THEN
	DO 5010 II=1,inbc
	IS=ISPBC(II)

	DO 5020 IL=1,L$
	  BVAL(IS,IL)=TDLAT(IS,IL)*BCSSINPUT(IS,IC)*1.E-9
          if (IS .eq. 20) 
     >        BVAL(IS,IL) = TDLAT(IS,IL)*BCSSINPUT(IS,IC)*1.E-6                 ! CO2 BC's are in PPMV
5020	CONTINUE

	print *,' ii=',ii,' ic=',ic,' iyearss=',iyearss,
     *  ' bval=',bval(is,9),' is=',is
c	print *,' bcssi=',bcssinput

5010	CONTINUE
	ENDIF

5000	CONTINUE
	ENDIF


c  read in SFC constituents from GMI simulation - gmi2dbcr(14,19,45,10)
C    - this was written out in IDL on Linux, so use "little_endian";  zzgmi(4) is .5, 1.5, 2.5, 3.5km
C
        OPEN (77, FILE = '/misc/kah07/fleming/gmic/gmi_2dbc.xdr',
     >        convert="little_endian", FORM = 'binary', 
     >        ACCESS='SEQUENTIAL', STATUS = 'OLD')
           READ (77) gmi2dbcr
           READ (77) zzgmi
        CLOSE (77)


C  interpolate to current latitude grid:  gmilat(45), bbin(45), bbout(L$), LAT4(L$)
C  -> gmi2dbc(14,19,L$,10) is in COMMON - everything is in ppp;  ik: 1=0.5 km;  2=1.5 km
C
        do 151 ij=1,45
 151	   gmilat(ij) = (ij-1)*4. - 88.


        do 205 ik=1,10
        do 205 iis=1,19
        do 205 iid=1,14

          do 210 ij=1,45
 210	       bbin(ij) = gmi2dbcr(iid,iis,ij,ik)

            CALL LINTERP(GMILAT, 45, BBIN, LAT4, L$, 0, BBOUT)

          do 220 ij=1,L$
 220	       gmi2dbc(iid,iis,ij,ik) = bbout(ij)

 205	 CONTINUE


C      define time index for 14 months:  timesb(14) in COMMON
C
        do 700 itj=1,14
 700	   timesb(itj)  = itj*30. - 15. - 30.


C
C  read in GMI lightning rates:  gmilt(14,45,20), xlatlt(45), zzlt(20) in COMMON
C         N production in #/cm^3/sec
C    - this was written out in IDL on Linux, so use "little_endian"
C
       OPEN (77, FILE = '/misc/kah07/fleming/gmic/gmi_2d_lightning.xdr',
     >       convert="little_endian", FORM = 'binary', 
     >       ACCESS='SEQUENTIAL', STATUS = 'OLD')
          READ (77) gmilt
          READ (77) xlatlt
          READ (77) zzlt
       CLOSE (77)


	RETURN
	END
