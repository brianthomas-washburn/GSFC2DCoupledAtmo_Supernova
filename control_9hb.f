	SUBROUTINE CONTROL

        include "com2d.h"

	SAVE

c ITRANS:  no of transported species
c INDAYS:  no of days for this calculation
C IOUTP:  interval in days of printout to unit 33
C ICOUP :  switch for MODEL DYNAMICAL FIELDS:
C     -1 => FIXED Model IAV;    0 => FIXED MODEL Clim;   1 => COUPLED MODEL
C
	READ(7,2111)ITRANS, INDAYS, IOUTP, ICOUP
2111       FORMAT(I5,I10,5I5)

	print *, itrans, indays, icoup

C READ IN LSTSTATE (LOGICAL DEFINING IF RUN IS STEADY-STATE OR TIME-DEPENDENT.  
C READ IN IYEARSS - YEAR of STEADY-STATE SCENARIO, or starting year of run (Jan 2012)
	READ(7,1900)LSTSTATE,IYEARSS
1900	FORMAT(1X,L1,1X,I4)

!BT 2018: Adding this separate startYear readin in order to allow for 
!         automating repeated runs (e.g. run 2010-2011, then restart 2012-2013, etc.)
!         Runscript updates the value in startYear.txt file
	open(unit=2551,file='startYear.txt',form='formatted')
	read(2551,'(1I4)') IYEARSS
	read(2551,'(1I4)') INDAYS !BT2023 adding length of run in days specified by runscript
	close(2551)
	print *, '------- start year (IYEARSS): ', IYEARSS
	print *, '------- number of days of run (INDAYS): ', INDAYS
!end BT

                             ! IFIXOH = 0/1 (model/TRANSCOM clim troposph OH)
	READ(7,1901) IFIXOH
1901	FORMAT(I2)

	print *, lststate, iyearss, ifixoh


	READ(7,1111)(INTTRA(I1),I1=1,ITRANS)
C   INTTRA IS THE ARRAY OF SPECIES WHICH WILL BE 
c	TRANSPORTED IN THIS RUN 

        print *,inttra

1111       FORMAT(16I5)

C  READ CONTROL LOGICALS
	READ(7,1114)LDAV,YSMSCAT,YSCORSR,YSLFPRNT,
     c		YSO3FIX,YSHSCT
1114	FORMAT(30L1)
C   LDAV=TRUE   THEN DIURNAL AVERAGE APPROXIMATION IS ON.
C   LDAV= FALSE   THEN DAYTIME (NOON)  CALCULATION IS DONE.
C   YSMSCAT= TRUE   THEN MULTIPLE SCATTERING IS USED TO CALCULATE
C     THE RADIATION FIELD.  BOTH RAYLEIGH SCATTERING AND A GROUND
C     ALBEDO OF 0.3 ARE USED.
C   YSCORSR=TRUE   THEN CORRECT SCHUMANN RUNGE BAND CROSS SECTIONS
C     IN THE MANNER SUGGESTED BY J. FREDERICK ON 5/26/83
C   YSLFPRNT Inital value should be false; set true when IYR changes
C   YSO3FIX= TRUE  THEN O3 IS FIXED
C   YSHSCT= TRUE THEN stratospheric planes fly!!
C
	print *,ldav,ysmscat,yscorsr,yslfprnt,yso3fix,yshsct

       do 1220 iij=1,10
 1220	  READ(7,1118)
 1118	FORMAT(1X)
C
C  now read in altitude grid information, Z$ is the total number of levels, defined in COMMON
C    deltaz1 is the 1st grid interval in km, done for inz1 grid points
C    deltaz2 is the 2nd grid interval in km, done for inz2 grid points
C    deltaz3 is the 3rd grid interval in km, done for inz3 grid points
C    any remaining grid points, Z$-nz2-nz1, are defined at intervals of deltaz4 
C      (usually the 2km standard)
C
C   DP=0.2844, the altitude and pressure grids are set up in INPUT
C
       READ(7,1113) deltaz1, inz1, deltaz2, inz2, deltaz3, inz3, deltaz4
1113       FORMAT(D9.4, I7, D7.2, I7, D7.2, I7, D7.2)


	RETURN
	END
