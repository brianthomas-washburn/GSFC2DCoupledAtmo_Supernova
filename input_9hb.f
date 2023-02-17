	SUBROUTINE INPUT
C
C  Base 8A  version  has P  and PRESS arrays (which are in COMMON) defined here
C  Base8sa has the PHIH array defined from -90 to +90 (instead of +90 to -90), 
C    which is used in NEWDIF only. this won't make a difference since the 
C    latitudes are symmetrical, but we'll change it to be consistent.

!       USE degree_trig

       include "com2d.h"


ccc  idaymn doesn't seem to be used...
ccc	DATA IDAYMN/5,35,65,95,125,155,185,215,245,275,
ccc     *		305,335/

       INTEGER iwwv(40)
       REAL lat45(45), zz76(76)
                                              ! integer values of photolysis wavelength bins
       data iwwv/  0000, 1695, 1724, 1739, 1754, 1770, 1786, 1802, 1818, 
     >             1835, 1852, 1869, 1887, 1905, 1923, 1942, 1961, 1980, 
     >             2000, 2020, 2105, 2198, 2299, 2410, 2532, 2667, 2817, 
     >             2857, 2985, 3030, 3077, 3125, 3175, 3225, 3375, 3575,
     >             3775, 3975, 5475, 7350/


       SAVE


C   set up LATITUDE, PRESSURE, and ZZ ARRAYS for INPUT FILES, 
C        LATIN(18), PRES46(46), ZZ46(46), ZZ30(30), all defined in COMMON, and the 37x59 arrays in STREAMF
c 
        do 121 ij=1,18
 121	   LATIN(IJ) = -85. + (ij-1)*10.

        do 122 ik=1,46
    	   PRES46(IK) = 1013.*exp(-.2844*(ik-.5))
 122	   zz46(ik) = 7.*alog(1013./pres46(ik))     

       do 109 ik=1,59 
          pres59(ik) = 1013.*DEXP(-.2844d0*(ik-1.)) 
 109	  zz59(ik) = 7.*dlog(1013.d0/pres59(ik))      

       do 111 ik=1,58 
          pres58(ik) = 1013.*exp(-.2844*(ik-.5))
 111      zz58(ik) = 7.*alog(1013./pres58(ik))

        do 133 ik=1,30
 133	   zz30(ik) = zz58(ik) 


       do 117 ij=1,37 
 117      xlat5(ij) = (ij-1)*5.-90.


C   set up LATITUDE ARRAY LAT(L$), number of latitudes, dth is the resolution (degrees)
C          dth and L$ are defined in COMMON 
C
c   LATST(L$S) are the Latitudes for the STREAMFUNCTION (the box edges), L$S=180./dthst - both are REAL*8
C     dthst defined in COMMON  -  also define LAT and LATST arrays as REAL*4 in  LAT4(L$) and LATST4(L$S)
C
C   LATEG(L$+1), COSEG(L$+1) are the latitudes and cosines of box edges (sides), REAL*8 in COMMON
C                              LATEG4(L$+1) is REAL*4
C
        do 101 ij=1,L$
    	   LAT(IJ) = -90.d0 + dth/2. + (ij-1)*dth
 101	   LAT4(IJ) = LAT(IJ)

        do 102 ij=1,L$S
    	   LATST(IJ) = -90.d0 + (ij-1)*dthst
 102	   LATST4(IJ) = LATST(IJ)


        do 107 ij=1,L$+1
    	   LATEG(IJ) = -90.d0 + (ij-1)*dth
           LATEG4(ij) = LATEG(IJ)
 107	   COSEG(ij) = DCOSD(lateg(ij))
           

C  compute AREA(L$) array, now based on % of total globe (not hemisphere as was done before)
C    for each specified latitude, first compute the total of the cosines
C         also set up COSC(L$) which are the cosines of the box centers (chemistry grid points)
 
        tcos = 0.0d0
        do 401 ij=1,L$
 401	  tcos = tcos + DCOSD(lat(ij)) 

        do 402 ij=1,L$
          AREA(ij) = DCOSD(lat(ij))/tcos
 402	  COSC(ij) = DCOSD(lat(ij))


C  For altitude and pressure arrays, since we have an irregular grid, 
C     first set up Box edges up to top of extended model,(~115 km), Z$X1=Z$+12+1,  ZALTE(Z$X1), PRESSE(Z$X1)
C
C       deltaz1, inz1, deltaz2, inz2, deltaz3, inz3, deltaz4 are read in CONTROL and are in COMMON
C          the delta's are REAL*8
C
        zalte(1) = 0.0 

        do 200 ii=2,inz1+1
 200	   zalte(ii) = zalte(ii-1) + SNGL(deltaz1)

      if (inz2 .gt. 0) then
        do 201 ii=inz1+2,inz1+inz2+1
 201	   zalte(ii) = zalte(ii-1) + SNGL(deltaz2)
      endif

      if (inz3 .gt. 0) then
        do 211 ii=inz1+inz2+2, inz1+inz2+inz3+1
 211	   zalte(ii) = zalte(ii-1) + SNGL(deltaz3)
      endif

        itz12 = inz1+inz2+inz3+1  
        if (itz12 .lt. Z$+1) then 
          do 202 ii=inz1+inz2+inz3+2,Z$+1 
 202	     zalte(ii) = zalte(ii-1) + SNGL(deltaz4)
        endif

c                                 add on the 12 levels for the extended model (to 115 km) at 2 km resolution
          do 203 ii=Z$+2, Z$X1
 203	     zalte(ii) = zalte(ii-1) + 2.0


C  now do the same thing for ZALTE8(Z$X1) in REAL*8 in COMMON

       zalte8(1) = 0.d0 

        do 700 ii=2,inz1+1
 700	   zalte8(ii) = zalte8(ii-1) + deltaz1

      if (inz2 .gt. 0) then
        do 701 ii=inz1+2,inz1+inz2+1
 701	   zalte8(ii) = zalte8(ii-1) + deltaz2
      endif

      if (inz3 .gt. 0) then
        do 711 ii=inz1+inz2+2, inz1+inz2+inz3+1
 711	   zalte8(ii) = zalte8(ii-1) + deltaz3
      endif

        itz12 = inz1+inz2+inz3+1  
        if (itz12 .lt. Z$+1) then 
          do 702 ii=inz1+inz2+inz3+2,Z$+1 
 702	     zalte8(ii) = zalte8(ii-1) + deltaz4
        endif
c                                 add on the 12 levels for the extended model (to 115 km) at 2 km resolution
          do 703 ii=Z$+2, Z$X1
 703	     zalte8(ii) = zalte8(ii-1) + 2.d0


C  and define the grid points as the center of each box,  ZALT(Z$X), 
c       and for the chemistry grid only (0-90 km), ZALT90(Z$),     
C    ZALT8(Z$X) is the same as ZALT but REAL*8, also PRESS8 is PRESS in REAL*8 (below)

       do 205 ii=1,Z$X
    	  zalt(ii) = (zalte(ii) + zalte(ii+1))/2.
 205	  zalt8(ii) = (zalte8(ii) + zalte8(ii+1))/2.d0

       do 206 ii=1,Z$
 206	  zalt90(ii) = zalt(ii)


C  to set up altitude grid points first, ZALT(Z$) - NOT USED
C
C       zalt(1) = 0.0 + deltaz1/2.
C       do 220 ii=2,inz1+1
C 220	   zalt(ii) = zalt(ii-1) + deltaz1
C       do 221 ii=inz1+2,inz1+inz2+1
C221	   zalt(ii) = zalt(ii-1) + deltaz2
C     itz12 = inz1+inz2+1
C     if (itz12 .lt. Z$) then 
C       do 222 ii=inz1+inz2+2,Z$
C222	   zalt(ii) = zalt(ii-1) + deltaz3
C     endif

C 
C  now set up pressure array, I don't think we need the P(Z$) array
C	DO 104 IK=1,Z$
C 104     P(IK) = -DP*(IK-0.5), where DP=0.2844 IS THE PRESSURE INTERVAL,  LOG(P(I+1))-LOG(P(I))
C
C   PRESS(Z$X) and PRESSE(Z$X1), the pressures of the box edges, PRESSE8(Z$X1) is REAL*8
C   DP(Z$X) array is the DELTAP between box edges, is REAL*8,  PRES90(Z$) is REAL*4 for 0-90 km

	DO 114 IK=1,Z$X1
      		PRESSE8(IK) = 1013.d0*DEXP(-zalte8(ik)/7.d0)
 114		PRESSE(IK) = 1013.*EXP(-zalte(ik)/7.)

	DO 104 IK=1,Z$X
     		PRESS(IK) = 1013.*EXP(-zalt(ik)/7.)
 104		PRESS8(IK) = 1013.*DEXP(-zalt8(ik)/7.)

	DO 105 IK=1,Z$X
 105		DP(IK) = DLOG(presse8(ik+1)/presse8(ik))

        DO 207 ii=1,Z$
 207	  PRES90(ii) = PRESS(ii)

C                                       !  TPRESS8 is REAL*8 in COMMON
        tpress8 = 0.0d0
        do 403 ik=1,Z$X
 403	  tpress8 = tpress8 + PRESS8(IK)

C
C	type *,dp,dth
C

C       INITIALIZE PHYSICAL CONSTANTS.

ccc	DP=0.2844
ccc	DZSTAR=7.E5*DP
	RD=2.87e6
	WTA=28.97e0
	RE=6.371e3
	RCGS=RE*1.e5
	RINV=1.0/RCGS  
	PI=4.*ATAN(1.)
	RPI=1./PI
	DTR=PI/180.
	RTD=1./DTR
	BK=1.380622e-16
	AMU=1.660531e-24
	CP=3.5
	H=6.6262e-27
	CL=3.e10
	OMEGAD=15.
	HBAR=6.79
	ZBAR=30.
	R0=6371.
	YRL=360.
	DT0=86400.
	A=R0*1.e5
c                                      compute area of globe in cm^2
        globe = 4.*pi*(637100000.**2)

C                         DY is latitudinal increment in radians, DTH*DTR
ccc	DY=.1745329252
	DY=DTH*DTR
	DY2=2.*DY
	DPHI=DTH*DTR
	DY1=A*DPHI
ccc	DP2=DP*DP
	DPHI2=DPHI*DPHI

        rr1 = 287.d0


c   define ZSTR altitude array (in KM), and pressures, PRST, for STREAMF and TEMPIN  -  ALL REAL*8

       zstr(1) = 0.d0
       prst(1) = 1013.
       rhost(1) = prst(1)*100./(rr1*239.27) 
       grst(1) = 9.8066*(6371.**2)/((6371. + zstr(1))**2)
c                                                          Note: rhost is the same as in GS83, top of p. 1381
       do 213 ik=2,Z$S
          zstr(ik) = zstr(ik-1) + delz/1000.d0
    	  prst(IK) = 1013.*EXP(-zstr(ik)/7.)
          rhost(ik) = prst(ik)*100./(rr1*239.27) 
 213      grst(ik) = 9.8066*(6371.**2)/((6371. + zstr(ik))**2)


C     also define ZSTR as REAL*4 in ZSTR4(Z$S)

       zstr4(1) = 0.0
       do 214 ik=2,Z$S
 214	  ZSTR4(ik) = ZSTR4(ik-1) + SNGL(delz)/1000.


c
c   ***************   BEGIN LOOP for PDF INDICIES  ****************************
C
c   for temperature/latitude PDFs, need to find proper lat/alt index if not 45x76 grid
C       first define 45x76 lat/hgt gridbox CENTERS,  lat45(45), zz76(76)

        do 1101 ij=1,45
 1101	   lat45(ij) = (ij-1)*4. - 88.

        do 1104 ik=1,60
 1104	   zz76(ik) = float(ik) - .5

        do 1105 ik=61,76
 1105	   zz76(ik) = float(ik-61)*2. + 61.


C  for each current grid latitude LAT4(L$) and altitude ZALT90(Z$), find the closest lat45(45) and zz76(76),
C       and load into ijl(L$), ikz(Z$) (INTEGER arrays in COMMON) - first initialize

       do 1701 ij=1,L$
 1701	   ijl(ij) = 1

       do 901 ij=1,L$
          ij45=1
1550  if (ABS(lat4(ij)-lat45(ij45)).lt.ABS(lat4(ij)-lat45(ij45+1))) then
              ijl(ij) = ij45
          else
              ij45 = ij45 + 1
c                                            ! need this loop if we get to the endpoint
                  if (ij45 .ge. 45) then
                     ijl(ij)=45
                     GO TO 901
                  endif

              GO TO 1550
        endif
 901   continue


        do 601 ik=1,Z$
 601	   ikz(ik) = 1

        do 801 ik=1,Z$
           ik76=1
 2550     if (ABS(zalt90(ik) - zz76(ik76)) .lt. 
     >        ABS(zalt90(ik) - zz76(ik76+1))) then
                 ikz(ik) = ik76
          else
                 ik76 = ik76 + 1
c                                            ! need this loop if we get to the endpoint
                  if (ik76 .ge. 76) then
                     ikz(ik)=76
                     GO TO 801
                  endif

              GO TO 2550
        endif
 801   continue
C
C
C
C  Now for insolation PDFs: define the 181 latitudes, 90S-90N;  LAT181(181) is REAL*8 in COMMON
C                                                               ijind(181) is INTEGER in COMMON
       do 321 ij=1,181
 321	  LAT181(ij) = -90.d0 + DBLE(ij-1)

c  for each lat181, find the appropriate index of LAT(1 -> L$-1) in IJIND(181)

      do 150 ij=1,181
       if (LAT181(ij) .le. LAT(1)) ijind(ij) = 1
       if (LAT181(ij) .ge. LAT(L$-1)) ijind(ij) = L$-1

       if (LAT181(ij) .gt. LAT(1) .and. LAT181(ij) .lt. LAT(L$-1)) then
             ijz = 1
 2000    if (LAT181(ij).ge.LAT(ijz) .and. LAT181(ij).lt.LAT(ijz+1)) then
                ijind(ij) = ijz
             else
                ijz = ijz + 1
                GO TO 2000
          endif
       endif
 150   CONTINUE

c
c   ***************   END LOOP for PDF INDICIES  ****************************
C


c  T is the number of the run in days
	T=0.	
C   FIX TIME OF DAY AT NOON.  USED IN ZENITH ANGLE CALCULATION ONLY.
	CLOCK=12.

C       POSITIONS OF O2 AND O3 IN THE COLUMN DENSITY ARRAY. BOTH ARE
C       ASSUMED TO BE DEFINED.
	ISCO2=1
	ISCO3=2
	JLRT1=1
	JURT1=IL$

	do 100 IJ=1,l$+1
           phih(IJ) = (-90. + 10.*(IJ-1))*pi/180.
100	CONTINUE

	DO 335 IJ=1,L$
c                               phih, phi, cosin, tangent, rcsinv are not used anywhere in the model
c
       		PHI(IJ)=LAT(IJ)*DTR
cccc       		COSIN=COS(PHI(IJ))
cccc       		TANGENT(IJ)=SIN(PHI(IJ))/COSIN
cccc       		RCSINV(IJ)=1./(RCGS*COSIN)
c  zero column and terms for lifetime calculations
	DO 335 IU=1,40
	       	COLUMN(IU,IJ)=0.0

          do 337 ilf=1,10
		CLOSS(IU,ilf,IJ)=0.0
	       	CLOSSST(IU,ilf,IJ)=0.0
 337		CLOSSTR(IU,ilf,IJ)=0.0
335       CONTINUE

          do 230 iu=1,40
             do 230 ij=1,l$
                do 230 ik=1,z$
                   rloss(iu,ij,ik)=0.0
 230      continue

C set surface area of NAT aerosols equal to zero
        do 300 ij=1,l$
           do 300 ik=1,z$
              nataer(ij,ik)=0.0e0
 300          continue

c  SET YSLFPRNT FOR LIFETIME CALCULATION - changed for SPARC lifetime assessment
	YSLFPRNT=.FALSE.

	LIFECHAR(1)='CFCL3     '
	LIFECHAR(2)='CF2CL2    '
	LIFECHAR(3)='CCL4      '
	LIFECHAR(4)='CH3CCL3   '
	LIFECHAR(5)='CHCLF2    '
	LIFECHAR(6)='N2O       '
	LIFECHAR(7)='CH4       '

	LIFECHAR(8)='Halon-1211'
	LIFECHAR(9)='Halon-1301'
	LIFECHAR(10)='C2CL3F3   '
	LIFECHAR(11)='C2CLF5    '
	LIFECHAR(12)='HFC-134a  '
	LIFECHAR(13)='HFC-143a  '
	LIFECHAR(14)='HFC-23    '

	LIFECHAR(15)='C2CL2F4   '
	LIFECHAR(16)='HCFC-141b '
	LIFECHAR(17)='HCFC-142b '
	LIFECHAR(18)='CH3CL     '
	LIFECHAR(19)='CH3BR     '
	LIFECHAR(20)='Halon-1202'
	LIFECHAR(21)='Halon-2402'
	LIFECHAR(22)='HFC-32    '
	LIFECHAR(23)='HFC-125   '
	LIFECHAR(24)='HFC-152a  '
	LIFECHAR(25)='HFC-227ea '
	LIFECHAR(26)='HFC-245fa '

C  Extras:
	LIFECHAR(27)='CO        '
	LIFECHAR(28)='H2        '
	LIFECHAR(29)='HCFC-123  '
	LIFECHAR(30)='CH2Br2    '
	LIFECHAR(31)='CHBr3     '
	LIFECHAR(32)='CO2       '


	DO 9200 II=33,40
		LIFECHAR(II)='XXX       '
9200   CONTINUE


C compute average of 1/photolysis wavelength bins, for heating rate calc.   INTEGER iwwv(40)
C   wvmid(IL$) (in nm) is in COMMON

        wvmid(1) = 1215.67/10.

        do 750 JL=2,IL$
           i1 = iwwv(JL)
           i2 = iwwv(JL+1)
           xnn = REAL(i2 - i1 + 1)
           aaa = 0.
c                                ! convert to nm
           do 755 iww = i1, i2
 755         aaa = aaa + 1./(REAL(iww)/10.)/xnn

           wvmid(JL) = 1./aaa
 750    CONTINUE


	RETURN
	END

	function dcosd(dgr_argument)
	real(4) dcosd
        real(8) dgr_argument
        real(16), parameter :: 
     &	   quadpi = 3.141592653589793238462643383279502884197Q0
	real(16), parameter :: dgr_to_rad = (quadpi/180Q0)
	dcosd = cos(dgr_to_rad * dgr_argument)
	end function
