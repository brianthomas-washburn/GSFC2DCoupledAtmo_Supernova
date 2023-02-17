
	SUBROUTINE COLDIUR(ICLOCK)

C
C   NOW, just compute DU in each grid box, using the species value at the grid point, 
C         and the DELTAZ(L$,Z$X) (in KM) value at the species grid points, which is the 
C         difference between the geopotential height of the top edge of the box minus 
C         the geop. hgt. of the bottom edge of the box - EF 2/21/02
C
C
C   this computes column amounts of Ozone and NO from the Diurnally varying fields 
C                                     using REAL*8 CNDC(S$,18+3,L$,Z$) in COMMON

        include "com2d.h"

C                                                 1=Ozone;   2=NO
        DIMENSION INSP(2), topfac(2)

        DATA INSP/4, 5/
        DATA TOPFAC/5.e5, 0.e0/


        SAVE 


C  FIRST DEFINE TOTAL COLUMN ABOVE TOP MODEL LEVEL, now use DELTAZ(L$,Z$X) (in KM) 
C                 use CORRECTED for integration below, NCOLD(2,18+3,L$,Z$) is in COMMON

      DO 1000 IP=1,2
        DO 2000 IJ=1,L$
          SPCOL = CNDC(INSP(ip),iclock,IJ,Z$)*topfac(ip)

c                 integrate top-down - SPCOL is the column at the CENTER of the box
C                 ie, DELTAZ/2., store 1/2 box in NCOLD, then add in the other 1/2 to SPCOL
C
          DO 4000 IK = Z$,1,-1
            colcen = CNDC(INSP(ip),iclock,IJ,IK)*DELTAZ(IJ,IK)*1.D5/2.
            SPCOL = SPCOL + colcen

            NCOLD(IP,iclock,IJ,IK) = SPCOL
            SPCOL = SPCOL + colcen
 4000	  CONTINUE

 2000	CONTINUE
 1000	CONTINUE

      RETURN
      END
