C  PUT IN NOX SOURCE OF GCRS (BOTH MIN AND MAX)  4/17/90

	subroutine GCRS


        include "com2d.h"


        DIMENSION gcrout(L$,Z$)


	OPEN(UNIT=96, FILE='F40GCRMM.DAT',STATUS='OLD',
     c		FORM='FORMATTED')
 
c  data is symmetric, so read in for northern hemisphere only,  GCRMIN(18,46),GCRMAX(18,46) in COMMON

	DO 2000 IJ=10,18
		READ(96,81)
		READ(96,80)(GCRMIN(IJ,IK),IK=1,30)
		READ(96,81)
2000	CONTINUE
80	FORMAT(6E11.3)
81	FORMAT(1X)
	DO 2015 IJ=10,18
		READ(96,81)
		READ(96,80)(GCRMAX(IJ,IK),IK=1,30)
		READ(96,81)
2015	CONTINUE
	CLOSE(UNIT=96)
c  set southern hemisphere values equal to northern hemisphere data

	DO 2001 IK=1,30
	DO 2001 IJ=1,9
		JJ=18-IJ+1
		GCRMIN(IJ,IK)=GCRMIN(JJ,IK)
		GCRMAX(IJ,IK)=GCRMAX(JJ,IK)
2001       CONTINUE


	DO 2002 IK=31,46
	DO 2002 IJ=1,18
		GCRMIN(IJ,IK)=0.0
 2002		GCRMAX(IJ,IK)=0.0


c  select solar max gcr's or solar min, interpolate to new model grid, GCR(L$,Z$), gcrout(L$,Z$)
C                            LATIN(18), ZZ46(46),  LAT4(L$), ZALT90(Z$) are REAL*4 in COMMON 
C 

            CALL BINTERP(LATIN, 18, zz46, 46, GCRMAX, 
     >                   LAT4, L$, ZALT90, Z$, 0, 0, GCROUT)


	DO 2003 IK=1,Z$
	DO 2003 IJ=1,L$
		GCR(IJ,IK)=GCROUT(IJ,IK)
2003	CONTINUE


c  zero any other values above 60 km

      DO 110 IK=1,Z$
         IF (zalt90(ik) .GT. 60.) THEN
                DO 111 IJ=1,L$
 111               GCR(ij,ik) = 0.0
         ENDIF
 110  CONTINUE
	
	RETURN
	END
