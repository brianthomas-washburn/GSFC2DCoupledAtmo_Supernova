	SUBROUTINE AEROSOLS

C  NOTE:  only difference between this and the original aerosols.f is the inclusion of the
C         base1 directory name in the read statement.   


        include "com2d.h"


	common/aeros/aerosol_all(2,L$,Z$)

	DIMENSION AEROREAD(2,18,30), aeroin(18,46), aerout(L$,Z$)


	OPEN (unit=21,file='aerosol.dat',status='OLD', form='FORMATTED')

	do 100 ii=1,2
	do 100 ij=1,18
		read(21,6663)(aeroread(ii,ij,ik),ik=1,30)
6663            format(10f6.2)
100	continue

	CLOSE(unit=21)

C 
C   now interpolate AEROREAD(2,18,30) array to new model grid, LATIN(18), ZZ46(46) are REAL*4 in COMMON 
C                                        LAT4(L$), ZALT90(Z$) are REAL*4 in COMMON

	do 500 ii=1,2

           do 601 ij=1,18
             do 602 ik=1,30
 602		aeroin(ij,ik) = aeroread(ii,ij,ik)
             do 603 ik=31,46
 603		aeroin(ij,ik) = 0.0
 601       CONTINUE


            CALL BINTERP(LATIN, 18, zz46, 46, aeroin, 
     >                   LAT4, L$, ZALT90, Z$, 0, 0, AEROUT)


	do 200 ij=1,L$
	  iu=L$+1-ij
	do 300 ik=1,Z$
 300	   aerosol_all(ii,iu,ik) = aerout(ij,ik)
 200	continue

 500	CONTINUE

c  reset to 0.0 above 58 km (original level 30), aerosol_all(2,L$,Z$), ZALT90(Z$)

      DO 110 IK=1,Z$
         IF (zalt90(ik) .GT. 58.) THEN
		DO 111 IJ=1,L$
		DO 111 II=1,2
 111		   AEROSOL_ALL(II,IJ,IK)=0.0
         ENDIF
 110		CONTINUE

	RETURN
	END
